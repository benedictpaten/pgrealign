"""Stage 2, `extract`: every walk's run through every region, in one streaming pass over a contig's graph.

A region is the interior between two anchor nodes L and R (regions.py). Anchors separate the graph, so
each walk relates to a region in one of these ways, read off the order of its visits to L and R:

| run | visits | interior |
|---|---|---|
| spanning, forward | >L ... >R | the steps between |
| spanning, reverse | <R ... <L | the steps between, reversed and flipped |
| prefix | the walk starts inside, then >R (or <L) | the steps before the anchor |
| suffix | >L (or <R), then the walk ends inside | the steps after the anchor |
| internal | the walk never visits an anchor but uses interior nodes | the whole walk, flipped if it runs backward |
| complex | any other order (a U-turn, a second entry before leaving) | - |

A spanning run can be empty: L and R are then adjacent, and it spells the empty allele. A prefix or
suffix cannot: a walk that starts or ends on an anchor has no step inside at that end, so it has no run
there.

Every interior is written in the region's forward orientation (the reference's). `orient` records
whether the walk ran the other way.

**Orientation.** Turned to the region's forward orientation, every run should visit each interior node
the same way. A node stored reversed, such as an insertion first added from a reverse-strand contig, is
visited `<` by every run, which is fine. A node that runs visit both ways is an inversion between
haplotypes: the reference crosses its nodes forward and an inverted haplotype crosses them backward.
Realign would spell the inverted allele as a plain forward sequence and lose the reversal, so the
region is made complex. An anchor-free walk has no anchor to orient it, so it takes its orientation from
the interior nodes the runs visit: one that visits them all the other way runs against the region and
is flipped, like a reverse run, and one that agrees with some and not others makes the region complex.

A fragment's **side** is the end of the region it touches, in that forward orientation, which depends
on both its type and its orient:
- `left` touches L: a forward suffix or a reverse prefix.
- `right` touches R: a forward prefix or a reverse suffix.
- `internal` touches neither.
Fragments are kept by side, which is what realign and replace need.

**Anchor-free walks.** The region an anchor-free walk belongs to is known only once the pass has found
every interior, so the walk is held until then. It is held whole if it spells at most
`max_internal_bp`, and otherwise only as the nodes it visits (one entry per node, however many such
walks there are). A longer walk is never recorded as a run, and replace stops at a walk that visits a
replaced interior outside a recorded run. So a longer walk that uses an interior makes the region
complex.

**Status.** Each region has one:
- `complex`: a walk's run is complex, runs visit an interior node both ways (above), an anchor-free walk
  uses the interiors of more than one region or runs both ways through one, or an anchor-free walk
  longer than `max_internal_bp` uses its interior.
- `no_runs`: no walk has a run in it.
- `invariant`: there are no fragments, and every spanning run spells one sequence.
- `ok`: every other region.

Only `ok` regions go on to realign and replace. The others keep their subgraph.

Outputs, in the output directory:
- `packages.jsonl.gz`: one JSON object per region, in reference order:
  - `id`, `contig`, `L`, `R`, `ref_start`, `ref_end`, `targets`, `classes`: its row in the regions
    file;
  - `status`, as above;
  - `n_runs`: {type: number of runs};
  - `complex`: the reasons it is complex, each as 'walk: reason', at most 20;
  - `alleles`: the distinct spanning sequences, each with an id (a1, a2, ...), `seq`, and a `weight`,
    its number of runs;
  - `fragments`: the distinct (side, sequence) pairs, each with an id (f1, f2, ...), `side`, `seq` and
    `weight`;
  - `ref`: the id of the allele that a spanning run of the reference sample spells (the last, if it
    has several), or null if it has none.
- `occurrences.tsv.gz`: one row per run: walk, region, type, orient, the run's interior as step
  indices [first, end) in the W line, and the allele or fragment it spells. `replace` rewrites these.
- `interiors.tsv.gz`: region, and the ids of the nodes its runs visit. These are the nodes `replace`
  removes.
- `summary.json`: the walks and regions read, the regions of each status, the anchor-free walks held
  whole (`internal_candidates`), the regions an anchor-free walk longer than `max_internal_bp`
  made complex (`internal_too_long`), the regions whose runs visit a node both ways (`both_ways`), and
  the anchor-free walks flipped to run with their region (`internal_flipped`).
"""
import collections
import gzip
import json
import os

from .graph import GfaStream, parse_w, spell, walk_key, STEP


def flip(steps):
    return [('<' if o == '>' else '>', n) for o, n in reversed(steps)]


FWD, REV = 1, 2         # the ways a run visits an interior node, in the region's forward orientation


class RegionRuns:
    """Everything one region collects during the pass."""
    __slots__ = ('alleles', 'fragments', 'interior', 'complex', 'n_runs', 'ref', 'both_ways')

    def __init__(self):
        self.ref = None         # the allele the reference sample's walk spells
        self.alleles = {}       # sequence -> [id, weight]
        self.fragments = {}     # (side, sequence) -> [id, weight]
        self.interior = {}      # node -> the ways runs visit it, FWD | REV
        self.complex = []
        self.n_runs = collections.Counter()
        self.both_ways = False  # some node is visited both ways

    def visit(self, key, steps):
        """Record the interior nodes of walk `key`'s run, `steps` in the region's forward orientation. A
        run that visits a node the other way from an earlier run, or from itself, makes the region
        complex."""
        turned = None
        for o, n in steps:
            way = FWD if o == '>' else REV
            had = self.interior.get(n, 0)
            if had and not had & way and turned is None:
                turned = n
            self.interior[n] = had | way
        if turned is not None:
            self.both_ways = True
            self.complex.append('%s: visits node %d both ways, as at an inversion' % (key, turned))


def _member(table, key, prefix):
    m = table.get(key)
    if m is None:
        m = table[key] = ['%s%d' % (prefix, len(table) + 1), 0]
    m[1] += 1
    return m[0]


def side_of(typ, orient):
    if typ == 'internal':
        return 'internal'
    return 'left' if (typ == 'suffix') == (orient == '+') else 'right'


def runs_of_walk(toks, occ):
    """The runs of one walk through one region, from its anchor visits `occ` = [(index, 'L'|'R', orient)].
    Yields (type, orient, first, end), [first, end) being the interior's step indices in the walk, or
    ('complex', reason, None, None). A prefix or suffix is never empty."""
    state, s, seen = None, None, False
    for i, side, o in occ:
        if state is None:
            if side == 'L' and o == '>':
                state, s = 'fwd', i
            elif side == 'R' and o == '<':
                state, s = 'rev', i
            elif not seen:
                # the walk reaches an anchor from inside without having entered: it started inside,
                # or on the anchor itself
                if i > 0:
                    yield ('prefix', '+' if side == 'R' else '-', 0, i)
            else:
                yield ('complex', 'leaves through %s%s without entering' % (o, side), None, None)
                return
        elif state == 'fwd':
            if side == 'R' and o == '>':
                yield ('spanning', '+', s + 1, i)
                state = None
            else:
                yield ('complex', 'entered by >L, then %s%s' % (o, side), None, None)
                return
        else:  # rev
            if side == 'L' and o == '<':
                yield ('spanning', '-', s + 1, i)
                state = None
            else:
                yield ('complex', 'entered by <R, then %s%s' % (o, side), None, None)
                return
        seen = True
    if state is not None and s + 1 < len(toks):
        yield ('suffix', '+' if state == 'fwd' else '-', s + 1, len(toks))


def extract(graph, regions, out_dir, vg=None, threads=4, max_internal_bp=200000, ref_sample='CHM13'):
    """Run the pass over `graph` (a GBZ or GFA) for the `regions` (regions.Region) with status 'ok'; the
    others get no package. A package's `ref` is the allele that a walk of `ref_sample` spells. An
    anchor-free walk longer than `max_internal_bp` that uses an interior makes its region complex."""
    os.makedirs(out_dir, exist_ok=True)
    regs = {r.id: r for r in regions if r.status == 'ok'}
    anchor_of = collections.defaultdict(list)     # node -> [(region id, 'L'|'R')]
    for r in regs.values():
        anchor_of[str(r.L)].append((r.id, 'L'))
        anchor_of[str(r.R)].append((r.id, 'R'))
    runs = {rid: RegionRuns() for rid in regs}
    seq = {}
    occ_out = gzip.open(os.path.join(out_dir, 'occurrences.tsv.gz'), 'wt')
    occ_out.write('#walk\tregion\ttype\torient\tfirst\tend\tmember\n')
    candidates = []       # anchor-free walks of at most max_internal_bp: possible internal fragments
    long_walks = []       # the keys of the longer anchor-free walks
    long_at = {}          # node -> index in long_walks of the first longer walk to visit it
    n_walks = 0
    with GfaStream(graph, vg, threads) as f:
        for line in f:
            t = line[0]
            if t == 'S':
                x = line.split('\t', 3)
                seq[int(x[1])] = x[2].rstrip('\n')
                continue
            if t != 'W':
                continue
            n_walks += 1
            w = parse_w(line)
            toks = STEP.findall(w.steps)
            key = walk_key(w)
            hits = collections.defaultdict(list)
            for i, (o, n) in enumerate(toks):
                a = anchor_of.get(n)
                if a is not None:
                    for rid, side in a:
                        hits[rid].append((i, side, o))
            if not hits:
                if regs and toks:
                    if sum(len(seq[int(n)]) for _, n in toks) <= max_internal_bp:
                        candidates.append((key, toks))
                    else:
                        for _, n in toks:
                            long_at.setdefault(int(n), len(long_walks))
                        long_walks.append(key)
                continue
            for rid, occ in hits.items():
                rr = runs[rid]
                for typ, orient, a, b in runs_of_walk(toks, occ):
                    if typ == 'complex':
                        rr.complex.append('%s: %s' % (key, orient))
                        continue
                    steps = [(o, int(n)) for o, n in toks[a:b]]
                    if orient == '-':
                        steps = flip(steps)
                    rr.visit(key, steps)
                    s = spell(steps, seq)
                    if typ == 'spanning':
                        mid = _member(rr.alleles, s, 'a')
                        if w.sample == ref_sample:
                            rr.ref = mid
                    else:
                        mid = _member(rr.fragments, (side_of(typ, orient), s), 'f')
                    rr.n_runs[typ] += 1
                    occ_out.write('%s\t%s\t%s\t%s\t%d\t%d\t%s\n' % (key, rid, typ, orient, a, b, mid))
    # internal fragments: the held anchor-free walks that use a region's interior nodes
    n_flipped = 0
    node_region = {}
    for rid, rr in runs.items():
        for n in rr.interior:
            node_region[n] = rid
    for key, toks in candidates:
        rids = {node_region.get(int(n)) for _, n in toks} - {None}
        if not rids:
            continue
        if len(rids) > 1:
            for rid in rids:
                runs[rid].complex.append('%s: an anchor-free walk through %d regions' % (key, len(rids)))
            continue
        rid = rids.pop()
        rr = runs[rid]
        steps = [(o, int(n)) for o, n in toks]
        # the ways the runs visit the interior nodes this walk visits, against the ways it does
        known = [(FWD if o == '>' else REV, rr.interior[n]) for o, n in steps if n in rr.interior]
        if all(had & way for way, had in known):
            orient = '+'
        elif all(had & ((FWD | REV) ^ way) for way, had in known):
            orient, steps = '-', flip(steps)
            n_flipped += 1
        else:
            rr.complex.append('%s: an anchor-free walk runs both ways through the region' % key)
            continue
        rr.visit(key, steps)
        mid = _member(rr.fragments, ('internal', spell(steps, seq)), 'f')
        rr.n_runs['internal'] += 1
        occ_out.write('%s\t%s\tinternal\t%s\t0\t%d\t%s\n' % (key, rid, orient, len(toks), mid))
    occ_out.close()
    # the longer anchor-free walks have no run, so a region whose interior one uses is left as it is
    n_long = 0
    if long_at:
        for rr in runs.values():
            ks = sorted({long_at[n] for n in rr.interior if n in long_at})
            for k in ks:
                rr.complex.append('%s: an anchor-free walk longer than max_internal_bp (%d)'
                                  % (long_walks[k], max_internal_bp))
            n_long += bool(ks)
    status = collections.Counter()
    with gzip.open(os.path.join(out_dir, 'packages.jsonl.gz'), 'wt') as pk, \
            gzip.open(os.path.join(out_dir, 'interiors.tsv.gz'), 'wt') as it:
        it.write('#region\tnodes\n')
        for rid in sorted(regs, key=lambda k: (regs[k].ref_start, regs[k].ref_end)):
            r, rr = regs[rid], runs[rid]
            if rr.complex:
                st = 'complex'
            elif not rr.alleles and not rr.fragments:
                st = 'no_runs'
            elif len(rr.alleles) <= 1 and not rr.fragments:
                st = 'invariant'
            else:
                st = 'ok'
            status[st] += 1
            pkg = {'id': rid, 'contig': r.contig, 'L': r.L, 'R': r.R, 'ref_start': r.ref_start,
                   'ref_end': r.ref_end, 'targets': r.targets, 'classes': r.classes, 'status': st,
                   'n_runs': dict(rr.n_runs), 'complex': rr.complex[:20], 'ref': rr.ref,
                   'alleles': [{'id': v[0], 'weight': v[1], 'seq': s} for s, v in rr.alleles.items()],
                   'fragments': [{'id': v[0], 'weight': v[1], 'side': k[0], 'seq': k[1]}
                                 for k, v in rr.fragments.items()]}
            pk.write(json.dumps(pkg, separators=(',', ':')) + '\n')
            it.write('%s\t%s\n' % (rid, ','.join(str(n) for n in sorted(rr.interior))))
    summary = {'walks': n_walks, 'regions': len(regs), 'status': dict(status),
               'internal_candidates': len(candidates), 'internal_too_long': n_long,
               'both_ways': sum(rr.both_ways for rr in runs.values()), 'internal_flipped': n_flipped}
    path = os.path.join(out_dir, 'summary.json')
    with open(path + '.tmp', 'w') as f:
        json.dump(summary, f, indent=1)
    os.replace(path + '.tmp', path)
    return summary


def read_packages(path, status=None):
    """Region packages from packages.jsonl.gz, optionally only those with one of the given statuses."""
    with gzip.open(path, 'rt') as f:
        for line in f:
            p = json.loads(line)
            if status is None or p['status'] in status:
                yield p
