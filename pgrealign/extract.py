"""Stage 2, `extract`: every walk's run through every region, in one streaming pass over a contig's graph.

A region is the interior between two anchor nodes L and R (regions.py). Anchors separate the graph, so
each walk relates to a region in one of these ways, read off the order of its visits to L and R:

| run | visits | interior |
|---|---|---|
| spanning, forward | >L ... >R | the steps between |
| spanning, reverse | <R ... <L | the steps between, reversed and flipped |
| prefix | the walk starts inside, then >R (or <L) | the steps before the anchor |
| suffix | >L (or <R), then the walk ends inside | the steps after the anchor |
| internal | the walk never visits an anchor but uses interior nodes | the whole walk |
| complex | any other order (a U-turn, a second entry before leaving) | - |

Every interior is written in the region's forward orientation (the reference's). `orient` records
whether the walk ran the other way.

A fragment's **side** is the end of the region it touches, in that forward orientation, which depends
on both its type and its orient:
- `left` touches L: a forward suffix or a reverse prefix.
- `right` touches R: a forward prefix or a reverse suffix.
- `internal` touches neither.

Fragments are kept by side, which is what realign and replace need.
- A region with a complex run is `complex` and is left as it is.
- A region whose spanning runs all spell one sequence, with no fragments, is `invariant`.
- Every other region is `ok`.

Outputs, in the output directory:
- `packages.jsonl.gz`: one JSON object per region:
  - its anchors and status;
  - `alleles`: the distinct spanning sequences, with weights (their number of runs);
  - `fragments`: the distinct fragment sequences, with side and weight.
- `occurrences.tsv.gz`: one row per run: walk, region, type, orient, the run's interior as step
  indices [first, end) in the W line, and the allele or fragment it spells. `replace` rewrites these.
- `interiors.tsv.gz`: region, interior node ids. These are the nodes `replace` removes.
- `summary.json`.
"""
import collections
import gzip
import json
import os

from .graph import GfaStream, parse_w, revcomp, walk_key, STEP


def flip(steps):
    return [('<' if o == '>' else '>', n) for o, n in reversed(steps)]


def spell(steps, seq):
    return ''.join(seq[n] if o == '>' else revcomp(seq[n]) for o, n in steps)


class RegionRuns:
    """Everything one region collects during the pass."""
    __slots__ = ('alleles', 'fragments', 'interior', 'complex', 'n_runs')

    def __init__(self):
        self.alleles = {}       # sequence -> [id, weight]
        self.fragments = {}     # (type, sequence) -> [id, weight]
        self.interior = set()
        self.complex = []
        self.n_runs = collections.Counter()


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
    ('complex', reason, None, None)."""
    state, s, seen = None, None, False
    for i, side, o in occ:
        if state is None:
            if side == 'L' and o == '>':
                state, s = 'fwd', i
            elif side == 'R' and o == '<':
                state, s = 'rev', i
            elif not seen:
                # the walk reaches an anchor from inside without having entered: it started inside
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
    if state is not None:
        yield ('suffix', '+' if state == 'fwd' else '-', s + 1, len(toks))


def extract(graph, regions, out_dir, vg=None, threads=4, max_internal_bp=200000):
    """Run the pass over `graph` (a GBZ or GFA) for `regions` (regions.Region with status 'ok')."""
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
    candidates = []       # short walks with no anchor visit: possible internal fragments
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
                if len(toks) and sum(len(seq[int(n)]) for _, n in toks) <= max_internal_bp:
                    candidates.append((key, toks))
                continue
            for rid, occ in hits.items():
                rr = runs[rid]
                for typ, orient, a, b in runs_of_walk(toks, occ):
                    if typ == 'complex':
                        rr.complex.append('%s: %s' % (key, orient))
                        continue
                    steps = [(o, int(n)) for o, n in toks[a:b]]
                    rr.interior.update(n for _, n in steps)
                    s = spell(steps if orient == '+' else flip(steps), seq)
                    if typ == 'spanning':
                        mid = _member(rr.alleles, s, 'a')
                    else:
                        mid = _member(rr.fragments, (side_of(typ, orient), s), 'f')
                    rr.n_runs[typ] += 1
                    occ_out.write('%s\t%s\t%s\t%s\t%d\t%d\t%s\n' % (key, rid, typ, orient, a, b, mid))
    # internal fragments: short anchor-free walks that use a region's interior nodes
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
        rr.interior.update(n for _, n in steps)
        mid = _member(rr.fragments, ('internal', spell(steps, seq)), 'f')
        rr.n_runs['internal'] += 1
        occ_out.write('%s\t%s\tinternal\t+\t0\t%d\t%s\n' % (key, rid, len(toks), mid))
    occ_out.close()
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
                   'n_runs': dict(rr.n_runs), 'complex': rr.complex[:20],
                   'alleles': [{'id': v[0], 'weight': v[1], 'seq': s} for s, v in rr.alleles.items()],
                   'fragments': [{'id': v[0], 'weight': v[1], 'side': k[0], 'seq': k[1]}
                                 for k, v in rr.fragments.items()]}
            pk.write(json.dumps(pkg, separators=(',', ':')) + '\n')
            it.write('%s\t%s\n' % (rid, ','.join(str(n) for n in sorted(rr.interior))))
    summary = {'walks': n_walks, 'regions': len(regs), 'status': dict(status),
               'internal_candidates': len(candidates)}
    with open(os.path.join(out_dir, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=1)
    return summary


def read_packages(path, status=None):
    """Region packages from packages.jsonl.gz, optionally only those with one of the given statuses."""
    with gzip.open(path, 'rt') as f:
        for line in f:
            p = json.loads(line)
            if status is None or p['status'] in status:
                yield p
