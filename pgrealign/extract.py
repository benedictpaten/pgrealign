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
With `allow_inversions`, neither makes the region complex: an inverted allele is kept as the forward
sequence it spells, and an anchor-free walk takes the orientation that agrees with more of the nodes.

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

**Exit edges** (`exit_edges`). Anchors that walks go around (regions' `max_bypass`) do not separate
the graph: a haplotype whose sequence of the region never aligned to the reference can leave it before L
and come back inside, or leave inside and come back after R. With exit edges, a run need not begin or end
at an anchor: each walk's steps near the region are classed by the reference index (prepare's
<contig>.ref.tsv.gz) as before L, L, between the anchors (interior), R, after R, or off the reference.
A stretch of steps off the reference belongs to the region when a step that bounds it is on the
interior's side: an interior reference step, or an anchor that the walk crosses into or out of the
interior. A run is a maximal stretch of a walk's steps in the region, and is bounded by reference steps
outside the interior, or by the walk's ends (exit_runs_of_walk):
- spanning: in on one side (at L or before it) and out on the other, forward or reverse. These are
  the alleles, whether a bound is an anchor or another reference node (an exit);
- prefix and suffix: the walk starts or ends in the region;
- uturn: in and out on the same side, a fragment on that side;
- internal: the walk lies in the region.
A walk that leaves the reference before L and comes back after R without a step in the region has no
run. The region's interior is every node its runs visit; those that some walk also visits outside every
run stay in the graph for it (`kept`), so replace removes only the others. A node in the interiors of
two regions makes both complex.

**Keeping fragments** (`keep_fragments`, with exit edges). Only the spanning runs are recorded. Every
other run, and every anchor-free walk, counts as outside every run: replace leaves its walk as it is,
on its nodes, which stay with the edges between them. The interior is then the nodes of the alleles, and
those of its nodes that fragments visit are `kept`. A region's fragments are not realigned, but it
gains no unaligned copy of them either.

**Trimming to the targets** (`trim_targets`, with keep_fragments). A region's anchors can lie far from
its targets, so its runs can carry much more than the targets. With trim_targets, the targets file that
regions read, each run is cut down to the haplotype's own copy of each of the region's targets, its
core, found by sequence (find_cores):
- A target's k-mers are those that occur at least twice in the reference's sequence of it, in no other
  target of the region, and `CORE_ENRICHMENT` times more densely in it than in the reference's sequence
  of the region outside its targets.
- A run's k-mers are sampled, and counted in windows. A window is dense for a target when enough of its
  samples are that target's k-mers, and dense windows close together make a stretch.
- A run's core for a target is its stretch with the most hits, out to the hits in the windows on either
  side, then by `trim_margin` on each side, within the run, and out to whole steps (cut_steps). A step
  that would add more than `trim_margin` again, such as a long node of a contig's unaligned sequence, is
  left out unless the core starts or ends in it. Overlapping cores of two targets are one core.
- A core with a step of its walk on either side is recorded as a spanning run, an allele, whatever the
  type of the run it is in. A core that reaches a walk's end may be a truncated copy, and is not.
- A core is read in the orientation in which its k-mers read as the reference's, so its allele is in
  the region's forward orientation. A core of a spanning run whose k-mers give no clear orientation
  takes its run's. One of another run is not recorded.
Everything else stays on its nodes, as a fragment does: the rest of a run with a core, and every run
without one. So a fragment that holds a whole copy of a target, such as a contig that begins between
the target and its anchor, is realigned too.

**Status.** Each region has one:
- `complex`: a walk's run is complex, runs visit an interior node both ways (above, unless
  `allow_inversions`), an anchor-free walk uses the interiors of more than one region or runs both ways
  through one (unless `allow_inversions`), or an anchor-free walk longer than `max_internal_bp` uses its
  interior.
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
  - with exit edges, `exit_edges` (true), `exit_runs` (the runs with a bound off the anchors),
    `kept_nodes` and `kept_bp` (the interior nodes that stay, as walks visit them outside every run);
  - with keep_fragments, `kept_fragments` and `kept_fragment_bp`: the runs left on their nodes, and the
    bp they spell, counted per run;
  - with trim_targets, `trim`: `margin`, the cores recorded (`cores`), those of them in runs other
    than spanning ones (`cores_in_fragments`), the spanning runs without a core (`no_core`), the cores
    not recorded as they reach a walk's end (`at_end`) or have no orientation (`no_orient`), and the bp
    of the runs with a recorded core that stay on their nodes (`trimmed_bp`). The runs without a core
    count in kept_fragments.
- `occurrences.tsv.gz`: one row per run: walk, region, type, orient, the run's interior as step
  indices [first, end) in the W line, and the allele or fragment it spells; with exit edges, also the
  steps before and after the run in its walk ('.' at a walk's end), which replace links it to.
  `replace` rewrites these.
- `interiors.tsv.gz`: region, and the ids of the nodes its runs visit. These are the nodes `replace`
  removes, except those of a third column, with exit edges, which stay.
- with trim_targets, `cores.tsv.gz`: one row per core found, and per spanning run without one: walk,
  region, targets, the type of the run it is in, status (`core`, `at_end`, `no_orient` or `no_core`),
  orient, the core as step indices [first, end), contig, and in the contig's coordinates the core's
  stretch of hits and the core.
- `summary.json`: the walks and regions read, the regions of each status, the anchor-free walks held
  whole (`internal_candidates`), the regions an anchor-free walk longer than `max_internal_bp`
  made complex (`internal_too_long`), the regions whose runs visit a node both ways (`both_ways`), and
  the anchor-free walks flipped to run with their region (`internal_flipped`).
"""
import bisect
import collections
import gzip
import itertools
import json
import os
from array import array

from .graph import GfaStream, parse_w, revcomp, spell, walk_key, STEP

# trimming to the targets (find_cores)
CORE_K = 21             # k-mer length
CORE_STRIDE = 32        # a run's k-mers are sampled one every this many bases
CORE_WINDOW = 10000     # in windows of this many bases
CORE_DENSITY = 0.3      # a window is dense for a target when at least this share of its samples are its k-mers
CORE_GAP = 200000       # dense windows of a target at most this far apart make one stretch
CORE_ENRICHMENT = 10    # a target's k-mer is at least this many times denser in it than in the rest of the region


class TargetKmers:
    """The k-mers by which find_cores finds a region's targets (see Trimming to the targets): `fwd` and
    `rev` map each k-mer, and its reverse complement, to its target's index in `ids`. `targets` is
    [(id, the reference's sequence of it)] and `outside` [the reference's sequences of the region outside
    its targets]. A k-mer that is at most `enrichment` times as dense in its target as in `outside` is
    left out, so that a satellite of the flanks that shares the target's k-mers is not taken for it, but
    a few copies of the target's repeat past its end, as the edge of an array can have, do not count."""

    def __init__(self, targets, outside, k=CORE_K, enrichment=CORE_ENRICHMENT):
        self.k, self.ids = k, [t for t, _ in targets]
        owner, bad = {}, set()      # k-mer -> (target, count in it)
        for ti, (_, s) in enumerate(targets):
            s = s.upper()
            counts = collections.Counter(s[i:i + k] for i in range(len(s) - k + 1))
            for km, c in counts.items():
                if c < 2 or 'N' in km:
                    continue
                if owner.setdefault(km, (ti, c))[0] != ti:
                    bad.add(km)
        out_counts = collections.Counter()
        out_bp = sum(len(s) for s in outside)
        for s in outside:
            s = s.upper()
            for i in range(len(s) - k + 1):
                km = s[i:i + k]
                if km in owner:
                    out_counts[km] += 1
        target_bp = [max(1, len(s)) for _, s in targets]
        for km, n in out_counts.items():
            ti, c = owner[km]
            if c / target_bp[ti] <= enrichment * n / out_bp:
                bad.add(km)
        fwd = {km: ti for km, (ti, _) in owner.items() if km not in bad}
        # a k-mer whose reverse complement is one too reads both ways, so says nothing of orientation
        self.fwd = {km: ti for km, ti in fwd.items() if revcomp(km) not in fwd}
        self.rev = {revcomp(km): ti for km, ti in self.fwd.items()}


def find_cores(seq, tk, stride=CORE_STRIDE, window=CORE_WINDOW, density=CORE_DENSITY, gap=CORE_GAP):
    """Each target's core in `seq` (upper case), from the k-mers of `tk` (TargetKmers), before any
    margin: [(target index, start, end, stretch start, stretch end, strand)], [start, end) out to the
    hits in the windows on either side of the stretch, and strand '+' or '-' when at least two thirds of
    those hits read that way, else None."""
    k = tk.k
    last = len(seq) - k
    if last < 0:
        return []
    fwd, rev = tk.fwd, tk.rev
    hits = collections.defaultdict(list)    # target -> [(position, +1 forward or -1 reverse)], in order
    for i in range(0, last + 1, stride):
        km = seq[i:i + k]
        t = fwd.get(km)
        if t is not None:
            hits[t].append((i, 1))
            continue
        t = rev.get(km)
        if t is not None:
            hits[t].append((i, -1))

    def samples(w):
        lo, hi = w * window, min((w + 1) * window, last + 1)
        return -(-hi // stride) - -(-lo // stride)

    out = []
    for t, hs in sorted(hits.items()):
        per = collections.Counter(p // window for p, _ in hs)
        dense = sorted(w for w, c in per.items() if c >= density * samples(w))
        if not dense:
            continue
        stretches = [[dense[0], dense[0]]]
        for w in dense[1:]:
            if (w - stretches[-1][1] - 1) * window <= gap:
                stretches[-1][1] = w
            else:
                stretches.append([w, w])
        pos = [p for p, _ in hs]

        def span(lo, hi):
            return bisect.bisect_left(pos, lo), bisect.bisect_left(pos, hi)
        best = None
        for s0, s1 in stretches:
            a, b = span(s0 * window, (s1 + 1) * window)
            if best is None or b - a > best[0]:
                best = (b - a, s0, s1)
        _, w0, w1 = best
        a, b = span(w0 * window, (w1 + 1) * window)
        sa, sb = pos[a], pos[b - 1] + k
        a, b = span((w0 - 1) * window, (w1 + 2) * window)
        f = sum(1 for _, s in hs[a:b] if s > 0)
        r = (b - a) - f
        strand = '+' if f >= 2 * r else '-' if r >= 2 * f else None
        out.append((t, pos[a], pos[b - 1] + k, sa, sb, strand))
    return out


def flip(steps):
    return [('<' if o == '>' else '>', n) for o, n in reversed(steps)]


FWD, REV = 1, 2         # the ways a run visits an interior node, in the region's forward orientation


class RegionRuns:
    """Everything one region collects during the pass."""
    __slots__ = ('alleles', 'fragments', 'interior', 'complex', 'n_runs', 'ref', 'both_ways', 'allow_inversions',
                 'exit_runs', 'kept_runs', 'kept_run_bp', 'trim')

    def __init__(self, allow_inversions=False):
        self.ref = None         # the allele the reference sample's walk spells
        self.alleles = {}       # sequence -> [id, weight]
        self.fragments = {}     # (side, sequence) -> [id, weight]
        self.interior = {}      # node -> the ways runs visit it, FWD | REV
        self.complex = []
        self.n_runs = collections.Counter()
        self.both_ways = False  # some node is visited both ways
        self.allow_inversions = allow_inversions
        self.exit_runs = 0      # exit edges: the runs with a bound off the anchors
        self.kept_runs = 0      # keep_fragments: the fragment runs left on their nodes, and their bp
        self.kept_run_bp = 0
        self.trim = collections.Counter()   # trim_targets: the counts of the package's `trim`

    def visit(self, key, steps):
        """Record the interior nodes of walk `key`'s run, `steps` in the region's forward orientation. A
        run that visits a node the other way from an earlier run, or from itself, makes the region
        complex, unless allow_inversions."""
        turned = None
        for o, n in steps:
            way = FWD if o == '>' else REV
            had = self.interior.get(n, 0)
            if had and not had & way and turned is None:
                turned = n
            self.interior[n] = had | way
        if turned is not None:
            self.both_ways = True
            if not self.allow_inversions:
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
    if typ == 'uturn':
        # exit edges: a run in and out on one side, '+' for the left (exit_runs_of_walk)
        return 'left' if orient == '+' else 'right'
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


IN, OUT = 1, 2      # which side of a region a stretch's bound leaves it on


def _bound(kind, o, leaving):
    """Where the bound of a stretch puts the stretch: IN (the region's interior side), OUT, None (a walk
    end) or 'x' (a reference node the reference visits both inside and outside the region). `kind` is
    the bound's class (exit_runs_of_walk) and `o` its orientation; `leaving` says that the walk leaves
    the bound into the stretch, rather than arriving at it from the stretch."""
    if kind == 'i':
        return IN
    if kind == 'L':
        return IN if (o == '>') == leaving else OUT
    if kind == 'R':
        return IN if (o == '<') == leaving else OUT
    if kind == 'x':
        return 'x'
    return OUT


def exit_runs_of_walk(toks, lo, hi, kinds):
    """The runs of one walk through one region with exit edges, from the steps [lo, hi) of the walk
    `toks`. `kinds[i - lo]` is the class of step i: None for a node not on the reference, 'L' or 'R' for
    an anchor, 'i' for a reference node between them, 'o' or 'O' for a reference node before L or after
    R, 'x' for a reference node the reference visits both inside and outside. Step lo is the walk's first
    or a reference step, and so is step hi - 1 its last; the steps outside [lo, hi) are outside every run.

    A stretch of steps off the reference belongs to the region when either of the steps that bound it
    is on the region's interior side: an interior reference node, or an anchor the walk crosses from or
    into the interior. A run is a maximal stretch of steps in the region: interior reference steps and
    the stretches that belong. Its bounds are reference steps outside the interior, or the walk's ends.
    Yields (type, orient, first, end, exit) with [first, end) the run's steps, `exit` whether a bound is
    a reference node other than an anchor, or ('complex', reason, None, None, None):
    - spanning: bounds on both sides, L or before it then R or after it ('+'), or the reverse ('-');
    - prefix / suffix: the walk starts / ends in the region, and the other bound gives the orientation;
    - uturn: both bounds on one side, left or right: a fragment on that side;
    - internal: the walk lies in the region, oriented as most of its interior reference steps."""
    n = len(toks)
    ref = [i for i in range(lo, hi) if kinds[i - lo] is not None]
    inx = bytearray(hi - lo)
    for i in ref:
        if kinds[i - lo] == 'i':
            inx[i - lo] = 1
    # the stretches off the reference, each between two reference steps or a reference step and a walk end
    bounds = [None] + ref + [None]
    for p, q in zip(bounds, bounds[1:]):
        a = lo if p is None else p + 1
        b = hi if q is None else q
        if a >= b:
            continue
        if p is None and lo > 0 or q is None and hi < n:
            continue        # bounded by a step outside [lo, hi), which is not in the region
        left = None if p is None else _bound(kinds[p - lo], toks[p][0], True)
        right = None if q is None else _bound(kinds[q - lo], toks[q][0], False)
        if IN in (left, right):
            if 'x' in (left, right):
                yield ('complex', 'a run reaches a node the reference visits inside and outside the region',
                       None, None, None)
                return
            for i in range(a, b):
                inx[i - lo] = 1
    i = lo
    while i < hi:
        if not inx[i - lo]:
            i += 1
            continue
        a = i
        while i < hi and inx[i - lo]:
            i += 1
        b = i
        before = kinds[a - 1 - lo] if a > lo else None
        after = kinds[b - lo] if b < hi else None
        if 'x' in (before, after):
            yield ('complex', 'a run ends at a node the reference visits inside and outside the region',
                   None, None, None)
            return
        sb = None if before is None else 'left' if before in ('L', 'o') else 'right'
        sa = None if after is None else 'left' if after in ('L', 'o') else 'right'
        exit = before in ('o', 'O') or after in ('o', 'O')
        if sb and sa:
            if sb != sa:
                yield ('spanning', '+' if sb == 'left' else '-', a, b, exit)
            else:
                yield ('uturn', '+' if sb == 'left' else '-', a, b, exit)
        elif sa:
            yield ('prefix', '+' if sa == 'right' else '-', a, b, exit)
        elif sb:
            yield ('suffix', '+' if sb == 'left' else '-', a, b, exit)
        else:
            fwd = sum(1 if toks[k][0] == '>' else -1 for k in range(a, b) if kinds[k - lo] == 'i')
            yield ('internal', '+' if fwd >= 0 else '-', a, b, exit)


class RefIndex:
    """Where the reference visits each node, for exit edges: `start[node]` is the start of its first
    visit, or -1 for a node the reference does not visit; `multi` {node: (first start, last start)} for
    the nodes it visits more than once."""

    def __init__(self, index, max_id):
        self.start = array('q', [-1]) * (max(max_id, max((n for n, _, _, _ in index), default=0)) + 1)
        first, last = {}, {}
        for n, _, st, _ in index:
            if n in first:
                last[n] = st
            else:
                first[n] = st
                self.start[n] = st
        self.multi = {n: (first[n], st) for n, st in last.items()}

    def kind(self, n, r):
        """The class of node `n` for region `r` (exit_runs_of_walk)."""
        st = self.start[n] if n < len(self.start) else -1
        if st < 0:
            return None
        if n == r.L:
            return 'L'
        if n == r.R:
            return 'R'
        ks = {self._side(st, r)}
        if n in self.multi:
            ks |= {self._side(s, r) for s in self.multi[n]}
        return ks.pop() if len(ks) == 1 else 'x'

    @staticmethod
    def _side(st, r):
        return 'o' if st < r.ref_start else 'O' if st >= r.ref_end else 'i'


TRIM_COUNTS = ('cores', 'cores_in_fragments', 'no_core', 'at_end', 'no_orient', 'trimmed_bp')


def target_kmers(regs, ref_indexes, seq, trim_targets, targets_of):
    """{region: TargetKmers} for the regions `regs`, from the targets file `trim_targets` and the
    reference's sequence, its indexes {contig: read_ref_index} and the node sequences `seq`."""
    out = {}
    for rid, r in regs.items():
        index = ref_indexes[r.contig]
        starts = [x[2] for x in index]
        ids = set(r.targets.split(','))
        ts = sorted((t.start, t.end, t.id) for t in targets_of(trim_targets, r.contig) if t.id in ids)
        if {t[2] for t in ts} != ids:
            raise ValueError('region %s: targets %s are not in %s' % (rid, ','.join(sorted(ids - {t[2] for t in ts})),
                                                                     trim_targets))
        targets = [(i, ref_spell(index, starts, seq, max(s, r.ref_start), min(e, r.ref_end))) for s, e, i in ts]
        outside, pos = [], r.ref_start
        for s, e, _ in ts:
            if s > pos:
                outside.append(ref_spell(index, starts, seq, pos, min(s, r.ref_end)))
            pos = max(pos, e)
        if pos < r.ref_end:
            outside.append(ref_spell(index, starts, seq, pos, r.ref_end))
        out[rid] = TargetKmers(targets, outside, CORE_K, CORE_ENRICHMENT)
    return out


def cut_steps(off, cs, ce, margin, start):
    """Where a core [cs, ce) of a run is cut, as an index into the run's step offsets `off` (off[j] the
    start of step j, off[-1] the run's length): at the start (`start`) or end, the step boundary `margin`
    out from the core, or past it; or, where the step there reaches more than `margin` further out, as an
    unaligned stretch of a contig can be one long node, the next boundary in, if that is still outside
    the core."""
    if start:
        x = max(0, cs - margin)
        j = bisect.bisect_right(off, x) - 1
        if x - off[j] > margin and off[j + 1] <= cs:
            j += 1
    else:
        x = min(off[-1], ce + margin)
        j = bisect.bisect_left(off, x)
        if off[j] - x > margin and off[j - 1] >= ce:
            j -= 1
    return j


def ref_spell(index, starts, seq, lo, hi):
    """The reference's sequence of [lo, hi), from its index (prepare's <contig>.ref.tsv.gz, read by
    read_ref_index; `starts` its starts) and the node sequences `seq`."""
    out = []
    for n, o, st, en in index[max(0, bisect.bisect_right(starts, lo) - 1):bisect.bisect_left(starts, hi)]:
        if en > lo:
            s = seq[n] if o == '+' else revcomp(seq[n])
            out.append(s[max(lo - st, 0):min(hi, en) - st])
    return ''.join(out)


def extract(graph, regions, out_dir, vg=None, threads=4, max_internal_bp=200000, ref_sample='CHM13',
            allow_inversions=False, exit_edges=False, ref_dir=None, keep_fragments=False, trim_targets=None,
            trim_margin=20000):
    """Run the pass over `graph` (a GBZ or GFA) for the `regions` (regions.Region) with status 'ok'; the
    others get no package. A package's `ref` is the allele that a walk of `ref_sample` spells. An
    anchor-free walk longer than `max_internal_bp` that uses an interior makes its region complex.
    `allow_inversions` keeps a region whose runs visit a node both ways (see Orientation). `exit_edges`
    finds the runs as exit_runs_of_walk does (see Exit edges), from the reference indexes
    <contig>.ref.tsv.gz that prepare wrote in `ref_dir`. `keep_fragments`, with exit edges, records only
    the spanning runs and leaves the fragments on their nodes (see Keeping fragments). `trim_targets`, the
    targets file regions read, with keep_fragments, cuts each run down to its cores, `trim_margin` bp
    either side of the haplotype's copy of each target (see Trimming to the targets)."""
    from .prepare import read_ref_index
    from .regions import targets_of
    if keep_fragments and not exit_edges:
        raise ValueError('keeping fragments needs exit edges')
    if trim_targets is not None and not keep_fragments:
        raise ValueError('trimming to the targets needs keep_fragments')
    os.makedirs(out_dir, exist_ok=True)
    regs = {r.id: r for r in regions if r.status == 'ok'}
    anchor_of = collections.defaultdict(list)     # node -> [(region id, 'L'|'R')]
    for r in regs.values():
        anchor_of[str(r.L)].append((r.id, 'L'))
        anchor_of[str(r.R)].append((r.id, 'R'))
    ref_indexes, ref_of = {}, {}
    if exit_edges:
        if ref_dir is None:
            raise ValueError('exit edges need the reference indexes (ref_dir)')
        # each region's anchors and the reference nodes between them, by which the pass finds its walks
        anchor_of = collections.defaultdict(list)
        for contig in sorted({r.contig for r in regs.values()}):
            index = ref_indexes[contig] = read_ref_index(os.path.join(ref_dir, contig + '.ref.tsv.gz'))
            starts = [x[2] for x in index]
            for r in regs.values():
                if r.contig != contig:
                    continue
                anchor_of[str(r.L)].append((r.id, 'L'))
                anchor_of[str(r.R)].append((r.id, 'R'))
                for n, _, _, _ in index[bisect.bisect_left(starts, r.ref_start):bisect.bisect_left(starts, r.ref_end)]:
                    if n != r.L and n != r.R:
                        anchor_of[str(n)].append((r.id, 'i'))
    runs = {rid: RegionRuns(allow_inversions) for rid in regs}
    seq = {}
    tks = None            # trim_targets: region -> TargetKmers, once the node sequences are read
    cores_out = None
    if trim_targets is not None:
        cores_out = gzip.open(os.path.join(out_dir, 'cores.tsv.gz'), 'wt')
        cores_out.write('#walk\tregion\ttargets\ttype\tstatus\torient\tfirst\tend\tcontig\tstretch_start\t'
                        'stretch_end\tcore_start\tcore_end\n')
    used = None           # exit edges: the nodes that walks visit outside every run
    occ_out = gzip.open(os.path.join(out_dir, 'occurrences.tsv.gz'), 'wt')
    occ_out.write('#walk\tregion\ttype\torient\tfirst\tend\tmember%s\n' % ('\tbefore\tafter' if exit_edges else ''))
    candidates = []       # anchor-free walks of at most max_internal_bp: possible internal fragments
    long_walks = []       # the keys of the longer anchor-free walks
    long_at = {}          # node -> index in long_walks of the first longer walk to visit it
    n_walks = 0

    def mark_used(toks, a, b):
        for _, n in toks[a:b]:
            used[int(n)] = 1

    def on_ref(n):
        return any(ri.start[n] >= 0 for ri in ref_of.values() if n < len(ri.start))

    def record(rr, rid, key, w, toks, typ, orient, a, b):
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
        row = '%s\t%s\t%s\t%s\t%d\t%d\t%s' % (key, rid, typ, orient, a, b, mid)
        if exit_edges:
            # the steps that bound the run, which replace links its new walk to
            row += '\t%s\t%s' % (''.join(toks[a - 1]) if a > 0 else '.', ''.join(toks[b]) if b < len(toks) else '.')
        occ_out.write(row + '\n')

    def trim_run(rid, rr, key, w, toks, typ, orient, a, b):
        """The cores of the run [a, b) of walk `key`, which has type `typ` and orient `orient`, as
        [(first, end, orient)]; the rest of the run stays on its nodes (see Trimming to the targets)."""
        steps = toks[a:b]
        off = list(itertools.accumulate((len(seq[int(n)]) for _, n in steps), initial=0))
        s = ''.join(seq[int(n)] if o == '>' else revcomp(seq[int(n)]) for o, n in steps).upper()
        tk = tks[rid]
        found = []
        for t, cs, ce, ss, se, strand in find_cores(s, tk, CORE_STRIDE, CORE_WINDOW, CORE_DENSITY, CORE_GAP):
            found.append([a + cut_steps(off, cs, ce, trim_margin, True), a + cut_steps(off, cs, ce, trim_margin, False),
                          {tk.ids[t]}, {strand}, ss, se])
        found.sort()
        merged = []
        for c in found:
            if merged and c[0] < merged[-1][1]:
                m = merged[-1]
                m[1], m[4], m[5] = max(m[1], c[1]), min(m[4], c[4]), max(m[5], c[5])
                m[2] |= c[2]
                m[3] |= c[3]
            else:
                merged.append(c)
        base = (int(w.start) if w.start.isdigit() else 0) + sum(len(seq[int(n)]) for _, n in toks[:a])
        cores = []
        for first, end, ids, strands, ss, se in merged:
            strand = strands.pop() if len(strands) == 1 else None
            if first == 0 or end == len(toks):
                status = 'at_end'
            elif strand is None and typ != 'spanning':
                status = 'no_orient'
            else:
                status = 'core'
                cores.append((first, end, strand or orient))
            rr.trim['cores' if status == 'core' else status] += 1
            cores_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n' % (
                key, rid, ','.join(sorted(ids)), typ, status, cores[-1][2] if status == 'core' else '.', first, end,
                w.contig, base + ss, base + se, base + off[first - a], base + off[end - a]))
        if typ == 'spanning' and not merged:
            rr.trim['no_core'] += 1
            cores_out.write('%s\t%s\t.\t%s\tno_core\t.\t%d\t%d\t%s\t.\t.\t%d\t%d\n' % (
                key, rid, typ, a, b, w.contig, base, base + off[-1]))
        if cores:
            rr.trim['cores_in_fragments'] += len(cores) if typ != 'spanning' else 0
            rr.trim['trimmed_bp'] += off[-1] - sum(off[e - a] - off[f - a] for f, e, _ in cores)
        else:
            rr.kept_runs += 1
            rr.kept_run_bp += off[-1]
        return cores

    with GfaStream(graph, vg, threads) as f:
        for line in f:
            t = line[0]
            if t == 'S':
                x = line.split('\t', 3)
                seq[int(x[1])] = x[2].rstrip('\n')
                continue
            if t != 'W':
                continue
            if exit_edges and used is None:
                max_id = max(seq) if seq else 0
                ref_of = {c: RefIndex(index, max_id) for c, index in ref_indexes.items()}
                used = bytearray(max(len(ri.start) for ri in ref_of.values()) if ref_of else max_id + 1)
                if trim_targets is not None:
                    tks = target_kmers(regs, ref_indexes, seq, trim_targets, targets_of)
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
                if exit_edges and toks and (keep_fragments or any(on_ref(int(n)) for _, n in toks)):
                    mark_used(toks, 0, len(toks))
                elif regs and toks:
                    if sum(len(seq[int(n)]) for _, n in toks) <= max_internal_bp:
                        candidates.append((key, toks))
                    elif exit_edges:
                        # outside every run, so its nodes stay for it, whichever regions they are in
                        mark_used(toks, 0, len(toks))
                    else:
                        for _, n in toks:
                            long_at.setdefault(int(n), len(long_walks))
                        long_walks.append(key)
                continue
            if exit_edges:
                spans = []
                for rid, occ in hits.items():
                    r, rr = regs[rid], runs[rid]
                    ri = ref_of[r.contig]
                    st = ri.start
                    lo = occ[0][0] - 1
                    while lo >= 0 and st[int(toks[lo][1])] < 0:
                        lo -= 1
                    lo = max(lo, 0)
                    hi = occ[-1][0] + 1
                    while hi < len(toks) and st[int(toks[hi][1])] < 0:
                        hi += 1
                    hi = min(hi + 1, len(toks))
                    kinds = [ri.kind(int(toks[i][1]), r) for i in range(lo, hi)]
                    for typ, orient, a, b, ex in exit_runs_of_walk(toks, lo, hi, kinds):
                        if typ == 'complex':
                            rr.complex.append('%s: %s' % (key, orient))
                            continue
                        if tks is not None:
                            # only the cores are runs: the rest stays on its nodes
                            cores = trim_run(rid, rr, key, w, toks, typ, orient, a, b)
                            for ca, cb, co in cores:
                                record(rr, rid, key, w, toks, 'spanning', co, ca, cb)
                                spans.append((ca, cb, rid))
                            rr.exit_runs += bool(ex and cores)
                            continue
                        if keep_fragments and typ != 'spanning':
                            # outside every run, so its nodes stay for it
                            rr.kept_runs += 1
                            rr.kept_run_bp += sum(len(seq[int(n)]) for _, n in toks[a:b])
                            continue
                        record(rr, rid, key, w, toks, typ, orient, a, b)
                        rr.exit_runs += bool(ex)
                        spans.append((a, b, rid))
                spans.sort()
                pos = 0
                for a, b, rid in spans:
                    if a < pos:
                        for x in {rid, prev}:
                            runs[x].complex.append('%s: runs of regions %s and %s overlap' % (key, prev, rid))
                    mark_used(toks, pos, a)
                    pos, prev = max(pos, b), rid
                mark_used(toks, pos, len(toks))
                continue
            for rid, occ in hits.items():
                rr = runs[rid]
                for typ, orient, a, b in runs_of_walk(toks, occ):
                    if typ == 'complex':
                        rr.complex.append('%s: %s' % (key, orient))
                        continue
                    record(rr, rid, key, w, toks, typ, orient, a, b)
    # internal fragments: the held anchor-free walks that use a region's interior nodes
    n_flipped = 0
    node_region = {}
    for rid, rr in runs.items():
        for n in rr.interior:
            other = node_region.setdefault(n, rid)
            if other != rid and exit_edges:
                # exit edges can take one stretch into two regions
                for x in (rid, other):
                    runs[x].complex.append('node %d is in the interiors of regions %s and %s' % (n, other, rid))
    for key, toks in candidates:
        rids = {node_region.get(int(n)) for _, n in toks} - {None}
        if not rids:
            if exit_edges:
                mark_used(toks, 0, len(toks))
            continue
        if len(rids) > 1:
            for rid in rids:
                runs[rid].complex.append('%s: an anchor-free walk through %d regions' % (key, len(rids)))
            if exit_edges:
                mark_used(toks, 0, len(toks))
            continue
        rid = rids.pop()
        rr = runs[rid]
        steps = [(o, int(n)) for o, n in toks]
        # the ways the runs visit the interior nodes this walk visits, against the ways it does
        known = [(FWD if o == '>' else REV, rr.interior[n]) for o, n in steps if n in rr.interior]
        agree = sum(bool(had & way) for way, had in known)
        against = sum(bool(had & ((FWD | REV) ^ way)) for way, had in known)
        if agree < len(known) and against < len(known) and not allow_inversions:
            rr.complex.append('%s: an anchor-free walk runs both ways through the region' % key)
            continue
        if agree == len(known) or agree >= against:
            orient = '+'
        else:
            orient, steps = '-', flip(steps)
            n_flipped += 1
        rr.visit(key, steps)
        mid = _member(rr.fragments, ('internal', spell(steps, seq)), 'f')
        rr.n_runs['internal'] += 1
        occ_out.write('%s\t%s\tinternal\t%s\t0\t%d\t%s%s\n' % (key, rid, orient, len(toks), mid,
                                                             '\t.\t.' if exit_edges else ''))
    occ_out.close()
    if cores_out is not None:
        cores_out.close()
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
        it.write('#region\tnodes%s\n' % ('\tkept' if exit_edges else ''))
        for rid in sorted(regs, key=lambda k: (regs[k].ref_start, regs[k].ref_end)):
            r, rr = regs[rid], runs[rid]
            # exit edges: the nodes of the runs that walks also visit outside every run stay for them
            kept = sorted(n for n in rr.interior if used is not None and used[n])
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
            if exit_edges:
                pkg.update({'exit_edges': True, 'exit_runs': rr.exit_runs, 'kept_nodes': len(kept),
                            'kept_bp': sum(len(seq[n]) for n in kept)})
            if keep_fragments:
                pkg.update({'kept_fragments': rr.kept_runs, 'kept_fragment_bp': rr.kept_run_bp})
            if trim_targets is not None:
                pkg['trim'] = dict({k: rr.trim[k] for k in TRIM_COUNTS}, margin=trim_margin)
            pk.write(json.dumps(pkg, separators=(',', ':')) + '\n')
            it.write('%s\t%s%s\n' % (rid, ','.join(str(n) for n in sorted(rr.interior)),
                                      '\t' + ','.join(str(n) for n in kept) if exit_edges else ''))
    summary = {'walks': n_walks, 'regions': len(regs), 'status': dict(status),
               'internal_candidates': len(candidates), 'internal_too_long': n_long,
               'both_ways': sum(rr.both_ways for rr in runs.values()), 'internal_flipped': n_flipped}
    if exit_edges:
        summary['exit_runs'] = sum(rr.exit_runs for rr in runs.values())
    if keep_fragments:
        summary['kept_fragments'] = sum(rr.kept_runs for rr in runs.values())
    if trim_targets is not None:
        summary['trim'] = {k: sum(rr.trim[k] for rr in runs.values()) for k in TRIM_COUNTS}
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
