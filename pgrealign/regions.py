"""Stage 1, `regions`: anchor every target interval in the graph's own snarl tree.

**Anchors.** A target's anchors are two boundary nodes L and R of one chain of the snarl tree. The
snarls come from `vg snarls -T -P <ref>`, trivial snarls included, so consecutive snarls of a chain
share a boundary node.
- Only boundary nodes the reference visits once, forward, are used.
- The target is padded by `pad` on each side, within the contig, and then clamped to lie between the
  end of the reference's first boundary node and the start of its last, so a target reaching into
  either is anchored on it. A target that lies wholly before the first node's end, or wholly after
  the last node's start, has no boundary nodes on both sides.
- L is the chain's last boundary node ending at or before the padded target start.
- R is its first boundary node starting at or after the padded target end, and is never L itself.
- Over every chain, top-level or nested, the pair with the smallest reference span wins.

Two boundary nodes of one chain separate the graph. So every path that crosses the region crosses both
anchors, and a path that does not cross both starts or ends inside the region (a fragment). The
anchors stay in the graph. Only the nodes strictly between them, the region's interior, are replaced.

**Anchors that walks go around** (`max_bypass`). Snarl anchors can lie far from a target: a few
haplotypes whose own sequence of the target never aligned to the reference leave it on one side and
come back past every nearby boundary node, so the nearest pair that separates the graph lies Mb away.
With `max_bypass`, a fraction, L and R are instead the nearest nodes on each side, among those the
reference visits once, forward, that at most that fraction of the haplotypes' walks go around: two
consecutive reference steps of the walk lie on either side of the node. The anchors then no longer
separate the graph, and extract must be run with exit edges, which take the walks that go around them.
A region's note gives the number of walks around each anchor.

**Merging.** Regions are taken in order of their left anchor's start, and one merges into the region
before it when its left anchor starts before that region's right anchor starts. That covers both
overlapping interiors and an anchor inside another region's interior, where replacing one region would
change the other. With `max_bypass`, a region also merges into an earlier one, and every region
between, when a walk goes straight from that region's interior into its own: two consecutive reference
steps of the walk lie in the two interiors, so its runs in them would overlap. The merged targets are
re-anchored as one, and the pass repeats until nothing merges. Otherwise two regions that only share a
boundary node, one's right anchor being the next one's left, stay separate, because each replaces only
its own interior.

**Size policy.** A region whose anchors cannot be within `max_span` of each other has status
`too_long` and is not realigned:
- A target with no chain enclosing it, or none within `max_span`, keeps its own interval. Its row has
  '.' anchors, and the target's start, end and length as ref_start, ref_end and span. Another region
  can still overlap it.
- Merged targets that cannot be re-anchored as one keep the outermost anchors of the regions merged.

**Output.** `regions.tsv`, one row per region, in reference order:

    id contig L R ref_start ref_end span targets classes status note

- `id` is the id of its first target.
- L and R are node ids, both forward on the reference.
- ref_start and ref_end are the 0-based interior on the reference: from L's end to R's start.
- `span` is the reference length from L's start to R's end.
- `targets` lists the merged target ids, and `classes` their classes.
- `note` says how many targets were merged, or why a region is `too_long`.
"""
import bisect
import collections
import gzip
import json
import os
import subprocess
from array import array

from . import catalog
from .graph import find_vg


class AnchorError(RuntimeError):
    pass


def snarl_pairs(snarls_pb, vg=None):
    """(start node, end node) of every snarl in a vg snarls file, via `vg view -Rj`."""
    with subprocess.Popen([find_vg(vg), 'view', '-Rj', snarls_pb], stdout=subprocess.PIPE, text=True,
                          bufsize=1 << 20) as p:
        for line in p.stdout:
            s = json.loads(line)
            yield int(s['start']['node_id']), int(s['end']['node_id'])
    if p.returncode:
        raise RuntimeError('vg view -Rj %s failed' % snarls_pb)


def chains_from_pairs(pairs, ref_walk):
    """Chains as sorted lists of (start, end, node) of their reference-forward, once-visited boundary
    nodes. A chain is a connected component of "shares a boundary node" over all snarls."""
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in pairs:
        for x in (a, b):
            parent.setdefault(x, x)
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    visits = collections.Counter(n for n, _, _, _ in ref_walk)
    chains = collections.defaultdict(list)
    for n, o, st, en in ref_walk:
        if n in parent and o == '+' and visits[n] == 1:
            chains[find(n)].append((st, en, n))
    return [sorted(v) for v in chains.values() if len(v) >= 2]


class Anchors:
    """The tightest pair of boundary nodes of one chain around a reference interval."""

    def __init__(self, chains):
        pts = []
        for ci, ch in enumerate(chains):
            for st, en, n in ch:
                pts.append((st, en, n, ci))
        pts.sort()
        if not pts:
            raise AnchorError('no chain boundary nodes on the reference')
        self.pts = pts
        self.starts = [p[0] for p in pts]
        self.ends_sorted = sorted((p[1], k) for k, p in enumerate(pts))
        self.ends = [e for e, _ in self.ends_sorted]
        self.first_end = min(p[1] for p in pts)
        self.last_start = max(p[0] for p in pts)

    def query(self, a0, b0, max_span):
        """(L, R) as (start, end, node) around [a0, b0), clamped to the outermost boundary nodes."""
        a0, b0 = max(a0, self.first_end), min(b0, self.last_start)
        if a0 >= b0:
            raise AnchorError('no boundary nodes on both sides of %d-%d' % (a0, b0))
        left = {}
        i = bisect.bisect_right(self.ends, a0)
        lo_lim = b0 - max_span
        while i > 0:
            i -= 1
            e, k = self.ends_sorted[i]
            st, en, n, ci = self.pts[k]
            if st < lo_lim:
                if e < lo_lim:
                    break
                continue
            if ci not in left or en > left[ci][1]:
                left[ci] = (st, en, n)
        right = {}
        j = bisect.bisect_left(self.starts, b0)
        hi_lim = a0 + max_span
        while j < len(self.pts):
            st, en, n, ci = self.pts[j]
            if st > hi_lim:
                break
            if en <= hi_lim and (ci not in right or st < right[ci][0]):
                right[ci] = (st, en, n)
            j += 1
        best = None
        for ci in set(left) & set(right):
            L, R = left[ci], right[ci]
            if R[0] < L[1]:
                continue       # one node, or anchors that overlap: nothing lies between them
            span = R[1] - L[0]
            if span <= max_span and (best is None or span < best[0]):
                best = (span, L, R)
        if best is None:
            raise AnchorError('no chain encloses %d-%d within %d bp' % (a0, b0, max_span))
        return best[1], best[2]


def ref_ranks(index):
    """({node: the rank of its first step on the reference}, [(start, end, node)] of the nodes the
    reference visits once, forward, in reference order) from a reference index (prepare.read_ref_index)."""
    rank, visits = {}, collections.Counter(n for n, _, _, _ in index)
    once = []
    for k, (n, o, st, en) in enumerate(index):
        rank.setdefault(n, k)
        if o == '+' and visits[n] == 1:
            once.append((st, en, n))
    return rank, once


def bypass_counts(graph, index, windows, ref_sample='CHM13', vg=None, threads=4, gaps=()):
    """How many walks go around each candidate anchor: the nodes the reference visits once, forward,
    that lie within one of `windows` [(start, end)] of the reference. A walk goes around a node when two
    of its consecutive reference steps lie on either side of it, so that it leaves the reference before
    the node and comes back after it without visiting it. The walks counted are every haplotype's: not
    the reference sample's, nor those of a sample whose name starts with '_' (minigraph's segment paths
    in a Minigraph-Cactus graph). Returns ([(start, end, node)] of the candidates in reference order,
    [the number of walks around each], the number of haplotypes, links): `links` are the distinct pairs
    of consecutive reference steps that go around a candidate within one of `gaps` [(start, end)], as
    the sorted [(start, start)] of the two steps on the reference."""
    from .graph import GfaStream
    rank, once = ref_ranks(index)
    windows = sorted(windows)
    cands = [c for c in once if any(a <= c[0] and c[1] <= b for a, b in windows)]
    crank = [rank[n] for _, _, n in cands]           # increasing, as the candidates are in reference order
    counts = [0] * len(cands)
    if not cands:
        return cands, counts, 0, []
    lo_rank, hi_rank = crank[0], crank[-1]
    diff = [0] * (len(cands) + 1)
    # in_gap[i]: the candidates before i that lie within a gap, so a jump over candidates [i, j) goes around
    # one when in_gap[j] > in_gap[i]
    in_gap = [0]
    for st, en, _ in cands:
        in_gap.append(in_gap[-1] + any(a <= st and en <= b for a, b in gaps))
    links = set()
    haps = set()
    rk = None
    max_id = max(rank)
    with GfaStream(graph, vg, threads) as f:
        for line in f:
            t = line[0]
            if t == 'S':
                max_id = max(max_id, int(line[2:line.index('\t', 2)]))
                continue
            if t != 'W':
                continue
            if rk is None:
                rk = array('l', [-1]) * (max_id + 1)
                for n, k in rank.items():
                    rk[n] = k
            x = line.split('\t', 7)
            if x[1] == ref_sample or x[1].startswith('_'):
                continue
            haps.add((x[1], x[2]))
            walk = x[6].rstrip('\n')
            r = [k for k in (rk[int(n)] for n in walk[1:].replace('<', '>').split('>')) if k >= 0]
            spans = []
            for a, b in zip(r, r[1:]):
                if a > b:
                    a, b = b, a
                # consecutive reference steps that are not neighbours on the reference, around a candidate
                if b - a > 1 and b > lo_rank and a < hi_rank:
                    i, j = bisect.bisect_right(crank, a), bisect.bisect_left(crank, b)
                    if i < j:
                        spans.append((i, j))
                        if in_gap[j] > in_gap[i]:
                            links.add((index[a][2], index[b][2]))
            if not spans:
                continue
            spans.sort()
            end = -1
            for i, j in spans:
                # each candidate once per walk, however many of its steps go around it
                i = max(i, end)
                if i < j:
                    diff[i] += 1
                    diff[j] -= 1
                    end = j
    k = 0
    for i in range(len(cands)):
        k += diff[i]
        counts[i] = k
    return cands, counts, len(haps), sorted(links)


class BypassAnchors:
    """The anchors of a target are the nearest candidate on each side that at most `max_walks` walks go
    around (bypass_counts). Exit edges in extract (`exit_edges`) take the walks that go around them."""

    def __init__(self, cands, counts, max_walks):
        ok = [c for c, k in zip(cands, counts) if k <= max_walks]
        self.pts = ok
        self.starts = [c[0] for c in ok]
        self.ends = [c[1] for c in ok]
        self.count = {c[2]: k for c, k in zip(cands, counts)}
        self.max_walks = max_walks

    def query(self, a0, b0, max_span):
        """(L, R) as (start, end, node) around [a0, b0)."""
        i = bisect.bisect_right(self.ends, a0) - 1
        j = bisect.bisect_left(self.starts, b0)
        if i < 0 or j >= len(self.pts):
            raise AnchorError('no node on %s side of %d-%d that at most %d walks go around' % (
                'the left' if i < 0 else 'the right', a0, b0, self.max_walks))
        L, R = self.pts[i], self.pts[j]
        if R[1] - L[0] > max_span:
            raise AnchorError('the nearest nodes around %d-%d that at most %d walks go around are %d bp apart, '
                              'more than %d' % (a0, b0, self.max_walks, R[1] - L[0], max_span))
        return L, R


Region = collections.namedtuple('Region', 'id contig L R ref_start ref_end span targets classes status note')
HEADER = '\t'.join(Region._fields)


def _linked(m, g, links):
    """Whether one of `links` [(start, start)], sorted, goes from region m's interior into g's, later on
    the reference."""
    lo, hi = bisect.bisect_left(links, (m['L'][1],)), bisect.bisect_left(links, (m['R'][0],))
    return any(g['L'][1] <= b < g['R'][0] for _, b in links[lo:hi])


def assign(targets, anchors, contig, pad=200, max_span=100000, contig_length=None, links=()):
    """Regions for the targets of one contig: anchor each, merge, and apply max_span (see the module
    docstring). `links` are bypass_counts' links: a region merges with an earlier one when one goes from
    that region's interior into its own."""
    links = sorted(links)
    groups, out = [], []
    for t in targets:
        a0 = max(0, t.start - pad)
        b0 = t.end + pad if contig_length is None else min(contig_length, t.end + pad)
        try:
            L, R = anchors.query(a0, b0, max_span)
        except AnchorError as e:
            out.append(Region(t.id, contig, '.', '.', t.start, t.end, t.end - t.start, t.id, t.cls,
                              'too_long', str(e)))
            continue
        groups.append({'targets': [t], 'a0': a0, 'b0': b0, 'L': L, 'R': R})
    while True:
        groups.sort(key=lambda g: (g['L'][0], g['R'][1]))
        merged, changed = [], False
        for g in groups:
            # this region's left anchor starts before the previous one's right anchor: they merge
            if merged and g['L'][0] < merged[-1]['R'][0]:
                i = len(merged) - 1
            else:
                # a walk goes straight from an earlier region's interior into this one's: it merges with
                # that region and every region between
                i = next((k for k, m in enumerate(merged) if _linked(m, g, links)), None)
            if i is None:
                merged.append(g)
                continue
            group = merged[i:] + [g]
            del merged[i + 1:]
            m = merged[i]
            for h in group[1:]:
                m['targets'] += h['targets']
                m['a0'], m['b0'] = min(m['a0'], h['a0']), max(m['b0'], h['b0'])
            changed = True
            try:
                m['L'], m['R'] = anchors.query(m['a0'], m['b0'], max_span)
            except AnchorError as e:
                m['error'] = str(e)
                m['L'] = min(h['L'] for h in group)
                m['R'] = max((h['R'] for h in group), key=lambda x: x[1])
        groups = merged
        if not changed:
            break
    for g in groups:
        ts = sorted(g['targets'], key=lambda t: t.start)
        rid = ts[0].id
        ids = ','.join(t.id for t in ts)
        cls = ','.join(sorted({t.cls for t in ts}))
        if g.get('error'):
            out.append(Region(rid, contig, g['L'][2], g['R'][2], g['L'][1], g['R'][0],
                              g['R'][1] - g['L'][0], ids, cls, 'too_long', g['error']))
            continue
        L, R = g['L'], g['R']
        out.append(Region(rid, contig, L[2], R[2], L[1], R[0], R[1] - L[0], ids, cls, 'ok',
                          'merged %d targets' % len(ts) if len(ts) > 1 else ''))
    out.sort(key=lambda r: (r.ref_start, r.ref_end))
    return out


def write_regions(regions, path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path + '.tmp', 'wt') as f:
        f.write('#' + HEADER + '\n')
        for r in regions:
            f.write('\t'.join(str(v) for v in r) + '\n')
    os.replace(path + '.tmp', path)


def read_regions(path):
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            x = line.rstrip('\n').split('\t')
            L = int(x[2]) if x[2] != '.' else None
            R = int(x[3]) if x[3] != '.' else None
            yield Region(x[0], x[1], L, R, int(x[4]), int(x[5]), int(x[6]), x[7], x[8], x[9], x[10])


def targets_of(catalog_path, contig):
    return [t for t in catalog.read_catalog(catalog_path) if t.chrom == contig]


def build(targets, contig, snarls, ref, out, vg=None, pad=200, max_span=100000, graph=None, max_bypass=None,
          ref_sample='CHM13', threads=4):
    """The regions stage: anchor `targets` (a catalogue or BED) on `contig`'s snarls and write `out`.
    `snarls` and `ref` are prepare's <contig>.snarls.pb and <contig>.ref.tsv.gz. They are read only if
    the contig has targets; with none, `out` has no rows and the summary's `chains` is null.

    With `max_bypass`, a fraction, the anchors are instead the nearest reference nodes on each side that
    at most that fraction of the haplotypes' walks in `graph` go around (bypass_counts, BypassAnchors),
    and `snarls` is not read. Such anchors do not separate the graph, so extract must take the walks that
    go around them through exit edges (extract's `exit_edges`). Returns a summary."""
    from .prepare import read_ref_index
    ts = targets_of(targets, contig)
    chains, rs, extra, links = None, [], {}, ()
    if ts:
        index = read_ref_index(ref)
        contig_length = index[-1][3] if index else None
        if max_bypass is None:
            chains = chains_from_pairs(snarl_pairs(snarls, vg), index)
            anchors = Anchors(chains)
        else:
            if graph is None:
                raise ValueError('max_bypass needs the graph, to count the walks around each node')
            padded = sorted((max(0, t.start - pad), t.end + pad) for t in ts)
            windows = [(a - max_span, b + max_span) for a, b in padded]
            # the gaps between the padded targets, where a link between two regions' interiors would
            # go around both regions' anchors
            gaps, hi = [], None
            for a, b in padded:
                if hi is not None and a > hi:
                    gaps.append((hi, a))
                hi = b if hi is None else max(hi, b)
            cands, counts, n_haps, links = bypass_counts(graph, index, windows, ref_sample, vg, threads, gaps)
            max_walks = int(max_bypass * n_haps)
            anchors = BypassAnchors(cands, counts, max_walks)
            extra = {'haplotypes': n_haps, 'max_bypass_walks': max_walks, 'candidates': len(cands),
                     'links': len(links)}
        rs = assign(ts, anchors, contig, pad=pad, max_span=max_span, contig_length=contig_length, links=links)
        if max_bypass is not None:
            rs = [r._replace(note=('%s; ' % r.note if r.note else '') + 'walks around L %d, R %d' % (
                anchors.count[r.L], anchors.count[r.R])) if r.status == 'ok' else r for r in rs]
    write_regions(rs, out)
    return {'targets': len(ts), 'chains': None if chains is None else len(chains), 'regions': len(rs),
            'status': dict(collections.Counter(r.status for r in rs)),
            'merged_targets': sum(len(r.targets.split(',')) for r in rs if r.status == 'ok'),
            'span_bp_ok': sum(r.span for r in rs if r.status == 'ok'), **extra}
