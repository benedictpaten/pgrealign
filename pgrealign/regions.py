"""Stage 1, `regions`: anchor every target interval in the graph's own snarl tree.

**Anchors.** A target's anchors are two boundary nodes L and R of one chain of the snarl tree. The
snarls come from `vg snarls -T -P <ref>`, trivial snarls included, so consecutive snarls of a chain
share a boundary node.
- L is the chain's last boundary node ending at or before the padded target start.
- R is its first boundary node starting at or after the padded target end.
- Over every chain, top-level or nested, the pair with the smallest reference span wins.
- Only boundary nodes the reference visits once, forward, are used.

Two boundary nodes of one chain separate the graph. So every path that crosses the region crosses both
anchors, and a path that does not cross both starts or ends inside the region (a fragment). The
anchors stay in the graph. Only the nodes strictly between them, the region's interior, are replaced.

**Merging.** Targets whose interiors overlap are merged, and the union is re-anchored, repeatedly. Two
regions that only share a boundary node are separate, because each replaces its own interior.

**Size policy.** A region whose anchor span exceeds `max_span` is kept as it is. It is listed with
status `too_long`.

**Output.** `regions.tsv`, one row per region:

    id contig L R ref_start ref_end span targets classes status note

- L and R are node ids, both forward on the reference.
- ref_start and ref_end are the 0-based interior on the reference: from L's end to R's start.
- `targets` lists the merged target ids, and `classes` their classes.
"""
import bisect
import collections
import gzip
import json
import subprocess

from . import catalog
from .graph import find_vg


class AnchorError(RuntimeError):
    pass


def snarl_pairs(snarls_pb, vg=None):
    """(start node, end node) of every snarl in a vg snarls file, via `vg view -Rj`."""
    p = subprocess.Popen([find_vg(vg), 'view', '-Rj', snarls_pb], stdout=subprocess.PIPE, text=True,
                         bufsize=1 << 20)
    for line in p.stdout:
        s = json.loads(line)
        yield int(s['start']['node_id']), int(s['end']['node_id'])
    if p.wait():
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
            span = R[1] - L[0]
            if span <= max_span and (best is None or span < best[0]):
                best = (span, L, R)
        if best is None:
            raise AnchorError('no chain encloses %d-%d within %d bp' % (a0, b0, max_span))
        return best[1], best[2]


Region = collections.namedtuple('Region', 'id contig L R ref_start ref_end span targets classes status note')
HEADER = '\t'.join(Region._fields)


def assign(targets, anchors, contig, pad=200, max_span=100000, contig_length=None):
    """Regions for the targets of one contig: anchor, merge overlapping interiors, apply max_span."""
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
            # interiors overlap: this region's left anchor starts before the previous one's right anchor
            if merged and g['L'][0] < merged[-1]['R'][0]:
                m = merged[-1]
                m['targets'] += g['targets']
                m['a0'], m['b0'] = min(m['a0'], g['a0']), max(m['b0'], g['b0'])
                changed = True
                try:
                    m['L'], m['R'] = anchors.query(m['a0'], m['b0'], max_span)
                except AnchorError as e:
                    m['error'] = str(e)
                    m['L'] = min(m['L'], g['L'])
                    m['R'] = max(m['R'], g['R'], key=lambda x: x[1])
            else:
                merged.append(g)
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
    with op(path, 'wt') as f:
        f.write('#' + HEADER + '\n')
        for r in regions:
            f.write('\t'.join(str(v) for v in r) + '\n')


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
