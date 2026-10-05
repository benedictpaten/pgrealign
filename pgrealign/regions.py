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

**Merging.** Regions are taken in order of their left anchor's start, and one merges into the region
before it when its left anchor starts before that region's right anchor starts. That covers both
overlapping interiors and an anchor inside another region's interior, where replacing one region would
change the other. The merged targets are re-anchored as one, and the pass repeats until nothing
merges. Two regions that only share a boundary node, one's right anchor being the next one's left,
stay separate, because each replaces only its own interior.

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


Region = collections.namedtuple('Region', 'id contig L R ref_start ref_end span targets classes status note')
HEADER = '\t'.join(Region._fields)


def assign(targets, anchors, contig, pad=200, max_span=100000, contig_length=None):
    """Regions for the targets of one contig: anchor each, merge, and apply max_span (see the module
    docstring)."""
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


def build(targets, contig, snarls, ref, out, vg=None, pad=200, max_span=100000):
    """The regions stage: anchor `targets` (a catalogue or BED) on `contig`'s snarls and write `out`.
    `snarls` and `ref` are prepare's <contig>.snarls.pb and <contig>.ref.tsv.gz. They are read only if
    the contig has targets; with none, `out` has no rows and the summary's `chains` is null. Returns a
    summary."""
    from .prepare import read_ref_index
    ts = targets_of(targets, contig)
    chains, rs = None, []
    if ts:
        index = read_ref_index(ref)
        chains = chains_from_pairs(snarl_pairs(snarls, vg), index)
        rs = assign(ts, Anchors(chains), contig, pad=pad, max_span=max_span,
                    contig_length=index[-1][3] if index else None)
    write_regions(rs, out)
    return {'targets': len(ts), 'chains': None if chains is None else len(chains), 'regions': len(rs),
            'status': dict(collections.Counter(r.status for r in rs)),
            'merged_targets': sum(len(r.targets.split(',')) for r in rs if r.status == 'ok'),
            'span_bp_ok': sum(r.span for r in rs if r.status == 'ok')}
