"""Stage 4, `induce`: the graph of one multiple alignment, with one walk per row.

**Column induction.** Each column gets one node per distinct base, and each row walks through its
bases' nodes. A gap column adds nothing to a row.

**Unchopping.** A node u is merged into v when v is its only successor, u is v's only predecessor, and
no row ends at u or starts at v.

**Chopping.** Nodes longer than `max_node` (1,024, vg's default) are cut, so the graph needs no id
translation when it becomes a GBZ.

The result is a forward-only DAG with local node ids 1..n in topological order, and every row's walk
spells the row's sequence. An empty row (an allele that deletes the whole interior) has an empty walk.
"""
import collections
import heapq

MAX_NODE = 1024


class Graph:
    def __init__(self):
        self.seq = {}            # node -> sequence
        self.edges = set()       # (from, to), all forward
        self.walks = {}          # member -> [node, ...]


def columns(rows, names):
    g = Graph()
    nid = 1
    prev = {n: None for n in names}
    walks = {n: [] for n in names}
    L = len(rows[names[0]]) if names else 0
    for c in range(L):
        col = {}
        for n in names:
            b = rows[n][c]
            if b == '-':
                continue
            v = col.get(b)
            if v is None:
                v = col[b] = nid
                g.seq[nid] = b
                nid += 1
            if prev[n] is not None:
                g.edges.add((prev[n], v))
            prev[n] = v
            walks[n].append(v)
    g.walks = walks
    return g


def unchop(g):
    out_n, in_n = collections.defaultdict(set), collections.defaultdict(set)
    for a, b in g.edges:
        out_n[a].add(b)
        in_n[b].add(a)
    ends = {w[-1] for w in g.walks.values() if w}
    starts = {w[0] for w in g.walks.values() if w}
    nxt, has_prev = {}, set()
    for u in g.seq:
        if len(out_n[u]) != 1 or u in ends:
            continue
        v = next(iter(out_n[u]))
        if len(in_n[v]) != 1 or v in starts or v == u:
            continue
        nxt[u] = v
        has_prev.add(v)
    ng, new_id, k = Graph(), {}, 0
    for u in sorted(g.seq):
        if u in has_prev:
            continue
        k += 1
        parts, x = [g.seq[u]], u
        new_id[u] = k
        while x in nxt:
            x = nxt[x]
            parts.append(g.seq[x])
            new_id[x] = k
        ng.seq[k] = ''.join(parts)
    for a, b in g.edges:
        if nxt.get(a) != b:
            ng.edges.add((new_id[a], new_id[b]))
    for n, w in g.walks.items():
        q = []
        for x in w:
            m = new_id[x]
            if not q or x not in has_prev:
                q.append(m)
        ng.walks[n] = q
    return topological(ng)


def topological(g):
    indeg, succ = collections.Counter(), collections.defaultdict(list)
    for a, b in g.edges:
        indeg[b] += 1
        succ[a].append(b)
    heap = [n for n in g.seq if indeg[n] == 0]
    heapq.heapify(heap)
    topo = []
    while heap:
        n = heapq.heappop(heap)
        topo.append(n)
        for m in succ[n]:
            indeg[m] -= 1
            if indeg[m] == 0:
                heapq.heappush(heap, m)
    if len(topo) != len(g.seq):
        raise ValueError('the induced graph has a cycle')
    ren = {n: i + 1 for i, n in enumerate(topo)}
    out = Graph()
    out.seq = {ren[n]: s for n, s in g.seq.items()}
    out.edges = {(ren[a], ren[b]) for a, b in g.edges}
    out.walks = {k: [ren[n] for n in w] for k, w in g.walks.items()}
    return out


def chop(g, max_node=MAX_NODE):
    """Cut nodes longer than max_node into pieces; ids stay in topological order."""
    if all(len(s) <= max_node for s in g.seq.values()):
        return g
    pieces, out, nid = {}, Graph(), 0
    for n in sorted(g.seq):
        s = g.seq[n]
        ids = []
        for i in range(0, max(1, len(s)), max_node):
            nid += 1
            out.seq[nid] = s[i:i + max_node]
            ids.append(nid)
        pieces[n] = ids
        for a, b in zip(ids, ids[1:]):
            out.edges.add((a, b))
    for a, b in g.edges:
        out.edges.add((pieces[a][-1], pieces[b][0]))
    out.walks = {k: [p for n in w for p in pieces[n]] for k, w in g.walks.items()}
    return out


def induce(rows, max_node=MAX_NODE):
    """Graph of an MSA {member: row}; checks that every member's walk spells its row's bases."""
    names = sorted(rows)
    g = chop(unchop(columns(rows, names)), max_node)
    for n in names:
        if ''.join(g.seq[x] for x in g.walks[n]) != rows[n].replace('-', ''):
            raise ValueError('walk of %s does not spell its row' % n)
    return g
