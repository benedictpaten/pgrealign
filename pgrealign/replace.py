"""Stage 5, `replace`: write a contig's graph with each realigned region's interior replaced, in one pass.

**Inputs.**
- The contig's graph: the GBZ from `prepare`, or a GFA. A GFA must have every S and L line before its
  first W line, as `vg convert -f` writes them.
- The `extract` outputs: `packages.jsonl.gz`, `occurrences.tsv.gz` and `interiors.tsv.gz`.
- The `realign` outputs: one `msas.jsonl.gz`, or several, such as one per shard. A region may appear in
  only one of them.

Only regions whose package is `ok` and whose MSA is `ok` are replaced. Every other region keeps its
subgraph.

**What changes, per replaced region (anchors L and R stay):**
- **Nodes.** The interior nodes are removed, and the MSA's induced graph (induce.py) is added.
  - By default ('dense') every node is renumbered: old nodes keep their order, each region's new nodes
    follow its left anchor, and ids run densely from `id_start`, or from the contig's smallest id. Every
    S, L, P and W line is renumbered, inside the regions or not, so the ids no longer match the input.
    merge needs every contig's ids to stay clear of every other contig's; the summary's `id_range` says
    where this contig's ids ended up.
  - With 'reuse', old nodes keep their ids, and new nodes take the region's freed ids, then ids from
    `id_start`. That range must be free in the whole genome: a new id that is a node of the contig, or
    that two new nodes would share, is an error. New nodes can land far from their regions, which slows
    tools that read the graph in windows of node ids, such as vg call.
- **Edges.** Edges touching a removed node are dropped. The induced graph's edges are added, plus
  the edges to the anchors: L->first for an allele or a left fragment, last->R for an allele or a
  right fragment, and L->R for an empty allele (one that deletes the whole interior) unless the input
  has it already.
- **Walks (W lines).** Every run recorded in `occurrences.tsv.gz` is rewritten to its row's walk:
  forward for orient '+', reversed and flipped for '-'.
  - A spanning run replaces the steps between the anchors.
  - A suffix run replaces the steps after the anchor to the walk's end.
  - A prefix run replaces the steps from the walk's start to the anchor.
  - An internal run replaces the whole walk.
- **Paths (P lines).** extract records runs of walks only, so a generic path cannot be rewritten
  through a region: a P line that visits a removed node is an error. Every other P line is kept.

**Other lines.** H lines are copied without vg's graph-name tags (NM, SG and TL), so the reference
samples (RS tag) are kept. The output is a different graph, and `vg gbwt -G` would otherwise record it
as a subgraph of the input. Lines of any other type are copied unchanged.

**Checks.**
- Each rewritten run spells the same sequence as before. The old steps are spelled from the input's
  node sequences, the new ones from the induced graph's.
- No walk visits a removed node outside a recorded run.
- A failure stops the stage and leaves no output: it means extract or realign is wrong, not the region.

The output is a GFA. `vg gbwt -G out.gfa --gbz-format -g out.gbz` turns it into a GBZ.
"""
import collections
import gzip
import json
import os

from . import induce as induce_mod
from .graph import GfaStream, parse_w, revcomp, walk_key, STEP
from .extract import read_packages

GRAPH_NAME_TAGS = ('NM:', 'SG:', 'TL:')     # H-line tags that name the graph and relate it to others
ORDER = 'replace needs every S and L line before the first W line, as vg convert -f writes them'
FLIP = {'+': '-', '-': '+'}


def load_plan(extract_dir, msas):
    """{region: (package, msa result)} for every region to replace, in package order. `msas` is the path
    of one msas.jsonl.gz, or a list of them; a region found twice is an error."""
    paths = [msas] if isinstance(msas, (str, os.PathLike)) else msas
    found, ok = {}, {}
    for path in paths:
        with gzip.open(path, 'rt') as f:
            for line in f:
                m = json.loads(line)
                if m['id'] in found:
                    raise ValueError('region %s appears twice in the msas files (%s, %s)' % (
                        m['id'], found[m['id']], path))
                found[m['id']] = path
                if m.get('status') == 'ok':
                    ok[m['id']] = m
    plan = {}
    for p in read_packages(os.path.join(extract_dir, 'packages.jsonl.gz'), status={'ok'}):
        if p['id'] in ok:
            plan[p['id']] = (p, ok[p['id']])
    return plan


def read_interiors(path, keep):
    """{region: its interior's node ids} for the regions in `keep`."""
    out = {}
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            rid, nodes = line.rstrip('\n').split('\t')
            if rid in keep:
                out[rid] = [int(x) for x in nodes.split(',')] if nodes else []
    return out


def read_occurrences(path, keep):
    """{walk key: [(first, end, region, type, orient, member)]} for the regions in `keep`, sorted."""
    by_walk = collections.defaultdict(list)
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            w, rid, typ, orient, a, b, mid = line.rstrip('\n').split('\t')
            if rid in keep:
                by_walk[w].append((int(a), int(b), rid, typ, orient, mid))
    for v in by_walk.values():
        v.sort()
    return by_walk


def old_node_ids(graph, vg=None, threads=4):
    """Every node id of the input graph, from its S lines (`vg convert -f -H` for a GBZ, so no haplotype
    walks are written). A GFA file is read only to its first W line, since replace needs every S line
    before it (ORDER), so that the walks, nearly all of the file, are read once, by rewrite. A node whose
    S line comes later is missed, and named as out of order by whichever check meets it first."""
    ids = []
    stream = GfaStream(graph, vg, threads, drop_haplotypes=True)
    with stream as f:
        for line in f:
            if line[0] == 'S':
                ids.append(int(line.split('\t', 2)[1]))
            elif line[0] == 'W' and stream.proc is None:
                break
    return ids


class RegionGraph:
    """The graph that replaces one region's interior (see region_graph)."""
    __slots__ = ('seq', 'edges', 'walks')

    def __init__(self, seq, edges, walks):
        self.seq, self.edges, self.walks = seq, edges, walks


def region_graph(pkg, rows):
    """The graph that replaces the interior of the region `pkg`, from its MSA `rows` {member: row}:
    - `seq`: {node: sequence} of the induced graph (induce.py). Its nodes are numbered -1, -2, ... in
      topological order, so they cannot be mistaken for the anchors, which keep the input's ids.
    - `edges`: {(from, from_orient, to, to_orient)}, all forward: the induced graph's edges, plus L->first
      for an allele or a left fragment, last->R for an allele or a right fragment, and L->R for an
      empty allele.
    - `walks`: {member: [node, ...]}, in the region's forward orientation. An empty member's is empty."""
    g = induce_mod.induce(rows)
    ids = {n: -(i + 1) for i, n in enumerate(sorted(g.seq))}
    L, R = pkg['L'], pkg['R']
    side = {f['id']: f['side'] for f in pkg['fragments']}
    edges = {(ids[a], '+', ids[b], '+') for a, b in g.edges}
    walks = {}
    for m, w in g.walks.items():
        gw = walks[m] = [ids[n] for n in w]
        sd = side.get(m, 'spanning')
        if not gw:
            if sd == 'spanning':
                edges.add((L, '+', R, '+'))
            continue
        if sd in ('spanning', 'left'):
            edges.add((L, '+', gw[0], '+'))
        if sd in ('spanning', 'right'):
            edges.add((gw[-1], '+', R, '+'))
    return RegionGraph({ids[n]: sq for n, sq in g.seq.items()}, edges, walks)


def number_dense(old_ids, removed, regions, id_start=None):
    """The 'dense' numbering. `old_ids` are the input's node ids, sorted; `removed` holds the replaced
    interiors' nodes; `regions` is [(region, left anchor, number of new nodes)], in plan order.

    The nodes that stay take consecutive ids in their old order, from `id_start` or else the smallest old
    id. Right after each left anchor come the new nodes of its regions, region by region, each region's
    in its own order. Returns ({old id: id} for the nodes that stay, {region: [id of its node -1, -2,
    ...]}, (first id, last id))."""
    after = collections.defaultdict(list)          # left anchor -> [(region, number of new nodes)]
    for rid, left, size in regions:
        after[left].append((rid, size))
    renumber, new_ids = {}, {}
    k = first = id_start if id_start is not None else (old_ids[0] if old_ids else 1)
    for n in old_ids:
        if n in removed:
            continue
        renumber[n] = k
        k += 1
        for rid, size in after.pop(n, ()):
            new_ids[rid] = list(range(k, k + size))
            k += size
    if after:
        raise ValueError('%d left anchors are not nodes of the graph, or their S lines come after a W line: '
                         '%s' % (len(after), ORDER))
    return renumber, new_ids, (first, k - 1)


def number_reuse(regions, interiors, id_start):
    """The 'reuse' numbering of the new nodes; the nodes that stay keep their ids. `regions` is
    [(region, number of new nodes)], in plan order, and `interiors` {region: its interior's node ids}.
    Each region's new nodes take its freed ids in increasing order, then the next ids from `id_start`.
    Returns ({region: [id of its node -1, -2, ...]}, the next unused id from `id_start`).

    Two new nodes given one id is an error: it happens when interiors overlap or when the ids from
    `id_start` run into freed ones. A new id that a node outside the regions holds is found by rewrite."""
    new_ids, next_id, owner = {}, id_start, {}
    for rid, size in regions:
        ids = sorted(interiors[rid])[:size]
        extra = size - len(ids)
        ids.extend(range(next_id, next_id + extra))
        next_id += extra
        for n in ids:
            if n in owner:
                raise ValueError('id %d would go to new nodes of region %s and of region %s: interiors must '
                                 'not overlap, and the ids from id_start (%d) on must be free' % (
                                     n, owner[n], rid, id_start))
            owner[n] = rid
        new_ids[rid] = ids
    return new_ids, next_id


def out_of_order(n):
    """The error for node `n`, which has no id under the dense numbering: its S line comes after a W line."""
    return ValueError('node %d has no S line before the first W line: %s' % (n, ORDER))


class Splice:
    """What `rewrite` applies to the input graph, in the output's ids:
    - `seq` {node: sequence} and `edges` {(from, from_orient, to, to_orient)}: the new nodes and edges;
    - `walks` {(region, member): [node, ...]} and `spelling` {(region, member): sequence}: each member's
      new walk, and the sequence it and every run recorded for it must spell, region-forward;
    - `removed`: {node: region} for the replaced interiors' nodes, in the input's ids;
    - `runs`: {walk key: [(first, end, region, type, orient, member)]}, from read_occurrences;
    - `renumber`: {input id: output id} for the nodes that stay, or None when they keep their ids."""
    __slots__ = ('seq', 'edges', 'walks', 'spelling', 'removed', 'runs', 'renumber')

    def __init__(self, removed, runs, renumber):
        self.seq, self.edges, self.walks, self.spelling = {}, set(), {}, {}
        self.removed, self.runs, self.renumber = removed, runs, renumber

    def add(self, rid, pkg, rg, ids):
        """Add region `rid`'s graph `rg` (region_graph), whose node -i gets the id ids[i - 1]."""
        def out_id(n):
            if n < 0:
                return ids[-n - 1]
            return n if self.renumber is None else self.renumbered(n)
        for n, sq in rg.seq.items():
            self.seq[out_id(n)] = sq
        for a, ao, b, bo in rg.edges:
            self.edges.add((out_id(a), ao, out_id(b), bo))
        for m, w in rg.walks.items():
            self.walks[(rid, m)] = [out_id(n) for n in w]
        for x in pkg['alleles'] + pkg['fragments']:
            self.spelling[(rid, x['id'])] = x['seq'].upper()

    def renumbered(self, n):
        """Input node `n`'s id in the output, under `renumber`. old_node_ids gave every node whose S line
        comes before the first W line an id, so a node without one is out of ORDER."""
        try:
            return self.renumber[n]
        except KeyError:
            raise out_of_order(n) from None

    def check_outside(self, key, toks):
        """Steps `toks` of walk `key` lie outside every run, so they must not visit a removed node."""
        if self.removed and any(int(n) in self.removed for _, n in toks):
            n = next(int(n) for _, n in toks if int(n) in self.removed)
            raise ValueError('walk %s visits node %d of region %s outside a recorded run' % (
                key, n, self.removed[n]))

    def copied(self, key, toks):
        """The steps `toks` of walk `key`, outside every run, in the output's ids."""
        self.check_outside(key, toks)
        if self.renumber is None:
            return ''.join(o + n for o, n in toks)
        # every step of every walk comes through here, so the lookup is the dict's own (see renumbered)
        renumber = self.renumber
        try:
            return ''.join(o + str(renumber[int(n)]) for o, n in toks)
        except KeyError as e:
            raise out_of_order(e.args[0]) from None

    def path_line(self, line):
        """P line `line`, renumbered. It must not visit a removed node: extract records no runs of P
        lines, so their steps through a region cannot be rewritten."""
        x = line.rstrip('\n').split('\t')
        steps = x[2].split(',')
        for s in steps:
            rid = self.removed.get(int(s[:-1]))
            if rid is not None:
                raise ValueError('path %s visits node %s of region %s: replace rewrites walks (W lines), '
                                 'not paths (P lines)' % (x[1], s[:-1], rid))
        if self.renumber is None:
            return line
        try:
            x[2] = ','.join('%d%s' % (self.renumber[int(s[:-1])], s[-1]) for s in steps)
        except KeyError as e:
            raise out_of_order(e.args[0]) from None
        return '\t'.join(x) + '\n'


def rewrite_walk(line, sp, old_seq, counts):
    """W line `line` with its recorded runs replaced by their members' new walks and its other steps
    renumbered. `old_seq` holds the removed nodes' sequences, to check what each run spelled."""
    w = parse_w(line)
    key = walk_key(w)
    runs = sp.runs.get(key)
    toks = STEP.findall(w.steps)
    if not runs:
        if sp.renumber is None:
            sp.check_outside(key, toks)
            return line
        x = line.rstrip('\n').split('\t')
        x[6] = sp.copied(key, toks)
        return '\t'.join(x) + '\n'
    parts, pos = [], 0
    for a, b, rid, typ, orient, mid in runs:
        if a < pos:
            raise ValueError('overlapping runs in walk %s' % key)
        parts.append(sp.copied(key, toks[pos:a]))
        try:
            old_spelling = ''.join(old_seq[int(n)] if o == '>' else revcomp(old_seq[int(n)])
                                   for o, n in toks[a:b])
        except KeyError as e:
            n = e.args[0]
            if n in sp.removed:
                raise ValueError('walk %s comes before the S line of node %d: %s' % (key, n, ORDER)) from None
            raise ValueError('run %s %d-%d of walk %s visits node %d, which no replaced interior holds: the '
                             'extract outputs are not from this graph' % (rid, a, b, key, n)) from None
        new_walk = sp.walks[(rid, mid)]
        new_spelling = ''.join(sp.seq[n] for n in new_walk)
        want = sp.spelling[(rid, mid)]
        if (old_spelling if orient == '+' else revcomp(old_spelling)) != want or new_spelling != want:
            raise ValueError('run %s %d-%d of walk %s does not spell member %s' % (rid, a, b, key, mid))
        if orient == '+':
            parts.append(''.join('>%d' % n for n in new_walk))
        else:
            parts.append(''.join('<%d' % n for n in reversed(new_walk)))
        counts['runs_rewritten'] += 1
        pos = b
    parts.append(sp.copied(key, toks[pos:]))
    x = line.rstrip('\n').split('\t')
    x[6] = ''.join(parts)
    counts['walks_rewritten'] += 1
    return '\t'.join(x) + '\n'


def rewrite(graph, out_gfa, sp, vg=None, threads=4):
    """Stream `graph` (a GBZ or GFA) to `out_gfa` with the Splice `sp` applied, as the module docstring
    describes: removed nodes and their edges dropped, the new nodes and edges written before the first W
    line, every W line rewritten by rewrite_walk, and every other line that names a node renumbered when
    `sp.renumber` is set. Returns the counts of edges removed, runs rewritten and walks rewritten. On a
    failure, nothing is left at `out_gfa`."""
    renumber, removed = sp.renumber, sp.removed
    # an empty allele's L->R edge is usually in the input already, from the walks that spelled the
    # allele; then the input's line is the one written
    new_edges = set(sp.edges)
    direct_ends = {n for a, _, b, _ in sp.edges if a not in sp.seq and b not in sp.seq for n in (a, b)}
    counts = collections.Counter()
    old_seq = {}            # sequences of removed nodes, to spell the old runs
    new_written = False
    tmp = out_gfa + '.tmp'
    try:
        with GfaStream(graph, vg, threads) as f, open(tmp, 'w') as out:
            for line in f:
                t = line[0]
                if new_written and t in ('S', 'L'):
                    raise ValueError('%s has an %s line after a W line: %s' % (graph, t, ORDER))
                if t == 'S':
                    x = line.split('\t', 3)
                    n = int(x[1])
                    if n in removed:
                        old_seq[n] = x[2].rstrip('\n').upper()
                        continue
                    if renumber is not None:
                        out.write('S\t%d\t%s' % (renumber[n], x[2] if len(x) == 3 else x[2] + '\t' + x[3]))
                    else:
                        if n in sp.seq:
                            raise ValueError('a new node was given id %d, which a node of the graph has: '
                                             'the ids from id_start on must be free' % n)
                        out.write(line)
                elif t == 'L':
                    x = line.split('\t', 5)
                    a, b = int(x[1]), int(x[3])
                    if a in removed or b in removed:
                        counts['edges_removed'] += 1
                        continue
                    if renumber is not None:
                        a, b = sp.renumbered(a), sp.renumbered(b)
                        x[1], x[3] = str(a), str(b)
                        line = '\t'.join(x)
                    if a in direct_ends:
                        new_edges.discard((a, x[2], b, x[4]))
                        new_edges.discard((b, FLIP[x[4]], a, FLIP[x[2]]))
                    out.write(line)
                elif t == 'W':
                    if not new_written:
                        write_new(out, sp.seq, new_edges)
                        new_written = True
                    out.write(rewrite_walk(line, sp, old_seq, counts))
                elif t == 'P':
                    out.write(sp.path_line(line))
                elif t == 'H':
                    tags = [x for x in line.rstrip('\n').split('\t')[1:] if not x.startswith(GRAPH_NAME_TAGS)]
                    if tags:
                        out.write('\t'.join(['H'] + tags) + '\n')
                else:
                    out.write(line)
            if not new_written:
                write_new(out, sp.seq, new_edges)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    os.replace(tmp, out_gfa)
    return counts


def write_new(out, seq, edges):
    """The new nodes' S lines and the new edges' L lines, each sorted."""
    for n in sorted(seq):
        out.write('S\t%d\t%s\n' % (n, seq[n]))
    for a, ao, b, bo in sorted(edges):
        out.write('L\t%d\t%s\t%d\t%s\t0M\n' % (a, ao, b, bo))


def replace(graph, extract_dir, msas, out_gfa, id_start=None, vg=None, threads=4, id_mode='dense'):
    """Write the contig `graph` to `out_gfa` with every planned region's interior replaced, and return the
    summary. `msas` is one msas.jsonl.gz or a list of them. `id_mode`:
    - 'dense' (default): every node is renumbered, in the order of the old ids, with each region's new
      nodes placed right after its left anchor, from `id_start` or else the contig's smallest id. Every
      walk and path is renumbered, not only those through the regions. A region's nodes stay next to its
      neighbours, which keeps id-windowed readers (vg call's read windows) local.
    - 'reuse': old nodes keep their ids; new ones take the region's freed ids, then ids from `id_start`,
      which must be free. Ids stay stable, but new nodes can land far from their regions, so readers that
      fetch by windows of node ids straddle windows and slow down."""
    if id_mode not in ('dense', 'reuse'):
        raise ValueError('id_mode must be dense or reuse')
    if id_mode == 'reuse' and id_start is None:
        raise ValueError("id_mode 'reuse' needs id_start")
    plan = load_plan(extract_dir, msas)
    interiors = read_interiors(os.path.join(extract_dir, 'interiors.tsv.gz'), plan)
    occ = read_occurrences(os.path.join(extract_dir, 'occurrences.tsv.gz'), plan)

    removed = {n: rid for rid in plan for n in interiors[rid]}
    for rid, (pkg, _) in plan.items():
        for anchor in (pkg['L'], pkg['R']):
            if anchor in removed:
                raise ValueError('anchor %d of region %s is in the interior of region %s' % (
                    anchor, rid, removed[anchor]))

    stats = collections.Counter()
    graphs = {}
    for rid, (pkg, msa) in plan.items():
        graphs[rid] = region_graph(pkg, msa['rows'])
        stats['regions'] += 1
        stats['nodes_removed'] += len(interiors[rid])
        stats['nodes_added'] += len(graphs[rid].seq)
    if id_mode == 'dense':
        renumber, new_ids, id_range = number_dense(
            sorted(old_node_ids(graph, vg, threads)), removed,
            [(rid, pkg['L'], len(graphs[rid].seq)) for rid, (pkg, _) in plan.items()], id_start)
        stats['id_range'] = list(id_range)
    else:
        renumber = None
        new_ids, next_id = number_reuse([(rid, len(graphs[rid].seq)) for rid in plan], interiors, id_start)
        if next_id > id_start:
            stats['ids_above_start'] = next_id - id_start
        stats['next_id'] = next_id
    sp = Splice(removed, occ, renumber)
    for rid, (pkg, _) in plan.items():
        sp.add(rid, pkg, graphs[rid], new_ids[rid])
    stats.update(rewrite(graph, out_gfa, sp, vg, threads))
    stats['id_mode'] = id_mode
    return dict(stats)


def to_gbz(gfa, gbz, vg=None, keep_gfa=False):
    """Build `gbz` from replace's GFA with `vg gbwt -G`, then remove the GFA unless `keep_gfa`."""
    import subprocess
    from .graph import find_vg
    subprocess.run([find_vg(vg), 'gbwt', '-G', gfa, '--gbz-format', '-g', gbz + '.tmp'], check=True)
    os.replace(gbz + '.tmp', gbz)
    if not keep_gfa:
        os.remove(gfa)
