"""Stage 5, `replace`: write a contig's graph with each realigned region's interior replaced, in one pass.

**Inputs.**
- The contig's graph: the GBZ from `prepare`, or a GFA.
- The `extract` outputs: `packages.jsonl.gz`, `occurrences.tsv.gz` and `interiors.tsv.gz`.
- The `realign` outputs: `msas.jsonl.gz`.

Only regions whose package is `ok` and whose MSA is `ok` are replaced. Every other region keeps its
subgraph.

**What changes, per replaced region R (anchors L and R stay):**
- **Nodes.** The interior nodes are removed, and the MSA's induced graph (induce.py) is added.
  - By default ('dense') every node is renumbered: old nodes keep their order, each region's new nodes
    follow its left anchor, and ids run densely from the contig's smallest id. The ids then no
    longer match the input, and contigs need a global renumbering before they are merged.
  - With 'reuse', old nodes keep their ids, and new nodes take the region's freed ids, then ids from
    `id_start`. That range must be free in the whole genome. On chr20 it put 182,106 nodes far from
    their regions, and vg call's id-windowed read fetches slowed about 8x.
- **Edges.** Edges touching a removed node are dropped. The induced graph's edges are added, plus
  the edges to the anchors: L->first for an allele or a left fragment, last->R for an allele or a
  right fragment, and L->R for an empty allele.
- **Walks.** Every run recorded in `occurrences.tsv.gz` is rewritten to its row's walk: forward for
  orient '+', reversed and flipped for '-'.
  - A spanning run replaces the steps between the anchors.
  - A suffix run replaces the steps after the anchor to the walk's end.
  - A prefix run replaces the steps from the walk's start to the anchor.
  - An internal run replaces the whole walk.

**Checks.**
- Each rewritten run spells the same sequence as before. The old steps are spelled from the input's
  node sequences, the new ones from the induced graph's.
- No walk visits a removed node outside a recorded run.
- A failure stops the stage: it means extract or realign is wrong, not the region.

The output is a GFA. `vg gbwt -G out.gfa --gbz-format -g out.gbz` turns it into a GBZ. The H line,
and with it the reference samples (RS tag), is copied.
"""
import collections
import gzip
import json
import os

from . import induce as induce_mod
from .graph import GfaStream, parse_w, revcomp, walk_key, STEP
from .extract import read_packages


def load_plan(extract_dir, msas_path):
    """{region: (package, msa result)} for every region to replace."""
    msas = {}
    with gzip.open(msas_path, 'rt') as f:
        for line in f:
            m = json.loads(line)
            if m.get('status') == 'ok':
                msas[m['id']] = m
    plan = {}
    for p in read_packages(os.path.join(extract_dir, 'packages.jsonl.gz'), status={'ok'}):
        if p['id'] in msas:
            plan[p['id']] = (p, msas[p['id']])
    return plan


def read_interiors(path, keep):
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
    walks are written)."""
    ids = []
    with GfaStream(graph, vg, threads, drop_haplotypes=True) as f:
        for line in f:
            if line[0] == 'S':
                ids.append(int(line.split('\t', 2)[1]))
    return ids


def replace(graph, extract_dir, msas_path, out_gfa, id_start=None, vg=None, threads=4, id_mode='dense'):
    """Write the contig with every planned region's interior replaced. `id_mode`:
    - 'dense' (default): every node is renumbered, in the order of the old ids, with each region's new
      nodes placed right after its left anchor, from the contig's smallest id. A region's nodes stay
      next to its neighbours, which keeps id-windowed readers (vg call's read windows) local.
    - 'reuse': old nodes keep their ids; new ones take the region's freed ids, then ids from `id_start`.
      Ids stay stable, but on chr20 182,106 new nodes landed far from their regions, which made vg call's
      read fetches straddle windows (a 54-minute short-read call instead of 6)."""
    if id_mode not in ('dense', 'reuse'):
        raise ValueError('id_mode must be dense or reuse')
    if id_mode == 'reuse' and id_start is None:
        raise ValueError("id_mode 'reuse' needs id_start")
    plan = load_plan(extract_dir, msas_path)
    interiors = read_interiors(os.path.join(extract_dir, 'interiors.tsv.gz'), plan)
    occ = read_occurrences(os.path.join(extract_dir, 'occurrences.tsv.gz'), plan)

    # induce every region's graph; new nodes get temporary negative ids until the numbering is fixed
    new_seq, new_edges = {}, set()      # temp id -> seq; (from, from_orient, to, to_orient)
    walks = {}                          # (region, member) -> [temp ids]
    member_seq = {}                     # (region, member) -> sequence (region-forward)
    removed = set()
    after_anchor = collections.defaultdict(list)   # left anchor -> temp ids, in topological order
    region_tmp = {}
    tmp = -1
    stats = collections.Counter()
    for rid, (pkg, msa) in plan.items():
        g = induce_mod.induce(msa['rows'])
        ids = {}
        for n in sorted(g.seq):
            ids[n] = tmp
            tmp -= 1
        region_tmp[rid] = [ids[n] for n in sorted(g.seq)]
        after_anchor[pkg['L']].extend(region_tmp[rid])
        removed.update(interiors[rid])
        for n, sq in g.seq.items():
            new_seq[ids[n]] = sq
        for a, b in g.edges:
            new_edges.add((ids[a], '+', ids[b], '+'))
        L, R = pkg['L'], pkg['R']
        side = {f['id']: f['side'] for f in pkg['fragments']}
        for m, w in g.walks.items():
            gw = [ids[n] for n in w]
            walks[(rid, m)] = gw
            sd = side.get(m, 'spanning')
            if not gw:
                if sd == 'spanning':
                    new_edges.add((L, '+', R, '+'))
                continue
            if sd in ('spanning', 'left'):
                new_edges.add((L, '+', gw[0], '+'))
            if sd in ('spanning', 'right'):
                new_edges.add((gw[-1], '+', R, '+'))
        for x in pkg['alleles'] + pkg['fragments']:
            member_seq[(rid, x['id'])] = x['seq'].upper()
        stats['regions'] += 1
        stats['nodes_removed'] += len(interiors[rid])
        stats['nodes_added'] += len(g.seq)

    # the final numbering
    M = {}
    if id_mode == 'dense':
        old = sorted(old_node_ids(graph, vg, threads))
        k = old[0] if old else 1
        kept = 0
        for n in old:
            if n in removed:
                continue
            M[n] = k
            k += 1
            kept += 1
            for t in after_anchor.pop(n, ()):
                M[t] = k
                k += 1
        if after_anchor:
            raise ValueError('%d left anchors are not nodes of the graph' % len(after_anchor))
        stats['id_range'] = [old[0] if old else 1, k - 1]
        mp = M.__getitem__
    else:
        next_id = id_start
        for rid, lst in region_tmp.items():
            freed = sorted(interiors[rid])
            for i, t in enumerate(lst):
                if i < len(freed):
                    M[t] = freed[i]
                else:
                    M[t] = next_id
                    next_id += 1
                    stats['ids_above_start'] += 1
        stats['next_id'] = next_id
        mp = lambda n: M.get(n, n)
    new_seq = {mp(t): sq for t, sq in new_seq.items()}
    new_edges = {(mp(a), ao, mp(b), bo) for a, ao, b, bo in new_edges}
    walks = {k2: [mp(t) for t in w] for k2, w in walks.items()}

    def steps_out(toks):
        if id_mode == 'dense':
            return ''.join(o + str(M[int(n)]) for o, n in toks)
        return ''.join(o + n for o, n in toks)

    old_seq = {}       # sequences of removed nodes, to spell the old runs
    new_written = False
    with GfaStream(graph, vg, threads) as f, open(out_gfa + '.tmp', 'w') as out:
        def write_new():
            for n in sorted(new_seq):
                out.write('S\t%d\t%s\n' % (n, new_seq[n]))
            for a, ao, b, bo in sorted(new_edges):
                out.write('L\t%d\t%s\t%d\t%s\t0M\n' % (a, ao, b, bo))
        for line in f:
            t = line[0]
            if t == 'H':
                out.write(line)
            elif t == 'S':
                x = line.split('\t', 3)
                n = int(x[1])
                if n in removed:
                    old_seq[n] = x[2].rstrip('\n').upper()
                    continue
                if id_mode == 'dense':
                    out.write('S\t%d\t%s' % (M[n], x[2] if len(x) == 3 else x[2] + '\t' + x[3]))
                else:
                    if n in new_seq:
                        raise ValueError('new node id %d is already in the graph' % n)
                    out.write(line)
            elif t == 'L':
                x = line.split('\t', 5)
                a, b = int(x[1]), int(x[3])
                if a in removed or b in removed:
                    stats['edges_removed'] += 1
                    continue
                if id_mode == 'dense':
                    x[1], x[3] = str(M[a]), str(M[b])
                    out.write('\t'.join(x))
                else:
                    out.write(line)
            elif t == 'W':
                if not new_written:
                    write_new()       # all S and L lines have been read
                    new_written = True
                w = parse_w(line)
                key = walk_key(w)
                runs = occ.get(key)
                toks = STEP.findall(w.steps)
                if not runs:
                    if removed and any(int(n) in removed for _, n in toks):
                        raise ValueError('walk %s visits a replaced region outside a recorded run' % key)
                    if id_mode == 'dense':
                        x = line.rstrip('\n').split('\t')
                        x[6] = steps_out(toks)
                        out.write('\t'.join(x) + '\n')
                    else:
                        out.write(line)
                    continue
                parts, pos = [], 0
                for a, b, rid, typ, orient, mid in runs:
                    if a < pos:
                        raise ValueError('overlapping runs in walk %s' % key)
                    if any(int(n) in removed for _, n in toks[pos:a]):
                        raise ValueError('walk %s visits a replaced region outside a recorded run' % key)
                    parts.append(steps_out(toks[pos:a]))
                    old = ''.join(old_seq[int(n)] if o == '>' else revcomp(old_seq[int(n)]) for o, n in toks[a:b])
                    gw = walks[(rid, mid)]
                    newf = ''.join(new_seq[n] for n in gw)
                    want = member_seq[(rid, mid)]
                    if (old if orient == '+' else revcomp(old)) != want or newf != want:
                        raise ValueError('run %s %d-%d of walk %s does not spell member %s' % (rid, a, b, key, mid))
                    if orient == '+':
                        parts.append(''.join('>%d' % n for n in gw))
                    else:
                        parts.append(''.join('<%d' % n for n in reversed(gw)))
                    stats['runs_rewritten'] += 1
                    pos = b
                if any(int(n) in removed for _, n in toks[pos:]):
                    raise ValueError('walk %s visits a replaced region outside a recorded run' % key)
                parts.append(steps_out(toks[pos:]))
                x = line.rstrip('\n').split('\t')
                x[6] = ''.join(parts)
                out.write('\t'.join(x) + '\n')
                stats['walks_rewritten'] += 1
            else:
                out.write(line)
        if not new_written:
            write_new()
    os.replace(out_gfa + '.tmp', out_gfa)
    stats['id_mode'] = id_mode
    return dict(stats)
