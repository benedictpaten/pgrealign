"""Stage 5, `replace`: write a contig's graph with each realigned region's interior replaced, in one pass.

**Inputs.**
- The contig's graph: the GBZ from `prepare`, or a GFA.
- The `extract` outputs: `packages.jsonl.gz`, `occurrences.tsv.gz` and `interiors.tsv.gz`.
- The `realign` outputs: `msas.jsonl.gz`.

Only regions whose package is `ok` and whose MSA is `ok` are replaced. Every other region keeps its
subgraph.

**What changes, per replaced region R (anchors L and R stay):**
- **Nodes.** The interior nodes are removed, and the MSA's induced graph (induce.py) is added.
  - New nodes reuse the removed nodes' ids, smallest first in topological order, so ids stay dense
    and near their neighbours.
  - When the new graph has more nodes than were removed, the rest take ids from `id_start` upwards.
    That range must be free in the whole genome.
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


def replace(graph, extract_dir, msas_path, out_gfa, id_start, vg=None, threads=4):
    plan = load_plan(extract_dir, msas_path)
    interiors = read_interiors(os.path.join(extract_dir, 'interiors.tsv.gz'), plan)
    occ = read_occurrences(os.path.join(extract_dir, 'occurrences.tsv.gz'), plan)

    # induce every region's graph and give its nodes global ids
    new_seq, new_edges = {}, set()      # global id -> seq; (from, from_orient, to, to_orient)
    walks = {}                          # (region, member) -> [global ids]
    member_seq = {}                     # (region, member) -> sequence (region-forward)
    removed = set()
    next_id = id_start
    stats = collections.Counter()
    for rid, (pkg, msa) in plan.items():
        g = induce_mod.induce(msa['rows'])
        freed = sorted(interiors[rid])
        ids = {}
        for k, n in enumerate(sorted(g.seq)):
            if k < len(freed):
                ids[n] = freed[k]
            else:
                ids[n] = next_id
                next_id += 1
                stats['ids_above_start'] += 1
        removed.update(freed)
        for n, s in g.seq.items():
            new_seq[ids[n]] = s
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
        stats['nodes_removed'] += len(freed)
        stats['nodes_added'] += len(g.seq)
    reused = removed & set(new_seq)

    old_seq = {}       # sequences of removed nodes, to spell the old runs
    with GfaStream(graph, vg, threads) as f, open(out_gfa + '.tmp', 'w') as out:
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
                if n in new_seq:
                    raise ValueError('new node id %d is already in the graph' % n)
                out.write(line)
            elif t == 'L':
                x = line.split('\t', 5)
                if int(x[1]) in removed or int(x[3]) in removed:
                    stats['edges_removed'] += 1
                    continue
                out.write(line)
            elif t == 'W':
                if new_seq is not None:
                    # all S and L lines have been read: write the new ones before the first walk
                    for n in sorted(new_seq):
                        out.write('S\t%d\t%s\n' % (n, new_seq[n]))
                    for a, ao, b, bo in sorted(new_edges):
                        out.write('L\t%d\t%s\t%d\t%s\t0M\n' % (a, ao, b, bo))
                    new_seq_written, new_seq = new_seq, None
                w = parse_w(line)
                key = walk_key(w)
                runs = occ.get(key)
                if not runs:
                    if removed and any(int(n) in removed for _, n in STEP.findall(w.steps)):
                        raise ValueError('walk %s visits a replaced region outside a recorded run' % key)
                    out.write(line)
                    continue
                toks = STEP.findall(w.steps)
                parts, pos = [], 0
                for a, b, rid, typ, orient, mid in runs:
                    if a < pos:
                        raise ValueError('overlapping runs in walk %s' % key)
                    if any(int(n) in removed for _, n in toks[pos:a]):
                        raise ValueError('walk %s visits a replaced region outside a recorded run' % key)
                    parts.append(''.join(o + n for o, n in toks[pos:a]))
                    old = ''.join(old_seq[int(n)] if o == '>' else revcomp(old_seq[int(n)]) for o, n in toks[a:b])
                    gw = walks[(rid, mid)]
                    newf = ''.join(new_seq_written[n] for n in gw)
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
                parts.append(''.join(o + n for o, n in toks[pos:]))
                x = line.rstrip('\n').split('\t')
                x[6] = ''.join(parts)
                out.write('\t'.join(x) + '\n')
                stats['walks_rewritten'] += 1
            else:
                out.write(line)
        if new_seq is not None:
            for n in sorted(new_seq):
                out.write('S\t%d\t%s\n' % (n, new_seq[n]))
            for a, ao, b, bo in sorted(new_edges):
                out.write('L\t%d\t%s\t%d\t%s\t0M\n' % (a, ao, b, bo))
    os.replace(out_gfa + '.tmp', out_gfa)
    stats['ids_reused'] = len(reused)
    stats['next_id'] = next_id
    return dict(stats)
