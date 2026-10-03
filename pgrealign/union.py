"""Realign on one panel and apply the result to a graph holding other haplotypes: `union` and `project`.

The protocol (vntr-eval results, 4u-4x) aligns the full panel, then applies the alignment to a graph
built from a subset or a sample of it. An example is the 34-haplotype graph `vg call` used, whose
recombinant haplotypes are walks in the full graph.

- **`union`** merges, region by region, the packages `extract` wrote for the panel's graph (`base`)
  and for the target graph (`add`), at the same anchors.
  - A target allele with the same sequence as a panel allele is that allele.
  - A target allele the panel lacks joins with weight 0, so it never pulls the centre.
  - Fragments of both are kept, deduplicated by side and sequence.
  - Each union package records `target_map`, from the target's member ids to the union's.
- **`project`** turns the union's MSAs into MSAs for the target. It keeps the target members' rows,
  under the target's member ids, and drops the columns that are gaps in all of them.
"""
import collections
import gzip
import json
import os

from .extract import read_packages
from .realign import drop_gap_columns


def union_packages(base_path, add_path, out_path):
    add = {p['id']: p for p in read_packages(add_path)}
    stats = collections.Counter()
    with gzip.open(out_path + '.tmp', 'wt') as out:
        for p in read_packages(base_path):
            q = add.pop(p['id'], None)
            if q is None:
                stats['base_only'] += 1
                continue
            if (p['L'], p['R']) != (q['L'], q['R']):
                raise ValueError('region %s has different anchors in the two extracts' % p['id'])
            alleles = [dict(a) for a in p['alleles']]
            by_seq = {a['seq']: a['id'] for a in alleles}
            frags = [dict(f) for f in p['fragments']]
            by_frag = {(f['side'], f['seq']): f['id'] for f in frags}
            tmap = {}
            for a in q['alleles']:
                if a['seq'] not in by_seq:
                    nid = 't' + a['id']
                    alleles.append({'id': nid, 'weight': 0, 'seq': a['seq']})
                    by_seq[a['seq']] = nid
                    stats['target_only_alleles'] += 1
                tmap[a['id']] = by_seq[a['seq']]
            for f in q['fragments']:
                k = (f['side'], f['seq'])
                if k not in by_frag:
                    nid = 't' + f['id']
                    frags.append({'id': nid, 'weight': 0, 'side': f['side'], 'seq': f['seq']})
                    by_frag[k] = nid
                tmap[f['id']] = by_frag[k]
            st = 'ok' if p['status'] == 'ok' and q['status'] in ('ok', 'invariant') else p['status']
            if q['status'] == 'complex':
                st = 'complex'
            u = dict(p, alleles=alleles, fragments=frags, status=st, target_map=tmap,
                     target_status=q['status'])
            out.write(json.dumps(u, separators=(',', ':')) + '\n')
            stats['regions'] += 1
    os.replace(out_path + '.tmp', out_path)
    stats['add_only'] = len(add)
    return dict(stats)


def project(union_path, msas_path, out_path):
    tmaps = {p['id']: p['target_map'] for p in read_packages(union_path)}
    stats = collections.Counter()
    with gzip.open(msas_path, 'rt') as f, gzip.open(out_path + '.tmp', 'wt') as out:
        for line in f:
            m = json.loads(line)
            if m.get('status') != 'ok' or m['id'] not in tmaps:
                stats['skipped'] += 1
                continue
            tmap = tmaps[m['id']]
            rows = drop_gap_columns({t: m['rows'][u] for t, u in tmap.items()})
            out.write(json.dumps(dict(m, rows=rows, projected_from=len(m['rows'])), separators=(',', ':')) + '\n')
            stats['projected'] += 1
    os.replace(out_path + '.tmp', out_path)
    return dict(stats)
