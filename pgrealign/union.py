"""Realign on one panel and apply the result to a graph holding other haplotypes: `union` and `project`.

`extract` is run at the same anchors on the panel's graph (`base`) and on the target graph (`add`),
for example a graph built from a subset or a sample of the panel, whose haplotypes are walks in the
panel's graph. The union of the two extracts is realigned, and its alignments are projected onto the
target, whose own `replace` splices them in.

- **`union`** merges the two extracts region by region (`union_package`).
  - Only regions in both are kept. A region in only one is dropped and counted as `base_only` or
    `add_only`. A region's two packages must have the same anchors.
  - A target allele with the same sequence as a panel allele is that allele, and the panel's weight
    stands. A target allele the panel lacks joins with weight 0, so the centre is chosen by the panel
    alone. Target fragments join in the same way, matched by side and sequence. A joining member's id
    is `t` followed by its id in the target.
  - `target_map` maps each of the target's member ids to the union's.
  - The status is the panel's, except that a region complex in the target is complex: replace cannot
    rewrite the target's runs there. So a region invariant in the panel stays invariant, and is not
    realigned, even when the target has alleles the panel lacks. `target_status` is the target's own.
- **`project`** turns the union's `ok` MSAs into MSAs for the target (`project_rows`). It keeps the
  rows of the target's members, under the target's member ids, and drops the columns that are gaps in
  all of them. `columns` is the projected width, `projected_from` the union's row count, and every
  other field is the union's. The target's `replace` then rewrites the regions that are `ok` in its own
  extract and have a projected MSA.
"""
import collections
import gzip
import json
import os

from .extract import read_packages
from .star import drop_gap_columns


def _join(panel, target, key, target_map):
    """The panel's members, then the target's members that no panel member matches by `key`, with
    weight 0. Records in `target_map` the member each target member became."""
    out = list(panel)
    by_key = {key(m): m['id'] for m in panel}
    for m in target:
        k = key(m)
        if k not in by_key:
            by_key[k] = 't' + m['id']
            out.append(dict(m, id=by_key[k], weight=0))
        target_map[m['id']] = by_key[k]
    return out


def union_package(p, q):
    """The union of one region's panel package `p` and target package `q` (see the module docstring)."""
    if (p['L'], p['R']) != (q['L'], q['R']):
        raise ValueError('region %s has different anchors in the two extracts' % p['id'])
    tmap = {}
    alleles = _join(p['alleles'], q['alleles'], lambda a: a['seq'], tmap)
    frags = _join(p['fragments'], q['fragments'], lambda f: (f['side'], f['seq']), tmap)
    status = 'complex' if q['status'] == 'complex' else p['status']
    return dict(p, alleles=alleles, fragments=frags, status=status, target_map=tmap, target_status=q['status'])


def union_packages(base_path, add_path, out_path):
    """Write the union of the packages at `base_path` (the panel's) and `add_path` (the target's) to
    `out_path`, in the panel's order. Returns counts."""
    add = {q['id']: q for q in read_packages(add_path)}
    stats = collections.Counter()
    with gzip.open(out_path + '.tmp', 'wt') as out:
        for p in read_packages(base_path):
            q = add.pop(p['id'], None)
            if q is None:
                stats['base_only'] += 1
                continue
            u = union_package(p, q)
            out.write(json.dumps(u, separators=(',', ':')) + '\n')
            stats['regions'] += 1
            stats['target_only_alleles'] += len(u['alleles']) - len(p['alleles'])
            stats['target_only_fragments'] += len(u['fragments']) - len(p['fragments'])
    os.replace(out_path + '.tmp', out_path)
    stats['add_only'] = len(add)
    return dict(stats)


def project_rows(rows, target_map):
    """A union MSA's rows {union member: row} as the target's {target member: row}, without the columns
    that are gaps in every target row."""
    return drop_gap_columns({t: rows[u] for t, u in target_map.items()})


def project(union_path, msas_path, out_path):
    """Write the union's MSAs at `msas_path`, projected onto the target, to `out_path`. Returns counts."""
    tmaps = {p['id']: p['target_map'] for p in read_packages(union_path)}
    stats = collections.Counter()
    with gzip.open(msas_path, 'rt') as f, gzip.open(out_path + '.tmp', 'wt') as out:
        for line in f:
            m = json.loads(line)
            if m.get('status') != 'ok' or m['id'] not in tmaps:
                stats['skipped'] += 1
                continue
            rows = project_rows(m['rows'], tmaps[m['id']])
            m = dict(m, rows=rows, columns=len(next(iter(rows.values()), '')), projected_from=len(m['rows']))
            out.write(json.dumps(m, separators=(',', ':')) + '\n')
            stats['projected'] += 1
    os.replace(out_path + '.tmp', out_path)
    return dict(stats)
