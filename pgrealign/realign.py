"""Stage 3, `realign`: one multiple alignment per region package, by the medoid star (star.py).

`realign_package` aligns one package from extract.py; `run` aligns a whole packages.jsonl.gz in
parallel and writes msas.jsonl.gz.

**Limits.** A region the star cannot be built for keeps its old alignment: replace leaves it as it is.
Its status says why.
- `too_big`: the unbanded pairwise alignment would exceed `max_cells` dynamic-programming cells, as
  centre length times the longest allele or unplaced fragment. Before the centre is known, the
  shortest candidate centre bounds it. abPOA needs about 15 bytes per cell, so the default of 450 M
  cells is about 6.8 GB, which is 21 kb against 21 kb. A region whose alleles hold more than
  MAX_TOTAL_BP bases in all is also too big: choosing the centre holds every allele's k-mer counts.
- `timeout`: the region ran longer than `timeout` seconds.
- `aligner_failed`: abPOA exited with an error or returned the wrong rows.
- `no_alleles`: no walk spans the region.
- `error`: any other failure, with its message; it means a bug.

**Output.** msas.jsonl.gz holds one JSON object per region, in package order:
- `id`, `method` ('medoid_star'), `status`, `seconds`;
- for `ok`: `centre` (the centre allele's id), `rows` {member id: row}, one per allele and fragment,
  all the same length, `columns`, `slot_alignments`, `fragments_exact` and `fragments_aligned`;
- otherwise `message`.
"""
import collections
import gzip
import json
import os
import time

from .star import Abpoa, RealignError, SCORES, TooBig, star_msa

MAX_CELLS = 450_000_000
MAX_TOTAL_BP = 40_000_000     # choosing the centre holds every allele's k-mer counts: up to about 80 bytes per base


def realign_package(pkg, abpoa='abpoa', timeout=900, max_cells=MAX_CELLS, workdir=None, scores=SCORES):
    """The medoid-star alignment of one region package, as the dict described in the module docstring."""
    t0 = time.time()
    deadline = t0 + timeout
    res = {'id': pkg['id'], 'method': 'medoid_star'}
    alleles = [(a['id'], a['seq'], a['weight']) for a in pkg['alleles']]
    fragments = [(f['id'], f['seq'], f['side']) for f in pkg['fragments']]
    if not alleles:
        res.update(status='no_alleles', message='no spanning allele to align', seconds=0.0)
        return res
    try:
        lens = [len(s) for _, s, _ in alleles]
        longest, total = max(lens), sum(lens)
        cand = [n for n, (_, _, w) in zip(lens, alleles) if w > 0]
        shortest = min(cand) if cand else 0
        # the centre is at least the shortest candidate, so this bounds the cells from below
        if max(1, shortest) * max(1, longest) > max_cells or total > MAX_TOTAL_BP:
            raise TooBig('shortest %d bp x longest %d bp > %d cells, or %d bp in all' % (
                shortest, longest, max_cells, total))
        with Abpoa(abpoa, scores, deadline, workdir) as align:
            msa = star_msa(alleles, fragments, align=align, deadline=deadline,
                           max_cells=max_cells)
        rows = msa.rows
        res.update(status='ok', centre=msa.centre, rows=rows, columns=len(next(iter(rows.values()), '')),
                   slot_alignments=msa.slot_alignments, fragments_exact=msa.fragments_exact,
                   fragments_aligned=msa.fragments_aligned)
    except RealignError as e:
        res.update(status=e.status, message=str(e))
    except Exception as e:        # a bug: record it, so that one region cannot stop a contig
        res.update(status='error', message='%s: %s' % (type(e).__name__, e))
    res['seconds'] = round(time.time() - t0, 3)
    return res


def write_msas(path, results):
    """Write realign results (dicts from realign_package) to a msas.jsonl.gz file."""
    with gzip.open(path, 'wt') as f:
        for r in results:
            f.write(json.dumps(r, separators=(',', ':')) + '\n')


def selected(pkg, select):
    """Whether a package's alleles differ in the sense `select` asks: 'any' (any two differ, or there are
    fragments) or 'length' (some two differ in length, or there are fragments)."""
    if pkg['fragments']:
        return True
    seqs = [a['seq'] for a in pkg['alleles']]
    if select == 'length':
        return len({len(x) for x in seqs}) > 1
    return len(set(seqs)) > 1


def run(packages_path, out_path, select='any', shard=(0, 1), jobs=1, **kw):
    """Realign the selected `ok` packages of one shard (index i of n, by package order) in parallel, and
    write them to out_path.

    Results go to out_path + '.partial', one gzip member per region, as each finishes. If that file
    exists, the regions it holds are not redone, so an interrupted run resumes. Once every region is
    done, the results are written to out_path in package order and the partial file is removed."""
    import concurrent.futures
    from .extract import read_packages
    partial = out_path + '.partial'
    # each region's result is held as its compressed line, which is also its member in the partial file
    done, brief = {}, {}
    for r in read_partial(partial) if os.path.exists(partial) else ():
        done[r['id']] = gzip.compress((json.dumps(r, separators=(',', ':')) + '\n').encode())
        brief[r['id']] = (r['status'], r['seconds'])
    todo, order, skipped = [], [], collections.Counter()
    for k, p in enumerate(read_packages(packages_path)):
        if k % shard[1] != shard[0]:
            continue
        if p['status'] != 'ok':
            skipped[p['status']] += 1
        elif not selected(p, select):
            skipped['not_selected'] += 1
        else:
            order.append(p['id'])
            if p['id'] not in done:
                todo.append(p)
    resumed = len(order) - len(todo)
    # rewritten from the records that read back whole, so that nothing follows a cut-off one
    with open(partial + '.tmp', 'wb') as f:
        for rid in order:
            if rid in done:
                f.write(done[rid])
    os.replace(partial + '.tmp', partial)
    with open(partial, 'ab') as f:
        with concurrent.futures.ProcessPoolExecutor(max(1, jobs)) as ex:
            futures = [ex.submit(realign_package, p, **kw) for p in todo]
            for fut in concurrent.futures.as_completed(futures):
                r = fut.result()
                done[r['id']] = gzip.compress((json.dumps(r, separators=(',', ':')) + '\n').encode())
                brief[r['id']] = (r['status'], r['seconds'])
                f.write(done[r['id']])
                f.flush()
    # concatenated gzip members are one valid gzip file
    with open(out_path + '.tmp', 'wb') as out:
        for rid in order:
            out.write(done[rid])
    os.replace(out_path + '.tmp', out_path)
    os.remove(partial)
    status = collections.Counter(brief[rid][0] for rid in order)
    secs = sum(brief[rid][1] for rid in order)
    return {'realigned': dict(status), 'skipped': dict(skipped), 'aligner_seconds': round(secs, 1),
            'resumed': resumed}


def read_partial(path):
    """The results in a partial file: one gzip member per region, the last possibly cut off."""
    out = []
    try:
        with gzip.open(path, 'rt') as f:
            for line in f:
                try:
                    out.append(json.loads(line))
                except ValueError:
                    break
    except (EOFError, OSError):
        pass
    return out
