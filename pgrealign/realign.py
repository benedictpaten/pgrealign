"""Stage 3, `realign`: one multiple alignment per region package, by the medoid star (star.py).

`realign_package` aligns one package from extract.py; `run` aligns a whole packages.jsonl.gz in
parallel and writes msas.jsonl.gz.

**Limits.** abPOA's alignment of two sequences works in a matrix of cells, one per pair of their bases.
It reserves 20 bytes of address space per cell, five matrices of 4-byte scores, banded or not, and an
unbanded alignment touches all of it, so `max_cells` cells take about 8 GB. An alignment to the
centre whose matrix would exceed `max_cells` cells is banded instead: abPOA then fills only a band
around its best path, adaptive in width (about 10 cells plus 1% of the sequence's length on either
side). Its memory is then set by the centre's length, not by the band: each base of the centre
touches at least one page of each of the five matrices, about 80 KB per base with 16 KB pages (macOS
on Apple silicon), and more where the two lengths differ by tens of kb (memory.py). A region the star
cannot be built for keeps its old alignment, and replace leaves it as it is. Its status says why:
- `too_big`: even a banded alignment of the region's longest allele or fragment would exceed
  `max_cells` cells (with the default, sequences longer than about 150 kb), or its alleles hold more
  than MAX_TOTAL_BP bases in all, or its memory estimate alone exceeds `max_memory` (below), which
  `run` records itself, without aligning the region.
- `timeout`: the region ran longer than `timeout` seconds.
- `aligner_failed`: abPOA exited with an error or returned the wrong rows, as when the system refuses
  its reservation of address space.
- `no_alleles`: no walk spans the region.
- `not_homologous`: an allele of the region shares less than `min_share` (MIN_SHARE) of its 15-mers
  with the centre (star.py), as an unrelated or reversed allele does at an inversion or a
  non-homologous structural variant. Aligned forward, it would only be threaded through the centre with
  chance matches, each a variant of its own. Decided before anything is aligned, for alleles and a
  centre of at least 100 bp (star.HOMOLOGY_BP).
- `error`: any other failure, with its message; it means a bug.

**Memory.** `peak_bytes` estimates, from a package alone, the most memory its region can take. `run`
starts regions in package order, in up to `jobs` worker processes, while the estimates of those
running stay within `max_memory` GB.

**Output.** msas.jsonl.gz holds one JSON object per region, in package order:
- `id`, `method` ('medoid_star'), `status`, `seconds`;
- for `ok`: `centre` (the centre allele's id), `rows` {member id: row}, one per allele and fragment,
  all the same length, `columns`, `slot_alignments`, `fragments_exact`, `fragments_aligned` and
  `banded_pairs`;
- otherwise `message`.
"""
import collections
import functools
import gzip
import json
import os
import time

from . import memory
from .star import Abpoa, RealignError, SCORES, TooBig, star_msa

MAX_CELLS = 450_000_000
MAX_TOTAL_BP = 40_000_000     # choosing the centre holds every allele's k-mer counts: KMER_BYTES per base
# Of the 1.12 million regions realigned in the HPRC CHM13 graph, MIN_SHARE 0.5 refuses 2,794 (0.6% of their
# span): every region of the chrX and chr15 loci where unrelated or inverted alleles had been threaded
# through the centre, and few others. On chrX that takes back 410,601 non-reference sample-records that
# realignment had added to the VCF, and gives up 2,904 it had removed (docs/method.md).
MIN_SHARE = 0.5
# Measured on a 1 Mbp allele: its k-mer counts hold 74 bytes per base, 85 at most while they grow, and
# sketching it takes 49 more per base, comparing it exactly with another 25.
KMER_BYTES, COMPARE_BYTES = 80, 50


def banded_cells(length):
    """About how many cells abPOA fills aligning a sequence of `length` bases to one of the same length
    with its adaptive band, which reaches 10 cells plus 1% of the length either side of the path, and
    never more than the whole matrix."""
    return min(length * length, length * 2 * (10 + length // 100))


def peak_bytes(pkg, max_cells=MAX_CELLS, page=memory.PAGE):
    """An upper bound in practice on the memory, in bytes, that realign_package takes for `pkg`: memory.BASE;
    the k-mer counts that choose the centre, which the worker may still hold while abPOA runs; and the
    larger of abPOA's worst pair of the region's alleles and fragments in either order, since the
    centre is not known before the region runs (memory.any_pair_bytes), and its longest sequence
    banded against one of the same length, the bound taken for the slot alignments. A region that
    realign_package refuses before aligning anything takes memory.BASE."""
    alleles = [len(a['seq']) for a in pkg['alleles']]
    lengths = alleles + [len(f['seq']) for f in pkg['fragments']]
    if not alleles or banded_cells(max(lengths)) > max_cells or sum(alleles) > MAX_TOTAL_BP:
        return memory.BASE
    longest = max(lengths)
    aligner = max(memory.any_pair_bytes(lengths, max_cells, page), memory.abpoa_bytes(longest, longest, True, page))
    return memory.BASE + KMER_BYTES * sum(alleles) + COMPARE_BYTES * max(alleles) + aligner


def realign_package(pkg, abpoa='abpoa', timeout=900, max_cells=MAX_CELLS, workdir=None, scores=SCORES,
                    min_share=MIN_SHARE):
    """The medoid-star alignment of one region package, as the dict described in the module docstring.
    `min_share` 0 or None aligns the region however little its alleles share with the centre."""
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
        longest = max(lens + [len(s) for _, s, _ in fragments])
        if banded_cells(longest) > max_cells or sum(lens) > MAX_TOTAL_BP:
            raise TooBig('longest %d bp needs %d cells banded, more than %d, or %d bp in all' % (
                longest, banded_cells(longest), max_cells, sum(lens)))
        with Abpoa(abpoa, scores, deadline, workdir) as align:
            msa = star_msa(alleles, fragments, align=align, deadline=deadline, max_cells=max_cells,
                           min_share=min_share)
        rows = msa.rows
        res.update(status='ok', centre=msa.centre, rows=rows, columns=len(next(iter(rows.values()), '')),
                   slot_alignments=msa.slot_alignments, fragments_exact=msa.fragments_exact,
                   fragments_aligned=msa.fragments_aligned, banded_pairs=msa.banded_pairs)
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


def run(packages_path, out_path, select='any', shard=(0, 1), jobs=1, max_memory=None, **kw):
    """Realign the selected `ok` packages of one shard (index i of n, by package order) in up to `jobs`
    worker processes, and write them to out_path. Regions start in package order while the sum of
    their peak_bytes stays within `max_memory` GB (None: memory.default_gb()); a region whose
    estimate alone exceeds it is recorded as too_big without being aligned.

    Results go to out_path + '.partial', one gzip member per region, as each finishes. If that file
    exists, the regions it holds are not redone, so an interrupted run resumes. Once every region is
    done, the results are written to out_path in package order and the partial file is removed."""
    from .extract import read_packages
    limit = memory.budget(max_memory)
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
    need = [peak_bytes(p, kw.get('max_cells', MAX_CELLS)) for p in todo]
    with open(partial, 'ab') as f:
        def take(r):
            done[r['id']] = gzip.compress((json.dumps(r, separators=(',', ':')) + '\n').encode())
            brief[r['id']] = (r['status'], r['seconds'])
            f.write(done[r['id']])
            f.flush()
        for p, b in zip(todo, need):
            if b > limit:
                take({'id': p['id'], 'method': 'medoid_star', 'status': 'too_big',
                      'message': 'needs about %.2f GB, more than the memory budget of %.2f GB' % (
                          b / memory.GB, limit / memory.GB), 'seconds': 0.0})
        fit = [k for k in range(len(todo)) if need[k] <= limit]
        peak = memory.run_within(functools.partial(realign_package, **kw), [todo[k] for k in fit],
                                 [need[k] for k in fit], jobs, limit, lambda i, r: take(r))
    # concatenated gzip members are one valid gzip file
    with open(out_path + '.tmp', 'wb') as out:
        for rid in order:
            out.write(done[rid])
    os.replace(out_path + '.tmp', out_path)
    os.remove(partial)
    status = collections.Counter(brief[rid][0] for rid in order)
    secs = sum(brief[rid][1] for rid in order)
    return {'realigned': dict(status), 'skipped': dict(skipped), 'aligner_seconds': round(secs, 1),
            'resumed': resumed, 'peak_estimate_gb': round(peak / memory.GB, 2)}


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
