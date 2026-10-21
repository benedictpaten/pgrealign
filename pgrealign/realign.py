"""Stage 3, `realign`: one multiple alignment per region package, by the star (star.py).

`realign_package` aligns one package from extract.py; `run` aligns a whole packages.jsonl.gz in
parallel and writes msas.jsonl.gz.

**Centre.** The star's centre is the reference allele, the package's `ref`, whenever it has one of
positive weight (`centre` 'ref'), and otherwise the medoid. Every allele is then aligned to the
reference on its own, so the alignment a genotyper reads between the reference and any allele is that
allele's optimal pairwise alignment, its indels left-normalised against the reference, as a variant
caller writes them. Around a medoid that is not the reference, it would instead be the composition of
two alignments to the medoid, which is neither optimal nor normalised against the reference, so that the
same difference from the reference can come out in different places in different alleles. `centre`
'medoid' keeps the medoid.

**Limits.** abPOA's alignment of two sequences works in a matrix of cells, one per pair of their bases.
It reserves 20 bytes of address space per cell, five matrices of 4-byte scores, banded or not, and an
unbanded alignment touches all of it, so `max_cells` cells take about 8 GB. An alignment to the
centre whose matrix would exceed `max_cells` cells is banded instead: abPOA then fills only a band
around its best path, adaptive in width (about 10 cells plus 1% of the sequence's length on either
side). Its memory is then set by the centre's length, not by the band: each base of the centre
touches at least one page of each of the five matrices, about 80 KB per base with 16 KB pages (macOS
on Apple silicon), and more where the two lengths differ by tens of kb (memory.py).

Whether a region is aligned is decided from its package alone, before anything runs, so the same input
gives the same graph on any machine, under any memory budget. A region the star is not built for keeps
its old alignment, and replace leaves it as it is. Its status says why:
- `too_big`, decided from the package:
  - even a banded alignment of the region's longest allele or fragment would exceed `max_cells` cells
    (with the default, sequences longer than about 150 kb);
  - its alleles hold more than MAX_TOTAL_BP bases in all;
  - its alignments would fill more than `max_work` cells of dynamic-programming matrices (`work`);
  - its memory estimate (`peak_bytes`) exceeds `max_region_memory` GB;
  - or abPOA would reserve more than `max_address_space` GB for one of its alignments (`reservation`).
    Linux refuses a reservation larger than the machine's memory and swap, so a region that needs one
    would otherwise be aligned on a large machine and not on a small one (memory.py).
- `timeout`: the region ran longer than `timeout` seconds. There is no time limit unless one is given,
  since a time limit makes a region's outcome depend on the machine it runs on.
- `aligner_failed`: abPOA exited with an error or returned the wrong rows.
- `no_alleles`: no walk spans the region.
- `not_homologous`: an allele of the region shares less than `min_share` (MIN_SHARE) of its 15-mers
  with the centre (star.py), as an unrelated or reversed allele does at an inversion or a
  non-homologous structural variant. Aligned forward, it would only be threaded through the centre with
  chance matches, each a variant of its own. Decided before anything is aligned, for alleles and a
  centre of at least 100 bp (star.HOMOLOGY_BP).
- `error`: any other failure, with its message; it means a bug.

Nor do the machine's failures decide anything. `run` stops before it aligns anything if a region it admits
would need abPOA to reserve more than this machine grants (memory.granted_bytes). A region the machine
fails all the same (MACHINE_FAILURES: abPOA's reservation refused, as abPOA reports with `posix_memalign
fail`; abPOA killed by a signal, as by the OOM killer; the worker out of memory) is not recorded: `run`
finishes the others, then stops with an error, and run again where the region can be aligned, it resumes
with it.

**Memory.** `peak_bytes` estimates, from a package alone, the most memory its region can take, with
pages of ESTIMATE_PAGE bytes whatever the machine's. `run` starts regions in package order, in up to
`jobs` worker processes, while the estimates of those running stay within `max_memory` GB. That budget
decides only how many regions run at once: a region admitted under `max_region_memory` whose estimate
alone exceeds the budget stops `run` before it starts anything, rather than being left out.

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
import math
import os
import time

from . import memory
from .star import Abpoa, RealignError, SCORES, TooBig, star_msa

MAX_CELLS = 450_000_000
MAX_TOTAL_BP = 40_000_000     # choosing the centre holds every allele's k-mer counts: KMER_BYTES per base
CENTRES = ('ref', 'medoid')
ESTIMATE_PAGE = 16 * 1024     # the page size peak_bytes assumes: macOS on Apple silicon's, 4x Linux's
# Of the 1.12 million regions realigned in the HPRC CHM13 graph, MIN_SHARE 0.5 refuses 2,794 (0.6% of their
# span): every region of the chrX and chr15 loci where unrelated or inverted alleles had been threaded
# through the centre, and few others. On chrX that takes back 410,601 non-reference sample-records that
# realignment had added to the VCF, and gives up 2,904 it had removed (docs/method.md). Measured around the medoid;
# around the reference allele it refuses nearly the same regions.
MIN_SHARE = 0.5
# The most a region may take, each decided from its package: cells of its alignments (`work`), GB of estimated
# memory (`peak_bytes`) and GB of abPOA's largest reservation (`reservation`). docs/method.md calibrates them.
MAX_WORK = 200_000_000_000
MAX_REGION_MEMORY = 64
MAX_ADDRESS_SPACE = 128
# Failures of the machine rather than of the region, which `run` never records (see the module docstring)
MACHINE_FAILURES = ('reservation_refused', 'aligner_killed', 'out_of_memory')
# Measured on a 1 Mbp allele: its k-mer counts hold 74 bytes per base, 85 at most while they grow, and
# sketching it takes 49 more per base, comparing it exactly with another 25.
KMER_BYTES, COMPARE_BYTES = 80, 50


def banded_cells(n, m=None):
    """About how many cells abPOA fills aligning a sequence of `m` bases (by default `n`) to one of `n` given
    before it, with its adaptive band: in each of the first's rows, the two lengths' difference and 10 cells
    plus 1% of m either side, and never more than the whole row, as for its memory (memory.abpoa_bytes)."""
    if m is None:
        m = n
    return n * min(m, abs(n - m) + 2 * (10 + m // 100))


def centre_of(pkg, centre='ref'):
    """The id of the allele the star of `pkg` is built around: its reference allele when `centre` is 'ref'
    and it has one of positive weight, or None, for the medoid."""
    if centre not in CENTRES:
        raise ValueError('centre must be one of %s, not %r' % (', '.join(CENTRES), centre))
    if centre == 'ref':
        for a in pkg['alleles']:
            if a['id'] == pkg.get('ref') and a['weight'] > 0:
                return a['id']
    return None


def _placed(frag, side, alleles):
    """Whether star_msa places fragment `frag` in one of `alleles` (upper-cased sequences) rather than aligning
    it: as a prefix of one if `side` is 'left', a suffix if 'right', and anywhere if 'internal'."""
    if side == 'left':
        return any(a.startswith(frag) for a in alleles)
    if side == 'right':
        return any(a.endswith(frag) for a in alleles)
    return any(frag in a for a in alleles)


def _centre_and_others(pkg, centre):
    """(the centre's length, the lengths of the members aligned to it, the centre's id or None): the alleles
    but the centre, and the fragments that occur in no allele. With the medoid, whose length is not known
    before the region runs, the longest allele's stands in for it."""
    cid = centre_of(pkg, centre)
    lens = [(a['id'], len(a['seq'])) for a in pkg['alleles']]
    if cid is None and lens:
        cid_or_longest = max(lens, key=lambda t: t[1])[0]
        c = dict(lens)[cid_or_longest]
    else:
        cid_or_longest = cid
        c = dict(lens)[cid] if cid is not None else 0
    seqs = [a['seq'].upper() for a in pkg['alleles']]
    others = [n for a, n in lens if a != cid_or_longest] + \
        [len(f['seq']) for f in pkg['fragments'] if not _placed(f['seq'].upper(), f['side'], seqs)]
    return c, others, cid


def _slots(c, others):
    """The insertion slots' sequences as the estimates take them, longest first: for each member aligned to
    a centre of c bases, the bases by which it is longer than the centre, all in one slot."""
    return sorted((n - c for n in others if n > c), reverse=True)


def work(pkg, centre='ref', max_cells=MAX_CELLS):
    """About how many cells of dynamic-programming matrices aligning `pkg` fills, the measure of its run time
    (docs/method.md): each allele and fragment aligned to the centre, the product of the two lengths, or
    banded_cells where that exceeds `max_cells` and the alignment is banded; and the insertion slots, each
    member's excess over the centre aligned, banded, to the longest such (_slots). Around a medoid, which is
    not known in advance, the longest allele stands in for it, and the slots, like choosing the medoid, are
    not counted."""
    c, others, cid = _centre_and_others(pkg, centre)
    total = sum(c * n if c * n <= max_cells else banded_cells(c, n) for n in others if c and n)
    if cid is not None:
        slots = _slots(c, others)
        total += sum(banded_cells(slots[0], e) for e in slots[1:])
    return total


def peak_bytes(pkg, max_cells=MAX_CELLS, page=ESTIMATE_PAGE, centre='ref'):
    """An upper bound in practice on the memory, in bytes, that realign_package takes for `pkg`: memory.BASE;
    the k-mer counts the region's alleles take, which the worker may still hold while abPOA runs; and abPOA's
    worst alignment. Around a reference centre, whose alignments are known in advance, that is the worst of
    its alignments to the centre and of the slots' as `work` takes them. Around a medoid, which is not known,
    it is the worst pair of the region's alleles and fragments in either order (memory.any_pair_bytes), or
    the longest sequence banded against one of the same length, the bound taken for the slots. A region
    that realign_package refuses for its size before aligning anything takes memory.BASE. The page size is
    fixed, ESTIMATE_PAGE, so the estimate is the same on every machine."""
    alleles = [len(a['seq']) for a in pkg['alleles']]
    lengths = alleles + [len(f['seq']) for f in pkg['fragments']]
    if not alleles or banded_cells(max(lengths)) > max_cells or sum(alleles) > MAX_TOTAL_BP:
        return memory.BASE
    c, others, cid = _centre_and_others(pkg, centre)
    if cid is None:
        aligned = [c] + others                  # every allele, and the fragments found in none
        longest = max(aligned)
        aligner = max(memory.any_pair_bytes(aligned, max_cells, page), memory.abpoa_bytes(longest, longest, True, page))
    else:
        slots = _slots(c, others)
        aligner = max([memory.abpoa_bytes(c, n, c * n > max_cells, page) for n in others if c and n] +
                      [memory.abpoa_bytes(slots[0], e, True, page) for e in slots[1:]], default=0)
    return memory.BASE + KMER_BYTES * sum(alleles) + COMPARE_BYTES * max(alleles) + aligner


def reservation(pkg, centre='ref'):
    """The most address space, in bytes, that abPOA reserves for one alignment of `pkg`
    (memory.abpoa_reservation): its alignments to the centre, and around a reference centre the slots' as
    `work` takes them, the two longest; around a medoid, the longest allele's length stands in for the
    centre's, and the longest sequence against one of the same length is the bound taken for the slots."""
    c, others, cid = _centre_and_others(pkg, centre)
    pairs = [memory.abpoa_reservation(c, n) for n in others if c and n]
    if cid is None:
        longest = max([c] + others)
        return max(pairs + [memory.abpoa_reservation(longest, longest)])
    slots = _slots(c, others)
    return max(pairs + ([memory.abpoa_reservation(slots[0], slots[1])] if len(slots) > 1 else []), default=0)


def _refusal(pkg, max_cells, centre, max_work, max_region_memory, max_address_space):
    """(the limit `pkg` exceeds, why), or None: see refusal."""
    lens = [len(a['seq']) for a in pkg['alleles']]
    if not lens:
        return None             # no_alleles, in realign_package
    longest = max(lens + [len(f['seq']) for f in pkg['fragments']])
    if banded_cells(longest) > max_cells or sum(lens) > MAX_TOTAL_BP:
        return 'max_cells', 'longest %d bp needs %d cells banded, more than %d, or %d bp in all' % (
            longest, banded_cells(longest), max_cells, sum(lens))
    w = work(pkg, centre, max_cells)
    if max_work and w > max_work:
        return 'max_work', 'aligning it takes about %.3g cells, more than max_work %.3g' % (w, max_work)
    need = peak_bytes(pkg, max_cells, centre=centre)
    if max_region_memory and need > max_region_memory * memory.GB:
        return 'max_region_memory', 'needs about %.2f GB, more than max_region_memory %g GB' % (
            need / memory.GB, max_region_memory)
    space = reservation(pkg, centre)
    if max_address_space and space > max_address_space * memory.GB:
        return 'max_address_space', 'abPOA would reserve %g GB for one alignment, more than max_address_space %g GB' % (
            space / memory.GB, max_address_space)
    return None


def refusal(pkg, max_cells=MAX_CELLS, centre='ref', max_work=MAX_WORK, max_region_memory=MAX_REGION_MEMORY,
            max_address_space=MAX_ADDRESS_SPACE):
    """Why `pkg` is too big to align, from the package alone (see the module docstring), or None."""
    r = _refusal(pkg, max_cells, centre, max_work, max_region_memory, max_address_space)
    return r and r[1]


def realign_package(pkg, abpoa='abpoa', timeout=None, max_cells=MAX_CELLS, workdir=None, scores=SCORES,
                    min_share=MIN_SHARE, centre='ref', max_work=MAX_WORK, max_region_memory=MAX_REGION_MEMORY,
                    max_address_space=MAX_ADDRESS_SPACE):
    """The star alignment of one region package, as the dict described in the module docstring.
    `min_share` 0 or None aligns the region however little its alleles share with the centre; `max_work`,
    `max_region_memory` and `max_address_space` 0 or None set no such limit; `timeout` None sets no time
    limit. A failure of the machine gives one of the statuses MACHINE_FAILURES, which `run` does not record."""
    t0 = time.time()
    deadline = t0 + timeout if timeout else None
    res = {'id': pkg['id'], 'method': 'medoid_star'}
    alleles = [(a['id'], a['seq'], a['weight']) for a in pkg['alleles']]
    fragments = [(f['id'], f['seq'], f['side']) for f in pkg['fragments']]
    if not alleles:
        res.update(status='no_alleles', message='no spanning allele to align', seconds=0.0)
        return res
    try:
        why = refusal(pkg, max_cells, centre, max_work, max_region_memory, max_address_space)
        if why:
            raise TooBig(why)
        with Abpoa(abpoa, scores, deadline, workdir) as align:
            msa = star_msa(alleles, fragments, align=align, centre=centre_of(pkg, centre), deadline=deadline,
                           max_cells=max_cells, min_share=min_share)
        rows = msa.rows
        res.update(status='ok', centre=msa.centre, rows=rows, columns=len(next(iter(rows.values()), '')),
                   slot_alignments=msa.slot_alignments, fragments_exact=msa.fragments_exact,
                   fragments_aligned=msa.fragments_aligned, banded_pairs=msa.banded_pairs)
    except RealignError as e:
        res.update(status=e.status, message=str(e))
    except MemoryError:
        res.update(status='out_of_memory', message='the worker ran out of memory')
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


def _shard(packages_path, select, shard, skipped):
    """The `ok` packages of shard i of n (`shard`, by package order) that `select` takes, counting the others
    in `skipped` by their status, or as not_selected."""
    from .extract import read_packages
    for k, p in enumerate(read_packages(packages_path)):
        if k % shard[1] != shard[0]:
            continue
        if p['status'] != 'ok':
            skipped[p['status']] += 1
        elif not selected(p, select):
            skipped['not_selected'] += 1
        else:
            yield p


def plan(packages_path, select='any', shard=(0, 1), max_cells=MAX_CELLS, centre='ref', max_work=MAX_WORK,
         max_region_memory=MAX_REGION_MEMORY, max_address_space=MAX_ADDRESS_SPACE):
    """What `run` would do with the same packages and limits, decided from the packages alone, without
    aligning anything: the number of regions it would align (`admitted`), and leave too_big, by the limit
    each exceeds (`too_big`); their `work` in all; `need_gb`, the largest admitted region's memory estimate
    in GB, the least max_memory with which `run` aligns them, and that region's id (`largest`); and
    `reservation_gb`, the largest reservation abPOA would make for them, in GB, against `granted_gb`, the
    largest this machine grants (memory.granted_bytes; None for no limit), for which `run` must be given a
    machine that grants it."""
    skipped, too_big = collections.Counter(), collections.Counter()
    admitted, cells, need, largest, space = 0, 0, 0, None, 0
    for p in _shard(packages_path, select, shard, skipped):
        r = _refusal(p, max_cells, centre, max_work, max_region_memory, max_address_space)
        if r:
            too_big[r[0]] += 1
            continue
        admitted += 1
        cells += work(p, centre, max_cells)
        b = peak_bytes(p, max_cells, centre=centre)
        if b > need:
            need, largest = b, p['id']
        space = max(space, reservation(p, centre))
    granted = memory.granted_bytes()
    return {'admitted': admitted, 'too_big': dict(too_big), 'skipped': dict(skipped), 'work': cells,
            'need_gb': math.ceil(need / memory.GB * 100) / 100, 'largest': largest,
            'reservation_gb': space / memory.GB, 'granted_gb': None if granted is None else round(granted / memory.GB, 2)}


def run(packages_path, out_path, select='any', shard=(0, 1), jobs=1, max_memory=None, **kw):
    """Realign the selected `ok` packages of one shard (index i of n, by package order) in up to `jobs`
    worker processes, and write them to out_path. `kw` goes to realign_package. A region that `refusal`
    finds too big is recorded as too_big without being aligned. The others start in package order while
    the sum of their peak_bytes stays within `max_memory` GB (None: memory.default_gb()); if one of them
    alone needs more, or would need abPOA to reserve more than this machine grants, nothing starts and a
    RuntimeError says so, so neither the budget nor the machine changes which regions are aligned.

    Results go to out_path + '.partial', one gzip member per region, as each finishes. If that file
    exists, the regions it holds are not redone, so an interrupted run resumes; the limits apply afresh,
    so a region the partial file holds is too_big if it is now and redone if it was then but is not now.
    A region that fails for the machine's sake (MACHINE_FAILURES) is not recorded: the others finish, then
    a RuntimeError names it, and a rerun redoes it. Once every region is done, the results are written to
    out_path in package order and the partial file is removed."""
    limit = memory.budget(max_memory)
    partial = out_path + '.partial'
    max_cells, centre = kw.get('max_cells', MAX_CELLS), kw.get('centre', 'ref')
    limits = (kw.get('max_work', MAX_WORK), kw.get('max_region_memory', MAX_REGION_MEMORY),
              kw.get('max_address_space', MAX_ADDRESS_SPACE))
    # each region's result is held as its compressed line, which is also its member in the partial file
    held = {r['id']: r for r in (read_partial(partial) if os.path.exists(partial) else ())}
    done, brief, refused = {}, {}, {}
    todo, order, skipped = [], [], collections.Counter()
    for p in _shard(packages_path, select, shard, skipped):
        order.append(p['id'])
        why = refusal(p, max_cells, centre, *limits)
        r = held.get(p['id'])
        if why is not None:
            refused[p['id']] = why
        elif r is not None and r['status'] != 'too_big':
            done[p['id']] = gzip.compress((json.dumps(r, separators=(',', ':')) + '\n').encode())
            brief[p['id']] = (r['status'], r['seconds'])
        else:
            todo.append(p)
    resumed = len(done)
    # rewritten from the records kept, so that nothing follows a cut-off one
    with open(partial + '.tmp', 'wb') as f:
        for rid in order:
            if rid in done:
                f.write(done[rid])
    os.replace(partial + '.tmp', partial)
    need = [peak_bytes(p, max_cells, centre=centre) for p in todo]
    over = [k for k in range(len(todo)) if need[k] > limit]
    if over:
        k = max(over, key=lambda k: need[k])
        raise RuntimeError('%d regions need more memory than the budget of %.2f GB (max_memory); the largest, %s, '
                           'needs about %.2f GB. Give realign that much, or lower max_region_memory below it, which '
                           'leaves such regions too_big.' % (len(over), limit / memory.GB, todo[k]['id'],
                                                             need[k] / memory.GB))
    granted = memory.granted_bytes()
    space = [reservation(p, centre) for p in todo] if granted is not None else []
    beyond = [k for k in range(len(space)) if space[k] > granted]
    if beyond:
        k = max(beyond, key=lambda k: space[k])
        raise RuntimeError('%d regions would need abPOA to reserve more address space than this machine grants, %.1f '
                           'GB (its memory and swap); the largest, %s, %g GB. Run realign on a machine with that much '
                           'memory and swap, or lower max_address_space below it, which leaves such regions too_big.'
                           % (len(beyond), granted / memory.GB, todo[k]['id'], space[k] / memory.GB))
    failed = []               # regions this machine failed to align: not recorded, so a rerun redoes them
    with open(partial, 'ab') as f:
        def take(r):
            if r['status'] in MACHINE_FAILURES:
                failed.append(r)
                return
            done[r['id']] = gzip.compress((json.dumps(r, separators=(',', ':')) + '\n').encode())
            brief[r['id']] = (r['status'], r['seconds'])
            f.write(done[r['id']])
            f.flush()
        for rid in order:
            if rid in refused:
                take({'id': rid, 'method': 'medoid_star', 'status': 'too_big', 'message': refused[rid], 'seconds': 0.0})
        peak = memory.run_within(functools.partial(realign_package, **kw), todo, need, jobs, limit,
                                 lambda i, r: take(r))
    if failed:
        raise RuntimeError('this machine failed to align %d regions (%s), among them %s: %s. Run realign again where '
                           'it can (it resumes from %s), or lower max_region_memory or max_address_space, which '
                           'leaves such regions too_big.' % (
                               len(failed), ', '.join(sorted({r['status'] for r in failed})), failed[0]['id'],
                               failed[0]['message'][-160:], partial))
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
