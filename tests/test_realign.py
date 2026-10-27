import contextlib
import gzip
import io
import json
import os
import random
import shutil
import tempfile
import unittest
from unittest import mock

from pgrealign import cli, memory, realign, star


def pkg(pid, alleles, fragments=(), ref=None):
    return {'id': pid, 'status': 'ok', 'ref': ref,
            'alleles': [{'id': a, 'weight': w, 'seq': s} for a, s, w in alleles],
            'fragments': [{'id': f, 'side': side, 'seq': s} for f, s, side in fragments]}


@unittest.skipUnless(shutil.which('abpoa'), 'abpoa not on PATH')
class TestRealignPackage(unittest.TestCase):
    def test_ok_fields(self):
        r = realign.realign_package(pkg('R', [('a1', 'ACGTACGT', 1), ('a2', 'ACGTTCGT', 1)],
                                        [('f1', 'GGGG', 'internal')]))
        self.assertEqual((r['status'], r['fragments_aligned'], r['fragments_exact']), ('ok', 1, 0))
        self.assertEqual(r['rows']['f1'].replace('-', ''), 'GGGG')
        self.assertEqual(r['columns'], len(r['rows']['a1']))
        self.assertEqual(r['method'], 'medoid_star')

    def test_reference_centre(self):
        # the heavier a2 is the medoid. The reference allele, a1, is the centre unless the medoid is asked for,
        # or too much of the haplotypes' bases lie beyond its length: here 42%, as a1 lacks two or three copies
        # of the unit that the other nine haplotypes carry
        unit = 'ACGTTGCAAGT'
        alleles = [('a1', 'GATTACA' + unit * 2 + 'CCGG', 1), ('a2', 'GATTACA' + unit * 4 + 'CCGG', 5),
                   ('a3', 'GATTACA' + unit * 5 + 'CCGG', 4)]
        self.assertAlmostEqual(realign.excess(pkg('R', alleles, ref='a1'), 'a1'), (5 * 22 + 4 * 33) / 572)
        r = realign.realign_package(pkg('R', alleles, ref='a1'), max_excess=0.5)
        self.assertEqual((r['status'], r['centre']), ('ok', 'a1'))
        self.assertEqual(realign.realign_package(pkg('R', alleles, ref='a1'), max_excess=None)['centre'], 'a1')
        self.assertEqual(realign.realign_package(pkg('R', alleles, ref='a1'))['centre'], 'a2')       # 0.42 > 0.1
        self.assertEqual(realign.realign_package(pkg('R', alleles, ref='a1'), centre='medoid')['centre'], 'a2')
        self.assertEqual(realign.realign_package(pkg('R', alleles))['centre'], 'a2')     # no reference allele
        # a reference allele longer than the others is the centre: none of their bases lie beyond it
        longer = [('a1', 'GATTACA' + unit * 6 + 'CCGG', 1)] + alleles[1:]
        self.assertEqual(realign.realign_package(pkg('R', longer, ref='a1'))['centre'], 'a1')

    def test_statuses(self):
        two = [('a1', 'ACGTACGT', 1), ('a2', 'ACGTTCGT', 1)]
        self.assertEqual(realign.realign_package(pkg('R', two), max_cells=10)['status'], 'too_big')
        # 400 bp: 160,000 cells unbanded. One substitution, so that they share their 15-mers (min_share)
        long_two = [('a1', 'ACGTACGTAC' * 40, 2), ('a2', 'ACGTACGTAC' * 20 + 'ACGTTCGTAC' + 'ACGTACGTAC' * 19, 1)]
        self.assertEqual(realign.realign_package(pkg('R', long_two), max_cells=5000)['status'], 'too_big')
        r = realign.realign_package(pkg('R', long_two), max_cells=50000)   # banded: about 11,200 cells
        self.assertEqual((r['status'], r['banded_pairs']), ('ok', 1))
        self.assertEqual(realign.realign_package(pkg('R', long_two))['banded_pairs'], 0)
        self.assertEqual(realign.realign_package(pkg('R', []))['status'], 'no_alleles')
        self.assertEqual(realign.realign_package(pkg('R', two), timeout=1e-9)['status'], 'timeout')
        r = realign.realign_package(pkg('R', long_two), max_work=1000)
        self.assertEqual(r['status'], 'too_big')
        self.assertIn('more than max_work', r['message'])
        r = realign.realign_package(pkg('R', two), abpoa='false')
        self.assertEqual(r['status'], 'aligner_failed', r)


class TestEstimates(unittest.TestCase):
    """work, peak_bytes and the centre, from packages alone."""

    def test_centre_of(self):
        alleles = [('a1', 'ACGT', 1), ('a2', 'ACGA', 3)]
        self.assertEqual(realign.centre_of(pkg('R', alleles, ref='a1')), 'a1')
        self.assertIsNone(realign.centre_of(pkg('R', alleles, ref='a1'), 'medoid'))
        self.assertIsNone(realign.centre_of(pkg('R', alleles)))                     # no reference allele
        self.assertIsNone(realign.centre_of(pkg('R', [('a1', 'ACGT', 0), ('a2', 'ACGA', 3)], ref='a1')))
        with self.assertRaises(ValueError):
            realign.centre_of(pkg('R', alleles), 'longest')

    def test_placed_fragments_are_not_aligned(self):
        # a fragment that occurs in an allele takes that allele's row and is never aligned, so it adds no work
        alleles = [('r', 'ACGT' * 2500, 1), ('x', 'ACGT' * 2400 + 'TTTT' * 100, 1)]
        bare = realign.work(pkg('R', alleles, ref='r'))
        self.assertEqual(realign.work(pkg('R', alleles, [('f', 'ACGT' * 1000, 'left')], ref='r')), bare)
        self.assertEqual(realign.work(pkg('R', alleles, [('f', 'GGCC' * 1000, 'left')], ref='r')), bare + 10000 * 4000)
        self.assertEqual(realign.work(pkg('R', alleles, [('f', 'TTTT' * 50, 'right')], ref='r')), bare)

    def test_work_and_memory_around_the_reference(self):
        # a 10 kb reference, and 100 kb and 30 kb alleles: around the reference, x's alignment is banded, its
        # band spanning the 90 kb the lengths differ by, y's is not, and the 90 kb and 20 kb by which x and y
        # are longer are aligned in a slot; around a medoid, not known in advance, the longest stands in. With
        # 83% of the haplotypes' bases beyond the reference's length, the medoid is the centre by default, so
        # the estimates around the reference are asked for with max_excess None
        p = pkg('R', [('r', 'A' * 10000, 1), ('x', 'C' * 100000, 2), ('y', 'G' * 30000, 1)], ref='r')
        self.assertAlmostEqual(realign.excess(p, 'r'), (2 * 90000 + 20000) / 240000)
        self.assertIsNone(realign.centre_of(p))
        self.assertEqual(realign.work(p), realign.work(p, 'medoid'))
        ref = dict(max_excess=None)
        self.assertEqual(realign.centre_of(p, **ref), 'r')
        self.assertEqual(realign.banded_cells(10000, 100000), 10000 * (90000 + 2 * (10 + 1000)))
        self.assertEqual(realign.banded_cells(90000, 20000), 90000 * 20000)       # the band would span the rows
        self.assertEqual(realign.banded_cells(100000), 100000 * 2 * (10 + 1000))  # equal lengths: the band alone
        self.assertEqual(realign.work(p, **ref), realign.long_pair_cells(10000, 100000) + 10000 * 30000 +
                         realign.banded_cells(90000, 20000))
        self.assertEqual(realign.work(p, 'medoid'),
                         realign.long_pair_cells(100000, 10000) + realign.long_pair_cells(100000, 30000))
        # x's: the band, abPOA's own since the lengths differ by more than a wider one allows within MAX_CELLS
        self.assertEqual(star.band_width(10000, 100000, realign.MAX_CELLS), star.default_band(100000))
        self.assertEqual(realign.long_pair_cells(10000, 100000), realign.banded_cells(10000, 100000))
        # where the lengths allow a wider band, about MAX_CELLS
        self.assertLessEqual(realign.long_pair_cells(30000, 31000), realign.MAX_CELLS)
        self.assertGreater(realign.long_pair_cells(30000, 31000), 0.99 * realign.MAX_CELLS)
        # its memory: the worst of x's and y's alignments to the reference and the slot's, the 20 kb that y
        # adds aligned, banded, to the 90 kb that x adds, whose band spans the 70 kb they differ by
        page = realign.ESTIMATE_PAGE
        aligner = max(memory.pair_bytes(10000, 100000, realign.MAX_CELLS, page),
                      memory.abpoa_bytes(10000, 30000, False, page), memory.abpoa_bytes(90000, 20000, True, page))
        self.assertEqual(aligner, memory.abpoa_bytes(90000, 20000, True, page))
        self.assertEqual(realign.peak_bytes(p, **ref), memory.BASE + realign.KMER_BYTES * 140000 +
                         realign.COMPARE_BYTES * 100000 + aligner)
        self.assertLess(realign.peak_bytes(p, **ref), realign.peak_bytes(p, centre='medoid'))
        # the estimates assume pages of ESTIMATE_PAGE bytes, whatever the machine's, and the page matters
        self.assertEqual(realign.ESTIMATE_PAGE, 16384)
        near = pkg('E', [('r', 'A' * 60000, 1), ('x', 'C' * 60100, 1)], ref='r')       # banded, page-bound
        self.assertEqual(realign.centre_of(near), 'r')
        self.assertGreater(realign.peak_bytes(near), realign.peak_bytes(near, page=4096))
        self.assertIn('max_work', realign.refusal(p, max_work=1e8, **ref))
        self.assertIn('max_region_memory', realign.refusal(p, max_region_memory=0.5, **ref))
        # abPOA reserves its matrices whole, rounded up to a power of two: 2^35 bytes to align the 100 kb allele
        # to the reference, and 2^36 to align the 20 kb that y adds to the 90 kb that x adds, in a slot
        self.assertEqual(memory.abpoa_reservation(10000, 100000), 2 ** 35)
        self.assertEqual(realign.reservation(p, **ref), 2 ** 36)
        self.assertIsNone(realign.refusal(p, **ref))
        self.assertIn('abPOA would reserve 64 GB for one alignment, more than max_address_space 32 GB',
                      realign.refusal(p, max_address_space=32, **ref))
        self.assertIsNone(realign.refusal(p, max_address_space=0, **ref))
        # around a medoid, the longest allele against one of the same length bounds the slots: 2^38 bytes
        self.assertEqual(realign.reservation(p, 'medoid'), 2 ** 38)
        # with one long allele only, no slot holds two long insertions
        q = pkg('Q', [('r', 'A' * 1000, 1), ('x', 'C' * 100000, 2), ('y', 'G' * 1200, 1)], ref='r')
        self.assertEqual(realign.reservation(q, **ref), memory.abpoa_reservation(1000, 100000))
        self.assertEqual(realign.peak_bytes(q, **ref), memory.BASE + realign.KMER_BYTES * 102200 +
                         realign.COMPARE_BYTES * 100000 + 20 * 1000 * 100000)    # x to r, unbanded


class TestHomology(unittest.TestCase):
    """A region with an allele unrelated to the centre is refused before anything is aligned, so these
    need no abPOA (abpoa='false' would fail any alignment)."""

    def setUp(self):
        rng = random.Random(9)
        rand = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        centre = rand(300)
        self.unrelated = pkg('U', [('a1', centre, 3), ('a2', centre[:150] + 'T' + centre[150:], 1),
                                   ('a3', rand(300), 1)])

    def test_realign_package(self):
        r = realign.realign_package(self.unrelated, abpoa='false')
        self.assertEqual(r['status'], 'not_homologous', r)
        self.assertIn('the least, a3 (weight 1), shares 0.', r['message'])
        self.assertNotIn('rows', r)
        self.assertEqual(realign.realign_package(self.unrelated, abpoa='false', min_share=0)['status'],
                         'aligner_failed')        # aligned, so abpoa='false' fails it

    def test_run(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        pk, out = os.path.join(d, 'packages.jsonl.gz'), os.path.join(d, 'msas.jsonl.gz')
        with gzip.open(pk, 'wt') as f:
            f.write(json.dumps(self.unrelated) + '\n')
        s = realign.run(pk, out, abpoa='false')
        self.assertEqual(s['realigned'], {'not_homologous': 1})
        self.assertEqual([r['status'] for r in map(json.loads, gzip.open(out, 'rt'))], ['not_homologous'])
        self.assertEqual(realign.run(pk, out, abpoa='false', min_share=0)['realigned'], {'aligner_failed': 1})


@unittest.skipUnless(shutil.which('abpoa'), 'abpoa not on PATH')
class TestRun(unittest.TestCase):
    def test_run_and_resume(self):
        d = tempfile.mkdtemp()
        try:
            pk = os.path.join(d, 'packages.jsonl.gz')
            pkgs = [pkg('R%d' % i, [('a1', 'ACGTACGT' * (i + 1), 1), ('a2', 'ACGTTCGT' * (i + 1), 1)])
                    for i in range(4)] + [pkg('S', [('a1', 'ACGT', 2)])]
            with gzip.open(pk, 'wt') as f:
                for p in pkgs:
                    f.write(json.dumps(p) + '\n')
            out = os.path.join(d, 'msas.jsonl.gz')
            s = realign.run(pk, out, jobs=2)
            self.assertEqual(s['realigned'], {'ok': 4})
            self.assertEqual(s['skipped'], {'not_selected': 1})
            first = [json.loads(l) for l in gzip.open(out, 'rt')]
            self.assertEqual([r['id'] for r in first], ['R0', 'R1', 'R2', 'R3'])
            # an interrupted run: two regions are in the partial file, and a third was cut off mid-write
            with open(out + '.partial', 'wb') as f:
                for r in first[:2]:
                    f.write(gzip.compress((json.dumps(r) + '\n').encode()))
                f.write(gzip.compress((json.dumps(first[2]) + '\n').encode())[:30])
            os.remove(out)
            s = realign.run(pk, out, jobs=1)
            self.assertEqual((s['resumed'], s['realigned']), (2, {'ok': 4}))
            # a second interruption after the first: the partial file was rewritten, so it reads whole
            with open(out + '.partial', 'wb') as f:
                f.write(gzip.compress((json.dumps(first[0]) + '\n').encode())[:30])
            os.remove(out)
            s = realign.run(pk, out, jobs=1)
            self.assertEqual((s['resumed'], s['realigned']), (0, {'ok': 4}))
            again = [json.loads(l) for l in gzip.open(out, 'rt')]
            self.assertEqual([r['rows'] for r in again], [r['rows'] for r in first])
            self.assertFalse(os.path.exists(out + '.partial'))
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_machine_failures(self):
        # a reservation the machine refuses, or an abPOA it kills, is a failure of the run, not of the region: the
        # run aligns the others, stops without writing its output, and run where abPOA succeeds, it resumes with
        # the region that failed
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        refusing, killed = os.path.join(d, 'refusing'), os.path.join(d, 'killed')
        with open(refusing, 'w') as f:
            f.write('#!/bin/sh\necho "[SIMDMalloc] posix_memalign fail! Size: 137438953472" >&2\nexit 1\n')
        with open(killed, 'w') as f:
            f.write('#!/bin/sh\nkill -9 $$\n')
        for exe in (refusing, killed):
            os.chmod(exe, 0o755)
        unit = 'ACGTTGCAAGTCAT'
        snp = pkg('SNP', [('a1', unit * 3, 1), ('a2', unit + 'G' + unit[1:] + unit, 1)], ref='a1')
        pk, out = os.path.join(d, 'packages.jsonl.gz'), os.path.join(d, 'msas.jsonl.gz')
        with gzip.open(pk, 'wt') as f:
            # the failing region first, then two that need no abPOA, an allele aligning to the centre as a deletion
            for p in (snp, pkg('DEL', [('a1', unit * 3, 1), ('a2', '', 1)], ref='a1'),
                      pkg('DEL2', [('a1', unit * 4, 1), ('a2', '', 1)], ref='a1')):
                f.write(json.dumps(p) + '\n')
        for exe, status in ((refusing, 'reservation_refused'), (killed, 'aligner_killed')):
            r = realign.realign_package(snp, abpoa=exe)
            self.assertEqual(r['status'], status, r)
            with self.assertRaisesRegex(RuntimeError, r'failed to align 1 regions \(%s\), among them SNP' % status):
                realign.run(pk, out, abpoa=exe, jobs=1)
            self.assertFalse(os.path.exists(out))
            self.assertEqual([(r['id'], r['status']) for r in realign.read_partial(out + '.partial')],
                             [('DEL', 'ok'), ('DEL2', 'ok')])
        s = realign.run(pk, out)
        self.assertEqual((s['resumed'], s['realigned']), (2, {'ok': 3}))
        self.assertEqual([r['id'] for r in map(json.loads, gzip.open(out, 'rt'))], ['SNP', 'DEL', 'DEL2'])

    def test_out_of_memory(self):
        with mock.patch.object(realign, 'star_msa', side_effect=MemoryError):
            r = realign.realign_package(pkg('R', [('a1', 'ACGT', 1), ('a2', 'ACGA', 1)], ref='a1'), abpoa='false')
        self.assertEqual(r['status'], 'out_of_memory')
        self.assertIn('out_of_memory', realign.MACHINE_FAILURES)

    def test_reservation_beyond_the_machine(self):
        # a region whose reservation this machine would refuse stops the run before anything is aligned; plan
        # gives the reservation, and what the machine grants
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        pk, out = os.path.join(d, 'packages.jsonl.gz'), os.path.join(d, 'msas.jsonl.gz')
        unit = 'ACGTTGCAAGTCAT'
        big = pkg('BIG', [('a1', unit * 200, 1), ('a2', unit * 100 + 'G' + unit * 99, 1)], ref='a1')
        with gzip.open(pk, 'wt') as f:
            f.write(json.dumps(big) + '\n')
        space = realign.reservation(big)
        self.assertEqual(space, memory.abpoa_reservation(2800, 2800))
        with mock.patch.object(memory, 'granted_bytes', return_value=space // 2):
            with self.assertRaisesRegex(RuntimeError, 'more address space than this machine grants.*the largest, BIG'):
                realign.run(pk, out, abpoa='false')
            pl = realign.plan(pk)
        self.assertEqual((pl['reservation_gb'], pl['granted_gb']), (space / memory.GB, round(space / 2 / memory.GB, 2)))
        self.assertFalse(os.path.exists(out))

    def test_resume_applies_the_limits_afresh(self):
        # a record in the partial file is kept only if the limits given now admit its region, and a region left
        # too_big then is aligned now if they do, so a resumed run gives what a fresh one would
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        pk = os.path.join(d, 'packages.jsonl.gz')
        unit = 'ACGTTGCAAGTCAT'
        small = pkg('SMALL', [('a1', unit * 3, 1), ('a2', '', 1)], ref='a1')        # no abPOA needed
        large = pkg('LARGE', [('a1', unit * 300, 1), ('a2', '', 1)], ref='a1')
        with gzip.open(pk, 'wt') as f:
            for p in (small, large):
                f.write(json.dumps(p) + '\n')
        low = (realign.peak_bytes(small) + realign.peak_bytes(large)) / 2 / memory.GB     # admits SMALL only
        pl = realign.plan(pk, max_region_memory=low)
        self.assertEqual((pl['admitted'], pl['too_big'], pl['largest']), (1, {'max_region_memory': 1}, 'SMALL'))

        def records(path):
            return [{k: v for k, v in r.items() if k != 'seconds'} for r in map(json.loads, gzip.open(path, 'rt'))]

        def resume(name, held, **kw):
            out = os.path.join(d, name)
            with open(out + '.partial', 'wb') as f:
                for r in held:
                    f.write(gzip.compress((json.dumps(r) + '\n').encode()))
            return realign.run(pk, out, abpoa='false', **kw), records(out)

        fresh = os.path.join(d, 'fresh.jsonl.gz')
        realign.run(pk, fresh, abpoa='false')
        both_ok = list(map(json.loads, gzip.open(fresh, 'rt')))
        self.assertEqual([r['status'] for r in both_ok], ['ok', 'ok'])
        lowered = os.path.join(d, 'lowered.jsonl.gz')
        realign.run(pk, lowered, abpoa='false', max_region_memory=low)
        s, got = resume('a.jsonl.gz', both_ok, max_region_memory=low)
        self.assertEqual((s['resumed'], s['realigned']), (1, {'ok': 1, 'too_big': 1}))
        self.assertEqual(got, records(lowered))
        s, got = resume('b.jsonl.gz', list(map(json.loads, gzip.open(lowered, 'rt'))))
        self.assertEqual((s['resumed'], s['realigned']), (1, {'ok': 2}))
        self.assertEqual(got, records(fresh))

    def test_plan_on_the_command_line(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        pk = os.path.join(d, 'packages.jsonl.gz')
        with gzip.open(pk, 'wt') as f:
            f.write(json.dumps(pkg('R', [('a1', 'ACGT' * 30, 1), ('a2', 'ACGT' * 29, 1)], ref='a1')) + '\n')
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            cli.main(['realign', '--plan', '--packages', pk, '--max-work', '1000'])
        self.assertEqual(json.loads(buf.getvalue())['too_big'], {'max_work': 1})
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            cli.main(['realign', '--packages', pk])                                    # no -o and no --plan

    def test_memory_decides_nothing(self):
        # whether a region is aligned is decided by max_region_memory, from the package; the budget,
        # max_memory, decides only how many regions run at once
        d = tempfile.mkdtemp()
        try:
            pk = os.path.join(d, 'packages.jsonl.gz')
            unit = 'ACGTTGCAAGTCAT'
            big = pkg('BIG', [('a1', unit * 215, 1), ('a2', unit * 100 + 'G' + unit * 115, 1)])   # 3 kb alleles
            pkgs = [pkg('R0', [('a1', unit * 3, 1), ('a2', unit * 4, 1)]),
                    pkg('R1', [('a1', unit * 2, 1), ('a2', unit, 1)]), big,
                    pkg('R2', [('a1', unit * 5, 1), ('a2', unit * 3, 1)])]
            with gzip.open(pk, 'wt') as f:
                for p in pkgs:
                    f.write(json.dumps(p) + '\n')
            limit = 0.2       # GB: two small regions fit at once, and BIG, about 0.27 by its unbanded pair, does not
            self.assertGreater(realign.peak_bytes(big), limit * memory.GB)
            self.assertLess(2 * max(realign.peak_bytes(p) for p in pkgs if p is not big), limit * memory.GB)
            out = os.path.join(d, 'msas.jsonl.gz')
            # BIG is admitted but needs more than the budget: nothing starts, and the error says how much
            with self.assertRaisesRegex(RuntimeError, 'the largest, BIG, needs about 0.2'):
                realign.run(pk, out, jobs=2, max_memory=limit)
            self.assertFalse(os.path.exists(out))
            # which plan says beforehand, from the packages alone; its need_gb is enough
            pl = realign.plan(pk)
            self.assertEqual((pl['admitted'], pl['too_big'], pl['largest']), (4, {}, 'BIG'))
            self.assertGreaterEqual(pl['need_gb'] * memory.GB, realign.peak_bytes(big))
            self.assertEqual(realign.run(pk, out, jobs=2, max_memory=pl['need_gb'])['realigned'], {'ok': 4})
            pl = realign.plan(pk, max_region_memory=0.25)
            self.assertEqual((pl['admitted'], pl['too_big']), (3, {'max_region_memory': 1}))
            self.assertNotEqual(pl['largest'], 'BIG')
            self.assertLess(pl['need_gb'] * memory.GB, realign.peak_bytes(big))
            # under a max_region_memory BIG exceeds, it is too_big whatever the budget, and the rest align
            outs = []
            for budget in (limit, 4):
                o = os.path.join(d, 'msas.%s.jsonl.gz' % budget)
                s = realign.run(pk, o, jobs=2, max_memory=budget, max_region_memory=0.25)
                self.assertEqual(s['realigned'], {'ok': 3, 'too_big': 1})
                self.assertLessEqual(s['peak_estimate_gb'], budget)
                outs.append([json.loads(l) for l in gzip.open(o, 'rt')])
            first = outs[0]
            self.assertEqual([r['id'] for r in first], ['R0', 'R1', 'BIG', 'R2'])
            r = first[2]
            self.assertEqual((r['status'], r['seconds']), ('too_big', 0.0))
            self.assertIn('more than max_region_memory 0.25 GB', r['message'])
            self.assertEqual([(r['id'], r['status'], r.get('rows')) for r in outs[1]],
                             [(r['id'], r['status'], r.get('rows')) for r in first])
            # a resumed run keeps the aligned records in its partial file, and decides too_big afresh
            o = os.path.join(d, 'msas.4.jsonl.gz')
            with open(o + '.partial', 'wb') as f:
                for r in first[1:3]:
                    f.write(gzip.compress((json.dumps(r) + '\n').encode()))
            os.remove(o)
            s = realign.run(pk, o, jobs=2, max_memory=4, max_region_memory=0.25)
            self.assertEqual((s['resumed'], s['realigned']), (1, {'ok': 3, 'too_big': 1}))
            again = [json.loads(l) for l in gzip.open(o, 'rt')]
            self.assertEqual(again[2], first[2])
            self.assertEqual([r['rows'] for r in again if r['status'] == 'ok'],
                             [r['rows'] for r in first if r['status'] == 'ok'])
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
