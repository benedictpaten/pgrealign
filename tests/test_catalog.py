import gzip
import hashlib
import os
import random
import shutil
import tempfile
import unittest
from unittest import mock

from pgrealign import catalog

RM = """chr1\t100\t140\t(TCCC)n\t33\t+\tSimple_repeat\tundefined\t32.5\t1
chr1\t150\t300\t(TGTGTTA)n\t11\t+\tSimple_repeat\tundefined\t24.1\t2
chr1\t1000\t1040\tGA-rich\t14\t+\tLow_complexity\tundefined\t27.3\t3
chr1\t5000\t5300\tL1M5\t233\t-\tLINE\tL1\t25.7\t4
chr2\t10\t20\tHSATII\t9\t+\tSatellite\tSatellite\t5.0\t5
chr2\t71\t90\t(A)n\t9\t+\tSimple_repeat\tundefined\t5.0\t6
chrY\t10\t50\t(CA)n\t9\t+\tSimple_repeat\tundefined\t5.0\t7
"""


HAVE_CC = bool(shutil.which('cc') or shutil.which('clang') or shutil.which('gcc'))


def short_repeats(seed):
    """A sequence of 4,500 bp holding perfect and mutated short tandem repeats of periods 1, 2, 4 and 6."""
    rng = random.Random(seed)
    rand = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    mutate = lambda s, rate: ''.join(rng.choice('ACGT'.replace(c, '')) if rng.random() < rate else c for c in s)
    return ''.join([rand(600), 'CA' * 30, rand(400), 'A' * 25, rand(400), mutate('AATG' * 20, 0.08), rand(400),
                    'ACGTTG' * 8, rand(300), mutate('TTAGGC' * 15, 0.12), rand(500)])


class TestCatalog(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d)
        self.rm = os.path.join(self.d, 'rm.bed')
        with open(self.rm, 'w') as f:
            f.write(RM)

    def read(self, text, name='t.bed'):
        path = os.path.join(self.d, name)
        with open(path, 'w') as f:
            f.write(text)
        return list(catalog.read_catalog(path))

    def test_classify(self):
        self.assertEqual(catalog.classify('Simple_repeat', '(TCCC)n'), ('STR', 4, 'TCCC'))
        self.assertEqual(catalog.classify('Simple_repeat', '(TGTGTTA)n'), ('VNTR', 7, 'TGTGTTA'))
        self.assertEqual(catalog.classify('Satellite', 'HSATII'), ('SAT', 0, 'HSATII'))
        self.assertEqual(catalog.classify('Low_complexity', 'GA-rich'), ('LC', 0, 'GA-rich'))
        self.assertIsNone(catalog.classify('LINE', 'L1M5'))
        # a Simple_repeat whose name is not (MOTIF)n is STR with period 0, however long the name
        self.assertEqual(catalog.classify('Simple_repeat', 'GA'), ('STR', 0, 'GA'))
        self.assertEqual(catalog.classify('Simple_repeat', 'polyA-rich'), ('STR', 0, 'polyA-rich'))

    def test_merge_and_classes(self):
        recs = list(catalog.read_repeatmasker(self.rm, catalog.DEFAULT_CONTIGS))
        self.assertEqual(len(recs), 5)  # the LINE and chrY dropped
        self.assertTrue(all(r.src == 'rm' for r in recs))
        regions = catalog.merge(recs, gap=50)
        self.assertEqual([(r.chrom, r.start, r.end) for r in regions],
                         [('chr1', 100, 300), ('chr1', 1000, 1040), ('chr2', 10, 20), ('chr2', 71, 90)])
        # the longest member decides: the 150 bp VNTR beats the 40 bp STR
        self.assertEqual((regions[0].cls, regions[0].period, regions[0].motif, regions[0].n_members),
                         ('VNTR', 7, 'TGTGTTA', 2))
        self.assertEqual(regions[0].member_classes, 'STR:1,VNTR:1')
        self.assertEqual(regions[0].sources, 'rm:2')
        self.assertEqual(regions[0].id, 'TR_chr1_100')
        # the SAT ends at 20 and the STR starts at 71, 51 bp apart: separate at gap 50, merged at 51
        self.assertEqual(len(catalog.merge(recs, gap=51)), 3)
        merged = catalog.merge(recs, gap=51)[2]
        self.assertEqual((merged.start, merged.end, merged.cls, merged.member_classes), (10, 90, 'STR', 'SAT:1,STR:1'))

    def test_roundtrip(self):
        out = os.path.join(self.d, 'cat.bed.gz')
        s = catalog.build(out, repeatmasker=self.rm, gap=50, contigs=['chr1', 'chr2'])
        self.assertEqual(s['regions'], 4)
        back = list(catalog.read_catalog(out))
        self.assertEqual(back, catalog.merge(list(catalog.read_repeatmasker(self.rm, ['chr1', 'chr2'])), 50))
        with gzip.open(out, 'rt') as f:
            self.assertTrue(f.readline().startswith('#chrom'))
        self.assertNotIn('trfind', s)  # no FASTA, so trfind did not run

    def test_read_bed3(self):
        got = self.read('track name=t\n# a comment\nchr1\t100\t200\nchr1\t300\t400\n\nchr2\t0\t5\n')
        self.assertEqual(got, [catalog.Region('chr1', 100, 200, 'R_chr1_100', '.', 0, '.', 1, '.', '.'),
                               catalog.Region('chr1', 300, 400, 'R_chr1_300', '.', 0, '.', 1, '.', '.'),
                               catalog.Region('chr2', 0, 5, 'R_chr2_0', '.', 0, '.', 1, '.', '.')])

    def test_read_bed6_with_strand(self):
        got = self.read('chr1\t100\t200\tA\t0\t-\nchr2\t5\t50\tB\t960\t+\nchr2\t60\t70\t.\t0\t.\n')
        self.assertEqual([(r.chrom, r.start, r.end, r.id, r.cls, r.period) for r in got],
                         [('chr1', 100, 200, 'A', '.', 0), ('chr2', 5, 50, 'B', '.', 0),
                          ('chr2', 60, 70, 'R_chr2_60', '.', 0)])

    def test_read_bed12(self):
        # twelve columns but no catalogue header: a BED, not a catalogue
        line = 'chr1\t1000\t5000\tgene\t0\t+\t1000\t5000\t0\t2\t500,400\t0,3600\n'
        got = self.read('#chrom\tstart\tend\tname\n' + line)
        self.assertEqual([(r.chrom, r.start, r.end, r.id, r.cls, r.n_members) for r in got],
                         [('chr1', 1000, 5000, 'gene', '.', 1)])

    def test_read_repeated_names(self):
        got = self.read('chr1\t100\t200\tX\n'
                        'chr2\t100\t200\tX\n'
                        'chr1\t300\t400\tY\n'
                        'chr1\t100\t150\tX\n'  # the same name, contig and start as the first
                        'chr1\t500\t600\tX_chr1_500\n'  # the id that the next X would take
                        'chr1\t500\t600\tX\n'
                        'chr1\t700\t800\n'
                        'chr1\t700\t900\n')
        ids = [r.id for r in got]
        self.assertEqual(ids, ['X_chr1_100', 'X_chr2_100', 'Y', 'X_chr1_100_2', 'X_chr1_500',
                               'X_chr1_500_2', 'R_chr1_700', 'R_chr1_700_2'])
        self.assertEqual([(r.chrom, r.start, r.end) for r in got][3], ('chr1', 100, 150))

    def test_read_catalogue(self):
        regions = catalog.merge(list(catalog.read_repeatmasker(self.rm, ['chr1', 'chr2'])), 50)
        out = os.path.join(self.d, 'cat.bed')
        catalog.write_catalog(regions, out)
        self.assertEqual(list(catalog.read_catalog(out)), regions)
        with open(out) as f:
            lines = f.readlines()
        # a repeated id in a catalogue is an error, not a silent collapse downstream
        with self.assertRaisesRegex(ValueError, 'TR_chr1_100 repeats'):
            self.read(''.join(lines + lines[1:2]))
        with self.assertRaisesRegex(ValueError, 'line 2: a catalogue line has 10 columns, not 3'):
            self.read(lines[0] + 'chr1\t1\t2\n')
        self.assertEqual(self.read(lines[0]), [])

    def test_read_bad_lines(self):
        for text, msg in [('chr1\t100\n', 'line 1: a BED line has at least 3 columns, not 2'),
                          ('chr1\t100\t200\nchr1\tx\t200\n', 'line 2: invalid literal'),
                          ('chr1\t200\t100\n', 'line 1: start 200 and end 100'),
                          ('chr1\t-1\t100\n', 'line 1: start -1')]:
            with self.assertRaisesRegex(ValueError, msg):
                self.read(text)

    @unittest.skipUnless(HAVE_CC, 'no C compiler')
    def test_tool_build_cache(self):
        for name in catalog.C_TOOLS:
            cache = os.path.join(self.d, 'cache_' + name)
            with mock.patch.dict(os.environ, {'PGREALIGN_CACHE': cache}):
                b = catalog.tool_build(name)
                path = catalog.tool_binary(name)
                self.assertEqual(path, os.path.join(cache, '%s-%s' % (name, b['key'][:16])))
                self.assertTrue(os.access(path, os.X_OK))
                self.assertEqual(os.listdir(cache), [os.path.basename(path)])  # no temporary left behind
                st = os.stat(path)
                self.assertEqual(catalog.tool_binary(name), path)  # cached: not compiled again
                self.assertEqual(os.stat(path).st_mtime_ns, st.st_mtime_ns)
            self.assertNotEqual(os.environ.get('PGREALIGN_CACHE'), cache)
            with open(os.path.join(os.path.dirname(catalog.__file__), 'c', name + '.c'), 'rb') as f:
                self.assertEqual(b['source_sha256'], hashlib.sha256(f.read()).hexdigest())
            # another platform is another key, so another binary
            with mock.patch('platform.machine', return_value='another'):
                self.assertNotEqual(catalog.tool_build(name)['key'], b['key'])
            self.assertEqual(catalog.tool_binary(name, cache), path)
        self.assertNotEqual(catalog.tool_build('trfind')['key'], catalog.tool_build('strf')['key'])

    def test_tool_compile_failure(self):
        for name in catalog.C_TOOLS:
            cache = os.path.join(self.d, 'cache_' + name)
            with mock.patch.dict(os.environ, {'CC': 'false', 'PGREALIGN_CACHE': cache}):
                with self.assertRaisesRegex(RuntimeError, 'compiling %s failed' % name):
                    catalog.tool_binary(name)
            self.assertEqual(os.listdir(cache), [])

    def test_class_by_bp(self):
        # a 30 bp VNTR inside a 200 bp STR span: the STR covers more bp and wins
        R = catalog.Record
        g = [R('c', 0, 200, 'STR', 2, 'CA', 'rm'), R('c', 50, 80, 'VNTR', 9, 'ACGTACGTA', 'trfind')]
        r = catalog.merge(g, 50)[0]
        self.assertEqual((r.cls, r.period, r.sources), ('STR', 2, 'rm:1,trfind:1'))
        # equal bp: VNTR beats STR
        g = [R('c', 0, 50, 'STR', 2, 'CA', 'rm'), R('c', 60, 110, 'VNTR', 9, 'ACGTACGTA', 'trfind')]
        self.assertEqual(catalog.merge(g, 50)[0].cls, 'VNTR')

    @unittest.skipUnless(HAVE_CC, 'no C compiler')
    def test_trfind_finds_a_vntr(self):
        rng = random.Random(1)
        rand = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        unit = rand(37)
        seq = rand(2000) + unit * 12 + rand(2000)
        fa = os.path.join(self.d, 'ref.fa')
        name = 'chrT_' + 'x' * 300  # a long name comes back whole
        with open(fa, 'w') as f:
            f.write('>%s description\n%s\n' % (name, seq))
        with mock.patch.dict(os.environ, {'PGREALIGN_CACHE': os.path.join(self.d, 'cache')}):
            recs = catalog.run_trfind([fa], contigs=None)
            s = catalog.build(os.path.join(self.d, 'cat.bed'), fastas=[fa])
        self.assertEqual(s['trfind'], catalog.tool_build('trfind'))
        hits = [r for r in recs if r.start < 2000 + 37 * 12 and r.end > 2000]
        self.assertTrue(hits, recs)
        best = max(hits, key=lambda r: r.end - r.start)
        self.assertEqual((best.chrom, best.cls, best.period), (name, 'VNTR', 37))
        self.assertLessEqual(best.start, 2000 + 37)
        self.assertGreaterEqual(best.end, 2000 + 37 * 11)
        self.assertEqual(len(best.motif), 37)

    @unittest.skipUnless(HAVE_CC, 'no C compiler')
    def test_prebuilt_trfind(self):
        rng = random.Random(2)
        rand = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        fa = os.path.join(self.d, 'ref.fa')
        with open(fa, 'w') as f:
            f.write('>chrT\n%s\n' % (rand(1000) + rand(23) * 10 + rand(1000)))
        with mock.patch.dict(os.environ, {'PGREALIGN_CACHE': os.path.join(self.d, 'cache')}):
            built = catalog.tool_binary('trfind')
            want = catalog.run_trfind([fa])
        # a binary built elsewhere, found by path or on PATH, finds the same repeats and compiles nothing
        binary = os.path.join(self.d, 'bin', 'trfind')
        os.makedirs(os.path.dirname(binary))
        shutil.copy(built, binary)
        nocc = {'PGREALIGN_CACHE': os.path.join(self.d, 'nocache'), 'CC': 'false'}
        with mock.patch.dict(os.environ, {**nocc, 'PATH': os.path.dirname(binary) + os.pathsep + os.environ['PATH']}):
            self.assertEqual(catalog.run_trfind([fa], binary=binary), want)
            s = catalog.build(os.path.join(self.d, 'cat.bed'), fastas=[fa], contigs=['chrT'], trfind='trfind')
        self.assertTrue(want)
        self.assertEqual(s['trfind'], {'binary': binary})
        self.assertFalse(os.path.exists(nocc['PGREALIGN_CACHE']))
        with self.assertRaisesRegex(RuntimeError, 'not found'):
            catalog.run_trfind([fa], binary=os.path.join(self.d, 'missing'))

    @unittest.skipUnless(HAVE_CC, 'no C compiler')
    def test_no_repeats_is_an_empty_catalogue(self):
        # a reference with no tandem repeat, such as a short contig, gives a catalogue with no entry
        fa = os.path.join(self.d, 'ref.fa')
        with open(fa, 'w') as f:
            f.write('>chrT\nACGTTGCAAGTCATGGA\n')
        out = os.path.join(self.d, 'cat.bed.gz')
        with mock.patch.dict(os.environ, {'PGREALIGN_CACHE': os.path.join(self.d, 'cache')}):
            s = catalog.build(out, fastas=[fa], contigs=['chrT'])
        self.assertEqual((s['records'], s['regions']), (0, 0))
        self.assertEqual(list(catalog.read_catalog(out)), [])
        with self.assertRaisesRegex(ValueError, 'no sources'):
            catalog.build(out)

    def test_trf_batches(self):
        seqs = {'a': 'x' * 25, 'b': 'y' * 10, 'c': 'z' * 4}
        got = [(c, o, len(s)) for c, o, s in catalog.trf_batches(seqs, size=10, overlap=3)]
        # windows of 10 every 7 bp, the last reaching the end; a sequence no longer than one window is one
        self.assertEqual(got, [('a', 0, 10), ('a', 7, 10), ('a', 14, 10), ('a', 21, 4), ('b', 0, 10), ('c', 0, 4)])

    def test_parse_trf(self):
        lines = ['@1 a description',
                 '11 40 2 15.0 2 100 0 60 50 50 0 0 1.00 AC ACAC . .',
                 '50 61 6 2.0 6 100 0 24 50 50 0 0 1.00 ACGTTG ACGTTGACGTTG . .',   # too few copies
                 '@0',
                 '1 300 9 33.3 9 95 1 500 25 25 25 25 2.00 ACGTTGCAA ACG . .']
        recs = catalog.parse_trf(lines, [('chrT', 0), ('chrT', 7000)], min_copies=4)
        R = catalog.Record
        self.assertEqual(recs, [R('chrT', 7010, 7040, 'STR', 2, 'AC', 'trf'),
                                R('chrT', 0, 300, 'VNTR', 9, 'ACGTTGCAA', 'trf')])

    @unittest.skipUnless(shutil.which('trf'), 'trf not on PATH')
    def test_run_trf(self):
        # a (CA)n of 40 copies across a batch boundary is found whole, on the contig's coordinates
        rng = random.Random(3)
        rand = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        fa = os.path.join(self.d, 'ref.fa')
        with open(fa, 'w') as f:
            f.write('>chrT description\n%s\n>chrU\n%s\n' % (rand(950) + 'CA' * 40 + rand(950), rand(500)))
        with mock.patch.object(catalog, 'TRF_BATCH', 1000), mock.patch.object(catalog, 'TRF_OVERLAP', 200):
            recs = catalog.run_trf([fa], 'trf', contigs=['chrT'], jobs=2)
        hits = [r for r in recs if r.period == 2]
        self.assertTrue(hits, recs)
        self.assertEqual({(r.chrom, r.cls, r.src) for r in hits}, {('chrT', 'STR', 'trf')})
        self.assertEqual(min(r.start for r in hits), 950)
        self.assertEqual(max(r.end for r in hits), 1030)
        if HAVE_CC:
            with mock.patch.dict(os.environ, {'PGREALIGN_CACHE': os.path.join(self.d, 'cache')}):
                s = catalog.build(os.path.join(self.d, 'cat.bed'), fastas=[fa], contigs=['chrT'], trf='trf')
            self.assertEqual(s['trf']['binary'], os.path.abspath(shutil.which('trf')))
            self.assertIn('trf/STR', s['records_by_source_class'])
        with self.assertRaisesRegex(ValueError, 'TRF needs reference FASTA'):
            catalog.build(os.path.join(self.d, 'cat.bed'), repeatmasker=self.rm, trf='trf')

    def short_repeats_fasta(self):
        fa = os.path.join(self.d, 'ref.fa')
        with open(fa, 'w') as f:
            f.write('>chrT description\n%s\n>chrU\n%s\n' % (short_repeats(4), short_repeats(5)))
        return fa

    @unittest.skipUnless(HAVE_CC, 'no C compiler')
    def test_strf(self):
        # the records TRF finds in chrT (test_strf_matches_trf), here on batches of 1,000 bp that cut the
        # (CA)n and the mutated (AATG)n
        fa = self.short_repeats_fasta()
        R = catalog.Record
        want = [R('chrT', 600, 660, 'STR', 2, 'CA', 'strf'), R('chrT', 1059, 1085, 'STR', 1, 'A', 'strf'),
                R('chrT', 1485, 1565, 'STR', 4, 'AATG', 'strf'), R('chrT', 1965, 2013, 'STR', 6, 'ACGTTG', 'strf'),
                R('chrT', 2313, 2404, 'STR', 6, 'TTAGGC', 'strf')]
        out = os.path.join(self.d, 'cat.bed')
        stages = {name: list(params) for name, params, _ in catalog.TRF_STAGES}
        with mock.patch.object(catalog, 'TRF_BATCH', 1000), mock.patch.object(catalog, 'TRF_OVERLAP', 200):
            with mock.patch.dict(os.environ, {'PGREALIGN_CACHE': os.path.join(self.d, 'cache')}):
                self.assertEqual(catalog.run_strf([fa], contigs=['chrT'], jobs=2), want)
                s = catalog.build(out, fastas=[fa], contigs=['chrT'], strf=True)
                built = catalog.tool_binary('strf')
            self.assertEqual(s['strf'], {**catalog.tool_build('strf'), 'stages': stages})
            self.assertEqual(s['records_by_source_class']['strf/STR'], 5)
            self.assertEqual([r.sources for r in catalog.read_catalog(out) if 'strf' in r.sources], ['strf:1'] * 5)
            # a strf built elsewhere, found by path, finds the same repeats and compiles nothing
            binary = os.path.join(self.d, 'bin', 'strf')
            os.makedirs(os.path.dirname(binary))
            shutil.copy(built, binary)
            nocache = os.path.join(self.d, 'nocache')
            with mock.patch.dict(os.environ, {'PGREALIGN_CACHE': nocache, 'CC': 'false'}):
                self.assertEqual(catalog.run_strf([fa], binary, contigs=['chrT']), want)
            self.assertFalse(os.path.exists(nocache))
            with mock.patch.dict(os.environ, {'PGREALIGN_CACHE': os.path.join(self.d, 'cache')}):
                s = catalog.build(out, fastas=[fa], contigs=['chrT'], strf=binary)
            self.assertEqual(s['strf'], {'binary': binary, 'stages': stages})
        with self.assertRaisesRegex(RuntimeError, 'strf binary .* not found'):
            catalog.run_strf([fa], os.path.join(self.d, 'missing'))
        with self.assertRaisesRegex(ValueError, 'not both'):
            catalog.build(out, fastas=[fa], trf='trf', strf=True)
        with self.assertRaisesRegex(ValueError, 'strf needs reference FASTA'):
            catalog.build(out, repeatmasker=self.rm, strf=True)

    @unittest.skipUnless(HAVE_CC and shutil.which('trf'), 'no C compiler, or trf not on PATH')
    def test_strf_matches_trf(self):
        fa = self.short_repeats_fasta()
        with mock.patch.object(catalog, 'TRF_BATCH', 1000), mock.patch.object(catalog, 'TRF_OVERLAP', 200), \
                mock.patch.dict(os.environ, {'PGREALIGN_CACHE': os.path.join(self.d, 'cache')}):
            trf = catalog.run_trf([fa], 'trf', jobs=2)
            strf = catalog.run_strf([fa], jobs=2)
        self.assertEqual(len(trf), 11)
        self.assertEqual([r._replace(src='trf') for r in strf], trf)


def mutated(seq, rate, rng):
    return ''.join(rng.choice('ACGT'.replace(c, '')) if rng.random() < rate else c for c in seq)


def hor_array(rng, copies, monomer=171, types=4):
    """A higher-order repeat array: `copies` copies, each 1% diverged, of a unit of `types` monomers that
    are each a copy of one monomer, 25% diverged."""
    base = ''.join(rng.choice('ACGT') for _ in range(monomer))
    unit = ''.join(mutated(base, 0.25, rng) for _ in range(types))
    return ''.join(mutated(unit, 0.01, rng) for _ in range(copies))


class TestArrays(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d)
        self.rng = random.Random(7)
        self.rand = lambda n: ''.join(self.rng.choice('ACGT') for _ in range(n))

    def test_monomer_multiple(self):
        for period, want in ((171, True), (170, True), (342, True), (2057, True), (2040, True), (153, False),
                             (150, False), (5, False), (256, False), (2100, False)):
            self.assertEqual(catalog.monomer_multiple(period, 171), want, period)

    def test_base_period(self):
        unit = self.rand(50)
        seq = (self.rand(300) + ''.join(mutated(unit, 0.03, self.rng) for _ in range(60))).encode()
        self.assertEqual(catalog.base_period(seq, 300, 1000, 351), 50)
        # two types of monomer that alternate have the period of the pair
        pair = unit + mutated(unit, 0.3, self.rng)
        seq = ''.join(mutated(pair, 0.03, self.rng) for _ in range(30)).encode()
        self.assertEqual(catalog.base_period(seq, 0, 1000, 351), 100)
        # too near the end to test
        self.assertIsNone(catalog.base_period(seq, len(seq) - 1350, 1000, 351))
        self.assertIsNotNone(catalog.base_period(seq, len(seq) - 1351, 1000, 351))

    def test_find_arrays(self):
        # an alpha-satellite array, a satellite of a 5 bp unit that trfind reported at two monomers, and
        # an array too short to keep
        alpha, short = hor_array(self.rng, 40), hor_array(self.rng, 4)
        hsat = ''.join(mutated('AATGG', 0.1, self.rng) for _ in range(6000))
        pad = self.rand(2000)
        seq = pad + alpha + pad + hsat + pad + short + pad
        at = lambda part: seq.index(part)
        records = [('c', at(alpha), at(alpha) + 12000, 684), ('c', at(alpha) + 12500, at(alpha) + len(alpha), 171),
                   ('c', at(alpha), at(alpha) + 5000, 37),              # no monomer's period: not part of it
                   ('c', at(hsat), at(hsat) + len(hsat), 342),
                   ('c', at(short), at(short) + len(short), 684),
                   ('other', 0, 50000, 171)]                            # a contig with no sequence
        seqs = {'c': seq.encode()}
        found = catalog.find_arrays(records, seqs, 171, 20000, 1000)
        self.assertEqual(found, [catalog.Array('c', at(alpha), at(alpha) + len(alpha), 'array_c_1', 3, 1.0)])
        # the two stretches are one array only when the gap allows, and neither is long enough alone
        self.assertEqual(catalog.find_arrays(records, seqs, 171, 20000, 400), [])
        both = catalog.find_arrays(records, seqs, 171, 10000, 400)
        self.assertEqual([(a.start, a.end, a.id) for a in both],
                         [(at(alpha), at(alpha) + 12000, 'array_c_1'),
                          (at(alpha) + 12500, at(alpha) + len(alpha), 'array_c_2')])
        # the short satellite is long enough, and its windows never have the monomer's period
        self.assertEqual(len(hsat), 30000)
        self.assertEqual(catalog.find_arrays(records[3:4], seqs, 171, 20000, 1000), [])
        self.assertEqual(len(catalog.find_arrays(records[3:4], seqs, 171, 20000, 1000, min_share=0)), 1)

    @unittest.skipUnless(HAVE_CC, 'no C compiler')
    def test_arrays_of_a_reference(self):
        from pgrealign import cli
        alpha = hor_array(self.rng, 40)
        fa, other = os.path.join(self.d, 'ref.fa'), os.path.join(self.d, 'other.fa')
        with open(fa, 'w') as f:
            f.write('>chrT\n%s\n' % (self.rand(3000) + alpha + self.rand(3000)).lower())
        with open(other, 'w') as f:
            f.write('>chrU\n%s\n' % self.rand(30000))
        out = os.path.join(self.d, 'arrays.bed')
        with mock.patch.dict(os.environ, {'PGREALIGN_CACHE': os.path.join(self.d, 'cache')}):
            s = catalog.build_arrays(out, [other, fa], min_length=20000, jobs=2)
            self.assertEqual(s['trfind'], catalog.tool_build('trfind'))
            self.assertEqual((s['arrays'], [a['id'] for a in s['by_array']]), (1, ['array_chrT_1']))
            a = s['by_array'][0]
            self.assertEqual((a['chrom'], a['share'], s['bp']), ('chrT', 1.0, a['end'] - a['start']))
            self.assertLessEqual(abs(a['start'] - 3000), 684)
            self.assertLessEqual(abs(a['end'] - 3000 - len(alpha)), 684)
            with open(out) as f:
                self.assertEqual(f.read(), 'chrT\t%d\t%d\tarray_chrT_1\n' % (a['start'], a['end']))
            # any BED of targets is read the same way
            self.assertEqual([r.id for r in catalog.read_catalog(out)], ['array_chrT_1'])
            # a contig left out, or a reference with no array, gives an empty BED
            gz = os.path.join(self.d, 'none.bed.gz')
            self.assertEqual(catalog.build_arrays(gz, [fa, other], min_length=20000, contigs=['chrU'])['arrays'], 0)
            with gzip.open(gz, 'rt') as f:
                self.assertEqual(f.read(), '')
            with self.assertRaisesRegex(ValueError, 'no reference'):
                catalog.build_arrays(out, [])
            # the command line
            with mock.patch('sys.stdout'):
                cli.main(['arrays', '--reference', fa, '--min-length', '20000', '-o', out])
            with open(out) as f:
                self.assertEqual(f.read(), 'chrT\t%d\t%d\tarray_chrT_1\n' % (a['start'], a['end']))
            self.assertTrue(os.path.exists(out + '.manifest.json'))


if __name__ == '__main__':
    unittest.main()
