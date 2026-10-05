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
    def test_trfind_build_cache(self):
        cache = os.path.join(self.d, 'cache')
        with mock.patch.dict(os.environ, {'PGREALIGN_CACHE': cache}):
            b = catalog.trfind_build()
            path = catalog.trfind_binary()
            self.assertEqual(path, os.path.join(cache, 'trfind-' + b['key'][:16]))
            self.assertTrue(os.access(path, os.X_OK))
            self.assertEqual(os.listdir(cache), [os.path.basename(path)])  # no temporary left behind
            st = os.stat(path)
            self.assertEqual(catalog.trfind_binary(), path)  # cached: not compiled again
            self.assertEqual(os.stat(path).st_mtime_ns, st.st_mtime_ns)
        self.assertNotEqual(os.environ.get('PGREALIGN_CACHE'), cache)
        with open(os.path.join(os.path.dirname(catalog.__file__), 'c', 'trfind.c'), 'rb') as f:
            self.assertEqual(b['source_sha256'], hashlib.sha256(f.read()).hexdigest())
        # another platform is another key, so another binary
        with mock.patch('platform.machine', return_value='another'):
            self.assertNotEqual(catalog.trfind_build()['key'], b['key'])
        self.assertEqual(catalog.trfind_binary(cache), path)

    def test_trfind_compile_failure(self):
        cache = os.path.join(self.d, 'cache')
        with mock.patch.dict(os.environ, {'CC': 'false', 'PGREALIGN_CACHE': cache}):
            with self.assertRaisesRegex(RuntimeError, 'compiling trfind failed'):
                catalog.trfind_binary()
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
        self.assertEqual(s['trfind'], catalog.trfind_build())
        hits = [r for r in recs if r.start < 2000 + 37 * 12 and r.end > 2000]
        self.assertTrue(hits, recs)
        best = max(hits, key=lambda r: r.end - r.start)
        self.assertEqual((best.chrom, best.cls, best.period), (name, 'VNTR', 37))
        self.assertLessEqual(best.start, 2000 + 37)
        self.assertGreaterEqual(best.end, 2000 + 37 * 11)
        self.assertEqual(len(best.motif), 37)


if __name__ == '__main__':
    unittest.main()
