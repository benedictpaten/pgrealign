import gzip
import os
import random
import shutil
import tempfile
import unittest

from pgrealign import catalog

RM = """chr1\t100\t140\t(TCCC)n\t33\t+\tSimple_repeat\tundefined\t32.5\t1
chr1\t150\t300\t(TGTGTTA)n\t11\t+\tSimple_repeat\tundefined\t24.1\t2
chr1\t1000\t1040\tGA-rich\t14\t+\tLow_complexity\tundefined\t27.3\t3
chr1\t5000\t5300\tL1M5\t233\t-\tLINE\tL1\t25.7\t4
chr2\t10\t20\tHSATII\t9\t+\tSatellite\tSatellite\t5.0\t5
chr2\t71\t90\t(A)n\t9\t+\tSimple_repeat\tundefined\t5.0\t6
chrY\t10\t50\t(CA)n\t9\t+\tSimple_repeat\tundefined\t5.0\t7
"""


class TestCatalog(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.rm = os.path.join(self.d, 'rm.bed')
        with open(self.rm, 'w') as f:
            f.write(RM)

    def test_classify(self):
        self.assertEqual(catalog.classify('Simple_repeat', '(TCCC)n'), ('STR', 4, 'TCCC'))
        self.assertEqual(catalog.classify('Simple_repeat', '(TGTGTTA)n'), ('VNTR', 7, 'TGTGTTA'))
        self.assertEqual(catalog.classify('Satellite', 'HSATII'), ('SAT', 0, 'HSATII'))
        self.assertEqual(catalog.classify('Low_complexity', 'GA-rich'), ('LC', 0, 'GA-rich'))
        self.assertIsNone(catalog.classify('LINE', 'L1M5'))

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

    def test_class_by_bp(self):
        # a 30 bp VNTR inside a 200 bp STR span: the STR covers more bp and wins
        R = catalog.Record
        g = [R('c', 0, 200, 'STR', 2, 'CA', 'rm'), R('c', 50, 80, 'VNTR', 9, 'ACGTACGTA', 'trfind')]
        r = catalog.merge(g, 50)[0]
        self.assertEqual((r.cls, r.period, r.sources), ('STR', 2, 'rm:1,trfind:1'))
        # equal bp: VNTR beats STR
        g = [R('c', 0, 50, 'STR', 2, 'CA', 'rm'), R('c', 60, 110, 'VNTR', 9, 'ACGTACGTA', 'trfind')]
        self.assertEqual(catalog.merge(g, 50)[0].cls, 'VNTR')

    @unittest.skipUnless(shutil.which('cc') or shutil.which('clang') or shutil.which('gcc'), 'no C compiler')
    def test_trfind_finds_a_vntr(self):
        rng = random.Random(1)
        rand = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        unit = rand(37)
        seq = rand(2000) + unit * 12 + rand(2000)
        fa = os.path.join(self.d, 'ref.fa')
        with open(fa, 'w') as f:
            f.write('>chrT\n' + seq + '\n')
        os.environ['PGREALIGN_CACHE'] = os.path.join(self.d, 'cache')
        recs = catalog.run_trfind([fa], contigs=None)
        hits = [r for r in recs if r.start < 2000 + 37 * 12 and r.end > 2000]
        self.assertTrue(hits, recs)
        best = max(hits, key=lambda r: r.end - r.start)
        self.assertEqual((best.cls, best.period), ('VNTR', 37))
        self.assertLessEqual(best.start, 2000 + 37)
        self.assertGreaterEqual(best.end, 2000 + 37 * 11)
        self.assertEqual(len(best.motif), 37)


if __name__ == '__main__':
    unittest.main()
