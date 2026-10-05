import random
import shutil
import unittest

from pgrealign import star


class TestPieces(unittest.TestCase):
    def test_kmer_counts_count_copies(self):
        a = star.kmer_counts('ACGTACGTACGTACGTACGT', k=4)
        b = star.kmer_counts('ACGTACGTACGT', k=4)
        self.assertTrue(all(a[x] >= n for x, n in b.items()) and a != b)   # a strict sub-multiset
        self.assertEqual(star.kmer_counts('ACGNACGT', k=4), star.kmer_counts('ACGT', k=4))

    def test_medoid(self):
        seqs = ['ACGTACGTAC' * 3, 'ACGTACGTAC' * 4, 'ACGTACGTAC' * 5, 'TTTTGGGGCC' * 3]
        self.assertEqual(star.medoid(seqs, [1, 1, 1, 1], k=5), 1)
        self.assertEqual(star.medoid(seqs, [1, 1, 10, 1], k=5), 2)  # weight pulls it to the heavy allele

    def test_left_normalise(self):
        # a CA deletion in a CA repeat moves to the left end of the repeat
        a, b = star.left_normalise_pair('GCACACAT', 'GCACA--T')
        self.assertEqual(b, 'G--CACAT')
        self.assertEqual(a, 'GCACACAT')

    def test_place_fragment(self):
        alleles = [('a1', 'AACCGGTT')]
        rows = {'a1': 'AA-CCGGTT'}
        self.assertEqual(star.place_fragment('AAC', 'left', alleles, rows), 'AA-C-----')
        self.assertEqual(star.place_fragment('GTT', 'right', alleles, rows), '------GTT')
        self.assertEqual(star.place_fragment('CCG', 'internal', alleles, rows), '---CCG---')
        self.assertIsNone(star.place_fragment('GGG', 'internal', alleles, rows))
        self.assertIsNone(star.place_fragment('ACC', 'left', alleles, rows))  # not a prefix


    def test_scores(self):
        self.assertEqual(star.SCORES.flags(), ['-M', '8', '-X', '16', '-O', '16,96', '-E', '8,1'])
        self.assertEqual(star.Scores.parse(''), star.Scores())
        self.assertEqual(star.Scores.parse('-X 6').mismatch, 6)
        self.assertEqual(star.Scores.parse('-M 8 -X 16 -O 16,96 -E 8,1'), star.SCORES)
        for bad in ('-X', '-Q 3', '-O 4'):
            with self.assertRaises(ValueError):
                star.Scores.parse(bad)

    def test_respell(self):
        self.assertEqual(star.respell('A-CN-G', 'ARYG'), 'A-RY-G')


class FakeAligner:
    """Pads every sequence with gaps on the right: a valid, if poor, global alignment."""
    def __call__(self, seqs, banded=False):
        w = max(len(s) for s in seqs)
        return [s + '-' * (w - len(s)) for s in seqs]


class TestStarMSAInputs(unittest.TestCase):
    def test_checks(self):
        al = FakeAligner()
        with self.assertRaises(ValueError):
            star.star_msa([('a', 'ACGT', 1)])                        # no aligner
        with self.assertRaises(star.RealignError):
            star.star_msa([], align=al)
        with self.assertRaises(ValueError):
            star.star_msa([('a', 'ACGT', 1)], [('a', 'AC', 'left')], align=al)   # duplicate id
        with self.assertRaises(ValueError):
            star.star_msa([('a', 'ACGT', 1)], [('f', 'AC', 'middle')], align=al)
        with self.assertRaises(star.TooBig):
            star.star_msa([('a', 'ACGT' * 10, 1), ('b', 'ACGT' * 20, 1)], align=al, max_cells=100)

    def test_centre_and_weight_zero(self):
        al = FakeAligner()
        alleles = [('a', 'ACGTACGTAC' * 3, 0), ('b', 'ACGTACGTAC' * 4, 0), ('c', 'TTTTGGGGCC' * 3, 1)]
        self.assertEqual(star.star_msa(alleles, align=al, k=5).centre, 'c')   # the only positive weight
        self.assertEqual(star.star_msa(alleles, align=al, centre='b', k=5).centre, 'b')


@unittest.skipUnless(shutil.which('abpoa'), 'abpoa not on PATH')
class TestStarMSA(unittest.TestCase):
    def test_region(self):
        unit = 'ACGTTGCAAGT'
        alleles = [('a1', 'GATTACA' + unit * 4 + 'CCGG', 5), ('a2', 'GATTACA' + unit * 6 + 'CCGG', 1),
                   ('a3', 'GATTACA' + unit * 2 + 'CCGG', 1), ('a4', '', 1)]
        fragments = [('f1', 'GATTACA' + unit * 5, 'left'), ('f2', unit + 'CCGG', 'right')]
        with star.Abpoa() as align:
            m = star.star_msa(alleles, fragments, align=align)
        self.assertEqual(m.centre, 'a1')
        self.assertEqual(len({len(x) for x in m.rows.values()}), 1)
        self.assertEqual(m.rows['a4'], '-' * len(m.rows['a4']))
        self.assertEqual((m.fragments_exact, m.fragments_aligned), (2, 0))
        # the fragment lies inside a2, the only allele long enough to hold it
        self.assertTrue(all(f == '-' or f == a for f, a in zip(m.rows['f1'], m.rows['a2'])))

    def test_slot_insertions_are_aligned_to_each_other(self):
        unit = 'ACGTTGCAAGT'
        alleles = [('c', 'GATTACA' + unit * 2 + 'CCGG', 5), ('x', 'GATTACA' + unit * 4 + 'CCGG', 1),
                   ('y', 'GATTACA' + unit * 3 + 'CCGG', 1)]
        with star.Abpoa() as align:
            m = star.star_msa(alleles, align=align)
        self.assertEqual(m.centre, 'c')
        self.assertEqual(m.slot_alignments, 1)
        # y's extra copy shares columns with one of x's two
        self.assertGreaterEqual(sum(1 for p, q in zip(m.rows['x'], m.rows['y']) if p != '-' and q != '-'),
                                len(alleles[2][1]))

    def test_unrelated_stretch_is_not_forced_into_columns(self):
        # two alleles share their flanks and differ by 150 unrelated random bases: SCORES deletes one
        # stretch and inserts the other, where abPOA's defaults thread chance matches through it
        rng = random.Random(1)
        x, y = (''.join(rng.choice('ACGT') for _ in range(150)) for _ in range(2))
        flank = 'ACGTTGCATCCGATTAGCCTAGGCTTACG'
        alleles = [('a1', flank + x + flank[::-1], 2), ('a2', flank + y + flank[::-1], 1)]
        def shared(scores):
            with star.Abpoa(scores=scores) as align:
                rows = star.star_msa(alleles, align=align).rows
            return sum(1 for p, q in zip(rows['a1'], rows['a2']) if p != '-' and q != '-') - 2 * len(flank)
        self.assertLess(shared(star.SCORES), 40)
        self.assertGreater(shared(star.Scores()), 100)

    def test_other_characters_are_restored(self):
        alleles = [('a1', 'ACGTRACGTACGTACGTYYACGT', 2), ('a2', 'ACGTRACGTACGTAAAACGTYYACGT', 1)]
        with star.Abpoa() as align:
            m = star.star_msa(alleles, align=align)
        for aid, s, _ in alleles:
            self.assertEqual(m.rows[aid].replace('-', ''), s)

    def test_deadline(self):
        alleles = [('a1', 'ACGTACGTAC', 1), ('a2', 'ACGTTCGTAC', 1)]
        with star.Abpoa(deadline=0) as align:
            with self.assertRaises(star.TimedOut):
                star.star_msa(alleles, align=align)


if __name__ == '__main__':
    unittest.main()
