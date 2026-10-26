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


    def test_sketch_estimates_the_exact_distance(self):
        rng = random.Random(2)
        unit = 'ACGTTGCAAGTCATGGA'
        seqs = [unit * rng.randint(20, 60) + ''.join(rng.choice('ACGT') for _ in range(80)) for _ in range(2)]
        c = [star.kmer_counts(s) for s in seqs]
        inter = sum(min(c[0][x], c[1][x]) for x in c[0].keys() & c[1].keys())
        exact = 1 - inter / (sum(c[0].values()) + sum(c[1].values()) - inter)
        est = star.sketch_distance(star.sketch(c[0]), star.sketch(c[1]))
        self.assertAlmostEqual(est, exact, delta=0.05)
        self.assertEqual(star.sketch(c[0]), star.sketch(c[0]))            # the hash is fixed
        self.assertEqual(star.sketch_distance(star.sketch(c[0]), star.sketch(c[0])), 0.0)

    def test_medoid_by_sketch(self):
        seqs = ['ACGTACGTAC' * 3, 'ACGTACGTAC' * 4, 'ACGTACGTAC' * 5, 'TTTTGGGGCC' * 3]
        self.assertEqual(star.medoid(seqs, [1, 1, 1, 1], k=5, max_work=0), 1)

    def test_score(self):
        s = star.SCORES                                          # match 8, mismatch 16, gaps min(16 + 8g, 96 + g)
        self.assertEqual(s.score('ACGT', 'ACGT'), 32)
        self.assertEqual(s.score('ACGT', 'ACCT'), 8)
        self.assertEqual(s.score('ANGT', 'ACGN'), 16)             # N scores nothing
        self.assertEqual(s.score('AC--GT', 'ACTTGT'), 32 - 32)
        self.assertEqual(s.score('A-C', 'AG-'), 8 - 24 - 24)       # an insertion then a deletion: two gaps
        self.assertEqual(s.score('A' + '-' * 100 + 'C', 'A' + 'G' * 100 + 'C'), 16 - 196)   # the long gap's cost
        self.assertEqual(s.score('A--C', 'A--C'), 16)              # a column with no base counts for nothing
        self.assertEqual(s.score('', ''), 0)
        self.assertEqual((s.gap(0), s.gap(1), s.gap(10), s.gap(11), s.gap(12)), (0, 24, 96, 104, 108))

    def test_band_width(self):
        self.assertEqual(star.default_band(30000), 310)
        # the widest band whose |n - m| + 2w cells a row keep within max_cells, never narrower than abPOA's own
        self.assertEqual(star.band_width(30000, 30000, 450_000_000), 7500)
        self.assertEqual(star.band_width(30000, 32000, 450_000_000), 6500)
        self.assertEqual(star.band_width(30000, 60000, 450_000_000), 610)
        n, m = 27709, 27752
        w = star.band_width(n, m, 450_000_000)
        self.assertLessEqual(n * (abs(n - m) + 2 * w), 450_000_000)
        self.assertGreater(n * (abs(n - m) + 2 * (w + 1)), 450_000_000)

    def test_abpoa_flags(self):
        al, s = star.Abpoa(), star.SCORES.flags()
        self.assertEqual(al.flags(['A', 'C'], False), s + ['-b', '-1'])          # unbanded
        self.assertEqual(al.flags(['A', 'C'], True), s)                         # abPOA's own band
        self.assertEqual(al.flags(['A', 'C'], 7500), s + ['-b', '7500', '-f', '0'])
        self.assertEqual(al.flags(['A', 'C'], True, seeded=True), s + ['-S'])

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
    def __init__(self):
        self.calls = []

    def __call__(self, seqs, banded=False, seeded=False):
        self.calls.append((len(seqs), banded, seeded))
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

    def test_large_pairs_are_banded(self):
        al = FakeAligner()
        m = star.star_msa([('a', 'ACGT' * 10, 2), ('b', 'ACGT' * 20, 1), ('c', 'ACGA' * 2, 1)], align=al,
                          centre='a', max_cells=100)
        # a x b is 3,200 cells and a x c 320, both over 100: each is aligned with the widest band within 100
        # cells, which is abPOA's own, 10, and seeded; the two are the same here, so the band's is kept
        self.assertEqual((m.banded_pairs, m.seeded_pairs, m.refined_windows), (2, 0, 0))
        self.assertEqual(al.calls, [(2, 10, False), (2, True, True)] * 2)
        al = FakeAligner()
        m = star.star_msa([('a', 'ACGT' * 10, 2), ('b', 'ACGT' * 20, 1)], align=al, centre='a')
        self.assertEqual((m.banded_pairs, al.calls), (0, [(2, False, False)]))

    def test_centre_and_weight_zero(self):
        al = FakeAligner()
        alleles = [('a', 'ACGTACGTAC' * 3, 0), ('b', 'ACGTACGTAC' * 4, 0), ('c', 'TTTTGGGGCC' * 3, 1)]
        self.assertEqual(star.star_msa(alleles, align=al, k=5).centre, 'c')   # the only positive weight
        self.assertEqual(star.star_msa(alleles, align=al, centre='b', k=5).centre, 'b')

    def test_share(self):
        a, b = star.kmer_counts('ACGTACGGTTCA', 4), star.kmer_counts('ACGTACGGTTCATTGA', 4)
        self.assertEqual(star.share(a, b), 1.0)                  # b is a plus an insertion
        self.assertEqual(star.share(b, a), 1.0)
        self.assertEqual(star.share(star.kmer_counts('AAAATTTT', 4), star.kmer_counts('AAAACCCC', 4)), 0.2)
        self.assertIsNone(star.share(star.kmer_counts('ACG', 4), b))   # too short to hold a 4-mer

    def test_homology(self):
        rng = random.Random(5)
        rand = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        rc = lambda s: s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        centre = rand(400)
        snps = ''.join(c if i % 40 else 'ACGT'[('ACGT'.index(c) + 1) % 4] for i, c in enumerate(centre))
        kept = [('c', centre, 5), ('snps', snps, 2), ('ins', centre[:200] + rand(600) + centre[200:], 1),
                ('del', centre[:100] + centre[300:], 1), ('empty', '', 1), ('short', rand(10), 1)]
        al = FakeAligner()
        m = star.star_msa(kept + [('unrelated', rand(400), 0)], align=al, min_share=0.5)
        self.assertEqual(m.centre, 'c')            # weight 0, so not judged
        for name, other, weight in (('unrelated', rand(400), 1), ('inverted', rc(centre), 2)):
            al = FakeAligner()
            with self.assertRaises(star.NotHomologous) as e:
                star.star_msa(kept + [(name, other, weight)], align=al, min_share=0.5)
            self.assertEqual(e.exception.status, 'not_homologous')
            self.assertEqual(al.calls, [])         # refused before anything is aligned
            msg = str(e.exception)
            self.assertIn('1 of 7 alleles (weight %d) share less than 0.5 of their 15-mers with the centre c'
                          % weight, msg)
            self.assertIn('the least, %s (weight %d)' % (name, weight), msg)
            self.assertTrue(msg.endswith('1.000 reverse-complemented') == (name == 'inverted'), msg)
            # without min_share, or with 0, the region is aligned
            self.assertEqual(star.star_msa(kept + [(name, other, weight)], align=FakeAligner()).centre, 'c')
            star.star_msa(kept + [(name, other, weight)], align=FakeAligner(), min_share=0)
        # a centre given by name is judged against too
        with self.assertRaises(star.NotHomologous):
            star.star_msa(kept, align=FakeAligner(), centre='snps', min_share=0.99)


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

    def guide(self):
        """A centre, an allele that is the centre with 12 bases inserted in its middle stretch, and an alignment
        of the two that matches their flanks but puts the insertion at the stretch's end, so that 15 of the
        stretch's bases are aligned out of register."""
        rng = random.Random(7)
        rand = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
        a, x, b, ins = rand(100), rand(30), rand(100), rand(12)
        c, s = a + x + b, a + x[:15] + ins + x[15:] + b
        return c, s, (a + x + '-' * 12 + b, a + x[:15] + ins + x[15:] + b)

    def test_refine_pair(self):
        c, s, (gc, gs) = self.guide()
        with star.Abpoa() as align:
            (rc, rs), better = star.refine_pair(gc, gs, align, 10 ** 6)
            self.assertEqual((rc.replace('-', ''), rs.replace('-', ''), better), (c, s, 1))
            # the stretch between the two flanks, realigned on its own: all matched, and one 12-base gap
            self.assertEqual(star.SCORES.score(rc, rs), 8 * 230 - star.SCORES.gap(12))
            self.assertGreater(star.SCORES.score(rc, rs), star.SCORES.score(gc, gs))
            self.assertEqual(star.align_pair(c, s, align), (rc, rs))
            # the window is at most the stretch's 15 centre bases and 27 allele bases out of register; over
            # max_cells, it is left as it is
            self.assertEqual(star.refine_pair(gc, gs, align, 15 * 27), ((rc, rs), 1))
            self.assertEqual(star.refine_pair(gc, gs, align, 100), (star.left_normalise_pair(gc, gs), 0))
            # with no run of `anchor` matched columns, the whole alignment is one window
            self.assertEqual(star.refine_pair(gc, gs, align, len(c) * len(s), anchor=200), ((rc, rs), 1))
            self.assertEqual(star.refine_pair(gc, gs, align, len(c) * len(s) - 1, anchor=200)[1], 0)

    def test_long_pair(self):
        c, s, (gc, gs) = self.guide()
        poor = (c + '-' * len(s), '-' * len(c) + s)              # nothing aligned
        for band, seeded, kept in ((poor, (gc, gs), True), ((gc, gs), poor, False)):
            calls = []
            with star.Abpoa() as abpoa:
                def align(seqs, banded=False, seeded_=False, **kw):
                    seeded_ = kw.get('seeded', seeded_)
                    calls.append((banded, seeded_))
                    if banded is True and seeded_:
                        return list(seeded)
                    return list(band) if banded else abpoa(seqs)
                rows, used, better = star.long_pair(c, s, align, len(c) * len(s) - 1)
            self.assertEqual(calls[:2], [(star.band_width(len(c), len(s), len(c) * len(s) - 1), False), (True, True)])
            self.assertEqual(used, kept)
            if kept:           # seeding's kept, and its out-of-register window realigned
                self.assertEqual((better, len(calls), calls[2]), (1, 3, (False, False)))
                with star.Abpoa() as abpoa:
                    self.assertEqual(rows, star.align_pair(c, s, abpoa))
            else:              # the band's kept as it is
                self.assertEqual((rows, better, len(calls)), (star.left_normalise_pair(gc, gs), 0, 2))


if __name__ == '__main__':
    unittest.main()
