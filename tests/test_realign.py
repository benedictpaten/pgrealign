import shutil
import unittest

from pgrealign import realign


class TestPieces(unittest.TestCase):
    def test_kmer_tokens_count_copies(self):
        a = realign.kmer_tokens('ACGTACGTACGTACGTACGT', k=4)
        b = realign.kmer_tokens('ACGTACGTACGT', k=4)
        self.assertTrue(b < a)                 # the shorter array's tokens are a strict subset
        self.assertEqual(realign.kmer_tokens('ACGNACGT', k=4), realign.kmer_tokens('ACGT', k=4))

    def test_medoid(self):
        seqs = ['ACGTACGTAC' * 3, 'ACGTACGTAC' * 4, 'ACGTACGTAC' * 5, 'TTTTGGGGCC' * 3]
        self.assertEqual(realign.medoid(seqs, [1, 1, 1, 1], k=5), 1)
        self.assertEqual(realign.medoid(seqs, [1, 1, 10, 1], k=5), 2)  # weight pulls it to the heavy allele

    def test_left_normalise(self):
        # a CA deletion in a CA repeat moves to the left end of the repeat
        a, b = realign.left_normalise_pair('GCACACAT', 'GCACA--T')
        self.assertEqual(b, 'G--CACAT')
        self.assertEqual(a, 'GCACACAT')

    def test_place_fragment(self):
        alleles = [('a1', 'AACCGGTT')]
        rows = {'a1': 'AA-CCGGTT'}
        self.assertEqual(realign.place_fragment('AAC', 'left', alleles, rows), 'AA-C-----')
        self.assertEqual(realign.place_fragment('GTT', 'right', alleles, rows), '------GTT')
        self.assertEqual(realign.place_fragment('CCG', 'internal', alleles, rows), '---CCG---')
        self.assertIsNone(realign.place_fragment('GGG', 'internal', alleles, rows))
        self.assertIsNone(realign.place_fragment('ACC', 'left', alleles, rows))  # not a prefix


@unittest.skipUnless(shutil.which('abpoa'), 'abpoa not on PATH')
class TestStar(unittest.TestCase):
    def test_region(self):
        unit = 'ACGTTGCAAGT'
        pkg = {'id': 'R', 'alleles': [
            {'id': 'a1', 'weight': 5, 'seq': 'GATTACA' + unit * 4 + 'CCGG'},
            {'id': 'a2', 'weight': 1, 'seq': 'GATTACA' + unit * 6 + 'CCGG'},
            {'id': 'a3', 'weight': 1, 'seq': 'GATTACA' + unit * 2 + 'CCGG'},
            {'id': 'a4', 'weight': 1, 'seq': ''}],
            'fragments': [{'id': 'f1', 'side': 'left', 'seq': 'GATTACA' + unit * 5},
                          {'id': 'f2', 'side': 'right', 'seq': unit + 'CCGG'}]}
        r = realign.realign_package(pkg)
        self.assertEqual(r['status'], 'ok', r)
        self.assertEqual(r['centre'], 'a1')
        rows = r['rows']
        self.assertEqual(len({len(x) for x in rows.values()}), 1)
        self.assertEqual(rows['a2'].replace('-', ''), pkg['alleles'][1]['seq'])
        self.assertEqual(rows['f1'].replace('-', ''), pkg['fragments'][0]['seq'])
        self.assertEqual(rows['a4'], '-' * len(rows['a4']))
        # the fragment lies inside a2, the only allele long enough to hold it
        self.assertTrue(all(f == '-' or f == a for f, a in zip(rows['f1'], rows['a2'])))

    def test_unplaced_and_too_big(self):
        pkg = {'id': 'R', 'alleles': [{'id': 'a1', 'weight': 1, 'seq': 'ACGTACGT'},
                                      {'id': 'a2', 'weight': 1, 'seq': 'ACGTTCGT'}],
               'fragments': [{'id': 'f1', 'side': 'internal', 'seq': 'GGGG'}]}
        self.assertEqual(realign.realign_package(pkg)['status'], 'unplaced_fragment')
        pkg['fragments'] = []
        self.assertEqual(realign.realign_package(pkg, max_cells=10)['status'], 'too_big')


class TestInduce(unittest.TestCase):
    def test_induce(self):
        from pgrealign import induce
        rows = {'a': 'ACG-T', 'b': 'ACGGT', 'c': 'A-G-T', 'd': '-----', 'f': '--G--'}
        g = induce.induce(rows)
        for n, r in rows.items():
            self.assertEqual(''.join(g.seq[x] for x in g.walks[n]), r.replace('-', ''))
        self.assertEqual(g.walks['d'], [])
        # ids increase along every edge
        self.assertTrue(all(a < b for a, b in g.edges))
        # A then C|- then G: the unchopped graph has A, C, G, G(insert), T -> 5 nodes
        self.assertEqual(len(g.seq), 5)

    def test_chop(self):
        from pgrealign import induce
        rows = {'a': 'A' * 2500, 'b': 'A' * 2500}
        g = induce.induce(rows, max_node=1024)
        self.assertEqual([len(g.seq[x]) for x in g.walks['a']], [1024, 1024, 452])


if __name__ == '__main__':
    unittest.main()
