import unittest

from pgrealign import induce


class TestInduce(unittest.TestCase):
    def test_induce(self):
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
        rows = {'a': 'A' * 2500, 'b': 'A' * 2500}
        g = induce.induce(rows, max_node=1024)
        self.assertEqual([len(g.seq[x]) for x in g.walks['a']], [1024, 1024, 452])


if __name__ == '__main__':
    unittest.main()
