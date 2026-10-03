import gzip
import json
import os
import tempfile
import unittest

from pgrealign import extract, regions
from pgrealign.graph import revcomp

# A small graph: 1 (anchor L) -> {2, 3} -> 4 -> 5 (anchor R) -> 6, plus 7 between 2 and 4 for a longer
# allele. Nodes:
#   1 AAAA  2 C  3 G  7 TTT  4 CC  5 GGGG  6 ACGT
S = {1: 'AAAA', 2: 'C', 3: 'G', 7: 'TTT', 4: 'CC', 5: 'GGGG', 6: 'ACGT'}


def gfa(walks):
    lines = ['H\tVN:Z:1.1']
    lines += ['S\t%d\t%s' % (n, s) for n, s in S.items()]
    for i, (sample, hap, walk) in enumerate(walks):
        lines.append('W\t%s\t%s\tchrT\t0\t0\t%s' % (sample, hap, walk))
    return '\n'.join(lines) + '\n'


WALKS = [
    ('CHM13', '0', '>1>2>4>5>6'),           # reference: allele C+CC
    ('s1', '1', '>1>3>4>5>6'),              # spanning forward: G+CC
    ('s1', '2', '<6<5<4<7<2<1'),            # spanning reverse: forward allele C TTT CC
    ('s2', '1', '>3>4>5>6'),                # prefix: starts inside (G CC), leaves by >R
    ('s2', '2', '>1>2>4'),                  # suffix: enters by >L, ends inside (C CC)
    ('s3', '1', '>7>4'),                    # internal: no anchor, interior nodes only
    ('s4', '1', '>1>2>4>5>6'),              # spanning forward again: C+CC (weight 2 with the reference)
]


class TestRunsOfWalk(unittest.TestCase):
    def runs(self, walk_string, L='1', R='5'):
        toks = [(w[0], w[1:]) for w in walk_string.replace('>', ' >').replace('<', ' <').split()]
        occ = [(i, 'L' if n == L else 'R', o) for i, (o, n) in enumerate(toks) if n in (L, R)]
        return list(extract.runs_of_walk(toks, occ))

    def test_side(self):
        self.assertEqual(extract.side_of('suffix', '+'), 'left')
        self.assertEqual(extract.side_of('prefix', '-'), 'left')
        self.assertEqual(extract.side_of('prefix', '+'), 'right')
        self.assertEqual(extract.side_of('suffix', '-'), 'right')

    def test_kinds(self):
        self.assertEqual(self.runs('>1>2>4>5>6'), [('spanning', '+', 1, 3)])
        self.assertEqual(self.runs('<6<5<4<2<1'), [('spanning', '-', 2, 4)])
        self.assertEqual(self.runs('>3>4>5>6'), [('prefix', '+', 0, 2)])
        self.assertEqual(self.runs('>1>2>4'), [('suffix', '+', 1, 3)])
        self.assertEqual(self.runs('<4<2<1'), [('prefix', '-', 0, 2)])
        self.assertEqual(self.runs('<6<5<4'), [('suffix', '-', 2, 3)])
        self.assertEqual(self.runs('>1>2>1')[0][0], 'complex')     # U-turn out through L
        # two crossings of the same region by one walk (a duplication) are two spanning runs
        self.assertEqual(self.runs('>1>2>4>5>1>3>4>5'), [('spanning', '+', 1, 3), ('spanning', '+', 5, 7)])


class TestExtract(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.g = os.path.join(self.d, 'g.gfa')
        with open(self.g, 'w') as f:
            f.write(gfa(WALKS))

    def test_extract(self):
        reg = regions.Region('TR_chrT_4', 'chrT', 1, 5, 4, 7, 11, 'TR_chrT_4', 'VNTR', 'ok', '')
        s = extract.extract(self.g, [reg], os.path.join(self.d, 'x'))
        self.assertEqual(s['status'], {'ok': 1})
        pk = list(extract.read_packages(os.path.join(self.d, 'x', 'packages.jsonl.gz')))[0]
        alleles = {a['seq']: a['weight'] for a in pk['alleles']}
        self.assertEqual(alleles, {'CCC': 2, 'GCC': 1, 'CTTTCC': 1})
        frags = {(f['side'], f['seq']) for f in pk['fragments']}
        self.assertEqual(frags, {('right', 'GCC'), ('left', 'CCC'), ('internal', 'TTTCC')})
        self.assertEqual(pk['n_runs'], {'spanning': 4, 'prefix': 1, 'suffix': 1, 'internal': 1})
        with gzip.open(os.path.join(self.d, 'x', 'interiors.tsv.gz'), 'rt') as f:
            rows = [l.rstrip('\n').split('\t') for l in f if not l.startswith('#')]
        self.assertEqual(rows, [['TR_chrT_4', '2,3,4,7']])
        with gzip.open(os.path.join(self.d, 'x', 'occurrences.tsv.gz'), 'rt') as f:
            occ = [l.rstrip('\n').split('\t') for l in f if not l.startswith('#')]
        rev = [o for o in occ if o[0] == 's1#2#chrT#0'][0]
        self.assertEqual(rev[2:6], ['spanning', '-', '2', '5'])

    def test_complex_and_invariant(self):
        g2 = os.path.join(self.d, 'g2.gfa')
        with open(g2, 'w') as f:
            f.write(gfa([('CHM13', '0', '>1>2>4>5>6'), ('s1', '1', '>1>2>4>5>6')]))
        reg = regions.Region('R', 'chrT', 1, 5, 4, 7, 11, 'R', 'VNTR', 'ok', '')
        self.assertEqual(extract.extract(g2, [reg], os.path.join(self.d, 'y'))['status'], {'invariant': 1})
        with open(g2, 'a') as f:
            f.write('W\ts9\t1\tchrT\t0\t0\t>1>2>1\n')
        self.assertEqual(extract.extract(g2, [reg], os.path.join(self.d, 'z'))['status'], {'complex': 1})


class TestRegions(unittest.TestCase):
    def test_anchor_and_merge(self):
        # reference: nodes 10..19, each 10 bp, at 0,10,...,90; snarls 10-12, 12-14, 14-16, 16-19
        ref = [(n, '+', (n - 10) * 10, (n - 10) * 10 + 10) for n in range(10, 20)]
        chains = regions.chains_from_pairs([(10, 12), (12, 14), (14, 16), (16, 19)], ref)
        self.assertEqual(len(chains), 1)
        A = regions.Anchors(chains)
        self.assertEqual(A.query(31, 39, 1000), ((20, 30, 12), (40, 50, 14)))
        # an interval overlapping boundary node 12 must take the anchor before it
        self.assertEqual(A.query(25, 35, 1000), ((0, 10, 10), (40, 50, 14)))
        T = lambda i, s, e: type('T', (), {'id': i, 'start': s, 'end': e, 'cls': 'STR'})
        # two targets in snarls 12-14 and 14-16: they share anchor 14 only, so they stay separate
        rs = regions.assign([T('a', 32, 33), T('b', 52, 53)], A, 'chrT', pad=0)
        self.assertEqual([(r.id, r.L, r.R, r.status) for r in rs], [('a', 12, 14, 'ok'), ('b', 14, 16, 'ok')])
        # a target crossing node 14 needs anchors 12 and 16; it overlaps both and they merge into one
        rs = regions.assign([T('a', 32, 33), T('c', 38, 55), T('b', 52, 53)], A, 'chrT', pad=0)
        self.assertEqual([(r.L, r.R, r.targets, r.status) for r in rs], [(12, 16, 'a,c,b', 'ok')])
        # too long for max_span
        rs = regions.assign([T('d', 25, 85)], A, 'chrT', pad=0, max_span=30)
        self.assertEqual(rs[0].status, 'too_long')


if __name__ == '__main__':
    unittest.main()
