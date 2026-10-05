import gzip
import os
import shutil
import tempfile
import unittest

from pgrealign import extract, regions

# A small graph: 0 -> 1 (anchor L) -> {2, 3} -> 4 -> 5 (anchor R) -> 6, plus 7 between 2 and 4 for a
# longer allele. Nodes:
#   0 TT  1 AAAA  2 C  3 G  7 TTT  4 CC  5 GGGG  6 ACGT
S = {0: 'TT', 1: 'AAAA', 2: 'C', 3: 'G', 7: 'TTT', 4: 'CC', 5: 'GGGG', 6: 'ACGT'}
REGION = regions.Region('R', 'chrT', 1, 5, 4, 7, 11, 'R', 'VNTR', 'ok', '')


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
        # adjacent anchors: an empty spanning run, which is the empty allele
        self.assertEqual(self.runs('>1>5'), [('spanning', '+', 1, 1)])

    def test_walk_ends_on_an_anchor(self):
        # a walk that starts or ends on an anchor has no step inside at that end, so no run there
        self.assertEqual(self.runs('>5>6'), [])             # starts on R, forward
        self.assertEqual(self.runs('<1<0'), [])             # starts on L, reverse
        self.assertEqual(self.runs('>0>1'), [])             # ends on L, forward
        self.assertEqual(self.runs('<6<5'), [])             # ends on R, reverse
        self.assertEqual(self.runs('>1>2>4>5'), [('spanning', '+', 1, 3)])
        self.assertEqual(self.runs('>5>6>0>1>2>4'), [('suffix', '+', 4, 6)])


class TestExtract(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
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
        # `ref` names the reference sample's allele, and is null for a sample with no spanning run
        self.assertEqual({a['id']: a['seq'] for a in pk['alleles']}[pk['ref']], 'CCC')
        extract.extract(self.g, [reg], os.path.join(self.d, 'x2'), ref_sample='s2')
        self.assertIsNone(next(extract.read_packages(os.path.join(self.d, 'x2', 'packages.jsonl.gz')))['ref'])

    def run_extract(self, walks, name, **kw):
        g = os.path.join(self.d, name + '.gfa')
        with open(g, 'w') as f:
            f.write(gfa(walks))
        out = os.path.join(self.d, name)
        s = extract.extract(g, [REGION], out, **kw)
        pk = list(extract.read_packages(os.path.join(out, 'packages.jsonl.gz')))
        with gzip.open(os.path.join(out, 'occurrences.tsv.gz'), 'rt') as f:
            occ = [l.rstrip('\n').split('\t') for l in f if not l.startswith('#')]
        return s, pk, occ

    def test_complex_and_invariant(self):
        walks = [('CHM13', '0', '>1>2>4>5>6'), ('s1', '1', '>1>2>4>5>6')]
        self.assertEqual(self.run_extract(walks, 'y')[0]['status'], {'invariant': 1})
        s, pk, _ = self.run_extract(walks + [('s9', '1', '>1>2>1')], 'z')
        self.assertEqual(s['status'], {'complex': 1})
        self.assertEqual(pk[0]['complex'], ['s9#1#chrT#0: entered by >L, then >L'])
        # no walk visits the anchors of a region that is not in the graph
        g = os.path.join(self.d, 'y.gfa')
        nowhere = REGION._replace(id='N', L=8, R=9)
        self.assertEqual(extract.extract(g, [nowhere], os.path.join(self.d, 'n'))['status'], {'no_runs': 1})

    def test_walk_ends_on_an_anchor(self):
        # a walk that starts or ends exactly on an anchor touches no interior node there: no fragment
        ref = ('CHM13', '0', '>0>1>2>4>5>6')
        for k, walk in enumerate(('>5>6', '<6<5', '>0>1', '<1<0')):
            with self.subTest(walk=walk):
                s, pk, occ = self.run_extract([ref, ('s1', '1', walk)], 'e%d' % k)
                self.assertEqual(s['status'], {'invariant': 1})
                self.assertEqual(pk[0]['fragments'], [])
                self.assertEqual(pk[0]['n_runs'], {'spanning': 1})
                self.assertEqual([o[0] for o in occ], ['CHM13#0#chrT#0'])

    def test_long_anchor_free_walk(self):
        # s3's anchor-free walk >7>4 spells TTTCC (5 bp) inside the interior
        walks = [('CHM13', '0', '>1>2>4>5>6'), ('s3', '1', '>7>4'), ('s4', '1', '>6')]
        s, pk, occ = self.run_extract(walks, 'l5', max_internal_bp=5)
        self.assertEqual(s['status'], {'ok': 1})
        self.assertEqual([(f['side'], f['seq']) for f in pk[0]['fragments']], [('internal', 'TTTCC')])
        # with a cap one bp shorter, s3's walk cannot be held whole, so the region is left as it is
        s, pk, occ = self.run_extract(walks, 'l4', max_internal_bp=4)
        self.assertEqual(s['status'], {'complex': 1})
        self.assertEqual(s['internal_too_long'], 1)
        self.assertEqual(pk[0]['complex'], ['s3#1#chrT#0: an anchor-free walk longer than max_internal_bp (4)'])
        self.assertEqual([o[0] for o in occ], ['CHM13#0#chrT#0'])
        # an anchor-free walk over the cap outside every interior (s4's, 4 bp) changes nothing
        s, pk, occ = self.run_extract(walks[:1] + walks[2:], 'l3', max_internal_bp=3)
        self.assertEqual((s['status'], s['internal_too_long']), ({'invariant': 1}, 0))


if __name__ == '__main__':
    unittest.main()
