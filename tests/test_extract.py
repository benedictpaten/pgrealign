import gzip
import os
import random
import shutil
import tempfile
import unittest
from unittest import mock

from pgrealign import extract, prepare, regions
from pgrealign.graph import revcomp

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

    def test_inversion(self):
        # s1 crosses the reference's interior nodes 2 and 4 backward between the anchors: an inversion
        walks = [('CHM13', '0', '>1>2>4>5>6'), ('s1', '1', '>1<4<2>5>6'), ('s2', '1', '>1>3>4>5>6')]
        s, pk, occ = self.run_extract(walks, 'inv')
        self.assertEqual((s['status'], s['both_ways']), ({'complex': 1}, 1))
        self.assertEqual(pk[0]['complex'], ['s1#1#chrT#0: visits node 4 both ways, as at an inversion'])
        # a reverse run that crosses the inversion forward meets the reference's way after its flip
        rev = [('CHM13', '0', '>1>2>4>5>6'), ('s1', '2', '<6<5>2>4<1')]
        s, pk, occ = self.run_extract(rev, 'inv2')
        self.assertEqual(pk[0]['complex'], ['s1#2#chrT#0: visits node 4 both ways, as at an inversion'])
        # one run that visits a node both ways, as in a hairpin, is complex on its own
        s, pk, occ = self.run_extract([('s1', '1', '>1>2<2>5>6')], 'hairpin')
        self.assertEqual(pk[0]['complex'], ['s1#1#chrT#0: visits node 2 both ways, as at an inversion'])
        # allowed: the region is kept, and the inverted allele is the forward sequence it spells
        s, pk, occ = self.run_extract(walks, 'inv_allowed', allow_inversions=True)
        self.assertEqual((s['status'], s['both_ways']), ({'ok': 1}, 1))
        self.assertEqual(pk[0]['complex'], [])
        self.assertEqual({a['seq'] for a in pk[0]['alleles']}, {'CCC', 'GGG', 'GCC'})

    def test_node_stored_reversed(self):
        # node 7 is stored reversed: every run visits it backward, forward runs and the reverse one alike
        walks = [('CHM13', '0', '>1>2>4>5>6'), ('s1', '1', '>1>2<7>4>5>6'), ('s2', '1', '<6<5<4>7<2<1'),
                 ('s3', '1', '>7<2')]       # anchor-free: runs against the region, so it is flipped
        s, pk, occ = self.run_extract(walks, 'stored')
        self.assertEqual((s['status'], s['both_ways'], s['internal_flipped']), ({'ok': 1}, 0, 1))
        self.assertEqual({a['seq']: a['weight'] for a in pk[0]['alleles']}, {'CCC': 1, 'CAAACC': 2})
        self.assertEqual([(f['side'], f['seq']) for f in pk[0]['fragments']], [('internal', 'CAAA')])
        self.assertEqual([o[2:6] for o in occ if o[0] == 's3#1#chrT#0'], [['internal', '-', '0', '2']])

    def test_anchor_free_walk_orientation(self):
        ref = [('CHM13', '0', '>1>2>4>5>6'), ('s1', '1', '>1>2>7>4>5>6')]
        # with the region: kept as it is
        s, pk, occ = self.run_extract(ref + [('s3', '1', '>7>4')], 'with')
        self.assertEqual((s['status'], s['internal_flipped']), ({'ok': 1}, 0))
        self.assertEqual([o[3] for o in occ if o[0] == 's3#1#chrT#0'], ['+'])
        # against it: flipped, to spell the same fragment as the walk with the region does
        s, pk, occ = self.run_extract(ref + [('s3', '1', '<4<7')], 'against')
        self.assertEqual((s['status'], s['internal_flipped']), ({'ok': 1}, 1))
        self.assertEqual([(f['side'], f['seq']) for f in pk[0]['fragments']], [('internal', 'TTTCC')])
        self.assertEqual([o[3] for o in occ if o[0] == 's3#1#chrT#0'], ['-'])
        # both ways: the region is complex
        s, pk, occ = self.run_extract(ref + [('s3', '1', '>7<4')], 'both')
        self.assertEqual(s['status'], {'complex': 1})
        self.assertEqual(pk[0]['complex'], ['s3#1#chrT#0: an anchor-free walk runs both ways through the region'])
        # both ways, allowed: it takes the orientation more of its nodes agree with, here against the region
        s, pk, occ = self.run_extract(ref + [('s3', '1', '<2>7<4')], 'both_allowed', allow_inversions=True)
        self.assertEqual((s['status'], s['internal_flipped']), ({'ok': 1}, 1))
        self.assertEqual([o[3] for o in occ if o[0] == 's3#1#chrT#0'], ['-'])
        self.assertEqual([(f['side'], f['seq']) for f in pk[0]['fragments']], [('internal', 'CCAAAC')])


# Exit edges: a reference 1..9 with anchors L=2 and R=8 (interior 3..7); 10..14 are off the reference.
EXIT_S = {1: 'GATTACA', 2: 'CCGG', 3: 'ACGTTGCAAGT', 4: 'TTAGGCATC', 5: 'AAAACCCC', 6: 'GGTTGG', 7: 'CATCAT',
          8: 'GGCC', 9: 'TGCATGCA', 10: 'ACACACAGT', 11: 'TTTTGGG', 12: 'CCCAAA', 13: 'AGAGAG', 14: 'GTGTGT'}
EXIT_WALKS = [
    ('CHM13', '0', '>1>2>3>4>5>6>7>8>9'),
    ('s1', '1', '>1>2>3>4>5>6>7>8>9'),
    ('s2', '1', '>1>10>4>5>6>7>8>9'),        # leaves before L, comes back inside through 10: an exit
    ('s3', '1', '<9<11<6<5<4<3<2<1'),        # reverse, leaves after R through 11
    ('s4', '1', '>1>12>9'),                  # misses the whole region: no run
    ('s5', '1', '>9>10'),                    # visits 10 outside every run, so 10 stays
    ('s6', '1', '>5>6>7>8>9'),               # starts inside: a prefix
    ('s7', '1', '>1>2>3>14>1'),              # enters by L and leaves before it: a U-turn
    ('s8', '1', '>13>4'),                    # lies inside, with an interior reference step
    ('s9', '1', '>1>4>5>6>7>8>9'),           # an edge from before L straight into the interior
    ('s10', '1', '>11'),                     # anchor-free, on a node of s3's run
]


def exit_region(S=EXIT_S):
    """The region between 2 and 8 of the exit-edge reference."""
    pos, at = 0, {}
    for n in range(1, 10):
        at[n] = (pos, pos + len(S[n]))
        pos += len(S[n])
    return regions.Region('X', 'chrT', 2, 8, at[2][1], at[8][0], at[8][1] - at[2][0], 'X', 'cen', 'ok', '')


class TestExitEdges(unittest.TestCase):
    def runs(self, walk, kinds):
        toks = [(w[0], w[1:]) for w in walk.replace('>', ' >').replace('<', ' <').split()]
        return list(extract.exit_runs_of_walk(toks, 0, len(toks), kinds))

    def test_runs_of_walk(self):
        # spanning through the anchors, forward and reverse
        self.assertEqual(self.runs('>1>2>3>8>9', ['o', 'L', 'i', 'R', 'O']), [('spanning', '+', 2, 3, False)])
        self.assertEqual(self.runs('<9<8<3<2<1', ['O', 'R', 'i', 'L', 'o']), [('spanning', '-', 2, 3, False)])
        # an exit: leaves before L and comes back inside through a node off the reference
        self.assertEqual(self.runs('>1>10>3>8', ['o', None, 'i', 'R']), [('spanning', '+', 1, 3, True)])
        # off the reference between L and R: the run is the stretch
        self.assertEqual(self.runs('>2>10>8', ['L', None, 'R']), [('spanning', '+', 1, 2, False)])
        # from before L to after R off the reference: no run
        self.assertEqual(self.runs('>1>10>9', ['o', None, 'O']), [])
        # leaving L on its outer side, or arriving at it from there, is outside the region
        self.assertEqual(self.runs('<2>10', ['L', None]), [])
        self.assertEqual(self.runs('>10>2', [None, 'L']), [])
        # entering by L and ending inside: a suffix; starting inside and leaving by R: a prefix
        self.assertEqual(self.runs('>2>10', ['L', None]), [('suffix', '+', 1, 2, False)])
        self.assertEqual(self.runs('>10>3>9', [None, 'i', 'O']), [('prefix', '+', 0, 2, True)])
        # in and out on one side: a U-turn
        self.assertEqual(self.runs('>2>3>10>1', ['L', 'i', None, 'o']), [('uturn', '+', 1, 3, True)])
        # inside throughout: internal, oriented by its interior reference steps
        self.assertEqual(self.runs('<4<10<3', ['i', None, 'i']), [('internal', '-', 0, 3, False)])
        # a node the reference visits inside and outside
        self.assertEqual(self.runs('>2>3>5>8', ['L', 'i', 'x', 'R'])[0][0], 'complex')

    def test_extract(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        g = os.path.join(d, 'g.gfa')
        with open(g, 'w') as f:
            f.write('H\tVN:Z:1.1\n')
            f.writelines('S\t%d\t%s\n' % kv for kv in EXIT_S.items())
            f.writelines('W\t%s\t%s\tchrT\t0\t0\t%s\n' % w for w in EXIT_WALKS)
        prepare.write_ref_index(g, 'CHM13', 'chrT', os.path.join(d, 'chrT.ref.tsv.gz'))
        s = extract.extract(g, [exit_region()], os.path.join(d, 'x'), exit_edges=True, ref_dir=d)
        self.assertEqual((s['status'], s['exit_runs']), ({'ok': 1}, 4))
        pk = next(extract.read_packages(os.path.join(d, 'x', 'packages.jsonl.gz')))
        S = EXIT_S
        alleles = {a['seq']: a['weight'] for a in pk['alleles']}
        self.assertEqual(alleles, {S[3] + S[4] + S[5] + S[6] + S[7]: 2,         # CHM13, s1
                                   S[10] + S[4] + S[5] + S[6] + S[7]: 1,        # s2
                                   S[3] + S[4] + S[5] + S[6] + S[11]: 1,        # s3, forward
                                   S[4] + S[5] + S[6] + S[7]: 1})               # s9
        frags = {(f['side'], f['seq']) for f in pk['fragments']}
        self.assertEqual(frags, {('right', S[5] + S[6] + S[7]), ('left', S[3] + S[14]),
                                 ('internal', S[13] + S[4]), ('internal', S[11])})
        self.assertEqual(pk['n_runs'], {'spanning': 5, 'prefix': 1, 'uturn': 1, 'internal': 2})
        self.assertEqual((pk['exit_edges'], pk['kept_nodes'], pk['kept_bp']), (True, 1, len(S[10])))
        with gzip.open(os.path.join(d, 'x', 'interiors.tsv.gz'), 'rt') as f:
            rows = [l.rstrip('\n').split('\t') for l in f if not l.startswith('#')]
        self.assertEqual(rows, [['X', '3,4,5,6,7,10,11,13,14', '10']])
        with gzip.open(os.path.join(d, 'x', 'occurrences.tsv.gz'), 'rt') as f:
            occ = {l.split('\t')[0]: l.rstrip('\n').split('\t')[2:] for l in f if not l.startswith('#')}
        self.assertEqual(occ['s2#1#chrT#0'], ['spanning', '+', '1', '6', 'a2', '>1', '>8'])
        self.assertEqual(occ['s3#1#chrT#0'], ['spanning', '-', '1', '6', 'a3', '<9', '<2'])
        self.assertEqual(occ['s7#1#chrT#0'][:2] + occ['s7#1#chrT#0'][5:], ['uturn', '+', '>2', '>1'])
        self.assertEqual(occ['s10#1#chrT#0'][:2], ['internal', '+'])
        self.assertNotIn('s4#1#chrT#0', occ)
        self.assertNotIn('s5#1#chrT#0', occ)

    def test_keep_fragments(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        g = os.path.join(d, 'g.gfa')
        with open(g, 'w') as f:
            f.write('H\tVN:Z:1.1\n')
            f.writelines('S\t%d\t%s\n' % kv for kv in EXIT_S.items())
            f.writelines('W\t%s\t%s\tchrT\t0\t0\t%s\n' % w for w in EXIT_WALKS)
        prepare.write_ref_index(g, 'CHM13', 'chrT', os.path.join(d, 'chrT.ref.tsv.gz'))
        with self.assertRaisesRegex(ValueError, 'exit edges'):
            extract.extract(g, [exit_region()], os.path.join(d, 'y'), keep_fragments=True)
        s = extract.extract(g, [exit_region()], os.path.join(d, 'x'), exit_edges=True, ref_dir=d,
                            keep_fragments=True)
        # s6's prefix, s7's U-turn and s8's internal run stay on their nodes; anchor-free s10 is not a run
        self.assertEqual((s['status'], s['exit_runs'], s['kept_fragments']), ({'ok': 1}, 3, 3))
        pk = next(extract.read_packages(os.path.join(d, 'x', 'packages.jsonl.gz')))
        S = EXIT_S
        self.assertEqual(len(pk['alleles']), 4)
        self.assertEqual((pk['fragments'], pk['n_runs']), ([], {'spanning': 5}))
        self.assertEqual(pk['kept_fragment_bp'], sum(len(S[n]) for n in (5, 6, 7, 3, 14, 13, 4)))
        # the alleles' nodes, every one of which a fragment or a walk outside every run visits
        with gzip.open(os.path.join(d, 'x', 'interiors.tsv.gz'), 'rt') as f:
            rows = [l.rstrip('\n').split('\t') for l in f if not l.startswith('#')]
        self.assertEqual(rows, [['X', '3,4,5,6,7,10,11', '3,4,5,6,7,10,11']])
        with gzip.open(os.path.join(d, 'x', 'occurrences.tsv.gz'), 'rt') as f:
            occ = {l.split('\t')[0] for l in f if not l.startswith('#')}
        self.assertEqual(occ, {'CHM13#0#chrT#0', 's1#1#chrT#0', 's2#1#chrT#0', 's3#1#chrT#0', 's9#1#chrT#0'})



# Trimming to the targets: a reference 1..7 with anchors L=2 and R=6, whose interior is two random flanks (3
# and 5) around an array (4) of a repeated unit, the target X; the haplotypes' arrays (8, 10, 11, 12) are
# off the reference.
_rng = random.Random(7)


def _rand(n):
    return ''.join(_rng.choice('ACGT') for _ in range(n))


TRIM_UNIT = _rand(60)


def _flank(n):
    """Random, but such that no k-mer across its border with an array is the array's."""
    return (_rng.choice([b for b in 'ACGT' if b != TRIM_UNIT[0]]) + _rand(n - 2) +
            _rng.choice([b for b in 'ACGT' if b not in (TRIM_UNIT[-1], TRIM_UNIT[0])]))


def _array(copies, mutate=None):
    s = TRIM_UNIT * copies
    if mutate is not None:
        s = s[:mutate] + ('A' if s[mutate] != 'A' else 'C') + s[mutate + 1:]
    return s


TRIM_S = {1: _rand(50), 2: _rand(20), 3: _flank(300), 4: _array(10), 5: _flank(300), 6: _rand(20), 7: _rand(50),
          8: _array(8, 100), 9: _flank(200), 10: _array(12), 11: _array(6), 12: _array(9, 200)}
TRIM_WALKS = [
    ('CHM13', '0', '>1>2>3>4>5>6>7'),
    ('h1', '1', '>1>2>3>8>5>6>7'),      # a shorter array
    ('h2', '1', '>9>10>5>6>7'),         # starts between its array and the anchor: a prefix with a whole array
    ('h3', '1', '>1>2>3>11'),           # ends in its array: a suffix whose core reaches the walk's end
    ('h4', '1', '<7<6<5<12<3<2<1'),     # reverse
    ('h5', '1', '>1>2>3>5>6>7'),        # no array: a spanning run without a core
]
STEP_RE = extract.STEP
FLIP = {'>': '<', '<': '>'}
# tuned to the small graph: every k-mer sampled, in windows of 30
TRIM_TUNING = {'CORE_K': 11, 'CORE_STRIDE': 1, 'CORE_WINDOW': 30, 'CORE_GAP': 60, 'CORE_MIN_STRETCH': 200}


def trim_fixture(d, walks=TRIM_WALKS, seqs=TRIM_S, targets=((4, 'X'),), n_ref=7):
    """The trimming graph in `d`, with its reference index and targets BED; returns (graph, region, BED).
    The reference walks nodes 1 to `n_ref`, the region's anchors are its second and second-last node, and
    each target is a whole reference node."""
    g = os.path.join(d, 'g.gfa')
    edges = set()
    for _, _, walk in walks:
        st = [(o, int(n)) for o, n in STEP_RE.findall(walk)]
        for (oa, a), (ob, b) in zip(st, st[1:]):
            e, r = (a, oa, b, ob), (b, FLIP[ob], a, FLIP[oa])
            edges.add(min(e, r))
    with open(g, 'w') as f:
        f.write('H\tVN:Z:1.1\n')
        f.writelines('S\t%d\t%s\n' % kv for kv in seqs.items())
        f.writelines('L\t%d\t%s\t%d\t%s\t0M\n' % (a, '+' if oa == '>' else '-', b, '+' if ob == '>' else '-')
                     for a, oa, b, ob in sorted(edges))
        f.writelines('W\t%s\t%s\tchrT\t0\t0\t%s\n' % w for w in walks)
    prepare.write_ref_index(g, 'CHM13', 'chrT', os.path.join(d, 'chrT.ref.tsv.gz'))
    pos, at = 0, {}
    for n in range(1, n_ref + 1):
        at[n] = (pos, pos + len(seqs[n]))
        pos += len(seqs[n])
    bed = os.path.join(d, 'targets.bed')
    with open(bed, 'w') as f:
        f.writelines('chrT\t%d\t%d\t%s\n' % (at[n] + (t,)) for n, t in targets)
    r = n_ref - 1
    region = regions.Region('X', 'chrT', 2, r, at[2][1], at[r][0], at[r][1] - at[2][0],
                            ','.join(t for _, t in targets), '.', 'ok', '')
    return g, region, bed


class TestTrim(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        p = mock.patch.multiple(extract, **TRIM_TUNING)
        p.start()
        self.addCleanup(p.stop)

    def test_find_cores(self):
        tk = extract.TargetKmers([('X', TRIM_S[4])], [TRIM_S[3], TRIM_S[5]], k=11)
        self.assertEqual(len(tk.fwd), 60)
        # a few copies of the unit past the target's end leave its k-mers, a flank as dense in them does not
        self.assertEqual(len(extract.TargetKmers([('X', TRIM_S[4])], [_flank(5000) + TRIM_UNIT * 2], k=11).fwd), 60)
        self.assertEqual(len(extract.TargetKmers([('X', TRIM_S[4])], [_array(10)], k=11).fwd), 0)
        kw = dict(stride=1, window=30, gap=60)
        flank = _flank(400)
        seq = flank + _array(8) + flank
        self.assertEqual(extract.find_cores(seq, tk, **kw), [(400, 880, 400, 880, '+')])
        self.assertEqual(extract.find_cores(revcomp(seq), tk, **kw), [(400, 880, 400, 880, '-')])
        self.assertEqual(extract.find_cores(flank, tk, **kw), [])
        # a short copy far from the array is a stretch of its own, and the array has more hits
        seq = _array(2) + _flank(300) + _array(8) + flank
        self.assertEqual(extract.find_cores(seq, tk, **kw), [(420, 900, 420, 900, '+')])
        # with sparser samples the core still ends within a stride of the array
        a, b, _, _, st = extract.find_cores(flank + _array(8) + flank, tk, stride=7, window=60, gap=60)[0]
        self.assertTrue(400 <= a < 407 and 873 < b <= 880 and st == '+', (a, b))

    def test_needs_keep_fragments(self):
        g, region, bed = trim_fixture(self.d)
        with self.assertRaisesRegex(ValueError, 'keep_fragments'):
            extract.extract(g, [region], os.path.join(self.d, 'x'), exit_edges=True, ref_dir=self.d,
                            trim_targets=bed)

    def test_extract(self):
        g, region, bed = trim_fixture(self.d)
        x = os.path.join(self.d, 'x')
        s = extract.extract(g, [region], x, exit_edges=True, ref_dir=self.d, keep_fragments=True,
                            allow_inversions=True, trim_targets=bed, trim_margin=0)
        S = TRIM_S
        self.assertEqual(s['status'], {'ok': 1})
        self.assertEqual(s['trim'], {'cores': 4, 'cores_in_fragments': 1, 'no_core': 1, 'at_end': 1,
                                     'no_orient': 0, 'trimmed_bp': 4 * 600 - 100})
        pk = next(extract.read_packages(os.path.join(x, 'packages.jsonl.gz')))
        # every allele is its array alone, forward; h3's array at its walk's end and h5's run stay
        self.assertEqual({a['seq']: a['weight'] for a in pk['alleles']}, {S[4]: 1, S[8]: 1, S[10]: 1, S[12]: 1})
        self.assertEqual(pk['ref'], 'a1')
        self.assertEqual((pk['n_runs'], pk['kept_fragments']), ({'spanning': 4}, 2))
        self.assertEqual(pk['kept_fragment_bp'], len(S[3] + S[11]) + len(S[3] + S[5]))
        self.assertEqual(pk['trim']['margin'], 0)
        with gzip.open(os.path.join(x, 'interiors.tsv.gz'), 'rt') as f:
            rows = [l.rstrip('\n').split('\t') for l in f if not l.startswith('#')]
        self.assertEqual(rows, [['X', '4,8,10,12', '']])
        with gzip.open(os.path.join(x, 'occurrences.tsv.gz'), 'rt') as f:
            occ = {l.split('\t')[0]: l.rstrip('\n').split('\t')[2:] for l in f if not l.startswith('#')}
        self.assertEqual(occ['CHM13#0#chrT#0'], ['spanning', '+', '3', '4', 'a1', '>3', '>5'])
        self.assertEqual(occ['h2#1#chrT#0'], ['spanning', '+', '1', '2', 'a3', '>9', '>5'])
        self.assertEqual(occ['h4#1#chrT#0'], ['spanning', '-', '3', '4', 'a4', '<5', '<3'])
        self.assertEqual(set(occ), {'CHM13#0#chrT#0', 'h1#1#chrT#0', 'h2#1#chrT#0', 'h4#1#chrT#0'})
        with gzip.open(os.path.join(x, 'cores.tsv.gz'), 'rt') as f:
            cores = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1:] for l in f if not l.startswith('#')}
        self.assertEqual(cores['CHM13#0#chrT#0'], ['X', 'X', 'spanning', 'core', '+', '3', '4', 'chrT',
                                                   '370', '970', '370', '970'])
        self.assertEqual(cores['h3#1#chrT#0'][2:4], ['suffix', 'at_end'])
        self.assertEqual(cores['h5#1#chrT#0'][2:4], ['spanning', 'no_core'])

    def test_several_targets(self):
        # two targets that are pieces of one array (4 and 6), whose k-mers they share, on the reference
        # with a flank between them. a target's own k-mers would be none
        S = {1: _rand(50), 2: _rand(20), 3: _flank(300), 4: _array(10), 5: _flank(500), 6: _array(8),
             7: _flank(300), 8: _rand(20), 9: _rand(50), 10: _array(9, 200), 11: _array(2)}
        walks = [
            ('CHM13', '0', '>1>2>3>4>5>6>7>8>9'),
            ('h1', '1', '>1>2>3>4>5>10>7>8>9'),      # a longer second array
            ('h2', '1', '>1>2>3>4>7>8>9'),           # only the first
            ('h3', '1', '>1>2>3>4>5>11>7>8>9'),      # the first, and a short copy of the unit, which stays
            ('h4', '1', '<9<8<7<10<5<4<3<2<1'),      # reverse
        ]
        tk = extract.TargetKmers([('X', S[4]), ('Y', S[6])], [S[3], S[5], S[7]], k=11)
        self.assertEqual((tk.ids, len(tk.fwd)), (['X', 'Y'], 60))
        kw = dict(stride=1, window=30, gap=60, min_stretch=200)
        seq = S[3] + S[4] + S[5] + S[10] + S[7]
        self.assertEqual(extract.find_cores(seq, tk, **kw),
                         [(300, 900, 300, 900, '+'), (1400, 1940, 1400, 1940, '+')])
        self.assertEqual(extract.find_cores(revcomp(seq), tk, **kw),
                         [(300, 840, 300, 840, '-'), (1340, 1940, 1340, 1940, '-')])
        # only the best stretch, without several targets, and for one that is too short
        one = extract.TargetKmers([('X', S[4])], [S[3], S[5], S[7]], k=11)
        self.assertEqual(extract.find_cores(seq, one, **kw), [(300, 900, 300, 900, '+')])
        self.assertEqual(extract.find_cores(S[3] + S[4] + S[5] + S[11] + S[7], tk, **kw),
                         [(300, 900, 300, 900, '+')])
        g, region, bed = trim_fixture(self.d, walks, S, ((4, 'X'), (6, 'Y')), 9)
        x = os.path.join(self.d, 'x')
        s = extract.extract(g, [region], x, exit_edges=True, ref_dir=self.d, keep_fragments=True,
                            allow_inversions=True, trim_targets=bed, trim_margin=0)
        self.assertEqual(s['status'], {'ok': 1})
        self.assertEqual((s['trim']['cores'], s['trim']['no_core'], s['trim']['at_end']), (8, 0, 0))
        pk = next(extract.read_packages(os.path.join(x, 'packages.jsonl.gz')))
        self.assertEqual({a['seq']: a['weight'] for a in pk['alleles']}, {S[4]: 5, S[6]: 1, S[10]: 2})
        self.assertEqual(pk['n_runs'], {'spanning': 8})
        with gzip.open(os.path.join(x, 'occurrences.tsv.gz'), 'rt') as f:
            occ = [l.rstrip('\n').split('\t') for l in f if not l.startswith('#')]
        self.assertEqual([r[4:6] for r in occ if r[0] == 'CHM13#0#chrT#0'], [['3', '4'], ['5', '6']])
        self.assertEqual([r[4:6] for r in occ if r[0] == 'h2#1#chrT#0'], [['3', '4']])
        with gzip.open(os.path.join(x, 'cores.tsv.gz'), 'rt') as f:
            cores = [l.rstrip('\n').split('\t') for l in f if not l.startswith('#')]
        self.assertEqual([c[2:5] for c in cores if c[0] == 'h3#1#chrT#0'], [['X,Y', 'spanning', 'core']])

    def test_cut_steps(self):
        off = [0, 100, 5000, 5100, 6000]       # steps of 100, 4900, 100 and 900 bp; a core at [5000, 5100)
        self.assertEqual([extract.cut_steps(off, 5000, 5100, m, True) for m in (0, 50, 1000, 3000, 6000)],
                         [2, 2, 2, 1, 0])
        self.assertEqual([extract.cut_steps(off, 5000, 5100, m, False) for m in (0, 50, 500, 900, 6000)],
                         [3, 3, 4, 4, 4])
        # a core that starts inside a long step takes it
        self.assertEqual(extract.cut_steps(off, 4000, 5100, 10, True), 1)

    def test_margin(self):
        # a margin shorter than the flanking nodes leaves them out; h2's array, after its 200 bp first node,
        # is still recorded
        g, region, bed = trim_fixture(self.d)
        S = TRIM_S
        x = os.path.join(self.d, 'x')
        extract.extract(g, [region], x, exit_edges=True, ref_dir=self.d, keep_fragments=True,
                        allow_inversions=True, trim_targets=bed, trim_margin=10)
        pk = next(extract.read_packages(os.path.join(x, 'packages.jsonl.gz')))
        self.assertEqual({a['seq'] for a in pk['alleles']}, {S[4], S[8], S[10], S[12]})
        # a margin longer than them takes their whole nodes, within the run; h2's core then takes the first
        # step of its walk, so it is not recorded
        x = os.path.join(self.d, 'y')
        extract.extract(g, [region], x, exit_edges=True, ref_dir=self.d, keep_fragments=True,
                        allow_inversions=True, trim_targets=bed, trim_margin=400)
        pk = next(extract.read_packages(os.path.join(x, 'packages.jsonl.gz')))
        self.assertEqual({a['seq'] for a in pk['alleles']},
                         {S[3] + S[4] + S[5], S[3] + S[8] + S[5], S[3] + S[12] + S[5]})
        self.assertEqual((pk['trim']['cores_in_fragments'], pk['trim']['at_end']), (0, 2))

if __name__ == '__main__':
    unittest.main()
