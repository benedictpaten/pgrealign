import collections
import os
import shutil
import tempfile
import unittest

from pgrealign import regions

Target = collections.namedtuple('Target', 'id start end cls')

# a reference of nodes 10..19, each 10 bp, at 0, 10, ..., 90
REF = [(n, '+', (n - 10) * 10, (n - 10) * 10 + 10) for n in range(10, 20)]


def anchors(pairs, ref=REF):
    return regions.Anchors(regions.chains_from_pairs(pairs, ref))


def rows(rs):
    return [(r.id, r.L, r.R, r.ref_start, r.ref_end, r.span, r.targets, r.status, r.note) for r in rs]


class FixedAnchors:
    """Anchors whose answers are given, {(a0, b0): (L, R)}: for merges no small chain can show."""

    def __init__(self, table):
        self.table = table

    def query(self, a0, b0, max_span):
        if (a0, b0) not in self.table:
            raise regions.AnchorError('no chain encloses %d-%d within %d bp' % (a0, b0, max_span))
        return self.table[(a0, b0)]


class TestChains(unittest.TestCase):
    def test_chains(self):
        # snarls that share a boundary node form one chain; the others another
        chains = regions.chains_from_pairs([(10, 12), (12, 14), (15, 17), (17, 19)], REF)
        self.assertEqual(sorted([n for _, _, n in c] for c in chains), [[10, 12, 14], [15, 17, 19]])
        # a boundary node the reference visits reversed or twice is not used
        ref = [(12, '-', 20, 30) if n == 12 else x for x, n in zip(REF, range(10, 20))] + [(14, '+', 100, 110)]
        chains = regions.chains_from_pairs([(10, 12), (12, 14), (14, 16)], ref)
        self.assertEqual([[n for _, _, n in c] for c in chains], [[10, 16]])

    def test_query(self):
        # snarls 10-12, 12-14, 14-16, 16-19
        A = anchors([(10, 12), (12, 14), (14, 16), (16, 19)])
        self.assertEqual(A.query(31, 39, 1000), ((20, 30, 12), (40, 50, 14)))
        # an interval overlapping boundary node 12 must take the anchor before it
        self.assertEqual(A.query(25, 35, 1000), ((0, 10, 10), (40, 50, 14)))
        # an interval reaching into the first boundary node is clamped, and anchored on it
        self.assertEqual(A.query(0, 5, 1000), ((0, 10, 10), (20, 30, 12)))
        # the tightest pair over every chain wins: a second chain, 11-13-15, is tighter around 41-49
        A = anchors([(10, 12), (12, 14), (14, 16), (16, 19), (11, 13), (13, 15)])
        self.assertEqual(A.query(31, 39, 1000), ((20, 30, 12), (40, 50, 14)))
        self.assertEqual(A.query(41, 49, 1000), ((30, 40, 13), (50, 60, 15)))
        self.assertRaises(regions.AnchorError, A.query, 41, 49, 29)


class TestAssign(unittest.TestCase):
    def test_shared_boundary_node(self):
        # targets in snarls 12-14 and 14-16 share anchor 14 only, so they stay separate
        A = anchors([(10, 12), (12, 14), (14, 16), (16, 19)])
        rs = regions.assign([Target('a', 32, 33, 'STR'), Target('b', 52, 53, 'STR')], A, 'chrT', pad=0)
        self.assertEqual(rows(rs), [('a', 12, 14, 30, 40, 30, 'a', 'ok', ''),
                                    ('b', 14, 16, 50, 60, 30, 'b', 'ok', '')])

    def test_overlapping_interiors_merge(self):
        # c crosses node 14, so needs anchors 12 and 16, and overlaps both a and b
        A = anchors([(10, 12), (12, 14), (14, 16), (16, 19)])
        ts = [Target('a', 32, 33, 'STR'), Target('c', 38, 55, 'VNTR'), Target('b', 52, 53, 'STR')]
        rs = regions.assign(ts, A, 'chrT', pad=0)
        self.assertEqual(rows(rs), [('a', 12, 16, 30, 60, 50, 'a,c,b', 'ok', 'merged 3 targets')])
        self.assertEqual(rs[0].classes, 'STR,VNTR')

    def test_anchor_inside_interior_merges(self):
        # every node a boundary node: a (in node 12) is 11-13, b (in node 13) is 12-14. Their interiors,
        # nodes 12 and 13, do not overlap, but each one's anchor is the other's interior, and b's left
        # anchor starts before a's right anchor.
        A = anchors([(n, n + 1) for n in range(10, 19)])
        ts = [Target('a', 22, 23, 'STR'), Target('b', 32, 33, 'STR')]
        self.assertEqual([(r.L, r.R) for r in regions.assign(ts[:1], A, 'chrT', pad=0)], [(11, 13)])
        self.assertEqual([(r.L, r.R) for r in regions.assign(ts[1:], A, 'chrT', pad=0)], [(12, 14)])
        rs = regions.assign(ts, A, 'chrT', pad=0)
        self.assertEqual(rows(rs), [('a', 11, 14, 20, 40, 40, 'a,b', 'ok', 'merged 2 targets')])
        # each fits in 30 bp but the merged region does not: too_long, with the outermost anchors
        rs = regions.assign(ts, A, 'chrT', pad=0, max_span=30)
        self.assertEqual(rows(rs), [('a', 11, 14, 20, 40, 40, 'a,b', 'too_long',
                                     'no chain encloses 22-33 within 30 bp')])

    def test_merging_repeats(self):
        # q and r merge, and re-anchored as one they start before p's right anchor, so p joins them on
        # the next pass
        A = FixedAnchors({(100, 110): ((90, 100, 1), (110, 120, 2)),
                          (130, 140): ((120, 130, 3), (140, 150, 4)),
                          (135, 160): ((130, 135, 5), (160, 170, 6)),
                          (130, 160): ((105, 108, 7), (160, 170, 6)),
                          (100, 160): ((90, 100, 1), (160, 170, 6))})
        ts = [Target('p', 100, 110, 'STR'), Target('q', 130, 140, 'STR'), Target('r', 135, 160, 'STR')]
        rs = regions.assign(ts, A, 'chrT', pad=0)
        self.assertEqual(rows(rs), [('p', 1, 6, 100, 160, 80, 'p,q,r', 'ok', 'merged 3 targets')])

    def test_too_long_target(self):
        # two chains, 10-12-14 and 15-17-19: nothing encloses a target across nodes 14 and 15
        A = anchors([(10, 12), (12, 14), (15, 17), (17, 19)])
        rs = regions.assign([Target('e', 45, 55, 'STR')], A, 'chrT', pad=0)
        self.assertEqual(rows(rs), [('e', '.', '.', 45, 55, 10, 'e', 'too_long',
                                     'no chain encloses 45-55 within 100000 bp')])
        # a target whose only enclosing pair is longer than max_span keeps its own interval too, and
        # a region that is anchored can overlap it
        A = anchors([(10, 12), (12, 14), (14, 16), (16, 19)])
        rs = regions.assign([Target('d', 25, 85, 'STR'), Target('a', 32, 33, 'STR')], A, 'chrT', pad=0,
                            max_span=30)
        self.assertEqual(rows(rs), [('d', '.', '.', 25, 85, 60, 'd', 'too_long',
                                     'no chain encloses 25-85 within 30 bp'),
                                    ('a', 12, 14, 30, 40, 30, 'a', 'ok', '')])

    def test_pad(self):
        # the padded target is anchored; its row still has the interior between the anchors
        A = anchors([(10, 12), (12, 14), (14, 16), (16, 19)])
        rs = regions.assign([Target('a', 32, 33, 'STR')], A, 'chrT', pad=10)
        self.assertEqual(rows(rs), [('a', 10, 16, 10, 60, 70, 'a', 'ok', '')])


class TestFiles(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)

    def test_round_trip(self):
        A = anchors([(10, 12), (12, 14), (14, 16), (16, 19)])
        rs = regions.assign([Target('d', 25, 85, 'STR'), Target('a', 32, 33, 'STR')], A, 'chrT', pad=0,
                            max_span=30)
        path = os.path.join(self.d, 'regions.tsv.gz')
        regions.write_regions(rs, path)
        back = list(regions.read_regions(path))
        self.assertEqual(back[1], rs[1])
        self.assertEqual((back[0].L, back[0].R), (None, None))
        self.assertEqual(back[0][4:], rs[0][4:])

    def test_build_without_targets(self):
        # a contig with no targets is written empty, without its snarls or reference index
        bed = os.path.join(self.d, 'targets.bed')
        with open(bed, 'w') as f:
            f.write('chrA\t100\t200\tTR1\n')
        out = os.path.join(self.d, 'regions.tsv.gz')
        missing = os.path.join(self.d, 'missing')
        s = regions.build(bed, 'chrB', missing + '.snarls.pb', missing + '.ref.tsv.gz', out, vg=missing)
        self.assertEqual(s, {'targets': 0, 'chains': None, 'regions': 0, 'status': {}, 'merged_targets': 0,
                             'span_bp_ok': 0})
        self.assertEqual(list(regions.read_regions(out)), [])


if __name__ == '__main__':
    unittest.main()
