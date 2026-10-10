"""`verify` and `prepare`, the stages that compare and cut whole GBZs, on small fixture graphs."""
import gzip
import os
import shutil
import subprocess
import tempfile
import unittest

from pgrealign import graph, prepare, verify


def write_gfa(path, segments, links, walks, ref='REF'):
    """`walks` is [(sample, hap, contig, [(orient, node)])]."""
    with open(path, 'w') as f:
        f.write('H\tVN:Z:1.1\tRS:Z:%s\n' % ref)
        for n, s in segments.items():
            f.write('S\t%d\t%s\n' % (n, s))
        for a, ao, b, bo in links:
            f.write('L\t%d\t%s\t%d\t%s\t0M\n' % (a, ao, b, bo))
        for sample, hap, contig, steps in walks:
            length = sum(len(segments[n]) for _, n in steps)
            f.write('W\t%s\t%s\t%s\t0\t%d\t%s\n' % (sample, hap, contig, length,
                                                     ''.join('%s%d' % s for s in steps)))


def to_gbz(gfa):
    gbz = gfa[:-4] + '.gbz'
    subprocess.run(['vg', 'gbwt', '-G', gfa, '--gbz-format', '-g', gbz], check=True, capture_output=True)
    return gbz


def fwd(*nodes):
    return [('>', n) for n in nodes]


def rev(*nodes):
    return [('<', n) for n in reversed(nodes)]


# Before: 1 ACGT -> {2 GG} -> 3 T -> 4 CA. After: the same paths, with every node renumbered and the
# first one split, as replace's dense numbering leaves them.
BEFORE = ({1: 'ACGT', 2: 'GG', 3: 'T', 4: 'CA'},
          [(1, '+', 2, '+'), (2, '+', 3, '+'), (1, '+', 3, '+'), (3, '+', 4, '+')],
          [('REF', 0, 'chrA', fwd(1, 2, 3, 4)), ('s1', 1, 'chrA', fwd(1, 3, 4)), ('s1', 2, 'chrA', rev(1, 2, 3, 4))])
AFTER = ({11: 'AC', 12: 'GT', 13: 'GG', 14: 'T', 15: 'CA'},
         [(11, '+', 12, '+'), (12, '+', 13, '+'), (13, '+', 14, '+'), (12, '+', 14, '+'), (14, '+', 15, '+')],
         [('REF', 0, 'chrA', fwd(11, 12, 13, 14, 15)), ('s1', 1, 'chrA', fwd(11, 12, 14, 15)),
          ('s1', 2, 'chrA', rev(11, 12, 13, 14, 15))])


@unittest.skipUnless(shutil.which('vg'), 'vg not on PATH')
class TestVerify(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def graph(self, name, segments, links, walks, gbz=True):
        gfa = os.path.join(self.d, name + '.gfa')
        write_gfa(gfa, segments, links, walks)
        return to_gbz(gfa) if gbz else gfa

    def test_renumbered(self):
        for gbz in (True, False):
            r = verify.verify(self.graph('before', *BEFORE, gbz=gbz), self.graph('after', *AFTER, gbz=gbz))
            self.assertTrue(r['ok'], r)
            self.assertEqual((r['paths'], r['bases'], r['metadata_identical']), (3, 25, True))
            self.assertEqual((r['paths_missing'], r['paths_extra'], r['paths_differ']), (0, 0, 0))
            self.assertEqual((r['before']['nodes'], r['after']['nodes']), (4, 5))

    def test_a_path_spells_differently(self):
        segments, links, walks = AFTER
        walks = walks[:1] + [('s1', 1, 'chrA', fwd(11, 12, 13, 14, 15))] + walks[2:]
        r = verify.verify(self.graph('before', *BEFORE), self.graph('after', segments, links, walks))
        self.assertFalse(r['ok'])
        self.assertTrue(r['metadata_identical'])
        self.assertEqual((r['paths_differ'], r['paths_differ_names']), (1, ['s1#1#chrA#0']))

    def test_counts_past_the_names_shown(self):
        segments, links, walks = BEFORE
        many = [('x%02d' % i, 1, 'chrA', fwd(1, 3, 4)) for i in range(verify.SHOW + 5)]
        before = self.graph('before', segments, links, walks + many)
        after = self.graph('after', *AFTER[:2], AFTER[2] + [('y', 1, 'chrA', fwd(11, 12, 14, 15))])
        r = verify.verify(before, after)
        self.assertFalse(r['ok'] or r['metadata_identical'])
        self.assertEqual((r['paths'], r['paths_missing'], r['paths_extra'], r['paths_differ']),
                         (verify.SHOW + 8, verify.SHOW + 5, 1, 0))
        self.assertEqual(r['paths_missing_names'], ['x%02d#1#chrA#0' % i for i in range(verify.SHOW)])
        self.assertEqual(r['paths_extra_names'], ['y#1#chrA#0'])


# A genome of two contigs. chrA's reference runs through node 4 in reverse, and the gref_CHM13 cover has
# a fragment on each contig.
GENOME = ({1: 'ACGTACGTAA', 2: 'G', 3: 'T', 4: 'GGTTAACCGG', 11: 'TTTTGGGGCC', 12: 'A', 13: 'C', 14: 'GGAATTCCAA'},
          [(1, '+', 2, '+'), (1, '+', 3, '+'), (2, '+', 4, '-'), (3, '+', 4, '-'),
           (11, '+', 12, '+'), (11, '+', 13, '+'), (12, '+', 14, '+'), (13, '+', 14, '+')],
          [('CHM13', 0, 'chrA', [('>', 1), ('>', 2), ('<', 4)]), ('gref_CHM13', 0, 'chrA', fwd(1, 3)),
           ('h1', 1, 'chrA', [('>', 1), ('>', 3), ('<', 4)]), ('h1', 2, 'chrA', [('>', 1), ('>', 2), ('<', 4)]),
           ('CHM13', 0, 'chrB', fwd(11, 12, 14)), ('gref_CHM13', 0, 'chrB', fwd(13, 14)),
           ('h1', 1, 'chrB', fwd(11, 13, 14))])
CHR_A_INDEX = [(1, '+', 0, 10), (2, '+', 10, 11), (4, '-', 11, 21)]


def path_names(gbz):
    out = subprocess.run(['vg', 'paths', '-L', '-x', gbz], check=True, capture_output=True, text=True).stdout
    return sorted(out.split())


def gbz_reference_samples(gbz):
    out = subprocess.run(['vg', 'gbwt', '-Z', gbz, '--tags'], check=True, capture_output=True, text=True).stdout
    return dict(l.split('\t', 1) for l in out.splitlines())['reference_samples'].split()


class TestRefIndex(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def index(self, lines):
        gfa, out = os.path.join(self.d, 'g.gfa.gz'), os.path.join(self.d, 'ref.tsv.gz')
        with gzip.open(gfa, 'wt') as f:
            f.write('H\tVN:Z:1.1\nS\t1\tACGTACGTAA\nS\t2\tG\nS\t4\tGGTTAACCGG\n' + ''.join(l + '\n' for l in lines))
        prepare.write_ref_index(gfa, 'CHM13', 'chrA', out)
        return prepare.read_ref_index(out)

    def test_w_and_p_lines(self):
        w = 'W\tCHM13\t0\tchrA\t100\t121\t>1>2<4'
        self.assertEqual(self.index(['W\th1\t1\tchrA\t0\t21\t>1>2<4', w]),
                         [(n, o, s + 100, e + 100) for n, o, s, e in CHR_A_INDEX])
        self.assertEqual(self.index(['P\tCHM13#0#chrA\t1+,2+,4-\t*']), CHR_A_INDEX)
        self.assertEqual(self.index(['P\tCHM13#chrA\t1+,2+,4-']), CHR_A_INDEX)
        for lines in ([], [w, w.replace('\t100\t121\t', '\t200\t221\t')], ['P\tCHM13#0#chrB\t1+']):
            with self.assertRaises(RuntimeError):
                self.index(lines)
        # a walk is named by its sample and contig alone: CHM13's walk on another contig, and another sample's
        # whose name starts like CHM13's, are passed over
        self.assertEqual(self.index(['W\tCHM13\t0\tchrAB\t0\t21\t>1>2<4', 'W\tCHM13x\t0\tchrA\t0\t21\t>1>2<4', w]),
                         [(n, o, s + 100, e + 100) for n, o, s, e in CHR_A_INDEX])
        self.assertEqual(graph.walk_name(w + '\n'), ('CHM13', '0', 'chrA'))

    def test_the_reference_cannot_be_dropped(self):
        with self.assertRaises(ValueError):
            prepare.prepare('genome.gbz', 'chrA', self.d, drop_samples=('gref_CHM13', 'CHM13'))


@unittest.skipUnless(shutil.which('vg'), 'vg not on PATH')
class TestPrepare(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        gfa = os.path.join(self.d, 'genome.gfa')
        write_gfa(gfa, *GENOME, ref='CHM13 gref_CHM13')
        self.gbz = to_gbz(gfa)

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def test_default_drops_the_gref_cover(self):
        out_dir = os.path.join(self.d, 'out')
        out = prepare.prepare(self.gbz, 'chrA', out_dir)
        self.assertEqual(out, {'gbz': os.path.join(out_dir, 'chrA.gbz'),
                               'snarls': os.path.join(out_dir, 'chrA.snarls.pb'),
                               'ref': os.path.join(out_dir, 'chrA.ref.tsv.gz')})
        self.assertEqual(path_names(out['gbz']), ['CHM13#0#chrA', 'h1#1#chrA#0', 'h1#2#chrA#0'])
        self.assertEqual(gbz_reference_samples(out['gbz']), ['CHM13'])
        self.assertGreater(os.path.getsize(out['snarls']), 0)
        self.assertEqual(prepare.read_ref_index(out['ref']), CHR_A_INDEX)
        self.assertEqual(sorted(os.listdir(out_dir)), ['chrA.gbz', 'chrA.ref.tsv.gz', 'chrA.snarls.pb'])
        # a second run finds every output and rewrites none
        mtimes = {k: os.path.getmtime(v) for k, v in out.items()}
        self.assertEqual(prepare.prepare(self.gbz, 'chrA', out_dir), out)
        self.assertEqual({k: os.path.getmtime(v) for k, v in out.items()}, mtimes)

    def test_drop_nothing(self):
        out = prepare.prepare(self.gbz, 'chrB', os.path.join(self.d, 'out'), drop_samples=())
        self.assertEqual(path_names(out['gbz']), ['CHM13#0#chrB', 'gref_CHM13#0#chrB', 'h1#1#chrB#0'])
        self.assertEqual(gbz_reference_samples(out['gbz']), ['CHM13', 'gref_CHM13'])
        self.assertEqual(prepare.read_ref_index(out['ref']), [(11, '+', 0, 10), (12, '+', 10, 11), (14, '+', 11, 21)])


if __name__ == '__main__':
    unittest.main()
