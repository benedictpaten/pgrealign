"""replace's pieces, and replace on a three-region graph with hand-made MSAs (no aligner needed)."""
import gzip
import json
import os
import shutil
import tempfile
import unittest

from pgrealign import extract, regions, replace, verify
from pgrealign.graph import revcomp

UNIT = 'ACGTTGCAAGT'
UNIT2 = 'TTAGGCATC'
# Three regions on one chain: R1 (1, 9) and R2 (9, 30) share the anchor 9, and R3 (41, 49) has an allele
# that deletes its whole interior (a2 walks >41>49).
S = {5: 'TTTT', 1: 'GATTACA', 11: UNIT, 12: UNIT, 13: UNIT, 14: UNIT, 9: 'CCGG', 31: UNIT2, 32: UNIT2,
     33: UNIT2, 34: 'A', 30: 'GGCC', 40: 'TGCA', 41: 'CATG', 51: 'ACACACAC', 52: 'AC', 49: 'GGTT', 10: 'AAAA'}
WALKS = [
    ('CHM13', '0', '>5>1>11>12>13>9>31>32>30>40>41>51>49>10'),
    ('a', '1', '>5>1>11>12>9>31>32>33>30>40>41>51>52>49>10'),
    ('a', '2', '>5>1>11>12>13>14>9>31>30>40>41>49>10'),
    ('b', '1', '<10<49<52<51<41<40<30<32<34<31<9<13<12<11<1<5'),
    ('b', '2', '>5>1>11>12>13'),                 # ends inside R1: a left fragment
    ('c', '1', '>13>9>31>32>30>40'),             # starts inside R1: a right fragment; spans R2
    ('c', '2', '<10<49<51'),                     # reverse, ends inside R3: a right fragment
    ('d', '1', '>12>13'),                        # internal fragment of R1
    ('e', '1', '>5>1>11>12>9>31>34>32>30>40>41>51>49>10'),
    ('f', '1', '>32>33>30>40>41>51>52>49>10'),   # starts inside R2
]
REGIONS = [('R1', 1, 9), ('R2', 9, 30), ('R3', 41, 49)]
PATHS = [('gen1', '5+,1+'), ('gen2', '10-,49-,41-,40-')]      # generic paths that avoid the interiors
FLIP = {'+': '-', '-': '+'}


def canonical(edge):
    """An edge in the spelling of the two that sorts first."""
    a, ao, b, bo = edge
    return min(edge, (b, FLIP[bo], a, FLIP[ao]))


def parse_steps(w):
    return [(x[0], int(x[1:])) for x in w.replace('>', ' >').replace('<', ' <').split()]


def write_gfa(path, walks=WALKS, paths=()):
    edges = set()
    for _, _, w in walks:
        st = parse_steps(w)
        for (oa, a), (ob, b) in zip(st, st[1:]):
            edges.add(canonical((a, '+' if oa == '>' else '-', b, '+' if ob == '>' else '-')))
    with open(path, 'w') as f:
        f.write('H\tVN:Z:1.1\tRS:Z:CHM13\n')
        for n, s in S.items():
            f.write('S\t%d\t%s\n' % (n, s))
        for a, ao, b, bo in sorted(edges):
            f.write('L\t%d\t%s\t%d\t%s\t0M\n' % (a, ao, b, bo))
        for name, steps in paths:
            f.write('P\t%s\t%s\t*\n' % (name, steps))
        for sm, hp, w in walks:
            f.write('W\t%s\t%s\tchrT\t0\t0\t%s\n' % (sm, hp, w))


def read_gfa(path):
    """Node sequences, edges, walks {sample#hap: steps} and paths {name: steps} of a GFA."""
    seq, edges, walks, paths = {}, [], {}, {}
    with open(path) as f:
        for line in f:
            x = line.rstrip('\n').split('\t')
            if x[0] == 'S':
                seq[int(x[1])] = x[2]
            elif x[0] == 'L':
                edges.append((int(x[1]), x[2], int(x[3]), x[4]))
            elif x[0] == 'W':
                walks['%s#%s' % (x[1], x[2])] = parse_steps(x[6])
            elif x[0] == 'P':
                paths[x[1]] = [('>' if s[-1] == '+' else '<', int(s[:-1])) for s in x[2].split(',')]
    return seq, edges, walks, paths


def spell(steps, seq):
    return ''.join(seq[n] if o == '>' else revcomp(seq[n]) for o, n in steps)


def has_edge(edges, a, b):
    (oa, na), (ob, nb) = a, b
    fa, fb = ('+' if oa == '>' else '-'), ('+' if ob == '>' else '-')
    return (na, fa, nb, fb) in edges or (nb, FLIP[fb], na, FLIP[fa]) in edges


def naive_msa(pkg):
    """A valid, if poor, MSA of a package: every member left-justified."""
    members = pkg['alleles'] + pkg['fragments']
    width = max(len(m['seq']) for m in members)
    return {'id': pkg['id'], 'status': 'ok', 'rows': {m['id']: m['seq'] + '-' * (width - len(m['seq']))
                                                      for m in members}}


def ref_interval(L, R):
    pos, at = 0, {}
    for _, n in parse_steps(WALKS[0][2]):
        at[n] = (pos, pos + len(S[n]))
        pos += len(S[n])
    return at[L][1], at[R][0]


class Fixture:
    """The three-region graph, extracted, with naive MSAs in one msas file per region."""

    def __init__(self, d, walks=WALKS, paths=()):
        self.dir = d
        self.graph = os.path.join(d, 'g.gfa')
        write_gfa(self.graph, walks, paths)
        regs = [regions.Region(rid, 'chrT', L, R, *ref_interval(L, R), 0, rid, 'VNTR', 'ok', '')
                for rid, L, R in REGIONS]
        self.extract = os.path.join(d, 'x')
        summary = extract.extract(self.graph, regs, self.extract)
        assert summary['status'] == {'ok': 3}, summary
        self.msas = []
        for pkg in extract.read_packages(os.path.join(self.extract, 'packages.jsonl.gz')):
            path = os.path.join(d, 'msas.%s.jsonl.gz' % pkg['id'])
            with gzip.open(path, 'wt') as f:
                f.write(json.dumps(naive_msa(pkg)) + '\n')
            self.msas.append(path)
        self.single_msas = os.path.join(d, 'msas.jsonl.gz')
        with gzip.open(self.single_msas, 'wt') as out:
            for path in self.msas:
                with gzip.open(path, 'rt') as f:
                    out.write(f.read())


class TestPieces(unittest.TestCase):
    PKG = {'id': 'R', 'L': 1, 'R': 9,
           'alleles': [{'id': 'a1', 'seq': 'ACGT'}, {'id': 'a2', 'seq': ''}],
           'fragments': [{'id': 'f1', 'side': 'left', 'seq': 'AC'}, {'id': 'f2', 'side': 'right', 'seq': 'GT'},
                         {'id': 'f3', 'side': 'internal', 'seq': 'CG'}]}
    ROWS = {'a1': 'ACGT', 'a2': '----', 'f1': 'AC--', 'f2': '--GT', 'f3': '-CG-'}

    def test_region_graph(self):
        rg = replace.region_graph(self.PKG, self.ROWS)
        self.assertEqual(rg.seq, {-1: 'A', -2: 'C', -3: 'G', -4: 'T'})
        self.assertEqual(rg.walks, {'a1': [-1, -2, -3, -4], 'a2': [], 'f1': [-1, -2], 'f2': [-3, -4],
                                    'f3': [-2, -3]})
        # L->first for a1 and the left fragment, last->R for a1 and the right fragment, L->R for the empty
        # allele; nothing for the internal fragment
        self.assertEqual({e for e in rg.edges if e[0] > 0 or e[2] > 0},
                         {(1, '+', -1, '+'), (-4, '+', 9, '+'), (1, '+', 9, '+')})
        self.assertEqual({e for e in rg.edges if e[0] < 0 and e[2] < 0},
                         {(-1, '+', -2, '+'), (-2, '+', -3, '+'), (-3, '+', -4, '+')})

    def test_number_dense(self):
        regs = [('R1', 1, 3), ('R2', 5, 1)]
        renumber, new_ids, id_range = replace.number_dense([1, 2, 3, 5, 9], {2, 3}, regs)
        self.assertEqual(renumber, {1: 1, 5: 5, 9: 7})
        self.assertEqual(new_ids, {'R1': [2, 3, 4], 'R2': [6]})
        self.assertEqual(id_range, (1, 7))
        renumber, new_ids, id_range = replace.number_dense([1, 2, 3, 5, 9], {2, 3}, regs, id_start=100)
        self.assertEqual(renumber, {1: 100, 5: 104, 9: 106})
        self.assertEqual(new_ids, {'R1': [101, 102, 103], 'R2': [105]})
        self.assertEqual(id_range, (100, 106))
        with self.assertRaises(ValueError):
            replace.number_dense([1, 2, 3, 9], {2, 3}, regs)        # no left anchor 5

    def test_number_reuse(self):
        new_ids, next_id = replace.number_reuse([('R1', 3), ('R2', 1)], {'R1': [3, 2], 'R2': [8, 7]}, 100)
        self.assertEqual(new_ids, {'R1': [2, 3, 100], 'R2': [7]})
        self.assertEqual(next_id, 101)
        # the ids from id_start run into the region's own freed ids
        with self.assertRaisesRegex(ValueError, 'id 6 would go'):
            replace.number_reuse([('R1', 3)], {'R1': [5, 6]}, 6)
        # two regions' interiors overlap
        with self.assertRaisesRegex(ValueError, 'id 6 would go'):
            replace.number_reuse([('R1', 2), ('R2', 2)], {'R1': [5, 6], 'R2': [6, 7]}, 100)


class TestReplace(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d)

    def edited(self, path, edit):
        """A copy of the GFA `path` with its lines passed through `edit`."""
        with open(path) as f:
            lines = f.readlines()
        out = path + '.edited.gfa'
        with open(out, 'w') as f:
            f.writelines(edit(lines))
        return out

    def assert_fails(self, regex, graph, fx, **kw):
        """replace fails with a ValueError matching `regex`, and leaves no output."""
        out = os.path.join(self.d, 'fails.gfa')
        with self.assertRaisesRegex(ValueError, regex):
            replace.replace(graph, fx.extract, fx.single_msas, out, **kw)
        self.assertFalse(os.path.exists(out))
        self.assertFalse(os.path.exists(out + '.tmp'))

    def check_roundtrip(self, before, after):
        """Every walk and path of `after` spells what it did in `before`, and uses only edges `after`
        has; no node is unused."""
        seq0, _, walks0, paths0 = read_gfa(before)
        seq1, edges1, walks1, paths1 = read_gfa(after)
        edges1 = set(edges1)
        self.assertEqual(set(walks0), set(walks1))
        self.assertEqual(set(paths0), set(paths1))
        for old, new in ((walks0, walks1), (paths0, paths1)):
            for k in old:
                self.assertEqual(spell(new[k], seq1), spell(old[k], seq0), k)
                for a, b in zip(new[k], new[k][1:]):
                    self.assertTrue(has_edge(edges1, a, b), (k, a, b))
        used = {n for w in list(walks1.values()) + list(paths1.values()) for _, n in w}
        self.assertEqual(used, set(seq1))
        for a, _, b, _ in edges1:
            self.assertIn(a, seq1)
            self.assertIn(b, seq1)
        return seq1, walks1

    def assert_edges_once(self, gfa):
        edges = read_gfa(gfa)[1]
        self.assertEqual(len({canonical(e) for e in edges}), len(edges))

    def test_modes(self):
        fx = Fixture(self.d, paths=PATHS)
        for mode, start in (('reuse', 1000), ('dense', None), ('dense', 500)):
            out = os.path.join(self.d, 'out.%s.%s.gfa' % (mode, start))
            st = replace.replace(fx.graph, fx.extract, fx.single_msas, out, id_start=start, id_mode=mode)
            self.assertEqual(st['regions'], 3)
            seq1, walks1 = self.check_roundtrip(fx.graph, out)
            # a2 deletes R3's whole interior: its L->R edge, 41->49, is the input's and is written once
            self.assert_edges_once(out)
            if mode == 'reuse':
                for n in (5, 1, 9, 30, 40, 41, 49, 10):
                    self.assertEqual(seq1[n], S[n])
            else:
                first = min(S) if start is None else start
                self.assertEqual(sorted(seq1), list(range(first, first + len(seq1))))
                self.assertEqual(st['id_range'], [first, first + len(seq1) - 1])
                # the nodes that stay come first, in input order; each region's new nodes take the ids
                # between its left anchor's and the next node's that stays, in topological order
                kept = dict(zip([5, 1, 9, 30, 40, 41, 49, 10], seq1))
                self.assertEqual({n: seq1[k] for n, k in kept.items()}, {n: S[n] for n in kept})
                ref = [n for _, n in walks1['CHM13#0']]
                for _, L, R in REGIONS:
                    a, b = kept[L], min(v for v in kept.values() if v > kept[L])
                    run = ref[ref.index(kept[L]) + 1:ref.index(kept[R])]
                    self.assertTrue(run)
                    self.assertEqual(run, sorted(run))
                    self.assertTrue(all(a < n < b for n in run), (L, a, b, run))

    def test_several_msas_files(self):
        fx = Fixture(self.d)
        one, several = os.path.join(self.d, 'one.gfa'), os.path.join(self.d, 'several.gfa')
        st1 = replace.replace(fx.graph, fx.extract, fx.single_msas, one)
        st2 = replace.replace(fx.graph, fx.extract, fx.msas, several)
        self.assertEqual(st1, st2)
        with open(one) as a, open(several) as b:
            self.assertEqual(a.read(), b.read())
        with self.assertRaisesRegex(ValueError, 'appears twice'):
            replace.replace(fx.graph, fx.extract, [fx.single_msas, fx.msas[1]], several)

    def test_anchor_free_walk_against_the_region(self):
        # d's anchor-free walk runs against R1: extract flips it, and replace writes it back reversed
        walks = [w if w[0] != 'd' else ('d', '1', '<13<12') for w in WALKS]
        fx = Fixture(self.d, walks=walks)
        with gzip.open(os.path.join(fx.extract, 'occurrences.tsv.gz'), 'rt') as f:
            self.assertEqual([l.split('\t')[2:4] for l in f if l.startswith('d#')], [['internal', '-']])
        for mode in ('reuse', 'dense'):
            out = os.path.join(self.d, 'out.%s.gfa' % mode)
            replace.replace(fx.graph, fx.extract, fx.single_msas, out, id_start=1000, id_mode=mode)
            seq1, walks1 = self.check_roundtrip(fx.graph, out)
            self.assertEqual({o for o, _ in walks1['d#1']}, {'<'})

    def test_path_through_a_region(self):
        fx = Fixture(self.d, paths=PATHS + [('gen3', '40+,41+,51+')])
        for mode, start in (('reuse', 1000), ('dense', None)):
            self.assert_fails('path gen3 visits node 51 of region R3', fx.graph, fx, id_start=start,
                              id_mode=mode)

    def test_header_drops_graph_name(self):
        fx = Fixture(self.d)
        g = self.edited(fx.graph, lambda lines: lines[:1] + ['H\tNM:Z:abc\n', 'H\tSG:Z:abc,def\n'] + lines[1:])
        out = os.path.join(self.d, 'out.gfa')
        replace.replace(g, fx.extract, fx.single_msas, out)
        with open(out) as f:
            self.assertEqual([line for line in f if line[0] == 'H'], ['H\tVN:Z:1.1\tRS:Z:CHM13\n'])

    def test_segment_after_walk(self):
        fx = Fixture(self.d)

        def last(t):
            return lambda lines: [x for x in lines if x[0] != t] + [x for x in lines if x[0] == t]
        # Keeping its ids ('reuse'), replace finds a removed node's S line after the walks missing from the
        # first run through it, and another node's, or any L line, when it is read
        late10 = lambda lines: [x for x in lines if not x.startswith('S\t10\t')] + ['S\t10\tAAAA\n']
        reuse = {'id_mode': 'reuse', 'id_start': 1000}
        self.assert_fails('walk CHM13#0#chrT#0 comes before the S line of node 11',
                          self.edited(fx.graph, last('S')), fx, **reuse)
        self.assert_fails('an S line after a W line', self.edited(fx.graph, late10), fx, **reuse)
        self.assert_fails('an L line after a W line', self.edited(fx.graph, last('L')), fx, **reuse)
        # Renumbering ('dense'), it reads a GFA's node ids only up to the first W line, so it misses those
        # S lines, and names them out of order when it needs their ids
        self.assert_fails('3 left anchors are not nodes of the graph, or their S lines come after a W line',
                          self.edited(fx.graph, last('S')), fx)
        self.assert_fails('node 10 has no S line before the first W line', self.edited(fx.graph, late10), fx)
        self.assert_fails('an L line after a W line', self.edited(fx.graph, last('L')), fx)

    def test_old_node_ids_stop_at_the_walks(self):
        # a GFA's node ids are read up to its first W line, which comes after every S line it may have
        fx = Fixture(self.d)
        with open(fx.graph) as f:
            want = sorted(int(x.split('\t')[1]) for x in f if x[0] == 'S')
        self.assertEqual(sorted(replace.old_node_ids(fx.graph)), want)
        g = self.edited(fx.graph, lambda lines: lines + ['S\tnot-an-id\tA\n'])
        self.assertEqual(sorted(replace.old_node_ids(g)), want)

    def test_reuse_id_start_in_use(self):
        fx = Fixture(self.d)
        out = os.path.join(self.d, 'out.gfa')
        st = replace.replace(fx.graph, fx.extract, fx.single_msas, out, id_start=1000, id_mode='reuse')
        self.assertGreater(st['ids_above_start'], 0)        # so ids from id_start are handed out
        # 40 and up: nodes 40, 41 and 49 stay; 10 and up: 11 is freed by R1
        self.assert_fails('id 40, which a node of the graph has', fx.graph, fx, id_start=40, id_mode='reuse')
        self.assert_fails('id 11 would go to new nodes of region R1 and of region R2', fx.graph, fx, id_start=10,
                          id_mode='reuse')

    def test_anchor_in_another_interior(self):
        fx = Fixture(self.d)
        path = os.path.join(fx.extract, 'interiors.tsv.gz')
        with gzip.open(path, 'rt') as f:
            lines = f.readlines()
        with gzip.open(path, 'wt') as f:
            f.writelines(x.rstrip('\n') + ',41\n' if x.startswith('R1\t') else x for x in lines)
        self.assert_fails('anchor 41 of region R3 is in the interior of region R1', fx.graph, fx)

    @unittest.skipUnless(shutil.which('vg'), 'vg not on PATH')
    def test_dense_output_loads(self):
        fx = Fixture(self.d, paths=PATHS)
        out, gbz = os.path.join(self.d, 'out.gfa'), os.path.join(self.d, 'out.gbz')
        replace.replace(fx.graph, fx.extract, fx.single_msas, out)
        replace.to_gbz(out, gbz)
        res = verify.verify(fx.graph, gbz)
        self.assertTrue(res['ok'], res)
        self.assertEqual(res['paths'], len(WALKS) + len(PATHS))


if __name__ == '__main__':
    unittest.main()
