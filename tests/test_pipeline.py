"""End to end on a small graph: extract -> realign -> replace, then check that every walk of the output
spells what it spelled before and uses only edges the output has."""
import os
import shutil
import tempfile
import unittest

from pgrealign import extract, realign, regions, replace
from pgrealign.graph import revcomp

UNIT = 'ACGTTGCAAGT'
# 1 = L anchor, 9 = R anchor, flanks 0 and 10. Interior: a tandem array of UNIT copies, one node per copy
# (11, 12, 13, 14), so alleles of 2, 3 and 4 copies are different walks.
S = {0: 'TTTT', 1: 'GATTACA', 11: UNIT, 12: UNIT, 13: UNIT, 14: UNIT, 9: 'CCGG', 10: 'AAAA', 20: 'G'}
WALKS = [
    ('CHM13', '0', '>0>1>11>12>13>9>10'),        # 3 copies
    ('a', '1', '>0>1>11>12>9>10'),               # 2 copies
    ('a', '2', '>0>1>11>12>13>14>9>10'),         # 4 copies
    ('b', '1', '<10<9<14<13<12<11<1<0'),         # 4 copies, reverse
    ('b', '2', '>0>1>11>12>13'),                 # ends inside: left fragment, 3 copies
    ('c', '1', '>13>9>10'),                      # starts inside: right fragment, 1 copy
    ('c', '2', '<10<9<12'),                      # reverse, ends inside: right fragment
    ('d', '1', '>12>13'),                        # internal fragment
    ('e', '1', '>0>1>11>20>12>9>10'),            # 2 copies with a G between
]


def write_gfa(path):
    edges = set()
    for _, _, w in WALKS:
        st = [(x[0], int(x[1:])) for x in w.replace('>', ' >').replace('<', ' <').split()]
        for (oa, a), (ob, b) in zip(st, st[1:]):
            edges.add((a, '+' if oa == '>' else '-', b, '+' if ob == '>' else '-'))
    with open(path, 'w') as f:
        f.write('H\tVN:Z:1.1\tRS:Z:CHM13\n')
        for n, s in S.items():
            f.write('S\t%d\t%s\n' % (n, s))
        for a, ao, b, bo in sorted(edges):
            f.write('L\t%d\t%s\t%d\t%s\t0M\n' % (a, ao, b, bo))
        for sm, hp, w in WALKS:
            f.write('W\t%s\t%s\tchrT\t0\t0\t%s\n' % (sm, hp, w))


def read_gfa(path):
    seq, edges, walks = {}, set(), {}
    with open(path) as f:
        for line in f:
            x = line.rstrip('\n').split('\t')
            if x[0] == 'S':
                seq[int(x[1])] = x[2]
            elif x[0] == 'L':
                edges.add((int(x[1]), x[2], int(x[3]), x[4]))
            elif x[0] == 'W':
                walks['%s#%s' % (x[1], x[2])] = [(s[0], int(s[1:])) for s in
                                                 x[6].replace('>', ' >').replace('<', ' <').split()]
    return seq, edges, walks


def spell(steps, seq):
    return ''.join(seq[n] if o == '>' else revcomp(seq[n]) for o, n in steps)


def has_edge(edges, a, b):
    (oa, na), (ob, nb) = a, b
    fa, fb = ('+' if oa == '>' else '-'), ('+' if ob == '>' else '-')
    flip = {'+': '-', '-': '+'}
    return (na, fa, nb, fb) in edges or (nb, flip[fb], na, flip[fa]) in edges


@unittest.skipUnless(shutil.which('abpoa'), 'abpoa not on PATH')
class TestPipeline(unittest.TestCase):
    def test_roundtrip(self):
        d = tempfile.mkdtemp()
        g = os.path.join(d, 'g.gfa')
        write_gfa(g)
        reg = regions.Region('R1', 'chrT', 1, 9, 11, 11 + 3 * len(UNIT), 0, 'R1', 'VNTR', 'ok', '')
        s = extract.extract(g, [reg], os.path.join(d, 'x'))
        self.assertEqual(s['status'], {'ok': 1})
        pk = next(extract.read_packages(os.path.join(d, 'x', 'packages.jsonl.gz')))
        # c1 (forward, starts inside) and c2 (reverse, ends inside) both spell one copy at the right end
        self.assertEqual(sorted((f['side'], f['weight']) for f in pk['fragments']),
                         [('internal', 1), ('left', 1), ('right', 2)])
        res = realign.realign_package(pk)
        self.assertEqual(res['status'], 'ok', res)
        realign.write_jsonl(os.path.join(d, 'msas.jsonl.gz'), [res])
        seq0, _, walks0 = read_gfa(g)
        for mode in ('reuse', 'dense'):
            out = os.path.join(d, 'out.%s.gfa' % mode)
            st = replace.replace(g, os.path.join(d, 'x'), os.path.join(d, 'msas.jsonl.gz'), out, id_start=1000,
                                 id_mode=mode)
            self.assertEqual(st['regions'], 1)
            seq1, edges1, walks1 = read_gfa(out)
            self.assertEqual(set(walks0), set(walks1))
            for k in walks0:
                self.assertEqual(spell(walks1[k], seq1), spell(walks0[k], seq0), (mode, k))
                for a, b in zip(walks1[k], walks1[k][1:]):
                    self.assertTrue(has_edge(edges1, a, b), (mode, k, a, b))
            used = {n for w in walks1.values() for _, n in w}
            self.assertEqual(used, set(seq1))  # no orphan nodes
            if mode == 'reuse':
                # anchors and flanks keep their ids and sequences
                for n in (0, 1, 9, 10):
                    self.assertEqual(seq1[n], S[n])
            else:
                # ids run densely from the smallest old id, and every edge goes forward except at the
                # reverse-walk flanks, which keep their order
                self.assertEqual(sorted(seq1), list(range(min(S), min(S) + len(seq1))))
                ref = walks1['CHM13#0']
                self.assertEqual([n for _, n in ref], sorted(n for _, n in ref))

if __name__ == '__main__':
    unittest.main()
