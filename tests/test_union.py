import copy
import gzip
import json
import os
import shutil
import tempfile
import unittest

from pgrealign import realign, union

STATUSES = ('ok', 'invariant', 'complex', 'no_runs')


def pkg(pid, alleles, fragments=(), status='ok', L=1, R=9):
    return {'id': pid, 'L': L, 'R': R, 'status': status,
            'alleles': [{'id': a, 'weight': w, 'seq': s} for a, s, w in alleles],
            'fragments': [{'id': f, 'weight': w, 'side': side, 'seq': s} for f, s, side, w in fragments]}


def write(path, records):
    with gzip.open(path, 'wt') as f:
        for r in records:
            f.write(json.dumps(r) + '\n')


def read(path):
    with gzip.open(path, 'rt') as f:
        return [json.loads(l) for l in f]


PANEL = pkg('R', [('a1', 'TACGT', 3), ('a2', 'ACGGT', 1)], [('f1', 'AC', 'left', 2)])
TARGET = pkg('R', [('a1', 'ACGGT', 1), ('a2', 'ACCT', 2)],
             [('f1', 'AC', 'left', 1), ('f2', 'AC', 'right', 1), ('f3', 'GG', 'internal', 1)])
# an MSA of the union of PANEL and TARGET
ROWS = {'a1': 'TAC-GT', 'a2': '-ACGGT', 'ta2': '-ACC-T', 'f1': '-AC---', 'tf2': '----AC', 'tf3': '---GG-'}


class TestUnionPackage(unittest.TestCase):
    def test_members(self):
        p, q = copy.deepcopy(PANEL), copy.deepcopy(TARGET)
        u = union.union_package(p, q)
        self.assertEqual((p, q), (PANEL, TARGET))                     # the inputs are not changed
        # a target allele the panel has is the panel's, at the panel's weight; one it lacks joins at 0
        self.assertEqual(u['alleles'], PANEL['alleles'] + [{'id': 'ta2', 'weight': 0, 'seq': 'ACCT'}])
        # fragments match by side and sequence: the right-side AC is not the left-side one
        self.assertEqual(u['fragments'], PANEL['fragments'] + [
            {'id': 'tf2', 'weight': 0, 'side': 'right', 'seq': 'AC'},
            {'id': 'tf3', 'weight': 0, 'side': 'internal', 'seq': 'GG'}])
        self.assertEqual(u['target_map'], {'a1': 'a2', 'a2': 'ta2', 'f1': 'f1', 'f2': 'tf2', 'f3': 'tf3'})
        self.assertEqual((u['status'], u['target_status'], u['L'], u['R']), ('ok', 'ok', 1, 9))

    def test_status(self):
        for ps in STATUSES:
            for qs in STATUSES:
                u = union.union_package(dict(PANEL, status=ps), dict(TARGET, status=qs))
                self.assertEqual(u['status'], 'complex' if qs == 'complex' else ps, (ps, qs))
                self.assertEqual(u['target_status'], qs)

    def test_invariant_panel_stays_invariant(self):
        p = pkg('R', [('a1', 'ACGT', 4)], status='invariant')
        u = union.union_package(p, pkg('R', [('a1', 'ACGT', 1), ('a2', 'AGGT', 1)]))
        self.assertEqual((u['status'], len(u['alleles'])), ('invariant', 2))

    def test_anchors_must_agree(self):
        with self.assertRaises(ValueError):
            union.union_package(PANEL, dict(TARGET, R=8))


class TestProjectRows(unittest.TestCase):
    def test_rows(self):
        self.assertEqual(union.project_rows(ROWS, {'x': 'a1', 'y': 'a2'}), {'x': 'TAC-GT', 'y': '-ACGGT'})
        # columns 0 and 4 are gaps in both of these rows
        self.assertEqual(union.project_rows(ROWS, {'x': 'ta2', 'y': 'f1'}), {'x': 'ACCT', 'y': 'AC--'})
        self.assertEqual(union.project_rows(ROWS, {}), {})


class TestFiles(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def path(self, name):
        return os.path.join(self.d, name)

    def test_union_and_project(self):
        write(self.path('base.jsonl.gz'), [PANEL, pkg('B', [('a1', 'A', 1)])])
        write(self.path('add.jsonl.gz'), [pkg('C', [('a1', 'A', 1)]), TARGET])
        s = union.union_packages(self.path('base.jsonl.gz'), self.path('add.jsonl.gz'), self.path('u.jsonl.gz'))
        self.assertEqual(s, {'regions': 1, 'base_only': 1, 'add_only': 1, 'target_only_alleles': 1,
                             'target_only_fragments': 2})
        self.assertEqual(read(self.path('u.jsonl.gz')), [union.union_package(PANEL, TARGET)])
        msas = [{'id': 'R', 'status': 'ok', 'centre': 'a1', 'columns': 6, 'rows': ROWS},
                {'id': 'X', 'status': 'ok', 'rows': {'a1': 'A'}},          # not in the union
                {'id': 'R2', 'status': 'too_big'}]
        write(self.path('msas.jsonl.gz'), msas)
        s = union.project(self.path('u.jsonl.gz'), self.path('msas.jsonl.gz'), self.path('p.jsonl.gz'))
        self.assertEqual(s, {'projected': 1, 'skipped': 2})
        (m,) = read(self.path('p.jsonl.gz'))
        # the target's a1 is the panel's a2, so the panel's a1 and its column 0 are left out
        self.assertEqual(m['rows'], {'a1': 'ACGGT', 'a2': 'ACC-T', 'f1': 'AC---', 'f2': '---AC', 'f3': '--GG-'})
        self.assertEqual((m['columns'], m['projected_from'], m['centre']), (5, 6, 'a1'))

    @unittest.skipUnless(shutil.which('abpoa'), 'abpoa not on PATH')
    def test_realigned_union_spells_the_target(self):
        """Union, realign and project: every target member's projected row spells its sequence."""
        p = pkg('R', [('a1', 'ACGTACGTACGT', 3), ('a2', 'ACGTACGT', 1)], [('f1', 'ACGTA', 'left', 1)])
        q = pkg('R', [('a1', 'ACGTACGTACGTACGT', 1), ('a2', 'ACGTACGT', 1)],
                [('f1', 'CGTACGT', 'right', 1)])
        write(self.path('base.jsonl.gz'), [p])
        write(self.path('add.jsonl.gz'), [q])
        union.union_packages(self.path('base.jsonl.gz'), self.path('add.jsonl.gz'), self.path('u.jsonl.gz'))
        s = realign.run(self.path('u.jsonl.gz'), self.path('msas.jsonl.gz'))
        self.assertEqual(s['realigned'], {'ok': 1})
        union.project(self.path('u.jsonl.gz'), self.path('msas.jsonl.gz'), self.path('p.jsonl.gz'))
        (m,) = read(self.path('p.jsonl.gz'))
        want = {x['id']: x['seq'] for x in q['alleles'] + q['fragments']}
        self.assertEqual({k: r.replace('-', '') for k, r in m['rows'].items()}, want)
        self.assertEqual({len(r) for r in m['rows'].values()}, {m['columns']})


if __name__ == '__main__':
    unittest.main()
