import gzip
import json
import os
import shutil
import tempfile
import unittest

from pgrealign import realign, star


def pkg(pid, alleles, fragments=()):
    return {'id': pid, 'status': 'ok',
            'alleles': [{'id': a, 'weight': w, 'seq': s} for a, s, w in alleles],
            'fragments': [{'id': f, 'side': side, 'seq': s} for f, s, side in fragments]}


@unittest.skipUnless(shutil.which('abpoa'), 'abpoa not on PATH')
class TestRealignPackage(unittest.TestCase):
    def test_ok_fields(self):
        r = realign.realign_package(pkg('R', [('a1', 'ACGTACGT', 1), ('a2', 'ACGTTCGT', 1)],
                                        [('f1', 'GGGG', 'internal')]))
        self.assertEqual((r['status'], r['fragments_aligned'], r['fragments_exact']), ('ok', 1, 0))
        self.assertEqual(r['rows']['f1'].replace('-', ''), 'GGGG')
        self.assertEqual(r['columns'], len(r['rows']['a1']))
        self.assertEqual(r['method'], 'medoid_star')

    def test_statuses(self):
        two = [('a1', 'ACGTACGT', 1), ('a2', 'ACGTTCGT', 1)]
        self.assertEqual(realign.realign_package(pkg('R', two), max_cells=10)['status'], 'too_big')
        long_two = [('a1', 'ACGTACGTAC' * 40, 2), ('a2', 'ACGTTCGTAC' * 40, 1)]   # 400 bp: 160,000 cells unbanded
        self.assertEqual(realign.realign_package(pkg('R', long_two), max_cells=5000)['status'], 'too_big')
        r = realign.realign_package(pkg('R', long_two), max_cells=50000)   # banded: about 11,200 cells
        self.assertEqual((r['status'], r['banded_pairs']), ('ok', 1))
        self.assertEqual(realign.realign_package(pkg('R', long_two))['banded_pairs'], 0)
        self.assertEqual(realign.realign_package(pkg('R', []))['status'], 'no_alleles')
        self.assertEqual(realign.realign_package(pkg('R', two), timeout=0)['status'], 'timeout')
        r = realign.realign_package(pkg('R', two), abpoa='false')
        self.assertEqual(r['status'], 'aligner_failed', r)


@unittest.skipUnless(shutil.which('abpoa'), 'abpoa not on PATH')
class TestRun(unittest.TestCase):
    def test_run_and_resume(self):
        d = tempfile.mkdtemp()
        try:
            pk = os.path.join(d, 'packages.jsonl.gz')
            pkgs = [pkg('R%d' % i, [('a1', 'ACGTACGT' * (i + 1), 1), ('a2', 'ACGTTCGT' * (i + 1), 1)])
                    for i in range(4)] + [pkg('S', [('a1', 'ACGT', 2)])]
            with gzip.open(pk, 'wt') as f:
                for p in pkgs:
                    f.write(json.dumps(p) + '\n')
            out = os.path.join(d, 'msas.jsonl.gz')
            s = realign.run(pk, out, jobs=2)
            self.assertEqual(s['realigned'], {'ok': 4})
            self.assertEqual(s['skipped'], {'not_selected': 1})
            first = [json.loads(l) for l in gzip.open(out, 'rt')]
            self.assertEqual([r['id'] for r in first], ['R0', 'R1', 'R2', 'R3'])
            # an interrupted run: two regions are in the partial file, and a third was cut off mid-write
            with open(out + '.partial', 'wb') as f:
                for r in first[:2]:
                    f.write(gzip.compress((json.dumps(r) + '\n').encode()))
                f.write(gzip.compress((json.dumps(first[2]) + '\n').encode())[:30])
            os.remove(out)
            s = realign.run(pk, out, jobs=1)
            self.assertEqual((s['resumed'], s['realigned']), (2, {'ok': 4}))
            # a second interruption after the first: the partial file was rewritten, so it reads whole
            with open(out + '.partial', 'wb') as f:
                f.write(gzip.compress((json.dumps(first[0]) + '\n').encode())[:30])
            os.remove(out)
            s = realign.run(pk, out, jobs=1)
            self.assertEqual((s['resumed'], s['realigned']), (0, {'ok': 4}))
            again = [json.loads(l) for l in gzip.open(out, 'rt')]
            self.assertEqual([r['rows'] for r in again], [r['rows'] for r in first])
            self.assertFalse(os.path.exists(out + '.partial'))
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
