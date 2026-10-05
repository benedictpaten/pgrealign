import os
import shutil
import subprocess
import tempfile
import unittest

from pgrealign import merge

CONTIG = ('H\tVN:Z:1.1\tRS:Z:REF\n'
          'S\t{a}\tACGT\nS\t{b}\tGG\nS\t{c}\tT\n'
          'L\t{a}\t+\t{b}\t+\t0M\nL\t{b}\t+\t{c}\t+\t0M\nL\t{a}\t+\t{c}\t+\t0M\n'
          'W\tREF\t0\t{name}\t0\t7\t>{a}>{b}>{c}\nW\ts1\t1\t{name}\t0\t5\t>{a}>{c}\n')


class TestCheck(unittest.TestCase):
    def test_disjoint(self):
        merge.check_disjoint({'a': (1, 3), 'b': (10, 12)})
        with self.assertRaises(ValueError):
            merge.check_disjoint({'a': (1, 10), 'b': (10, 12)})


@unittest.skipUnless(shutil.which('vg'), 'vg not on PATH')
class TestMerge(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def gbz(self, name, a, b, c):
        gfa, gbz = os.path.join(self.d, name + '.gfa'), os.path.join(self.d, name + '.gbz')
        with open(gfa, 'w') as f:
            f.write(CONTIG.format(name=name, a=a, b=b, c=c))
        subprocess.run(['vg', 'gbwt', '-G', gfa, '--gbz-format', '-g', gbz], check=True, capture_output=True)
        return gbz

    def test_merge(self):
        out = os.path.join(self.d, 'genome.gbz')
        s = merge.merge([self.gbz('chrB', 10, 11, 12), self.gbz('chrA', 1, 2, 3)], out)
        self.assertEqual((s['nodes'], s['edges'], s['length'], s['paths']), (6, 6, 14, 4))
        self.assertEqual(s['reference_samples'], ['REF'])
        ref = subprocess.run(['vg', 'paths', '-F', '-x', out, '-Q', 'REF'], check=True, capture_output=True,
                             text=True).stdout
        self.assertEqual(ref.split(), ['>REF#0#chrA', 'ACGTGGT', '>REF#0#chrB', 'ACGTGGT'])

    def test_overlap(self):
        with self.assertRaises(ValueError):
            merge.merge([self.gbz('chrA', 1, 2, 3), self.gbz('chrB', 3, 4, 5)], os.path.join(self.d, 'g.gbz'))


if __name__ == '__main__':
    unittest.main()
