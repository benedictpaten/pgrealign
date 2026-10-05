"""End to end, with vg: a two-contig genome, each contig through `run_contig`, then `merge`."""
import gzip
import json
import os
import shutil
import subprocess
import tempfile
import unittest

from pgrealign import merge, pipeline

FLANK_L = 'CCTAGGCTTAGGACCTGATCGGATCCAGTA'
FLANK_R = 'GGATCCTAGCATTCGGATCCAAGTTCCAGA'
UNIT = 'ACGTTGCAAGTCAT'


def genome_gfa(path):
    """chrT1: a 3-copy tandem repeat with 1-, 2- and 4-copy alleles. chrT2: a SNP, no repeat target."""
    s = {1: FLANK_L, 2: UNIT, 3: UNIT, 4: UNIT, 5: UNIT[:5] + 'G' + UNIT[6:], 9: FLANK_R,
         21: FLANK_L, 22: 'A', 23: 'C', 29: FLANK_R}
    e = [(1, 2), (2, 3), (3, 4), (4, 9), (4, 5), (5, 9), (2, 9), (3, 9), (21, 22), (21, 23), (22, 29), (23, 29)]
    walks = [('CHM13', 0, 'chrT1', [1, 2, 3, 4, 9]), ('h1', 1, 'chrT1', [1, 2, 9]), ('h1', 2, 'chrT1', [1, 2, 3, 9]),
             ('h2', 1, 'chrT1', [1, 2, 3, 4, 5, 9]), ('h2', 2, 'chrT1', [1, 2, 3, 9]),
             ('CHM13', 0, 'chrT2', [21, 22, 29]), ('h1', 1, 'chrT2', [21, 23, 29]), ('h1', 2, 'chrT2', [21, 22, 29]),
             ('h2', 1, 'chrT2', [21, 23, 29]), ('h2', 2, 'chrT2', [21, 23, 29])]
    with open(path, 'w') as f:
        f.write('H\tVN:Z:1.1\tRS:Z:CHM13\n')
        for n, q in s.items():
            f.write('S\t%d\t%s\n' % (n, q))
        for a, b in e:
            f.write('L\t%d\t+\t%d\t+\t0M\n' % (a, b))
        for sample, hap, contig, w in walks:
            length = sum(len(s[n]) for n in w)
            f.write('W\t%s\t%d\t%s\t0\t%d\t%s\n' % (sample, hap, contig, length, ''.join('>%d' % n for n in w)))


def path_seqs(gbz):
    out, name = {}, None
    for line in subprocess.run(['vg', 'paths', '-F', '-x', gbz], check=True, capture_output=True,
                               text=True).stdout.splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            out[name] = ''
        else:
            out[name] += line.strip()
    return out


@unittest.skipUnless(shutil.which('vg') and shutil.which('abpoa'), 'vg and abpoa needed')
class TestRun(unittest.TestCase):
    def test_two_contigs(self):
        d = tempfile.mkdtemp()
        try:
            gfa, gbz = os.path.join(d, 'genome.gfa'), os.path.join(d, 'genome.gbz')
            genome_gfa(gfa)
            subprocess.run(['vg', 'gbwt', '-G', gfa, '--gbz-format', '-g', gbz], check=True, capture_output=True)
            targets = os.path.join(d, 'targets.bed')
            with open(targets, 'w') as f:
                f.write('chrT1\t%d\t%d\tTR1\n' % (len(FLANK_L), len(FLANK_L) + 3 * len(UNIT)))
            out = os.path.join(d, 'out')
            quiet = lambda m: None
            s1 = pipeline.run_contig(gbz, 'chrT1', targets, out, drop_samples=(), log=quiet)
            self.assertEqual(s1['realign']['realigned'], {'ok': 1})
            self.assertTrue(s1['verify']['ok'])
            self.assertEqual(s1['qc']['regions'], 1)
            self.assertIsNotNone(s1['qc']['graph_pair_f1'] if s1['qc']['pair_opt'] else s1['qc']['graph_ref_f1'])
            pkg = json.loads(gzip.open(os.path.join(out, 'chrT1', 'extract', 'packages.jsonl.gz'), 'rt').readline())
            self.assertIsNotNone(pkg['ref'])
            again = pipeline.run_contig(gbz, 'chrT1', targets, out, drop_samples=(), log=quiet)
            self.assertEqual(set(again), {'prepare', 'id_range'})          # every later stage resumed
            ver = os.path.join(out, 'chrT1', 'verify.json')
            good = open(ver).read()
            with open(ver, 'w') as f:
                json.dump(dict(json.loads(good), ok=False), f)
            with self.assertRaises(RuntimeError):                          # a recorded failure stops every run
                pipeline.run_contig(gbz, 'chrT1', targets, out, drop_samples=(), log=quiet)
            with open(ver, 'w') as f:
                f.write(good)
            s2 = pipeline.run_contig(gbz, 'chrT2', targets, out, drop_samples=(), log=quiet)
            self.assertEqual(s2['regions']['regions'], 0)
            g = os.path.join(d, 'merged.gbz')
            merge.merge([os.path.join(out, c, c + '.realigned.gbz') for c in ('chrT1', 'chrT2')], g)
            self.assertEqual(path_seqs(g), path_seqs(gbz))
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
