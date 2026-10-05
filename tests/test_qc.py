import shutil
import unittest

from pgrealign import qc, star


class TestPieces(unittest.TestCase):
    def test_induced_and_pairs(self):
        a, b = qc.induced('AC-GT-', 'A--GTT')
        self.assertEqual((a, b), ('ACGT-', 'A-GTT'))
        self.assertEqual(qc.aligned_pairs(a, b), {(0, 0), (2, 1), (3, 2)})

    def test_long_gap_runs(self):
        self.assertEqual(qc.long_gap_runs('A' + '-' * 50 + 'C', 'A' + 'G' * 50 + 'C'), 1)
        self.assertEqual(qc.long_gap_runs('A' + '-' * 49 + 'C', 'A' + 'G' * 49 + 'C'), 0)
        self.assertEqual(qc.long_gap_runs('G' * 60, '-' * 60), 1)

    def test_forced_columns(self):
        same = 'ACGT' * 20
        self.assertEqual(qc.forced_columns(same, same), 0)
        other = ''.join({'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[c] for c in same)
        self.assertEqual(qc.forced_columns(same, other), len(same))
        self.assertEqual(qc.forced_columns(same[:30], other[:30]), 0)      # shorter than a window

    def test_sequence_bp(self):
        self.assertEqual(qc.sequence_bp({'a': 'AC-T', 'b': 'AGGT'}), 5)

    def test_summarise(self):
        r = dict.fromkeys(qc.FIELDS, 0)
        r.update(id='R', longest=10, sequence=12, pair_tp=8, pair_opt=10, pair_msa=8)
        s = qc.summarise([r])
        self.assertEqual((s['pair_precision'], s['pair_recall'], s['growth']), (1.0, 0.8, 1.2))
        self.assertIsNone(s['ref_f1'])


@unittest.skipUnless(shutil.which('abpoa'), 'abpoa not on PATH')
class TestRegion(unittest.TestCase):
    def test_star_agrees_with_itself_through_the_centre(self):
        unit = 'ACGTTGCAAGT'
        pkg = {'id': 'R', 'ref': 'a1', 'alleles': [
            {'id': 'a1', 'weight': 5, 'seq': 'GATTACA' + unit * 4 + 'CCGG'},
            {'id': 'a2', 'weight': 1, 'seq': 'GATTACA' + unit * 6 + 'CCGG'},
            {'id': 'a3', 'weight': 1, 'seq': 'GATTACA' + unit * 2 + 'CCGG'}], 'fragments': []}
        with star.Abpoa() as align:
            rows = star.star_msa([(a['id'], a['seq'], a['weight']) for a in pkg['alleles']], align=align).rows
            r = qc.region_qc(pkg, rows, align)
        # the reference is the centre here, so its induced alignments are the optimal ones
        self.assertEqual((r['ref_tp'], r['ref_opt'], r['ref_msa']), (r['ref_opt'],) * 3)
        self.assertGreater(r['ref_opt'], 0)
        self.assertEqual(r['alleles'], 3)
        self.assertGreaterEqual(r['sequence'], r['longest'])


if __name__ == '__main__':
    unittest.main()
