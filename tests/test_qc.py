import gzip
import json
import os
import shutil
import tempfile
import unittest

from pgrealign import memory, qc, star


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

    def test_graph_pairs(self):
        length = {1: 2, 2: 3, 3: 1}
        a = qc.graph_positions([('>', 1), ('>', 2)], length)               # 5 bases
        b = qc.graph_positions([('>', 1), ('>', 3), ('>', 2)], length)     # 6 bases: an insertion
        self.assertEqual(qc.graph_pairs(a, b), {(0, 0), (1, 1), (2, 3), (3, 4), (4, 5)})
        c = qc.graph_positions([('<', 2), ('<', 1)], length)               # the reverse strand: no pairs
        self.assertEqual(qc.graph_pairs(a, c), set())
        d = qc.graph_positions([('>', 1), ('>', 1)], length)               # a cycle: every copy pairs
        self.assertEqual(qc.graph_pairs(a, d), {(0, 0), (1, 1), (0, 2), (1, 3)})

    def test_sequence_bp(self):
        self.assertEqual(qc.sequence_bp({'a': 'AC-T', 'b': 'AGGT'}), 5)

    def test_summarise(self):
        r = dict.fromkeys(qc.FIELDS, 0)
        r.update(id='R', longest=10, sequence=12, pair_tp=8, pair_opt=10, pair_msa=8)
        s = qc.summarise([r])
        self.assertEqual((s['pair_precision'], s['pair_recall'], s['growth']), (1.0, 0.8, 1.2))
        self.assertIsNone(s['ref_f1'])

    def test_pairs_over_the_budget(self):
        seq = {'a1': 'ACGTACGTAC', 'a2': 'ACGTTCGTAC', 'a3': 'GATTACA' * 400}       # 10, 10 and 2,800 bp
        pkg = {'id': 'R', 'ref': 'a1', 'alleles': [{'id': a, 'weight': 1, 'seq': x} for a, x in seq.items()]}
        width = sum(len(x) for x in seq.values())
        rows, k = {}, 0
        for a, x in seq.items():                     # an alignment that aligns nothing
            rows[a] = '-' * k + x + '-' * (width - k - len(x))
            k += len(x)
        calls = []

        def align(seqs, banded=False):                # no base aligned either: one deletion, one insertion
            calls.append((len(seqs[0]), len(seqs[1]), banded))
            return [seqs[0] + '-' * len(seqs[1]), '-' * len(seqs[0]) + seqs[1]]
        self.assertEqual(qc.compared_pairs(pkg, rows), [('pair', 'a2', 'a3'), ('ref', 'a1', 'a2'), ('ref', 'a1', 'a3')])
        limit = qc.pair_bytes(10, 10) + 1000                  # the pairs with a3 need more
        r = qc.region_qc(pkg, rows, align, max_cells=50, max_bytes=limit)
        # a1 and a2 are compared, banded since 100 cells exceed 50; the pairs with a3, over max_cells too, are not
        self.assertEqual((r['skipped'], r['banded'], calls), (2, 1, [(10, 10, True)]))
        self.assertEqual(qc.peak_bytes(pkg, rows, max_cells=50, max_bytes=limit), qc.pair_bytes(10, 10, 50))
        self.assertEqual(qc.peak_bytes(pkg, rows, max_cells=50), qc.pair_bytes(10, 2800, 50))
        self.assertEqual(qc.peak_bytes(pkg, rows, max_bytes=1), memory.BASE)
        r = qc.region_qc(pkg, rows, align, max_cells=50)
        self.assertEqual((r['skipped'], r['banded']), (0, 3))
        # with the graph, every pair's estimate holds the alleles' graph positions too
        held = qc.GRAPH_BYTES * width
        graph = ({a: [('>', k)] for k, a in enumerate(seq)}, {k: len(x) for k, x in enumerate(seq.values())}, 0)
        self.assertEqual(qc.peak_bytes(pkg, rows, max_cells=50, graph=True), held + qc.pair_bytes(10, 2800, 50))
        self.assertEqual(qc.peak_bytes(pkg, rows, max_cells=50, max_bytes=limit, graph=True), held + memory.BASE)
        for budget, skipped in ((limit, 3), (limit + held, 2)):
            r = qc.region_qc(pkg, rows, align, max_cells=50, graph=graph, max_bytes=budget)
            self.assertEqual(r['skipped'], skipped)


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

    def test_graph_baseline(self):
        # a graph that aligns the two alleles' flanks and leaves their middles unaligned
        pkg = {'id': 'R', 'ref': 'a1', 'alleles': [
            {'id': 'a1', 'weight': 1, 'seq': 'GATTACA' + 'ACGTACGTAC' + 'CCGG'},
            {'id': 'a2', 'weight': 1, 'seq': 'GATTACA' + 'ACGTACGTTC' + 'CCGG'}], 'fragments': []}
        length = {1: 7, 2: 10, 3: 10, 4: 4}
        graph = ({'a1': [('>', 1), ('>', 2), ('>', 4)], 'a2': [('>', 1), ('>', 3), ('>', 4)]}, length, 31)
        with star.Abpoa() as align:
            rows = star.star_msa([(a['id'], a['seq'], a['weight']) for a in pkg['alleles']], align=align).rows
            r = qc.region_qc(pkg, rows, align, graph=graph)
        self.assertEqual((r['ref_opt'], r['ref_msa'], r['ref_tp']), (21, 21, 21))
        self.assertEqual((r['graph_ref_msa'], r['graph_ref_tp'], r['graph_sequence']), (11, 11, 31))
        s = qc.summarise([r])
        self.assertEqual((s['ref_recall'], s['graph_ref_recall'], s['graph_ref_precision']), (1.0, 0.5238, 1.0))



@unittest.skipUnless(shutil.which('abpoa'), 'abpoa not on PATH')
class TestRun(unittest.TestCase):
    def test_memory_budget(self):
        d = tempfile.mkdtemp()
        try:
            unit = 'ACGTTGCAAGTCAT'
            regions = {'S': [unit * 4, unit * 2 + 'G' + unit[1:] + unit, unit * 3 + 'T' + unit[1:]],
                       'B': [unit * 215, unit * 100 + 'G' + unit[1:] + unit * 114, unit * 214 + 'ACGTTGCAAGTCAA']}
            os.makedirs(os.path.join(d, 'x'))
            with gzip.open(os.path.join(d, 'x', 'packages.jsonl.gz'), 'wt') as p, \
                    gzip.open(os.path.join(d, 'msas.jsonl.gz'), 'wt') as m:
                for rid in ('B', 'S'):
                    seqs = regions[rid]
                    self.assertEqual(len({len(x) for x in seqs}), 1)      # so the sequences are their own rows
                    p.write(json.dumps({'id': rid, 'status': 'ok', 'ref': 'a0', 'fragments': [], 'alleles': [
                        {'id': 'a%d' % i, 'weight': 1, 'seq': x} for i, x in enumerate(seqs)]}) + '\n')
                    m.write(json.dumps({'id': rid, 'status': 'ok', 'rows': {'a%d' % i: x for i, x in enumerate(seqs)}})
                            + '\n')
            limit = 0.2               # GB: S's pairs fit, and B's, about 0.27 each unbanded, do not
            self.assertGreater(qc.pair_bytes(len(unit) * 215, len(unit) * 215), limit * memory.GB)
            out = os.path.join(d, 'qc.tsv')
            args = (os.path.join(d, 'x'), os.path.join(d, 'msas.jsonl.gz'), out)
            # the budget decides only how many regions run at once: B needs more than 0.2 GB, so nothing starts
            with self.assertRaisesRegex(RuntimeError, 'the largest, B, needs about 0.2'):
                qc.run(*args, jobs=2, max_memory=limit)
            # max_pair_memory decides which pairs are compared, the same under any budget
            results = []
            for budget in (limit, 4):
                s = qc.run(*args, jobs=2, max_memory=budget, max_pair_memory=limit)
                lines = [l.rstrip('\n').split('\t') for l in open(out)]
                got = {x[0]: dict(zip(lines[0], x)) for x in lines[1:]}
                self.assertEqual([x[0] for x in lines[1:]], ['B', 'S'])
                self.assertEqual((got['B']['skipped'], got['B']['ref_opt'], got['B']['pair_opt']), ('3', '0', '0'))
                self.assertEqual(got['S']['skipped'], '0')
                self.assertGreater(int(got['S']['ref_opt']), 0)
                self.assertEqual(s['skipped'], 3)
                self.assertLessEqual(s['peak_estimate_gb'], budget)
                results.append(lines)
            self.assertEqual(results[0], results[1])
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
