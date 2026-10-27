import argparse
import random
import time
import unittest

from pgrealign import cli, memory, qc, realign, star
from pgrealign.memory import BASE, GB, abpoa_bytes, any_pair_bytes, pair_bytes


def pkg(alleles, fragments=()):
    return {'id': 'R', 'status': 'ok', 'alleles': [{'id': 'a%d' % i, 'weight': 1, 'seq': 'A' * n}
                                                   for i, n in enumerate(alleles)],
            'fragments': [{'id': 'f%d' % i, 'side': 'left', 'seq': 'A' * n} for i, n in enumerate(fragments)]}


class TestEstimate(unittest.TestCase):
    def test_bounds_the_measured_footprints(self):
        # memory.py's table: macOS footprints, 16 KB pages; abpoa_bytes with BASE must reach each one
        for n, m, banded, gb in [(7500, 7500, False, 1.05), (15000, 15000, False, 4.20), (30000, 30000, False, 16.8),
                                 (7500, 7500, True, 0.60), (15000, 15000, True, 1.24), (30000, 30000, True, 2.65),
                                 (58700, 58700, True, 5.85), (58700, 40100, True, 19.2), (58700, 103200, True, 47.1)]:
            self.assertGreaterEqual(BASE + abpoa_bytes(n, m, banded, page=16384), gb * GB, (n, m, banded))
        # star.long_pair's band, the widest within MAX_CELLS (or abPOA's own where the lengths differ too much for
        # that): peak resident memory on Linux, 4 KB pages, abPOA 1.5.7
        for n, m, gb in [(27709, 27752, 7.64), (27709, 37174, 6.91), (32286, 34412, 8.13), (27709, 48811, 10.19)]:
            self.assertGreaterEqual(BASE + pair_bytes(n, m, realign.MAX_CELLS, page=4096), gb * GB, (n, m))

    def test_monotone_in_length(self):
        # unbanded, then banded: every pair's matrix is under 450,000,000 cells in the first, over it in the second
        for lengths in ([100, 1000, 5000, 10000, 20000], [25000, 40000, 60000, 80000, 100000]):
            est = [realign.peak_bytes(pkg([n, n - 7, 3 * n // 4])) for n in lengths]
            self.assertEqual(est, sorted(est))
            self.assertLess(est[0], est[-1])
        # another allele or fragment never lowers it
        self.assertLessEqual(realign.peak_bytes(pkg([30000, 20000])), realign.peak_bytes(pkg([30000, 20000, 8000])))
        self.assertLessEqual(realign.peak_bytes(pkg([30000, 20000])), realign.peak_bytes(pkg([30000, 20000], [9000])))

    def test_over_max_cells(self):
        n = 20000
        self.assertEqual(any_pair_bytes([n, n], n * n), abpoa_bytes(n, n, False))
        self.assertEqual(any_pair_bytes([n, n], n * n - 1), pair_bytes(n, n, n * n - 1))
        # just over max_cells, a pair takes about what max_cells cells unbanded do; far over it, much less than
        # its own matrix unbanded, and never less than refinement's windows
        self.assertLess(abs(pair_bytes(n, n, n * n - 1, 0) - abpoa_bytes(n, n, False)), abpoa_bytes(n, n, False) / 20)
        big = 100000
        self.assertLess(pair_bytes(big, big, realign.MAX_CELLS), abpoa_bytes(big, big, False) / 10)
        self.assertGreaterEqual(pair_bytes(big, big, realign.MAX_CELLS), 20 * realign.MAX_CELLS)
        self.assertEqual(pair_bytes(big, big, realign.MAX_CELLS, 0),
                         abpoa_bytes(big, big, star.band_width(big, big, realign.MAX_CELLS), 0))
        self.assertEqual(qc.pair_bytes(n, n, n * n) - BASE, abpoa_bytes(n, n, False))
        # qc's estimate, like realign's, assumes ESTIMATE_PAGE whatever the machine's page size
        self.assertEqual(qc.pair_bytes(n, n, n * n - 1) - BASE, pair_bytes(n, n, n * n - 1, memory.ESTIMATE_PAGE))
        # a region realign refuses before aligning takes only BASE
        self.assertEqual(realign.peak_bytes(pkg([n, n]), max_cells=1000), BASE)
        self.assertEqual(realign.peak_bytes(pkg([], [n])), BASE)

    def test_page_size(self):
        n = 50000
        for page in (4096, 16384):
            self.assertEqual(abpoa_bytes(n, n, True, page) - abpoa_bytes(n, n, True, 0), 5 * n * page)
            self.assertEqual(abpoa_bytes(n, n, False, page), 20 * n * n)
        small, large = (realign.peak_bytes(pkg([n, n - 100]), page=p) for p in (4096, 16384))
        self.assertEqual(large - small, 5 * n * (16384 - 4096))

    def test_any_pair_is_the_worst_pair(self):
        rng = random.Random(1)
        for _ in range(300):
            lengths = [rng.choice([0, 50, 3000, 9000, 20000, 41000, 90000]) for _ in range(rng.randint(2, 6))]
            max_cells = rng.choice([10 ** 6, 10 ** 8, 450_000_000])
            worst = max(pair_bytes(lengths[i], lengths[j], max_cells)
                        for i in range(len(lengths)) for j in range(len(lengths)) if i != j)
            self.assertEqual(any_pair_bytes(lengths, max_cells), worst, (lengths, max_cells))

    def test_budget(self):
        self.assertEqual(memory.budget(2), 2 * GB)
        self.assertEqual(memory.budget(None), int(memory.default_gb() * GB))
        self.assertGreater(memory.default_gb(), 0)
        for bad in (0, -1, float('nan'), float('inf')):
            with self.assertRaises(ValueError):
                memory.budget(bad)
        self.assertEqual(cli.gigabytes('2.5'), 2.5)
        for bad in ('0', '-3', 'nan', 'many'):
            with self.assertRaises(argparse.ArgumentTypeError):
                cli.gigabytes(bad)


def nap(item):
    """A fake worker: sleeps item[1] seconds, and returns item[0] with when it ran."""
    start = time.time()
    time.sleep(item[1])
    return item[0], start, time.time()


def fail(item):
    if item == 'bad':
        raise KeyError(item)
    return item


class TestRunWithin(unittest.TestCase):
    def check(self, costs, jobs, limit):
        items = [(i, 0.05 + 0.03 * (i % 4)) for i in range(len(costs))]
        got = {}

        def take(i, r):
            self.assertNotIn(i, got)
            got[i] = r
        peak = memory.run_within(nap, items, costs, jobs, limit, take)
        self.assertEqual(sorted(got), list(range(len(costs))))
        self.assertTrue(all(got[i][0] == i for i in got))
        span = {i: (s, e) for i, s, e in got.values()}
        for i, (s, _) in span.items():                  # what runs when each item starts
            now = [j for j, (s2, e2) in span.items() if s2 <= s < e2]
            self.assertLessEqual(len(now), jobs)
            if costs[i] > limit:
                self.assertEqual(now, [i])
            else:
                self.assertLessEqual(sum(costs[j] for j in now), limit, now)
        return peak

    def test_within_the_budget(self):
        costs = [4, 4, 4, 9, 2, 2, 2, 1, 3, 3]
        self.assertLessEqual(self.check(costs, 3, 10), 10)

    def test_oversized_runs_alone(self):
        self.assertEqual(self.check([2, 3, 15, 2, 2], 3, 10), 15)

    def test_jobs(self):
        self.assertLessEqual(self.check([0] * 8, 2, 10), 0)
        self.check([1] * 7, 1, 100)

    def test_schedule_seconds(self):
        # in order, at most `jobs` at once, while the costs of those running fit the limit
        self.assertEqual(memory.schedule_seconds([1, 1, 1], [5, 3, 4], 3, 10), 5)       # all at once
        self.assertEqual(memory.schedule_seconds([1, 1, 1], [5, 3, 4], 1, 10), 12)      # one at a time
        self.assertEqual(memory.schedule_seconds([10, 10, 10], [10, 10, 10], 8, 20), 20)  # two fit the limit
        # the second waits for the first, which leaves room, and the third, small, waits behind it: in order
        self.assertEqual(memory.schedule_seconds([15, 10, 1], [4, 6, 1], 8, 20), 10)
        self.assertEqual(memory.schedule_seconds([30, 1], [2, 2], 8, 20), 4)           # one over the limit, alone
        self.assertEqual(memory.schedule_seconds([], [], 4, 20), 0)

    def test_exception(self):
        got = []
        with self.assertRaises(KeyError):
            memory.run_within(fail, ['a', 'bad', 'c'], [1, 1, 1], 1, 10, lambda i, r: got.append(r))
        self.assertEqual(got, ['a'])
        self.assertEqual(memory.run_within(fail, [], [], 2, 10, None), 0)


if __name__ == '__main__':
    unittest.main()
