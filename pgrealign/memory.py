"""Memory: what one alignment takes, and running regions in parallel within a memory budget.

realign and qc work on regions in worker processes, and each worker runs abPOA once per alignment.
Before a region starts, each estimates from its package the most memory the region can take
(realign.peak_bytes, qc.peak_bytes), and `run_within` starts regions in order while the estimates of
the regions running, the next one's included, stay within a budget: by default half the machine's
physical memory (`default_gb`).

**abPOA's memory.** abPOA aligns a sequence of m bases to a graph, here the sequence of n bases given
before it, in five matrices of 4-byte scores, with one row per base of the graph and one column per
base of the sequence. It reserves all of them, banded or not, and only the pages it touches take
memory. Unbanded, it touches every cell: 20 bytes a cell. Banded, it touches in each row of each
matrix a run of cells from the alignment's best position so far to the diagonal that ends at the
matrix's corner, widened by abPOA's adaptive band, 10 cells plus 1% of m, on either side. The run is
about |n - m| cells long where the two lengths differ, and it starts at an arbitrary offset in a
page, so a row of a matrix touches about one page plus 4 bytes for each of its cells, and never
more than the whole matrix. `abpoa_bytes` computes this, with the machine's page size. Measured as
macOS's peak memory footprint, which counts memory compressed or swapped out as well as resident
(abPOA 1.5.7, 16 KB pages, scores -M 8 -X 16 -O 16,96 -E 8,1), against `abpoa_bytes` in parentheses:

| alignment | footprint, GB of 2^30 bytes |
|---|---|
| unbanded, two sequences of 7.5, 15, 30 kb | 1.05, 4.20, 16.8 (1.05, 4.19, 16.8) |
| banded, two near-equal sequences of 7.5, 15, 30, 58.7 kb | 0.60, 1.24, 2.65, 5.85 (0.60, 1.23, 2.64, 5.80) |
| banded, 58.7 kb then 40.1 kb | 19.2 (25.7) |
| banded, 58.7 kb then 103.2 kb | 47.1 (55.4) |

The band can reach further than |n - m|, as where indels in opposite directions cancel, so this is
an upper bound in practice, not a guarantee. The memory abPOA takes beyond its matrices, up to about
50 MB in these runs, and the worker's own are in BASE. With 4 KB pages, as on most Linux systems, the
per-page part is a quarter as large; that has not been measured.
"""
import collections
import concurrent.futures
import os
import resource

GB = 2 ** 30
PAGE = resource.getpagesize()
CELL_BYTES = 20                 # five matrices of 4-byte scores
BASE = 100 * 2 ** 20            # per region: a worker's interpreter, about 20 MB, and abPOA's own tables


def abpoa_bytes(n, m, banded, page=PAGE):
    """The memory, in bytes, that abPOA touches aligning a sequence of m bases to one of n bases given
    before it, with its adaptive band if `banded` (see the module docstring)."""
    full = CELL_BYTES * n * m
    if not banded:
        return full
    return min(full, 5 * n * (page + 4 * (abs(n - m) + 2 * (10 + m // 100))))


def any_pair_bytes(lengths, max_cells, page=PAGE):
    """The most abpoa_bytes can be for two of the sequences whose lengths are `lengths`, given in either
    order, banded when the product of their lengths exceeds `max_cells`."""
    count = collections.Counter(lengths)
    distinct = sorted(count)
    worst = 0
    for n in distinct:
        for m in reversed(distinct):
            if m == n and count[n] == 1:
                continue
            if n * m <= max_cells:          # the longest sequence aligned to this one unbanded
                worst = max(worst, abpoa_bytes(n, m, False))
                break
            worst = max(worst, abpoa_bytes(n, m, True, page))
    return worst


def default_gb():
    """Half the machine's physical memory, in GB: the default budget."""
    return os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / 2 / GB


def budget(max_memory=None):
    """`max_memory` GB as bytes; None means default_gb()."""
    if max_memory is None:
        max_memory = default_gb()
    if not 0 < max_memory < float('inf'):
        raise ValueError('the memory budget must be more than 0 GB, not %r' % (max_memory,))
    return int(max_memory * GB)


def run_within(fn, items, costs, jobs, limit, take):
    """Call fn(item) for each of `items` in up to `jobs` worker processes, and take(i, result) for each
    as it finishes, i being the item's index. Items start in order, each once the costs of those
    running, with its own, are at most `limit`; an item whose cost alone exceeds `limit` starts when
    no other runs, and runs alone. Returns the largest sum of the costs of items running at once.
    An exception that fn raises, or take, stops any further item from starting, and is raised once
    the items running have finished."""
    jobs = max(1, jobs)
    running, used, peak, nxt = {}, 0, 0, 0
    with concurrent.futures.ProcessPoolExecutor(jobs) as ex:
        while nxt < len(items) or running:
            while nxt < len(items) and len(running) < jobs and (not running or used + costs[nxt] <= limit):
                running[ex.submit(fn, items[nxt])] = nxt
                used += costs[nxt]
                nxt += 1
            peak = max(peak, used)
            finished, _ = concurrent.futures.wait(running, return_when=concurrent.futures.FIRST_COMPLETED)
            for fut in sorted(finished, key=running.get):
                i = running.pop(fut)
                used -= costs[i]
                take(i, fut.result())
    return peak
