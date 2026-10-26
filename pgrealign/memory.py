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
matrix's corner, widened by the adaptive band on either side: abPOA's own, 10 cells plus 1% of m
(star.default_band), or as many cells as it is given. The run is about |n - m| cells long where the
two lengths differ, and it starts at an arbitrary offset in a page, so a row of a matrix touches about
one page plus 4 bytes for each of its cells, and never more than the whole matrix. `abpoa_bytes`
computes this, with the machine's page size. Measured as
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

**A pair too large to align unbanded** (star.long_pair) takes the most of its alignments (`pair_bytes`):
the band of star.band_width, about `max_cells` cells, or more where the two lengths differ by more than
that allows; and refinement's windows, unbanded, at most `max_cells` cells each. Seeding's alignment,
window by window under abPOA's own band, is not counted: it took at most 1.3 GB resident in pairs of up to
71 kb (Linux, abPOA 1.5.7), against the band's 8 GB.

**abPOA's reservation.** The matrices are reserved whole, 20 bytes a cell, rounded up to a power of two
(`abpoa_reservation`): two 62,146 bp sequences reserve 2^37 bytes, banded or not. macOS grants any
reservation. Linux, under its default overcommit policy, refuses one larger than the machine's memory
and swap, and abPOA then fails ("posix_memalign fail"), so whether such a region can be aligned at all
depends on the machine. realign therefore refuses, from the package, a region that would need a larger
reservation than a fixed limit, and treats a refused reservation as a failure of the run, not of the
region.
"""
import collections
import concurrent.futures
import os
import resource

from .star import band_width, default_band

GB = 2 ** 30
PAGE = resource.getpagesize()
ESTIMATE_PAGE = 16 * 1024       # the page size the estimates that decide anything assume: macOS on Apple silicon's
CELL_BYTES = 20                 # five matrices of 4-byte scores
BASE = 100 * 2 ** 20            # per region: a worker's interpreter, about 20 MB, and abPOA's own tables


def abpoa_bytes(n, m, banded, page=PAGE):
    """The memory, in bytes, that abPOA touches aligning a sequence of m bases to one of n bases given
    before it, with its adaptive band if `banded`: abPOA's own if it is True, `banded` cells either side if
    it is a number (see the module docstring)."""
    full = CELL_BYTES * n * m
    if not banded:
        return full
    band = default_band(m) if banded is True else banded
    return min(full, 5 * n * (page + 4 * (abs(n - m) + 2 * band)))


def pair_bytes(n, m, max_cells, page=PAGE):
    """The memory, in bytes, that abPOA touches when star.align_pair aligns a sequence of m bases to one of
    n: unbanded if n * m is at most `max_cells`, or max_cells is None, and otherwise the most of
    star.long_pair's alignments (see the module docstring)."""
    if max_cells is None or n * m <= max_cells:
        return abpoa_bytes(n, m, False)
    return max(abpoa_bytes(n, m, band_width(n, m, max_cells), page), CELL_BYTES * max_cells)


def abpoa_reservation(n, m):
    """The address space, in bytes, that abPOA reserves to align a sequence of m bases to one of n: its
    five matrices whole, CELL_BYTES a cell, rounded up to a power of two."""
    full = CELL_BYTES * (n + 1) * (m + 1)
    return 1 << (full - 1).bit_length()


def granted_bytes():
    """The largest reservation this machine grants one process, in bytes, as far as can be told, or None for
    no limit: under Linux's default overcommit policy its memory and swap (MemTotal and SwapTotal in
    /proc/meminfo), under its strict policy its commit limit, and under 'always overcommit', or elsewhere,
    none. It decides nothing about a region, only whether realign can align here what it admits."""
    try:
        with open('/proc/sys/vm/overcommit_memory') as f:
            mode = int(f.read())
        info = {}
        with open('/proc/meminfo') as f:
            for line in f:
                key, value = line.split(':', 1)
                info[key] = int(value.split()[0]) * 1024
    except (OSError, ValueError, IndexError):
        return None
    if mode == 1:
        return None
    if mode == 2:
        return info.get('CommitLimit')
    return info.get('MemTotal', 0) + info.get('SwapTotal', 0)


def any_pair_bytes(lengths, max_cells, page=PAGE):
    """The most pair_bytes can be for two of the sequences whose lengths are `lengths`, given in either
    order."""
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
            worst = max(worst, pair_bytes(n, m, max_cells, page))
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
