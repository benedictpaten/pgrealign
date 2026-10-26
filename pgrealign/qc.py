"""Stage `qc`: how good each region's multiple alignment is, measured against pairwise alignments.

A multiple alignment implies an alignment of every pair of its rows: the columns where both rows have a
base. For a star alignment, the pair through the centre is optimal, and every other pair is only
induced. `qc` measures, region by region, how far the induced pairs fall from the pairs' own optimal
alignments, under the same scores. It also measures how much sequence the alignment adds to the graph,
and how much of it forces unrelated sequence into shared columns. Given the graph extract read, it
measures that graph's own alignment of the same pairs too, so each region has a before and an after.

**Metrics, per region.** All counts are summed over the sampled pairs.
- `sequence`: the bases of the induced graph, i.e. the sum over columns of the distinct bases in the
  column. `longest` is the longest allele, so `sequence / longest` is the alignment's growth.
- `pair_tp`, `pair_opt`, `pair_msa`: up to `pairs` random pairs of alleles, the reference allele
  excluded. `pair_opt` counts the base pairs the optimal alignment aligns, `pair_msa` those the
  induced alignment aligns, and `pair_tp` those both align. Precision is tp/msa, recall tp/opt, and
  F1 2tp/(opt+msa).
- `ref_tp`, `ref_opt`, `ref_msa`: the same for the reference allele against up to `refs` of the
  heaviest other alleles. These alignments decide how a genotyper's records are written against the
  reference.
- `ref_sv` and `ref_sv_opt`: runs of 50 or more gap columns in those induced and optimal alignments.
  An excess of induced runs means indels the reference comparison does not need.
- `forced` and `forced_opt`: columns of the induced (resp. optimal) alignments that lie in a
  forced-homology window: 40 columns, of which at least 20 have both bases, with more than 30% of
  those mismatched. Homologous tandem-repeat alleles differ mostly by gaps, so such a window shows
  unrelated sequence aligned base against base.
- `banded`: the pairs compared with star.long_pair's alignment rather than an optimal one (below).
- `skipped`: the pairs not compared, because aligning them would exceed `max_pair_memory` (below).
- With the graph: `graph_pair_tp` and `graph_pair_msa`, `graph_ref_tp` and `graph_ref_msa`, the same
  counts for the alignment the graph implies, against the same optimal alignments; and
  `graph_sequence`, the bases of the nodes the region's runs visit, which `replace` removes. The graph
  aligns two bases when both alleles' walks read them from one node, at one offset, in one
  direction. Each allele takes the walk of its first spanning run. A walk that visits a node twice
  aligns each of its bases there with each base the other walk reads there, so a graph that collapses
  copies of a repeat into a cycle aligns more pairs than any alignment can, and loses precision.

The optimal alignment of a pair is global, unbanded and left-normalised, under the scores given, which
should be those the alignment was built with. A pair whose unbanded alignment would exceed `max_cells`
cells is compared with the alignment star.long_pair makes of it instead, which may fall short of the
optimum. For a pair through the centre that is too large for `realign` too, that is the alignment
`realign` used. Pairs are drawn with a seed fixed per region, so a rerun measures the same pairs.

**Memory.** `pair_bytes` estimates the memory that aligning one pair takes (memory.py), and with
the graph a region also holds GRAPH_BYTES per base of its alleles, their graph positions, while it
aligns. A pair's estimate is the two together, with pages of memory.ESTIMATE_PAGE bytes whatever the
machine's. A pair whose estimate exceeds `max_pair_memory` GB is not compared, and counts in `skipped`
and in no other field, the same on every machine. `run` measures regions in up to `jobs` worker
processes, starting them in order while the estimates of those running, each that of its largest
compared pair (`peak_bytes`), stay within `max_memory` GB. That budget decides only how many run at
once: if a region needs more on its own, `run` stops before it measures anything.

**Output.** A TSV with one row per region that realign made `ok` and that has at least two alleles,
and a summary (written next to it as .json) with the sums and the F1 values over all regions.
"""
import collections
import gzip
import json
import os
import random
import zlib

from . import memory
from .extract import flip, read_packages
from .graph import STEP, GfaStream, parse_w, walk_key
from .star import Abpoa, SCORES, align_pair

MAX_PAIR_MEMORY = 64    # GB: a pair whose estimate exceeds it is not compared, on any machine

FIELDS = ['id', 'alleles', 'longest', 'sequence', 'pair_tp', 'pair_opt', 'pair_msa', 'ref_tp', 'ref_opt',
          'ref_msa', 'ref_sv', 'ref_sv_opt', 'forced', 'forced_opt', 'banded', 'skipped']
GRAPH_FIELDS = ['graph_pair_tp', 'graph_pair_msa', 'graph_ref_tp', 'graph_ref_msa', 'graph_sequence']
SV_RUN = 50
WINDOW, MIN_SHARED, MAX_MISMATCH = 40, 20, 0.30
GRAPH_BYTES = 120       # per base of the alleles, their graph positions: measured 72-112, more with longer nodes


def induced(ra, rb):
    """The alignment of two rows that a multiple alignment implies: its columns where either has a base."""
    keep = [i for i in range(len(ra)) if ra[i] != '-' or rb[i] != '-']
    return ''.join(ra[i] for i in keep), ''.join(rb[i] for i in keep)


def aligned_pairs(ra, rb):
    """The set of (i, j): base i of the first sequence aligned to base j of the second."""
    out, i, j = set(), 0, 0
    for x, y in zip(ra, rb):
        if x != '-' and y != '-':
            out.add((i, j))
        i += x != '-'
        j += y != '-'
    return out


def long_gap_runs(ra, rb, run=SV_RUN):
    """The number of runs of `run` or more consecutive gaps, in either row, opposite bases."""
    n = 0
    for row, other in ((ra, rb), (rb, ra)):
        length = 0
        for x, y in zip(row, other):
            if x == '-' and y != '-':
                length += 1
            else:
                n += length >= run
                length = 0
        n += length >= run
    return n


def forced_columns(ra, rb):
    """Columns of a two-row alignment that lie in a forced-homology window (see the module docstring)."""
    st = [0 if x == y else 1 if (x != '-' and y != '-') else 2 for x, y in zip(ra, rb)
          if not (x == '-' and y == '-')]
    n = len(st)
    if n < WINDOW:
        return 0
    flag = bytearray(n)
    shared = sum(1 for s in st[:WINDOW] if s < 2)
    mism = sum(1 for s in st[:WINDOW] if s == 1)
    for i in range(0, n - WINDOW + 1):
        if i:
            a, b = st[i - 1], st[i + WINDOW - 1]
            shared += (b < 2) - (a < 2)
            mism += (b == 1) - (a == 1)
        if shared >= MIN_SHARED and mism > MAX_MISMATCH * shared:
            flag[i:i + WINDOW] = b'\x01' * WINDOW
    return sum(flag)


def sequence_bp(rows):
    """The sum over the alignment's columns of the distinct bases in each."""
    names = list(rows)
    if not names:
        return 0
    return sum(len({rows[m][c] for m in names} - {'-'}) for c in range(len(rows[names[0]])))


def graph_positions(path, length):
    """(node, offset, reverse) of each base a walk's steps read, `path` being [(orient, node)]."""
    out = []
    for o, n in path:
        k = length[n]
        if o == '>':
            out.extend((n, i, False) for i in range(k))
        else:
            out.extend((n, i, True) for i in range(k - 1, -1, -1))
    return out


def graph_pairs(pa, pb):
    """The set of (i, j): base i of one walk and base j of another read at one graph position (see
    graph_positions)."""
    at = collections.defaultdict(list)
    for j, key in enumerate(pb):
        at[key].append(j)
    return {(i, j) for i, key in enumerate(pa) for j in at.get(key, ())}


def compared_pairs(pkg, rows, pairs=20, refs=20):
    """The pairs of alleles region_qc compares, in order, as (kind, a, b): up to `pairs` random pairs of
    alleles other than the reference allele ('pair'), then, when `rows` has the reference allele, it
    with up to `refs` of the heaviest others ('ref')."""
    alleles = [(a['id'], a['weight']) for a in pkg['alleles']]
    ref = pkg.get('ref')
    others = [a for a, _ in alleles if a != ref]
    rng = random.Random(zlib.crc32(pkg['id'].encode()))
    cand = [(i, j) for i in range(len(others)) for j in range(i + 1, len(others))]
    out = [('pair', others[i], others[j]) for i, j in sorted(rng.sample(cand, min(pairs, len(cand))))]
    if ref is not None and ref in rows:
        heaviest = sorted((a for a in alleles if a[0] != ref), key=lambda a: (-a[1], a[0]))[:refs]
        out += [('ref', ref, b) for b, _ in heaviest]
    return out


def pair_bytes(n, m, max_cells=None, page=memory.ESTIMATE_PAGE):
    """The memory, in bytes, that aligning an allele of m bases to one of n bases takes: memory.BASE and
    abPOA's matrices (memory.pair_bytes), with pages of `page` bytes."""
    return memory.BASE + memory.pair_bytes(n, m, max_cells, page)


def peak_bytes(pkg, rows, pairs=20, refs=20, max_cells=None, max_bytes=None, graph=False):
    """The most memory, in bytes, that region_qc takes for a region: what it holds, GRAPH_BYTES per base
    of the alleles with `graph`, plus the largest pair_bytes of the pairs it compares, or plus
    memory.BASE when it compares none. A pair whose sum exceeds `max_bytes` is not compared."""
    length = {a['id']: len(a['seq']) for a in pkg['alleles']}
    held = GRAPH_BYTES * sum(length.values()) if graph else 0
    need = [held + pair_bytes(length[a], length[b], max_cells)
            for _, a, b in compared_pairs(pkg, rows, pairs, refs)]
    return max([b for b in need if max_bytes is None or b <= max_bytes], default=held + memory.BASE)


def region_qc(pkg, rows, align, pairs=20, refs=20, max_cells=None, graph=None, max_bytes=None, scores=SCORES):
    """The metrics of one region (see the module docstring), as a dict keyed by FIELDS, and by
    GRAPH_FIELDS too when `graph` is given: ({allele: [(orient, node)]}, {node: length}, interior bp).
    A pair whose pair_bytes, with GRAPH_BYTES per base of the alleles when `graph` is given, exceeds
    `max_bytes` is not compared, and counts in `skipped`. `scores` are the aligner's."""
    alleles = [(a['id'], a['seq'].upper(), a['weight']) for a in pkg['alleles']]
    seqs = {a: s for a, s, _ in alleles}
    out = dict.fromkeys(FIELDS + (GRAPH_FIELDS if graph else []), 0)
    out.update(id=pkg['id'], alleles=len(alleles), longest=max(len(s) for s in seqs.values()),
               sequence=sequence_bp(rows))
    held = GRAPH_BYTES * sum(len(s) for s in seqs.values()) if graph else 0
    if graph:
        paths, length, out['graph_sequence'] = graph
        position = {a: graph_positions(paths[a], length) for a in seqs}
        for a in seqs:
            if len(position[a]) != len(seqs[a]):
                raise ValueError('%s: allele %s spells %d bp, its walk %d' % (pkg['id'], a, len(seqs[a]),
                                                                            len(position[a])))

    def compare(kind, a, b):
        if max_bytes is not None and held + pair_bytes(len(seqs[a]), len(seqs[b]), max_cells) > max_bytes:
            out['skipped'] += 1
            return
        banded = max_cells is not None and len(seqs[a]) * len(seqs[b]) > max_cells
        out['banded'] += banded
        opt = align_pair(seqs[a], seqs[b], align, max_cells, scores)
        ind = induced(rows[a], rows[b])
        o, m = aligned_pairs(*opt), aligned_pairs(*ind)
        out[kind + '_tp'] += len(o & m)
        out[kind + '_opt'] += len(o)
        out[kind + '_msa'] += len(m)
        out['forced'] += forced_columns(*ind)
        out['forced_opt'] += forced_columns(*opt)
        if kind == 'ref':
            out['ref_sv'] += long_gap_runs(*ind)
            out['ref_sv_opt'] += long_gap_runs(*opt)
        if graph:
            g = graph_pairs(position[a], position[b])
            out['graph_%s_tp' % kind] += len(o & g)
            out['graph_%s_msa' % kind] += len(g)

    for kind, a, b in compared_pairs(pkg, rows, pairs, refs):
        compare(kind, a, b)
    return out


def summarise(rows):
    """Sums over regions, with precision, recall and F1 of the pair and reference agreement, for the
    alignment and, when the rows have them, for the graph."""
    fields = FIELDS + (GRAPH_FIELDS if rows and 'graph_sequence' in rows[0] else [])
    s = collections.Counter()
    for r in rows:
        for k in fields[1:]:
            s[k] += r[k]
    out = {'regions': len(rows), **{k: s[k] for k in fields[1:]}}
    for g in ('', 'graph_')[:1 + ('graph_sequence' in s)]:
        for p in ('pair', 'ref'):
            tp, opt, msa = s[g + p + '_tp'], s[p + '_opt'], s[g + p + '_msa']
            out[g + p + '_precision'] = round(tp / msa, 4) if msa else None
            out[g + p + '_recall'] = round(tp / opt, 4) if opt else None
            out[g + p + '_f1'] = round(2 * tp / (opt + msa), 4) if opt + msa else None
        seq = s[g + 'sequence']
        out[g + 'growth'] = round(seq / s['longest'], 4) if s['longest'] else None
    return out


def _one(args):
    pkg, rows, graph, kw = args
    with Abpoa(kw['abpoa'], kw['scores'], workdir=kw['workdir']) as align:
        return region_qc(pkg, rows, align, kw['pairs'], kw['refs'], kw['max_cells'], graph, kw['max_bytes'],
                         kw['scores'])


def graph_walks(extract_dir, graph, regions, vg=None, threads=4):
    """{region: ({allele: [(orient, node)]}, {node: length}, interior bp)} for `regions`, from the graph
    extract read: each allele's walk is that of its first spanning run in occurrences.tsv.gz, in the
    region's forward orientation."""
    want = collections.defaultdict(list)          # walk -> [(first, end, orient, region, allele)]
    seen = set()
    with gzip.open(os.path.join(extract_dir, 'occurrences.tsv.gz'), 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            walk, region, typ, orient, first, end, member = line.rstrip('\n').split('\t')
            if typ == 'spanning' and region in regions and (region, member) not in seen:
                seen.add((region, member))
                want[walk].append((int(first), int(end), orient, region, member))
    interior = {r: [] for r in regions}
    with gzip.open(os.path.join(extract_dir, 'interiors.tsv.gz'), 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            region, nodes = line.rstrip('\n').split('\t')
            if region in regions:
                interior[region] = [int(n) for n in nodes.split(',') if n]
    needed = {n for ns in interior.values() for n in ns}
    length, paths = {}, collections.defaultdict(dict)
    with GfaStream(graph, vg, threads) as f:
        for line in f:
            if line[0] == 'S':
                x = line.split('\t', 3)
                n = int(x[1])
                if n in needed:
                    length[n] = len(x[2].rstrip('\n'))
            elif line[0] == 'W':
                w = parse_w(line)
                runs = want.get(walk_key(w))
                if runs:
                    toks = STEP.findall(w.steps)
                    for a, b, orient, region, member in runs:
                        p = [(o, int(n)) for o, n in toks[a:b]]
                        paths[region][member] = p if orient == '+' else flip(p)
    return {r: (paths[r], {n: length[n] for n in interior[r]}, sum(length[n] for n in interior[r]))
            for r in regions}


def reference_alleles(occurrences_path, ref_sample):
    """{region: allele id} for the spanning runs of `ref_sample`'s walks, read from extract's
    occurrences.tsv.gz. For packages written before extract recorded `ref`."""
    out = {}
    prefix = ref_sample + '#'
    with gzip.open(occurrences_path, 'rt') as f:
        for line in f:
            if line.startswith(prefix):
                walk, region, typ, _, _, _, member = line.rstrip('\n').split('\t')
                if typ == 'spanning':
                    out[region] = member
    return out


def run(extract_dir, msas_path, out_path, pairs=20, refs=20, sample=None, jobs=1, abpoa='abpoa',
        scores=SCORES, workdir=None, ref_sample='CHM13', max_cells=450_000_000, graph=None, vg=None, threads=4,
        max_memory=None, max_pair_memory=MAX_PAIR_MEMORY):
    """Measure every region realign made `ok` (or a fixed random `sample` of them) and write the TSV and
    its summary. `extract_dir` is extract's output directory, and `graph`, if given, the graph extract
    read (GBZ or GFA), whose own alignment is then measured too. A pair whose estimate exceeds
    `max_pair_memory` GB (0 or None: no limit) is not compared; `max_memory` is the budget in GB that
    decides how many regions run at once (see the module docstring; None: memory.default_gb()), and a
    region that needs more alone raises RuntimeError before anything is measured. Returns the summary,
    with `peak_estimate_gb`, the largest sum of the estimates of regions measured at once."""
    limit = memory.budget(max_memory)
    max_bytes = int(max_pair_memory * memory.GB) if max_pair_memory else None
    packages_path = os.path.join(extract_dir, 'packages.jsonl.gz')
    msas = {}
    with gzip.open(msas_path, 'rt') as f:
        for line in f:
            m = json.loads(line)
            if m.get('status') == 'ok' and len(m['rows']) > 1:
                msas[m['id']] = m['rows']
    if sample is not None and sample < len(msas):
        keep = set(random.Random(7).sample(sorted(msas), sample))
        msas = {k: v for k, v in msas.items() if k in keep}
    todo = [p for p in read_packages(packages_path) if p['id'] in msas and len(p['alleles']) > 1]
    if any('ref' not in p for p in todo):
        ref = reference_alleles(os.path.join(extract_dir, 'occurrences.tsv.gz'), ref_sample)
        for p in todo:
            p.setdefault('ref', ref.get(p['id']))
    walks = graph_walks(extract_dir, graph, {p['id'] for p in todo}, vg, threads) if graph else {}
    kw = {'pairs': pairs, 'refs': refs, 'abpoa': abpoa, 'scores': scores, 'workdir': workdir, 'max_cells': max_cells,
          'max_bytes': max_bytes}
    need = [peak_bytes(p, msas[p['id']], pairs, refs, max_cells, max_bytes, bool(graph)) for p in todo]
    over = [k for k in range(len(todo)) if need[k] > limit]
    if over:
        k = max(over, key=lambda k: need[k])
        raise RuntimeError('%d regions need more memory than the budget of %.2f GB (max_memory); the largest, %s, '
                           'needs about %.2f GB. Give qc that much, or lower max_pair_memory below it, which leaves '
                           'such pairs uncompared.' % (len(over), limit / memory.GB, todo[k]['id'],
                                                       need[k] / memory.GB))
    out = [None] * len(todo)
    peak = memory.run_within(_one, [(p, msas[p['id']], walks.get(p['id']), kw) for p in todo], need, jobs, limit,
                             out.__setitem__)
    fields = FIELDS + (GRAPH_FIELDS if graph else [])
    summary = summarise(out)
    summary['peak_estimate_gb'] = round(peak / memory.GB, 2)
    with open(out_path + '.json.tmp', 'w') as f:
        json.dump(summary, f, indent=1)
        f.write('\n')
    with open(out_path + '.tmp', 'w') as f:
        f.write('\t'.join(fields) + '\n')
        for r in out:
            f.write('\t'.join(str(r[k]) for k in fields) + '\n')
    os.replace(out_path + '.json.tmp', out_path + '.json')
    os.replace(out_path + '.tmp', out_path)
    return summary
