"""Stage 0, `catalog`: a tandem-repeat catalogue of the reference.

This is the only stage that knows the regions are tandem repeats. Any BED of target intervals can take
the catalogue's place in the stages after it. The catalogue takes no sample or truth data as input.

**Sources.**
- **RepeatMasker** (a BED of ten columns: chrom start end name score strand class family divergence
  id). Simple_repeat records become STR when the motif in "(MOTIF)n" is 1-6 bp and VNTR when it is
  longer. Satellite records become SAT and Low_complexity records LC. RepeatMasker reports no
  Simple_repeat motif longer than 30 bp, so on its own it misses most VNTRs: 511 of the 671 VNTRs
  that the earlier chr6 test found.
- **trfind** (c/trfind.c, run on the reference FASTA), a seed-based tandem-repeat finder for periods
  7-2,000 bp.
  - Every exact 10-mer proposes the distance to its previous copy as a period.
  - Seeds of about the same period chain into a repeat.
  - A repeat is reported when it has at least about two copies.
  - The reported period is the smallest divisor of the chained period that repeats as well.
  - A repeat becomes VNTR at a period of 7 or more, and STR below.

**Merging.** Records on a contig that overlap or lie within `gap` bp of each other merge into one
region.
- **Class.** The class covering the most bp of the region, with ties broken VNTR, SAT, STR, LC.
- **Period and motif.** Those of the longest record of that class. trfind's motif is the reference
  sequence of its first period.

**Output.** A BED with a header line, gzipped when the name ends .gz, one region per line:

    chrom start end id class period motif n_members member_classes sources

- `id` is TR_<chrom>_<start>, stable as long as the region's start does not move.
- `period` is 0 when the motif is not a "(...)n" repeat, as for most satellites and low-complexity
  records.
- `member_classes` and `sources` count the members, e.g. "STR:2,VNTR:1" and "rm:2,trfind:1".
"""
import collections
import concurrent.futures
import gzip
import os
import re
import shutil
import subprocess
import tempfile

KEEP = {'Simple_repeat', 'Satellite', 'Low_complexity'}
MOTIF = re.compile(r'^\(([ACGTN]+)\)n$', re.IGNORECASE)
DEFAULT_CONTIGS = ['chr%d' % i for i in range(1, 23)] + ['chrX']
PRIORITY = {'VNTR': 0, 'SAT': 1, 'STR': 2, 'LC': 3}

Record = collections.namedtuple('Record', 'chrom start end cls period motif src')
Region = collections.namedtuple('Region',
                                'chrom start end id cls period motif n_members member_classes sources')


# ----------------------------------------------------------------------------- RepeatMasker

def classify(rm_class, name):
    """(class, period, motif) of one RepeatMasker record, or None for a class we do not keep."""
    if rm_class not in KEEP:
        return None
    m = MOTIF.match(name)
    motif = m.group(1).upper() if m else name
    period = len(motif) if m else 0
    if rm_class == 'Simple_repeat':
        if not m:
            return ('STR' if len(name) <= 6 else 'VNTR'), 0, name
        return ('STR' if period <= 6 else 'VNTR'), period, motif
    return ('SAT' if rm_class == 'Satellite' else 'LC'), period, motif


def _open(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


def read_repeatmasker(path, contigs=None):
    """Kept records of a RepeatMasker BED, as Records, in file order."""
    contigs = set(contigs) if contigs else None
    with _open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            x = line.rstrip('\n').split('\t')
            if len(x) < 7:
                raise ValueError('RepeatMasker BED needs at least 7 columns: %r' % line[:200])
            if contigs is not None and x[0] not in contigs:
                continue
            c = classify(x[6], x[3])
            if c is None:
                continue
            yield Record(x[0], int(x[1]), int(x[2]), *c, 'rm')


# ----------------------------------------------------------------------------- trfind

def trfind_binary(cache=None):
    """Compile c/trfind.c once into a cache directory and return the binary's path."""
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'c', 'trfind.c')
    cache = cache or os.environ.get('PGREALIGN_CACHE', os.path.expanduser('~/.cache/pgrealign'))
    os.makedirs(cache, exist_ok=True)
    out = os.path.join(cache, 'trfind')
    if os.path.exists(out) and os.path.getmtime(out) >= os.path.getmtime(src):
        return out
    cc = os.environ.get('CC') or shutil.which('cc') or shutil.which('clang') or shutil.which('gcc')
    if cc is None:
        raise RuntimeError('trfind needs a C compiler (set CC)')
    subprocess.run([cc, '-O2', '-o', out, src, '-lm'], check=True)
    return out


def fasta_records(path):
    """(name, sequence) of each record of a FASTA file (plain or gzipped)."""
    name, chunks = None, []
    with _open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, ''.join(chunks)


def _trfind_one(args):
    """Run trfind on one FASTA file and return its Records (the motif is cut from the sequence)."""
    binary, fasta, minp, maxp, contigs = args
    seqs = {n: s.upper() for n, s in fasta_records(fasta) if contigs is None or n in contigs}
    if not seqs:
        return []
    with tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False) as t:
        for n, s in seqs.items():
            t.write('>%s\n%s\n' % (n, s))
        tmp = t.name
    try:
        p = subprocess.run([binary, tmp, str(minp), str(maxp)], capture_output=True, text=True, check=True)
    finally:
        os.unlink(tmp)
    out = []
    for line in p.stdout.splitlines():
        x = line.split('\t')
        chrom, a, b, q = x[0], int(x[1]), int(x[2]), int(x[4])
        out.append(Record(chrom, a, b, 'VNTR' if q >= 7 else 'STR', q, seqs[chrom][a:a + q], 'trfind'))
    return out


def run_trfind(fastas, minp=7, maxp=2000, contigs=None, jobs=1):
    """trfind Records from one or more FASTA files (one contig per file runs them in parallel)."""
    binary = trfind_binary()
    contigs = set(contigs) if contigs else None
    recs = []
    with concurrent.futures.ThreadPoolExecutor(max(1, jobs)) as ex:
        for r in ex.map(_trfind_one, [(binary, f, minp, maxp, contigs) for f in fastas]):
            recs.extend(r)
    return recs


# ----------------------------------------------------------------------------- merge and write

def merge(records, gap=50):
    """Merge records of one contig that overlap or lie within `gap` bp, in coordinate order."""
    by_contig = collections.defaultdict(list)
    for r in records:
        by_contig[r.chrom].append(r)
    out = []
    for chrom in sorted(by_contig, key=_contig_key):
        rs = sorted(by_contig[chrom], key=lambda r: (r.start, r.end))
        group = [rs[0]]
        end = rs[0].end
        for r in rs[1:]:
            if r.start <= end + gap:
                group.append(r)
                end = max(end, r.end)
            else:
                out.append(_region(group, end))
                group, end = [r], r.end
        out.append(_region(group, end))
    return out


def _union_bp(ivs):
    tot, cs, ce = 0, None, None
    for s, e in sorted(ivs):
        if cs is None or s > ce:
            if cs is not None:
                tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + (ce - cs if cs is not None else 0)


def _region(group, end):
    bp = {k: _union_bp([(r.start, r.end) for r in group if r.cls == k]) for k in {r.cls for r in group}}
    cls = min(bp, key=lambda k: (-bp[k], PRIORITY[k]))
    longest = max((r for r in group if r.cls == cls), key=lambda r: r.end - r.start)  # first of equals
    counts = collections.Counter(r.cls for r in group)
    srcs = collections.Counter(r.src for r in group)
    start = min(r.start for r in group)
    return Region(group[0].chrom, start, end, 'TR_%s_%d' % (group[0].chrom, start), cls,
                  longest.period, longest.motif, len(group),
                  ','.join('%s:%d' % (k, counts[k]) for k in sorted(counts)),
                  ','.join('%s:%d' % (k, srcs[k]) for k in sorted(srcs)))


def _contig_key(c):
    m = re.match(r'^chr(\d+)$', c)
    return (0, int(m.group(1)), '') if m else (1, 0, c)


HEADER = '#chrom\tstart\tend\tid\tclass\tperiod\tmotif\tn_members\tmember_classes\tsources\n'


def write_catalog(regions, path):
    with _open(path, 'wt') as f:
        f.write(HEADER)
        for r in regions:
            f.write('\t'.join(str(v) for v in r) + '\n')


def read_catalog(path):
    """Regions of a catalogue written by write_catalog, or of a plain BED (chrom start end [id ...])."""
    with _open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            x = line.rstrip('\n').split('\t')
            chrom, start, end = x[0], int(x[1]), int(x[2])
            rid = x[3] if len(x) > 3 and x[3] else 'R_%s_%d' % (chrom, start)
            if len(x) >= 10:
                yield Region(chrom, start, end, rid, x[4], int(x[5]), x[6], int(x[7]), x[8], x[9])
            else:
                yield Region(chrom, start, end, rid, '.', 0, '.', 1, '.', '.')


def build(out, repeatmasker=None, fastas=(), gap=50, contigs=None, minp=7, maxp=2000, jobs=1):
    """Read the sources, merge and write; returns a summary dict."""
    recs = []
    if repeatmasker:
        recs.extend(read_repeatmasker(repeatmasker, contigs))
    if fastas:
        recs.extend(run_trfind(fastas, minp, maxp, contigs, jobs))
    if not recs:
        raise ValueError('no records: give a RepeatMasker BED and/or reference FASTA files')
    regions = merge(recs, gap)
    write_catalog(regions, out)
    return {'records': len(recs), 'regions': len(regions),
            'records_by_source_class': {'%s/%s' % k: v for k, v in
                                        sorted(collections.Counter((r.src, r.cls) for r in recs).items())},
            'regions_by_class': dict(collections.Counter(r.cls for r in regions)),
            'bp': sum(r.end - r.start for r in regions)}
