"""Stage 0, `catalog`: a tandem-repeat catalogue of the reference.

This is the only stage that knows the regions are tandem repeats. Any BED of target intervals can take
the catalogue's place in the stages after it (see Reading). The catalogue takes no sample or truth data
as input.

**Sources.**
- **RepeatMasker** (a BED of ten columns: chrom start end name score strand class family divergence
  id). Simple_repeat records become STR when the motif in "(MOTIF)n" is 1-6 bp and VNTR when it is
  longer. Satellite records become SAT and Low_complexity records LC. A record whose name is not
  "(MOTIF)n" has period 0 and its name as the motif, and is STR if it is a Simple_repeat.
  RepeatMasker reports no Simple_repeat motif longer than 30 bp, so on its own it misses the VNTRs of
  longer periods.
- **trfind** (c/trfind.c, run on the reference FASTA), a seed-based tandem-repeat finder for periods
  7-2,000 bp.
  - Every exact 10-mer proposes the distance to its previous copy as a period.
  - Seeds of about the same period chain into a repeat.
  - A repeat is reported when it has at least 3 seeds spanning about two copies. The header of
    c/trfind.c states the rules in full.
  - The reported period is the smallest divisor of the chained period that repeats as well.
  - A repeat becomes VNTR at a period of 7 or more, and STR below.
  - trfind is compiled on first use into the cache directory ($PGREALIGN_CACHE, else
    ~/.cache/pgrealign), under a key that hashes its source, the compiler flags and the platform
    (trfind_build).

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

**Reading.** read_catalog reads a catalogue, which it recognises by that header line ("#chrom" and
the columns above), or any other BED of target intervals (BED3, BED6, BED12, ...). Every region it
reads has a unique id.
"""
import collections
import concurrent.futures
import gzip
import hashlib
import itertools
import os
import platform
import re
import shlex
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
    """(class, period, motif) of one RepeatMasker record, or None for a class we do not keep. A name
    "(MOTIF)n" gives the period and motif, and any other name period 0 and the name as the motif. A
    Simple_repeat is VNTR at a period above 6 and STR otherwise, period 0 included."""
    if rm_class not in KEEP:
        return None
    m = MOTIF.match(name)
    period, motif = (len(m.group(1)), m.group(1).upper()) if m else (0, name)
    if rm_class == 'Simple_repeat':
        return ('VNTR' if period > 6 else 'STR'), period, motif
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

CC_FLAGS = ('-O2', '-lm')


def _trfind_source():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'c', 'trfind.c')


def trfind_build():
    """How trfind is built here: the sha256 of c/trfind.c, the compiler flags, the platform, and `key`,
    a sha256 of all of them that names the cached binary."""
    with open(_trfind_source(), 'rb') as f:
        b = {'source_sha256': hashlib.sha256(f.read()).hexdigest(), 'cc_flags': ' '.join(CC_FLAGS),
             'system': platform.system(), 'machine': platform.machine()}
    b['key'] = hashlib.sha256('\t'.join(b[k] for k in sorted(b)).encode()).hexdigest()
    return b


def trfind_binary(cache=None):
    """Path of the trfind binary for trfind_build()'s key, compiled on first use into the cache directory
    (`cache`, else $PGREALIGN_CACHE, else ~/.cache/pgrealign) with $CC, else cc, clang or gcc. It is
    compiled under a temporary name and renamed into place, so no run ever sees a partial binary."""
    cache = cache or os.environ.get('PGREALIGN_CACHE') or os.path.expanduser('~/.cache/pgrealign')
    out = os.path.join(cache, 'trfind-' + trfind_build()['key'][:16])
    if os.path.exists(out):
        return out
    cc = shlex.split(os.environ.get('CC', '')) or [shutil.which('cc') or shutil.which('clang') or
                                                   shutil.which('gcc')]
    if cc[0] is None:
        raise RuntimeError('trfind needs a C compiler (set CC)')
    os.makedirs(cache, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix='.trfind-', dir=cache)
    try:
        p = subprocess.run([*cc, '-o', os.path.join(tmp, 'trfind'), _trfind_source(), *CC_FLAGS],
                           capture_output=True, text=True)
        if p.returncode:
            raise RuntimeError('compiling trfind failed: %s\n%s' % (' '.join(p.args), p.stderr))
        os.replace(os.path.join(tmp, 'trfind'), out)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
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
    """Regions of a catalogue, or of any BED of target intervals, in file order.

    A file is a catalogue when a header line before its first data line is HEADER. Its lines are read
    as write_catalog wrote them, and its ids must be unique. Any other file is a BED: chrom, start, end
    and an optional name, with further columns (strand, BED12 blocks) ignored. A BED region has class
    '.', period 0, motif '.' and one member, and its id is
    - its name, when no other line has that name;
    - otherwise <name>_<chrom>_<start>, or R_<chrom>_<start> for a line without a name (none, empty or
      '.'), with _2, _3, ... appended, in file order, where that would repeat an id already given."""
    with _open(path) as f:
        lines = enumerate(f, 1)
        catalogue = False
        for n, line in lines:
            if line.startswith('#'):
                catalogue = catalogue or line.rstrip('\r\n') == HEADER.rstrip('\n')
            elif line.strip() and not line.startswith(('track', 'browser')):
                break
        else:
            return
        lines = itertools.chain([(n, line)], lines)
        if catalogue:
            yield from _unique(path, _parse(path, lines, _catalogue_line))
        else:
            yield from _bed_ids(list(_parse(path, lines, _bed_line)))


def _parse(path, lines, parse):
    """parse(columns) of each data line of (line number, line) pairs. A line that does not parse, or
    whose start and end are not 0 <= start <= end, raises a ValueError naming it."""
    for n, line in lines:
        if not line.strip() or line.startswith(('#', 'track', 'browser')):
            continue
        try:
            r = parse(line.rstrip('\r\n').split('\t'))
            if not 0 <= r.start <= r.end:
                raise ValueError('start %d and end %d are not 0 <= start <= end' % (r.start, r.end))
        except ValueError as e:
            raise ValueError('%s line %d: %s' % (path, n, e)) from None
        yield r


def _catalogue_line(x):
    if len(x) < 10:
        raise ValueError('a catalogue line has 10 columns, not %d' % len(x))
    return Region(x[0], int(x[1]), int(x[2]), x[3], x[4], int(x[5]), x[6], int(x[7]), x[8], x[9])


def _bed_line(x):
    """A BED line as a Region whose id is its name, or None when it has none."""
    if len(x) < 3:
        raise ValueError('a BED line has at least 3 columns, not %d' % len(x))
    name = x[3] if len(x) > 3 and x[3] not in ('', '.') else None
    return Region(x[0], int(x[1]), int(x[2]), name, '.', 0, '.', 1, '.', '.')


def _unique(path, regions):
    seen = set()
    for r in regions:
        if r.id in seen:
            raise ValueError('%s: the catalogue id %s repeats' % (path, r.id))
        seen.add(r.id)
        yield r


def _bed_ids(regions):
    """The BED regions with the ids that read_catalog defines."""
    names = collections.Counter(r.id for r in regions if r.id is not None)
    given = {k for k, c in names.items() if c == 1}
    for r in regions:
        if r.id is not None and names[r.id] == 1:
            yield r
            continue
        base = '%s_%s_%d' % (r.id or 'R', r.chrom, r.start)
        rid, k = base, 1
        while rid in given:
            k += 1
            rid = '%s_%d' % (base, k)
        given.add(rid)
        yield r._replace(id=rid)


def build(out, repeatmasker=None, fastas=(), gap=50, contigs=None, minp=7, maxp=2000, jobs=1):
    """Read the sources, merge and write; returns a summary dict, with trfind_build() as `trfind` when
    trfind ran."""
    recs = []
    if repeatmasker:
        recs.extend(read_repeatmasker(repeatmasker, contigs))
    if fastas:
        recs.extend(run_trfind(fastas, minp, maxp, contigs, jobs))
    if not recs:
        raise ValueError('no records: give a RepeatMasker BED and/or reference FASTA files')
    regions = merge(recs, gap)
    write_catalog(regions, out)
    summary = {'records': len(recs), 'regions': len(regions),
               'records_by_source_class': {'%s/%s' % k: v for k, v in
                                           sorted(collections.Counter((r.src, r.cls) for r in recs).items())},
               'regions_by_class': dict(collections.Counter(r.cls for r in regions)),
               'bp': sum(r.end - r.start for r in regions)}
    if fastas:
        summary['trfind'] = trfind_build()
    return summary
