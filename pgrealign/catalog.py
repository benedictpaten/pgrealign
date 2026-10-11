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
    (tool_build), unless `--trfind` names a binary built from c/trfind.c elsewhere, as where there is
    no C compiler (prebuilt_tool).
- **TRF** (Tandem Repeats Finder, `--trf`, run on the reference FASTA), the way RepeatMasker 4.1 finds
  its Simple_repeat records, for a reference with no RepeatMasker annotation. Each sequence is cut into
  RepeatMasker's batches of 60 kb that overlap by 2 kb, which also keeps TRF from spending hours on a
  centromere's satellites, and each batch goes through RepeatMasker's two TRF stages (TRF_STAGES): one
  for young repeats and one for old. A repeat of more copies than the stage's minimum becomes a
  Simple_repeat record named by TRF's consensus motif, so STR or VNTR by the rule above. RepeatMasker
  then rescores each repeat's alignment and masks interspersed repeats first; catalog does neither, so
  TRF finds about every Simple_repeat RepeatMasker would, and more besides.
- **strf** (c/strf.c, `--strf`), pgrealign's own reimplementation of the part of TRF that catalog
  uses, run in TRF's place on the same batches and stages. It takes TRF's command line and writes
  TRF's -ngs output, for periods up to 20 and PM 75 or 80 only, and finds nearly the same records
  (docs/method.md). It is compiled on first use as trfind is, unless `--strf` names a binary built from
  c/strf.c elsewhere.

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


# ----------------------------------------------------------------------------- the C tools

CC_FLAGS = ('-O2', '-lm')
C_TOOLS = ('trfind', 'strf')   # c/trfind.c and c/strf.c


def _c_source(name):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'c', name + '.c')


def tool_build(name):
    """How the C tool `name` (C_TOOLS) is built here: the sha256 of c/<name>.c, the compiler flags, the
    platform, and `key`, a sha256 of all of them that names the cached binary."""
    with open(_c_source(name), 'rb') as f:
        b = {'source_sha256': hashlib.sha256(f.read()).hexdigest(), 'cc_flags': ' '.join(CC_FLAGS),
             'system': platform.system(), 'machine': platform.machine()}
    b['key'] = hashlib.sha256('\t'.join(b[k] for k in sorted(b)).encode()).hexdigest()
    return b


def tool_binary(name, cache=None):
    """Path of the C tool `name`'s binary for tool_build(name)'s key, compiled on first use into the cache
    directory (`cache`, else $PGREALIGN_CACHE, else ~/.cache/pgrealign) with $CC, else cc, clang or gcc.
    It is compiled under a temporary name and renamed into place, so no run ever sees a partial binary."""
    cache = cache or os.environ.get('PGREALIGN_CACHE') or os.path.expanduser('~/.cache/pgrealign')
    out = os.path.join(cache, '%s-%s' % (name, tool_build(name)['key'][:16]))
    if os.path.exists(out):
        return out
    cc = shlex.split(os.environ.get('CC', '')) or [shutil.which('cc') or shutil.which('clang') or
                                                   shutil.which('gcc')]
    if cc[0] is None:
        raise RuntimeError('%s needs a C compiler (set CC)' % name)
    os.makedirs(cache, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix='.%s-' % name, dir=cache)
    try:
        p = subprocess.run([*cc, '-o', os.path.join(tmp, name), _c_source(name), *CC_FLAGS],
                           capture_output=True, text=True)
        if p.returncode:
            raise RuntimeError('compiling %s failed: %s\n%s' % (name, ' '.join(p.args), p.stderr))
        os.replace(os.path.join(tmp, name), out)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return out


def prebuilt_tool(name, binary):
    """The path of a binary of the C tool `name` built elsewhere from c/<name>.c, such as by a package
    that has no C compiler where it runs. `binary` is a path, or a name to find on PATH."""
    path = shutil.which(binary)
    if path is None:
        raise RuntimeError('%s binary %s not found, or not executable' % (name, binary))
    return os.path.abspath(path)


# ----------------------------------------------------------------------------- trfind


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


def run_trfind(fastas, minp=7, maxp=2000, contigs=None, jobs=1, binary=None):
    """trfind Records from one or more FASTA files (one contig per file runs them in parallel). `binary` is
    a prebuilt trfind (prebuilt_tool), or None to compile one on first use (tool_binary)."""
    binary = prebuilt_tool('trfind', binary) if binary else tool_binary('trfind')
    contigs = set(contigs) if contigs else None
    recs = []
    with concurrent.futures.ThreadPoolExecutor(max(1, jobs)) as ex:
        for r in ex.map(_trfind_one, [(binary, f, minp, maxp, contigs) for f in fastas]):
            recs.extend(r)
    return recs


# ----------------------------------------------------------------------------- TRF and strf

# RepeatMasker 4.1's two TRF stages (runTRFStage in its RepeatMasker script): the parameters it gives
# TRF, as match, mismatch, delta, PM, PI, minscore and maxperiod, and the copy number a repeat must
# exceed to be kept
TRF_STAGES = (('perfect', (2, 7, 7, 80, 10, 50, 10), 4),
              ('old', (2, 3, 5, 75, 20, 33, 7), 5))
TRF_BATCH, TRF_OVERLAP = 60000, 2000     # RepeatMasker's batches: their length, and their overlap


def trf_batches(seqs, size=None, overlap=None):
    """(contig, offset, sequence) of the batches that cover each of `seqs` {contig: sequence}: windows of
    `size` bp (TRF_BATCH), each starting `size - overlap` bp (TRF_OVERLAP) after the one before, the last
    reaching the end."""
    size = size or TRF_BATCH
    overlap = TRF_OVERLAP if overlap is None else overlap
    for contig, s in seqs.items():
        start = 0
        while True:
            yield contig, start, s[start:start + size]
            if start + size >= len(s):
                break
            start += size - overlap


def parse_trf(lines, batches, min_copies, src='trf'):
    """Records from TRF's -ngs output, `lines`: an '@<name>' line for each sequence with repeats, then a
    line per repeat: start and end (1-based, inclusive), period, copy number, consensus size, percent
    matches, percent indels, score, A, C, G, T, entropy, consensus, sequence, and the two flanks. A
    sequence's name is its index in `batches`, [(contig, offset)], which places it on its contig. A
    repeat of `min_copies` copies or fewer is dropped. `src` is the Records' source, 'trf' or 'strf'."""
    out, where = [], None
    for line in lines:
        if line.startswith('@'):
            where = batches[int(line[1:].split()[0])]
            continue
        x = line.split()
        if not x or float(x[3]) <= min_copies:
            continue
        cls, period, motif = classify('Simple_repeat', '(%s)n' % x[13])
        out.append(Record(where[0], where[1] + int(x[0]) - 1, where[1] + int(x[1]), cls, period, motif, src))
    return out


def _trf_one(args):
    """Run one TRF stage on one FASTA file of batches, and return its Records."""
    binary, fasta, params, min_copies, batches, src = args
    p = subprocess.run([binary, fasta, *map(str, params), '-h', '-ngs'], capture_output=True, text=True,
                       cwd=os.path.dirname(fasta))
    if p.returncode:
        raise RuntimeError('%s failed (%d) on %s: %s' % (binary, p.returncode, fasta, p.stderr[-2000:]))
    return parse_trf(p.stdout.splitlines(), batches, min_copies, src)


def run_trf(fastas, binary, contigs=None, jobs=1, src='trf'):
    """TRF Records from one or more FASTA files (see the module docstring), with the batches of every file
    shared out over `jobs` TRF processes. A repeat in the overlap of two batches can be found
    in each, as two records that merge. `binary` is a path, or a name to find on PATH, of TRF, or of strf
    with `src` 'strf'."""
    exe = shutil.which(binary)
    if exe is None:
        raise RuntimeError('TRF binary %s not found, or not executable' % binary)
    contigs = set(contigs) if contigs else None
    batches = []
    for f in fastas:
        seqs = {n: s.upper() for n, s in fasta_records(f) if contigs is None or n in contigs}
        batches.extend(trf_batches(seqs))
    n = max(1, min(jobs, len(batches)))
    places = [b[:2] for b in batches]
    tmp = tempfile.mkdtemp(prefix='pgrealign-trf.')
    try:
        work = []
        for k in range(n):
            # every n-th batch, so the slow ones, which cluster in a centromere, are shared out
            path = os.path.join(tmp, 'batches_%d.fa' % k)
            with open(path, 'w') as t:
                for i in range(k, len(batches), n):
                    t.write('>%d\n%s\n' % (i, batches[i][2]))
            work.extend((exe, path, params, min_copies, places, src) for _, params, min_copies in TRF_STAGES)
        recs = []
        with concurrent.futures.ThreadPoolExecutor(n) as ex:
            for r in ex.map(_trf_one, work):
                recs.extend(r)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return sorted(set(recs))


def run_strf(fastas, binary=None, contigs=None, jobs=1):
    """strf Records from one or more FASTA files: run_trf with strf in TRF's place. `binary` is a prebuilt
    strf (prebuilt_tool), or None to compile one on first use (tool_binary)."""
    binary = prebuilt_tool('strf', binary) if binary else tool_binary('strf')
    return run_trf(fastas, binary, contigs, jobs, src='strf')


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


def build(out, repeatmasker=None, fastas=(), gap=50, contigs=None, minp=7, maxp=2000, jobs=1, trfind=None,
          trf=None, strf=None):
    """Read the sources, merge and write; returns a summary dict, with `trfind` when trfind ran: its
    tool_build(), or {'binary': its path} for a prebuilt `trfind` (run_trfind). `trf`, a TRF binary,
    runs TRF on the `fastas` too (run_trf), and the summary's `trf` gives its path. `strf`, True or a
    prebuilt strf binary, runs strf in TRF's place (run_strf), and the summary's `strf` gives its build or
    path; each also gives the stages' parameters. Sources that hold no record on the selected contigs,
    such as a reference with no tandem repeat, give an empty catalogue."""
    if not repeatmasker and not fastas:
        raise ValueError('no sources: give a RepeatMasker BED and/or reference FASTA files')
    if trf and strf:
        raise ValueError('give TRF or strf, not both: strf finds what TRF does')
    if (trf or strf) and not fastas:
        raise ValueError('%s needs reference FASTA files to run on' % ('TRF' if trf else 'strf'))
    recs = []
    if repeatmasker:
        recs.extend(read_repeatmasker(repeatmasker, contigs))
    if fastas:
        recs.extend(run_trfind(fastas, minp, maxp, contigs, jobs, binary=trfind))
    if trf:
        recs.extend(run_trf(fastas, trf, contigs, jobs))
    if strf:
        recs.extend(run_strf(fastas, None if strf is True else strf, contigs, jobs))
    regions = merge(recs, gap)
    write_catalog(regions, out)
    summary = {'records': len(recs), 'regions': len(regions),
               'records_by_source_class': {'%s/%s' % k: v for k, v in
                                           sorted(collections.Counter((r.src, r.cls) for r in recs).items())},
               'regions_by_class': dict(collections.Counter(r.cls for r in regions)),
               'bp': sum(r.end - r.start for r in regions)}
    stages = {name: list(params) for name, params, _ in TRF_STAGES}
    if fastas:
        summary['trfind'] = {'binary': prebuilt_tool('trfind', trfind)} if trfind else tool_build('trfind')
    if trf:
        summary['trf'] = {'binary': os.path.abspath(shutil.which(trf)), 'stages': stages}
    if strf:
        summary['strf'] = {**(tool_build('strf') if strf is True else {'binary': prebuilt_tool('strf', strf)}),
                           'stages': stages}
    return summary
