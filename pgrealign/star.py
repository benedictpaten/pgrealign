"""The medoid star: one multiple alignment of a region's alleles and fragments.

This module is the alignment method on its own. It takes sequences and returns rows, and knows nothing
of graphs, files or the pipeline's stages, so other tools can call `star_msa` directly. realign.py
applies it to the region packages that extract.py writes.

**Input.**
- Alleles are a region's distinct whole sequences, each with a weight, the number of haplotypes that
  carry it. Only an allele of positive weight can be the centre.
- Fragments are sequences that cover only part of the region, each with a side. A `left` fragment
  starts at the region's left end, a `right` fragment ends at its right end, and an `internal`
  fragment touches neither end.

**Method.**
1. **Centre.** The centre is the allele of positive weight with the least weighted sum of k-mer
   distances to all the alleles. The distance between two alleles is one minus the Jaccard index of
   their k-mer multisets: the copies of each k-mer the two share, over the copies in either. Ties go
   to the earlier allele.
2. **Pairwise alignment.** Every other allele is aligned to the centre with a global, unbanded
   alignment. An allele identical to the centre takes the centre's row; an empty allele, or an empty
   centre, is aligned as one deletion or insertion.
3. **Normalisation.** Each pair's indels are shifted left wherever the alignment's score stays the
   same.
4. **Merge.** The alignment has one column per centre base, and an insertion slot between each two
   centre bases and at both ends. The distinct insertions that fall in one slot are aligned to each
   other, longest first, with a banded alignment.
5. **Fragments.** A fragment that occurs in an allele, at the end its side names or anywhere if it is
   internal, takes that allele's row, cut to the matching bases. The first such allele in input order
   is used. Any other fragment is aligned to the centre like an allele and joins the merge.

Every row spells its input exactly, and every row has the same length; `star_msa` checks both.

**Alphabet.** The aligner sees A, C, G, T and N only. Any other character is aligned as N and restored
in the returned rows.
"""
import collections
import os
import shutil
import subprocess
import tempfile
import time

K = 15
_CODE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
_ALIGNABLE = frozenset('ACGTN')


class RealignError(Exception):
    """A region that could not be aligned. `status` names the reason in realign's output."""
    status = 'error'

    def __init__(self, message):
        super().__init__(message)


class TooBig(RealignError):
    status = 'too_big'


class TimedOut(RealignError):
    status = 'timeout'


class AlignerFailed(RealignError):
    status = 'aligner_failed'


# ----------------------------------------------------------------------------- scores

class Scores:
    """abPOA's scores: match and mismatch, and two-piece gap open and extension penalties (o1, o2) and
    (e1, e2). A gap of length n costs min(o1 + n * e1, o2 + n * e2)."""
    __slots__ = ('match', 'mismatch', 'gap_open', 'gap_ext')

    def __init__(self, match=2, mismatch=4, gap_open=(4, 24), gap_ext=(2, 1)):
        self.match, self.mismatch = int(match), int(mismatch)
        self.gap_open = tuple(int(x) for x in gap_open)
        self.gap_ext = tuple(int(x) for x in gap_ext)

    def flags(self):
        """abPOA's command-line flags for these scores."""
        return ['-M', str(self.match), '-X', str(self.mismatch),
                '-O', '%d,%d' % self.gap_open, '-E', '%d,%d' % self.gap_ext]

    @classmethod
    def parse(cls, text):
        """Scores from abPOA flags, e.g. '-M 8 -X 16 -O 16,96 -E 8,1'. A flag left out takes abPOA's
        default, so '' gives abPOA's defaults."""
        s, toks = cls(), text.split()
        if len(toks) % 2:
            raise ValueError('scores need a value after every flag: %r' % text)
        for flag, value in zip(toks[0::2], toks[1::2]):
            if flag == '-M':
                s.match = int(value)
            elif flag == '-X':
                s.mismatch = int(value)
            elif flag in ('-O', '-E'):
                pair = tuple(int(x) for x in value.split(','))
                if len(pair) != 2:
                    raise ValueError('%s takes two values, as in 4,24: %r' % (flag, value))
                if flag == '-O':
                    s.gap_open = pair
                else:
                    s.gap_ext = pair
            else:
                raise ValueError('unknown scoring flag %r (use -M, -X, -O, -E)' % flag)
        return s

    def __eq__(self, other):
        return isinstance(other, Scores) and self.flags() == other.flags()

    def __repr__(self):
        return 'Scores(match=%d, mismatch=%d, gap_open=%r, gap_ext=%r)' % (
            self.match, self.mismatch, self.gap_open, self.gap_ext)


# abPOA's defaults times four, except the long-gap extension, which stays 1: a long indel then costs a
# quarter as much per base, relative to a match. Under the defaults, two unrelated stretches are
# cheaper to align base against base, threading short gaps through chance matches, than to delete one
# and insert the other, so one replacement becomes many small bubbles.
SCORES = Scores(match=8, mismatch=16, gap_open=(16, 96), gap_ext=(8, 1))


# ----------------------------------------------------------------------------- the aligner

def read_fasta_rows(path):
    out, name, chunks = collections.OrderedDict(), None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    out[name] = ''.join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        out[name] = ''.join(chunks)
    return out


class Abpoa:
    """The default aligner: the abpoa program, run once per alignment.

    Use it as a context manager, which owns a temporary directory. Calling it with a list of sequences
    returns their global multiple alignment, one upper-case row per sequence in input order:
    `abpoa -m 0 -r 1 <scores>`, with `-b -1` (no band) unless `banded` is true, in which case abPOA's
    adaptive band applies. Past `deadline` (a time.time() value), a call raises TimedOut."""

    def __init__(self, exe='abpoa', scores=SCORES, deadline=None, workdir=None):
        self.exe, self.scores, self.deadline, self.workdir = exe, scores, deadline, workdir
        self.dir = None

    def __enter__(self):
        self.dir = tempfile.mkdtemp(prefix='realign.', dir=self.workdir)
        return self

    def __exit__(self, *exc):
        shutil.rmtree(self.dir, ignore_errors=True)
        self.dir = None
        return False

    def flags(self, seqs, banded):
        """abPOA's flags for aligning `seqs`, after `-m 0 -r 1`."""
        return self.scores.flags() + ([] if banded else ['-b', '-1'])

    def __call__(self, seqs, banded=False):
        timeout = None
        if self.deadline is not None:
            timeout = self.deadline - time.time()
            if timeout <= 0:
                raise TimedOut('time limit reached')
        fa, out = os.path.join(self.dir, 'in.fa'), os.path.join(self.dir, 'out.fa')
        with open(fa, 'w') as f:
            for i, s in enumerate(seqs):
                f.write('>s%d\n%s\n' % (i, s))
        flags = self.flags(seqs, banded)
        try:
            with open(out, 'w') as o:
                subprocess.run([self.exe, '-m', '0', '-r', '1', *flags, fa], stdout=o,
                               stderr=subprocess.PIPE, check=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            raise TimedOut('time limit reached in abPOA')
        except subprocess.CalledProcessError as e:
            raise AlignerFailed('abPOA exited with %d: %s' % (
                e.returncode, (e.stderr or b'').decode(errors='replace').strip()[-300:]))
        rows = read_fasta_rows(out)
        try:
            return [rows['s%d' % i].upper() for i in range(len(seqs))]
        except KeyError:
            raise AlignerFailed('abPOA returned %d rows for %d sequences' % (len(rows), len(seqs)))


# ----------------------------------------------------------------------------- centre

def kmer_counts(seq, k=K):
    """The multiset of k-mers, {k-mer code: copies}. k-mers holding a base other than A, C, G or T are
    skipped."""
    out = collections.Counter()
    mask = (1 << (2 * k)) - 1
    code, valid = 0, 0
    for c in seq:
        v = _CODE.get(c)
        if v is None:
            code, valid = 0, 0
            continue
        code = ((code << 2) | v) & mask
        valid += 1
        if valid >= k:
            out[code] += 1
    return out


def medoid(seqs, weights, k=K, deadline=None):
    """Index of the sequence minimising sum_j w_j (1 - Jaccard(i, j)) over all the others; ties go to the
    lower index. Jaccard is over k-mer multisets: the copies of each k-mer the two share, over the
    copies in either. Only sequences of positive weight are candidates, when there are any. Past
    `deadline`, raises TimedOut."""
    n = len(seqs)
    counts = [kmer_counts(s, k) for s in seqs]
    sizes = [sum(c.values()) for c in counts]
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        if deadline is not None and time.time() > deadline:
            raise TimedOut('time limit reached choosing the centre')
        for j in range(i + 1, n):
            a, b = counts[i], counts[j]
            inter = sum(min(a[x], b[x]) for x in a.keys() & b.keys())
            uni = sizes[i] + sizes[j] - inter
            dist[i][j] = dist[j][i] = 1.0 - (inter / uni if uni else 1.0)
    cands = [i for i in range(n) if weights[i] > 0] or list(range(n))
    best = None
    for i in cands:
        d = 0.0
        for j in range(n):
            if i != j:
                d += weights[j] * dist[i][j]
        if best is None or d < best[0]:
            best = (d, i)
    return best[1]


# ----------------------------------------------------------------------------- pieces of the star

def left_normalise_pair(a, b):
    """Shift the indels of a two-row alignment left wherever its score stays the same, until nothing
    moves."""
    a, b = list(a), list(b)
    n = len(a)
    moved = True
    while moved:
        moved = False
        for g, o in ((a, b), (b, a)):
            c = 0
            while c < n:
                if g[c] != '-' or o[c] == '-':
                    c += 1
                    continue
                e = c
                while e < n and g[e] == '-' and o[e] != '-':
                    e += 1
                end, s = e, c
                while s > 0 and g[s - 1] != '-' and o[s - 1] != '-' and o[s - 1] == o[e - 1]:
                    g[e - 1], g[s - 1] = g[s - 1], '-'
                    s, e = s - 1, e - 1
                    moved = True
                c = end
    return ''.join(a), ''.join(b)


def align_pair(centre, seq, align):
    """The left-normalised global alignment of `seq` to `centre`, as (centre row, row)."""
    if seq == centre:
        return centre, seq
    if not centre or not seq:
        return centre + '-' * len(seq), '-' * len(centre) + seq
    ca, sa = align([centre, seq], banded=False)
    keep = [i for i in range(len(ca)) if ca[i] != '-' or sa[i] != '-']
    ca, sa = ''.join(ca[i] for i in keep), ''.join(sa[i] for i in keep)
    if ca.replace('-', '') != centre or sa.replace('-', '') != seq:
        raise AlignerFailed('the aligner\'s rows do not spell their inputs')
    return left_normalise_pair(ca, sa)


def star_merge(centre, pairs, align):
    """Merge pairwise alignments to the centre, {name: (centre row, row)}, into one alignment
    {name: row}, and the number of slot alignments it ran."""
    n = len(centre)
    cols, ins = {}, {}
    for k, (ca, sa) in pairs.items():
        col, slot, p = ['-'] * n, collections.defaultdict(list), 0
        for x, y in zip(ca, sa):
            if x == '-':
                if y != '-':
                    slot[p].append(y)
            else:
                col[p] = y
                p += 1
        if p != n:
            raise ValueError('pairwise centre row of %s spells %d of %d bases' % (k, p, n))
        cols[k], ins[k] = col, {q: ''.join(v) for q, v in slot.items()}
    byslot = collections.defaultdict(set)
    for d in ins.values():
        for q, s in d.items():
            byslot[q].add(s)
    slot_rows, slot_alignments = {}, 0
    for q, strs in byslot.items():
        if len(strs) == 1:
            s = next(iter(strs))
            slot_rows[q] = (len(s), {s: s})
            continue
        recs = sorted(strs, key=lambda x: (-len(x), x))
        got = align(recs, banded=True)
        al = dict(zip(recs, got))
        slot_rows[q] = (len(al[recs[0]]), al)
        slot_alignments += 1
    out = {}
    for k in pairs:
        parts = []
        for q in range(n + 1):
            if q in slot_rows:
                w, al = slot_rows[q]
                s = ins[k].get(q)
                parts.append(al[s] if s else '-' * w)
            if q < n:
                parts.append(cols[k][q])
        out[k] = ''.join(parts)
    return out, slot_alignments


def drop_gap_columns(rows):
    """The alignment without its all-gap columns."""
    names = list(rows)
    if not names:
        return rows
    L = len(rows[names[0]])
    keep = [c for c in range(L) if any(rows[n][c] != '-' for n in names)]
    return {n: ''.join(rows[n][c] for c in keep) for n in names}


def place_fragment(frag, side, alleles, rows):
    """A row for a fragment, cut from the row of the first allele that contains it, or None.
    A left fragment must be a prefix of the allele, a right fragment a suffix of it, and an internal
    fragment anywhere in it, at its leftmost occurrence. `alleles` is [(id, sequence)]."""
    for aid, seq in alleles:
        if side == 'left':
            pos = 0 if seq.startswith(frag) else -1
        elif side == 'right':
            pos = len(seq) - len(frag) if seq.endswith(frag) else -1
        else:
            pos = seq.find(frag)
        if pos < 0:
            continue
        row, out, k = rows[aid], [], 0
        for ch in row:
            if ch == '-':
                out.append('-')
            else:
                out.append(ch if pos <= k < pos + len(frag) else '-')
                k += 1
        return ''.join(out)
    return None


def respell(row, seq):
    """`row` with its bases replaced, in order, by those of `seq`."""
    it = iter(seq)
    return ''.join('-' if c == '-' else next(it) for c in row)


# ----------------------------------------------------------------------------- the star

class StarMSA:
    """The result of star_msa: `rows` {id: row}, every row the same length; `centre`, the centre allele's
    id; `slot_alignments`, the number of insertion slots that needed an alignment; `fragments_exact`
    and `fragments_aligned`, the fragments placed in an allele and aligned to the centre."""
    __slots__ = ('rows', 'centre', 'slot_alignments', 'fragments_exact', 'fragments_aligned')

    def __init__(self, rows, centre, slot_alignments, fragments_exact, fragments_aligned):
        self.rows, self.centre = rows, centre
        self.slot_alignments = slot_alignments
        self.fragments_exact, self.fragments_aligned = fragments_exact, fragments_aligned


def star_msa(alleles, fragments=(), align=None, centre=None, k=K, deadline=None, max_cells=None):
    """The medoid star of a region (see the module docstring).

    `alleles` is a sequence of (id, sequence, weight), with distinct ids; `fragments` is a sequence of
    (id, sequence, side), side being 'left', 'right' or 'internal'. `align` is the aligner: called with
    a list of sequences and `banded`, it returns their global alignment, one row per sequence in order.
    Pass an Abpoa inside its `with` block. `centre` names the centre allele instead of the medoid.
    Past `deadline` (a time.time() value), raises TimedOut. If the centre's length times that of the
    longest sequence aligned to it exceeds `max_cells`, raises TooBig before aligning anything.
    Returns a StarMSA."""
    if align is None:
        raise ValueError('star_msa needs an aligner, e.g. `with Abpoa() as align`')
    if not alleles:
        raise RealignError('no allele to align')
    ids = [a[0] for a in alleles] + [f[0] for f in fragments]
    if len(set(ids)) != len(ids):
        raise ValueError('allele and fragment ids must be distinct')
    for f in fragments:
        if f[2] not in ('left', 'right', 'internal'):
            raise ValueError('fragment %s has side %r' % (f[0], f[2]))
    orig = {m: s.upper() for m, s in [(a[0], a[1]) for a in alleles] + [(f[0], f[1]) for f in fragments]}
    seq = {m: s if set(s) <= _ALIGNABLE else ''.join(c if c in _ALIGNABLE else 'N' for c in s)
           for m, s in orig.items()}
    names = [a[0] for a in alleles]
    if centre is None:
        ci = medoid([seq[m] for m in names], [a[2] for a in alleles], k, deadline)
        cid = names[ci]
    elif centre in seq and centre in names:
        cid = centre
    else:
        raise ValueError('centre %r is not an allele' % (centre,))
    cseq = seq[cid]
    found = [(m, orig[m]) for m in names]
    exact, unplaced = [], []
    for fid, _, side in fragments:
        if place_fragment(orig[fid], side, found, orig) is not None:
            exact.append((fid, side))
        else:
            unplaced.append(fid)
    longest = max(len(seq[m]) for m in names + unplaced)
    if max_cells is not None and max(1, len(cseq)) * max(1, longest) > max_cells:
        raise TooBig('centre %d bp x longest %d bp > %d cells' % (len(cseq), longest, max_cells))
    pairs = collections.OrderedDict((m, align_pair(cseq, seq[m], align)) for m in names + unplaced)
    rows, slot_alignments = star_merge(cseq, pairs, align)
    for fid, side in exact:
        rows[fid] = place_fragment(orig[fid], side, found, rows)
    rows = drop_gap_columns(rows)
    rows = {m: respell(r, orig[m]) if seq[m] is not orig[m] else r for m, r in rows.items()}
    widths = {len(r) for r in rows.values()}
    if len(widths) != 1:
        raise ValueError('rows of unequal length')
    for m, r in rows.items():
        if r.replace('-', '') != orig[m]:
            raise ValueError('row %s does not spell its sequence' % m)
    return StarMSA(rows, cid, slot_alignments, len(exact), len(unplaced))
