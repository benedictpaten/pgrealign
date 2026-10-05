"""Stage 3, `realign`: one multiple alignment per region, from its package (extract.py).

The default method is the **medoid star**, chosen on chr20 (vntr-eval results, section 4x):

1. **Centre.** The allele with the least weighted sum of multiset 15-mer Jaccard distances to all the
   others. The weight is the allele's number of paths. The multiset counts the j-th copy of a k-mer
   as its own token, so copy number counts. k-mers are encoded as integers, so the centre does not
   depend on Python's per-process string hashing.
2. **Pairwise alignment.** Every other allele is aligned to the centre with
   `abpoa -m 0 -r 1 -b -1 SCORES`: global and unbanded, with the scores in SCORES.
3. **Normalisation.** Each pair's indels are shifted left wherever the cost stays the same.
4. **Merge.** One column per centre base, with an insertion slot between each pair of centre bases and
   at both ends. The distinct insertions of a slot are aligned to each other by abPOA, longest first,
   with the same scores.

**Fragments.** A fragment is a walk that starts or ends inside the region.
- Its side (extract.py) says which part of an allele it is: `left` touches anchor L, `right` touches
  R, and `internal` touches neither.
- A fragment is placed at the leftmost occurrence of its sequence, in the right position, in an
  allele that contains it. Its row is that allele's row, cut to the matching bases.
- A fragment found in no allele (on chr6, 433 of 593) is aligned to the centre like an allele, and
  joins the merge. abPOA's global alignment keeps every base, and the walk's missing end becomes one
  end gap.
- Fragments do not take part in choosing the centre.

**Limits.**
- A region is `too_big` when its centre times its longest allele exceeds `max_cells`. Unbanded abPOA
  needs about 15 bytes per cell; 22 kb x 22 kb fits 12 GB.
- A region is `timeout` when it runs longer than `timeout` seconds.
- Either way it is left as it is.

Output: `msas.jsonl.gz`, one object per region:
- `id`, `status`, `centre` (the allele id), `seconds`;
- `rows`: {member id: MSA row}, one per allele and fragment, every row the same length.
"""
import collections
import gzip
import json
import os
import shutil
import subprocess
import tempfile
import time

K = 15
MAX_CELLS = 450_000_000
MAX_TOTAL_BP = 40_000_000     # the medoid holds every allele's k-mer counts: up to about 50 bytes per base
# abPOA's scores, for the pairwise and the slot alignments alike: its defaults (match 2, mismatch 4, gap
# open 4 and 24, extension 2 and 1) times four, except the long-gap extension, which stays 1. A long
# indel then costs a quarter as much per base, relative to a match. Under the defaults, unrelated
# sequence is cheaper to align base against base, threading short gaps through chance matches, than to
# delete and insert; one replacement then becomes many small bubbles.
SCORES = ('-M', '8', '-X', '16', '-O', '16,96', '-E', '8,1')
_CODE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


class Capped(Exception):
    def __init__(self, status, message):
        super().__init__(message)
        self.status = status


# ----------------------------------------------------------------------------- centre

def kmer_counts(seq, k=K):
    """The multiset of k-mers, {k-mer code: copies}. k-mers holding a base other than ACGT are skipped."""
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


def medoid(seqs, weights, k=K):
    """Index of the sequence minimising sum_j w_j (1 - Jaccard(i, j)) over all the others; ties go to the
    lower index. Jaccard is over multisets: the shared copies of each k-mer over the copies in either.
    Only sequences of positive weight (carried by some path of the panel) are candidates, when there
    are any."""
    n = len(seqs)
    counts = [kmer_counts(s, k) for s in seqs]
    sizes = [sum(c.values()) for c in counts]
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
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


# ----------------------------------------------------------------------------- abPOA

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


def abpoa_rows(recs, wd, flags, deadline, abpoa='abpoa', scores=SCORES):
    """`abpoa -m 0 -r 1 SCORES FLAGS` on [(name, seq)] in that order -> {name: upper-case row}."""
    left = deadline - time.time()
    if left <= 0:
        raise Capped('timeout', 'time limit reached')
    fa, out = os.path.join(wd, 'in.fa'), os.path.join(wd, 'out.fa')
    with open(fa, 'w') as f:
        for n, s in recs:
            f.write('>%s\n%s\n' % (n, s))
    try:
        with open(out, 'w') as o:
            subprocess.run([abpoa, '-m', '0', '-r', '1', *scores, *flags, fa], stdout=o,
                           stderr=subprocess.DEVNULL, check=True, timeout=left)
    except subprocess.TimeoutExpired:
        raise Capped('timeout', 'time limit reached in abPOA')
    rows = read_fasta_rows(out)
    return {n: r.upper() for n, r in rows.items()}


# ----------------------------------------------------------------------------- star

def left_normalise_pair(a, b):
    """Left-normalise the indels of a two-row alignment without changing its cost (to a fixed point)."""
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


def star_merge(centre, pairs, wd, deadline, stats, abpoa='abpoa', scores=SCORES):
    """Merge pairwise alignments to the centre, {name: (centre row, row)}, into one MSA {name: row}."""
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
    slot_rows = {}
    for q, strs in byslot.items():
        if len(strs) == 1:
            s = next(iter(strs))
            slot_rows[q] = (len(s), {s: s})
            continue
        recs = sorted(strs, key=lambda x: (-len(x), x))
        got = abpoa_rows([('i%d' % i, s) for i, s in enumerate(recs)], wd, [], deadline, abpoa, scores)
        al = {s: got['i%d' % i] for i, s in enumerate(recs)}
        slot_rows[q] = (len(al[recs[0]]), al)
        stats['slot_alignments'] += 1
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
    return out


def drop_gap_columns(rows):
    names = list(rows)
    if not names:
        return rows
    L = len(rows[names[0]])
    keep = [c for c in range(L) if any(rows[n][c] != '-' for n in names)]
    return {n: ''.join(rows[n][c] for c in keep) for n in names}


# ----------------------------------------------------------------------------- fragments

def place_fragment(frag, side, alleles, rows):
    """A row for a fragment cut from the row of an allele that contains it, or None.
    A left fragment must be a prefix of the allele, a right fragment a suffix of it, and an internal
    fragment anywhere in it. Alleles are tried in the given order."""
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


# ----------------------------------------------------------------------------- one region

def realign_package(pkg, abpoa='abpoa', timeout=900, max_cells=MAX_CELLS, workdir=None, scores=SCORES):
    """The medoid-star MSA of one region package -> {'id', 'status', 'centre', 'rows', 'seconds', ...}."""
    t0 = time.time()
    deadline = t0 + timeout
    res = {'id': pkg['id'], 'method': 'medoid_star'}
    alleles = [(a['id'], a['seq'].upper()) for a in pkg['alleles']]
    weights = [a['weight'] for a in pkg['alleles']]
    stats = collections.Counter()
    wd = tempfile.mkdtemp(prefix='realign.', dir=workdir)
    try:
        if not alleles:
            raise Capped('no_alleles', 'no spanning allele to align')
        longest = max(len(s) for _, s in alleles)
        shortest = min(len(s) for (_, s), w in zip(alleles, weights) if w > 0) if any(weights) else 0
        total = sum(len(s) for _, s in alleles)
        # the centre is at least the shortest candidate, so this bounds the cells from below
        if max(1, shortest) * max(1, longest) > max_cells or total > MAX_TOTAL_BP:
            raise Capped('too_big', 'shortest %d bp x longest %d bp > %d cells, or %d bp in all' % (
                shortest, longest, max_cells, total))
        ci = medoid([s for _, s in alleles], weights)
        cid, cseq = alleles[ci]
        if max(1, len(cseq)) * max(1, longest) > max_cells:
            raise Capped('too_big', 'centre %d bp x longest %d bp > %d cells' % (len(cseq), longest, max_cells))
        pairs = collections.OrderedDict()
        for aid, s in alleles:
            if s == cseq:
                pairs[aid] = (cseq, s)
                continue
            if not cseq or not s:
                # an empty allele (a deletion of the whole interior) or an empty centre: all gaps
                pairs[aid] = (cseq + '-' * len(s), '-' * len(cseq) + s)
                continue
            got = abpoa_rows([('c', cseq), ('q', s)], wd, ['-b', '-1'], deadline, abpoa, scores)
            ca, sa = got['c'], got['q']
            keep = [i for i in range(len(ca)) if ca[i] != '-' or sa[i] != '-']
            ca, sa = ''.join(ca[i] for i in keep), ''.join(sa[i] for i in keep)
            if ca.replace('-', '') != cseq or sa.replace('-', '') != s:
                raise ValueError('abPOA rows do not spell their inputs (%s)' % aid)
            pairs[aid] = left_normalise_pair(ca, sa)
        exact, aligned = [], 0
        for f in pkg['fragments']:
            fs = f['seq'].upper()
            if place_fragment(fs, f['side'], alleles, {aid: s for aid, s in alleles}) is not None:
                exact.append(f)
                continue
            if not cseq or not fs:
                pairs[f['id']] = (cseq + '-' * len(fs), '-' * len(cseq) + fs)
                continue
            got = abpoa_rows([('c', cseq), ('q', fs)], wd, ['-b', '-1'], deadline, abpoa, scores)
            ca, sa = got['c'], got['q']
            keep = [i for i in range(len(ca)) if ca[i] != '-' or sa[i] != '-']
            ca, sa = ''.join(ca[i] for i in keep), ''.join(sa[i] for i in keep)
            if ca.replace('-', '') != cseq or sa.replace('-', '') != fs:
                raise ValueError('abPOA rows do not spell their inputs (%s)' % f['id'])
            pairs[f['id']] = left_normalise_pair(ca, sa)
            aligned += 1
        rows = star_merge(cseq, pairs, wd, deadline, stats, abpoa, scores)
        for f in exact:
            rows[f['id']] = place_fragment(f['seq'].upper(), f['side'], alleles, rows)
        stats['fragments_exact'], stats['fragments_aligned'] = len(exact), aligned
        rows = drop_gap_columns(rows)
        for mid, s in alleles + [(f['id'], f['seq'].upper()) for f in pkg['fragments']]:
            if rows[mid].replace('-', '') != s:
                raise ValueError('row %s does not spell its sequence' % mid)
        res.update(status='ok', centre=cid, rows=rows, columns=len(next(iter(rows.values()), '')),
                   slot_alignments=stats['slot_alignments'], fragments_exact=stats['fragments_exact'],
                   fragments_aligned=stats['fragments_aligned'])
    except Capped as e:
        res.update(status=e.status, message=str(e))
    finally:
        shutil.rmtree(wd, ignore_errors=True)
    res['seconds'] = round(time.time() - t0, 3)
    return res


def write_jsonl(path, objs):
    with gzip.open(path, 'wt') as f:
        for o in objs:
            f.write(json.dumps(o, separators=(',', ':')) + '\n')


def selected(pkg, select):
    """Whether a package's alleles differ in the sense `select` asks: 'any' (any two differ, or there are
    fragments) or 'length' (some two differ in length, or there are fragments)."""
    if pkg['fragments']:
        return True
    seqs = [a['seq'] for a in pkg['alleles']]
    if select == 'length':
        return len({len(x) for x in seqs}) > 1
    return len(set(seqs)) > 1


def _one(args):
    pkg, kw = args
    return realign_package(pkg, **kw)


def run(packages_path, out_path, select='any', shard=(0, 1), jobs=1, **kw):
    """Realign the selected `ok` packages of one shard (index i of n, by package order), in parallel."""
    import concurrent.futures
    from .extract import read_packages
    todo, skipped = [], collections.Counter()
    for k, p in enumerate(read_packages(packages_path)):
        if k % shard[1] != shard[0]:
            continue
        if p['status'] != 'ok':
            skipped[p['status']] += 1
        elif not selected(p, select):
            skipped['not_selected'] += 1
        else:
            todo.append(p)
    status, secs = collections.Counter(), 0.0
    with gzip.open(out_path + '.tmp', 'wt') as f, concurrent.futures.ProcessPoolExecutor(max(1, jobs)) as ex:
        for r in ex.map(_one, [(p, kw) for p in todo], chunksize=8):
            status[r['status']] += 1
            secs += r['seconds']
            f.write(json.dumps(r, separators=(',', ':')) + '\n')
    os.replace(out_path + '.tmp', out_path)
    return {'realigned': dict(status), 'skipped': dict(skipped), 'aligner_seconds': round(secs, 1)}
