// strf: a short-period tandem-repeat finder that reproduces the output of Tandem Repeats Finder (TRF)
// for the part of TRF that pgrealign's catalog uses.
//
// Clean-room provenance. Written from (1) the published algorithm, Benson G. (1999) "Tandem repeats
// finder: a program to analyze DNA sequences", Nucleic Acids Research 27(2):573-580; (2) TRF's public
// documentation; (3) black-box runs of the TRF binary, comparing outputs. A behavioural description of
// TRF written by someone who had read its source was also used; every point of it was checked against
// the binary before use. No TRF source code was read in writing this file.
//
// Usage: strf File Match Mismatch Delta PM PI Minscore MaxPeriod [options]
//   File is a FASTA file, or '-' for stdin. Options: -r (no redundancy elimination); -l N or -l=N (the
//   longest repeat expected, in millions of bp; default 2); -h, -ngs, -m, -f, -d and -noprune are
//   accepted and ignored. Any other option is an error, as in TRF.
//
// Limits.
// - PM must be 80 or 75 (TRF also accepts only these), and MaxPeriod 1..20 (TRF allows up to 2000).
// - Output is always TRF's -ngs format, on stdout, with or without -ngs: '@' and the header line (up to
//   199 characters) for each record with repeats, then a line per repeat: start end period copies
//   consensus_size %match %indel score %A %C %G %T entropy consensus sequence left_flank right_flank
//   (1-based, inclusive; flanks of 50 bp, '.' when empty). The order is that in which repeats are found.
// - PI is accepted and ignored: TRF uses it only for periods above 20.
//
// Input, as TRF reads it. The file must start with '>' (an empty file, or one starting with anything
// else, is an error, exit status 1). A header runs from '>' to the end of the line ('\n' or '\r'). In
// sequence, letters are kept, upper-cased; every other character (digits, spaces, '-', '*', ...) is
// skipped; a '>' anywhere starts a new header. A record with no sequence gives no output.
// Non-ACGT letters. N, IUPAC codes (R, Y, K, M, S, W, B, D, H, V) and any other letter (U, X, ...) are
// treated alike. They never take part in a k-tuple match, so a run of them is never searched; they get
// no vote in a consensus, which is therefore made of A, C, G and T; against a consensus they score as
// mismatches. A starting pattern may hold one, and there letters compare as letters (N matches N), as
// they do in %match and in the period. They are printed as read in the sequence and flanks, and count in
// a repeat's length but not in its %A %C %G %T.
//
// Method.
// 1. Detection (Benson 1999: k-tuple matches, distance lists, statistical criteria). K = 4 for PM 80
//    and 3 for PM 75 (the paper's Table 1, periods up to 29). For each distance d = 1..min(MaxPeriod+3,
//    20), a position i that ends a run of >= K bases with s[x] == s[x-d] (A, C, G or T) is a hit when,
//    in the window of the last W = max(d, 20) bases (cut at the record's start),
//    (a) at least 5 bases lie in such runs ("sum of heads"), and
//    (b) the first K-tuple match in the window ends fewer than WT bases after the window start
//        (a waiting-time form of the paper's "apparent size" test). WT is the 99th percentile of where
//        the first run of K heads ends in a Bernoulli(PM) sequence of length W that meets (a), plus 1,
//        computed exactly in waiting_time(): 18 for PM 80, 15 for PM 75.
//    Hits at d inside the last alignment made at d (the realignment of 4, if it was made, else the
//    first alignment) are not analysed again.
// 2. Alignment (wraparound dynamic programming; match +Match, mismatch -Mismatch, indel -Delta). The
//    pattern is the d bases ending at the hit. A backward pass from the hit, local within the detection
//    window and thereafter alive only while a cell scores > 0, gives the left end: the leftmost row with
//    a live cell. A forward pass from there is local up to the hit; past it a cell that falls to <= 0
//    dies and the pass stops at the first all-dead row. The best cell (earliest of equals) and its
//    traceback are the alignment. Each pass reaches at most the -l length from the hit, so a repeat
//    longer than that is found in pieces, in time linear in its length.
// 3. Tests: the alignment holds >= 1.9 copies; and d is among the 3 largest residuals of the histogram
//    of distances (up to 500) between identical 2-mers of A, C, G and T in the aligned region, after a
//    least-squares linear trend is removed (so a repeat of period p is analysed at p, 2p and 3p only),
//    or, for d = 1, one base makes up >= 80% of the region. A region longer than 65,536 bp (HREGION) is
//    histogrammed over its first 65,536 bp, which bounds the cost on whole chromosomes; pgrealign's 60 kb
//    batches never reach that.
// 4. Consensus by majority per column: ties between bases go to the alphabetically first (A, C, G, T),
//    and a gap wins a tie with a base (deleting the column); a base is inserted after a column when
//    >= 50% of the copies that span it have an insertion there. The sequence is realigned to the
//    consensus as in 2 (from the same hit); again >= 1.9 copies; kept when score >= Minscore.
// 5. Reported statistics: copies = consensus columns consumed / consensus size; period = the most
//    common distance between matching bases at the same column of adjacent copies; %match and %indel
//    over the alignments between adjacent copies that the consensus alignment induces; composition
//    (truncated percentages) and entropy over the repeat. Repeats with period > MaxPeriod are dropped.
// 6. Redundancy: a repeat is dropped when >= 90% of it lies in another kept repeat whose period equals
//    its own and whose score is at least as high, or whose period divides its own and whose score is at
//    least 1/1.1 of its own. Of identical repeats (found from different hits), the last found is kept,
//    which places it in the output.
//
// Differences from TRF that remain. On pgrealign's batches of CHM13 chr20, GRCh38 chr21 and S288C, 99.9%
// to 100% of the 'perfect' stage's repeats (2 7 7 80 10 50 10) and 99.1% to 99.3% of the 'old' stage's
// (2 3 5 75 20 33 7) are identical to TRF's in position and period, and almost all of those in every
// field. The rest are marginal repeats, mostly scoring near Minscore, whose detection, consensus or
// multiples test comes out differently; ties in the alignment path that move %match or %indel by a
// point; and, in a few batches, the order of two repeats found at nearly the same place. Malformed input
// that crashes TRF (a header at the end of a file with no newline after it) is read as an empty record.
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXC 64            // longest consensus
#define FLANK 50
#define WIN 20             // detection window, and the backward pass's local range, for d <= 20
#define HEADS 5            // sum-of-heads minimum
#define NPEAKS 3           // pattern sizes analysed per repeat
#define HMAX 500           // longest distance in the 2-mer histogram
#define HREGION 65536      // longest region histogrammed (> pgrealign's 60 kb batches)
#define MINCOPIES 1.9
#define NAMEMAX 199        // TRF keeps this much of a header line

static int MA, MI, DE, PM, MINSCORE, MAXP, K, DMAX, WT, NOREDUND;
static long MAXLEN = 2000000;               // -l: how far an alignment pass may reach from its hit

static unsigned char *seq; static long slen, scap;   // the record, upper-cased letters
static char name[NAMEMAX + 1]; static int namelen;
static unsigned char CODE[256];                       // A C G T -> 0..3, anything else 4
static const char *ACGT = "ACGT";

static void *xrealloc(void *p, size_t n) {
    void *q = realloc(p, n ? n : 1);
    if (!q) { fprintf(stderr, "strf: out of memory\n"); exit(2); }
    return q;
}

typedef struct {
    long s, e;              // 0-based, inclusive
    long order;             // when it was found
    int per, csz, score, pm, pi, acgt[4], dead;
    double copies, ent;
    char cons[MAXC + 1];
} Rec;
static Rec *recs; static long nrec, caprec;

// ---------------------------------------------------------------- wraparound alignment from a hit
typedef struct { char op; long x; int c; } Op;   // 'M' x aligned to column c, 'I' x inserted, 'D' c deleted
static Op *path; static long npath, cappath;
static unsigned char *tb; static size_t tbcap;
static const int NEG = -1000000000;

// Deletions within a row: left to right, then around the wrap for as long as it improves a cell.
static inline void deletions(int *N, unsigned char *t, int p) {
    for (int c = 1; c < p; c++) if (N[c - 1] - DE > N[c]) { N[c] = N[c - 1] - DE; if (t) t[c] = 3; }
    if (N[p - 1] - DE > N[0]) {
        N[0] = N[p - 1] - DE; if (t) t[0] = 3;
        for (int c = 1; c < p && N[c - 1] - DE > N[c]; c++) { N[c] = N[c - 1] - DE; if (t) t[c] = 3; }
    }
}

// Letters compare as letters: a consensus holds only A, C, G and T, so against one any other letter is a
// mismatch; a starting pattern may hold other letters, and there N matches N.
static long back_pass(long i, const unsigned char *pat, int p) {
    unsigned char rp[MAXC];
    for (int c = 0; c < p; c++) rp[c] = pat[p - 1 - c];
    int H[MAXC], N[MAXC];
    long wlo = i - ((WIN > p ? WIN : p) - 1), left = i, stop = i - MAXLEN;
    for (long x = i; x >= 0 && x > stop; x--) {
        unsigned char ch = seq[x];
        int fresh = x >= wlo, rowbest = 0;
        for (int c = 0; c < p; c++) {
            int sc = ch == rp[c] ? MA : -MI, cp = c == 0 ? p - 1 : c - 1;
            int v = (x < i && H[cp] > 0) ? H[cp] + sc : (fresh ? sc : NEG);
            int w = (x < i && H[c] > 0) ? H[c] - DE : NEG;
            N[c] = v > w ? v : w;
        }
        deletions(N, NULL, p);
        for (int c = 0; c < p; c++) { if (N[c] > rowbest) rowbest = N[c]; H[c] = N[c] > 0 ? N[c] : 0; }
        if (rowbest > 0) left = x;
        else if (!fresh) break;
    }
    return left;
}

static int fwd_pass(long lo, long i, const unsigned char *pat, int p, long *a, long *b) {
    int best = 0, bc = -1, H[MAXC], N[MAXC];
    long bx = -1, rows = 0, end = i + MAXLEN < slen ? i + MAXLEN : slen;   // up to -l past the hit
    for (int c = 0; c < p; c++) H[c] = 0;
    for (long x = lo; x < end; x++, rows++) {
        if ((size_t)(rows + 1) * p > tbcap) { tbcap = (size_t)(rows + 1) * p * 2 + 4096; tb = xrealloc(tb, tbcap); }
        unsigned char *t = tb + (size_t)rows * p, ch = seq[x];
        int local = x <= i, first = x == lo;
        for (int c = 0; c < p; c++) {
            int sc = ch == pat[c] ? MA : -MI, cp = c == 0 ? p - 1 : c - 1, v;
            unsigned char tt = 0;
            if (!first && H[cp] > 0) { v = H[cp] + sc; tt = 1; }       // tt 0: a fresh start
            else v = local ? sc : NEG;
            if (!first && H[c] > 0 && H[c] - DE > v) { v = H[c] - DE; tt = 2; }
            N[c] = v; t[c] = tt;
        }
        deletions(N, t, p);
        int rowbest = NEG, rc = 0;
        for (int c = 0; c < p; c++) { if (N[c] > rowbest) { rowbest = N[c]; rc = c; } H[c] = N[c]; }
        if (rowbest > best) { best = rowbest; bx = x; bc = rc; }
        if (rowbest <= 0 && !local) break;
    }
    npath = 0;
    if (bx < 0) return 0;
    for (long x = bx, c = bc;;) {
        if (npath + 1 > cappath) { cappath = cappath * 2 + 1024; path = xrealloc(path, cappath * sizeof(Op)); }
        unsigned char tt = tb[(size_t)(x - lo) * p + c];
        if (tt == 0) { path[npath++] = (Op){'M', x, (int)c}; break; }
        if (tt == 1) { path[npath++] = (Op){'M', x, (int)c}; x--; c = c == 0 ? p - 1 : c - 1; }
        else if (tt == 2) { path[npath++] = (Op){'I', x, (int)c}; x--; }
        else { path[npath++] = (Op){'D', x, (int)c}; c = c == 0 ? p - 1 : c - 1; }
    }
    for (long u = 0, v = npath - 1; u < v; u++, v--) { Op t = path[u]; path[u] = path[v]; path[v] = t; }
    *a = path[0].x; *b = bx;
    return best;
}

static int align_hit(long i, const unsigned char *pat, int p, long *a, long *b) {
    return fwd_pass(back_pass(i, pat, p), i, pat, p, a, b);
}

static double path_copies(int p) {
    long used = 0;
    for (long k = 0; k < npath; k++) used += path[k].op != 'I';
    return (double)used / p;
}

// ---------------------------------------------------------------- multiples test
static int multiples_ok(long a, long b, int d) {
    long L = b - a + 1;
    if (d == 1) {
        long comp[5] = {0, 0, 0, 0, 0};
        for (long x = a; x <= b; x++) comp[CODE[seq[x]]]++;
        for (int k = 0; k < 4; k++) if (comp[k] * 5 >= L * 4) return 1;
        return 0;
    }
    if (L > HREGION) L = HREGION;
    long D = L - 2 < HMAX ? L - 2 : HMAX;
    if (D < d) return 0;
    static double hist[HMAX + 1];
    static unsigned char c2[HREGION];
    long n2 = L - 1;
    for (long x = 0; x < n2; x++) {
        unsigned char u = CODE[seq[a + x]], v = CODE[seq[a + x + 1]];
        c2[x] = (u < 4 && v < 4) ? (unsigned char)(u * 4 + v) : 255;
    }
    for (long dd = 1; dd <= D; dd++) {
        long n = 0, m = n2 - dd;
        for (long x = 0; x < m; x++) n += (c2[x] == c2[x + dd]) & (c2[x] != 255);
        hist[dd] = n;
    }
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (long dd = 1; dd <= D; dd++) { sx += dd; sy += hist[dd]; sxx += (double)dd * dd; sxy += dd * hist[dd]; }
    double den = D * sxx - sx * sx, slope = den != 0 ? (D * sxy - sx * sy) / den : 0, icpt = (sy - slope * sx) / D;
    double rd = hist[d] - (icpt + slope * d);
    int above = 0;                                   // distances ranked above d (ties: the smaller first)
    for (long dd = 1; dd <= D; dd++) {
        double r = hist[dd] - (icpt + slope * dd);
        if ((r > rd + 1e-9 || (fabs(r - rd) <= 1e-9 && dd < d)) && ++above >= NPEAKS) return 0;
    }
    return 1;
}

// ---------------------------------------------------------------- consensus by majority rule
static int consensus(int p, unsigned char *out) {
    int votes[MAXC][5], insv[MAXC][4], span[MAXC], insyes[MAXC];
    memset(votes, 0, sizeof votes); memset(insv, 0, sizeof insv);
    memset(span, 0, sizeof span); memset(insyes, 0, sizeof insyes);
    int lastc = -1, ninsert = 0, firstins = 4;
    for (long i = 0; i < npath; i++) {
        Op *o = &path[i];
        if (o->op == 'I') { if (ninsert++ == 0) firstins = CODE[seq[o->x]]; continue; }
        if (lastc >= 0) {
            span[lastc]++;
            if (ninsert > 0) { insyes[lastc]++; if (firstins < 4) insv[lastc][firstins]++; }
        }
        ninsert = 0;
        if (o->op == 'M') { if (CODE[seq[o->x]] < 4) votes[o->c][CODE[seq[o->x]]]++; }
        else votes[o->c][4]++;
        lastc = o->c;
    }
    int n = 0;
    for (int c = 0; c < p; c++) {
        int bk = 0;
        for (int k = 1; k < 4; k++) if (votes[c][k] > votes[c][bk]) bk = k;
        if (votes[c][4] < votes[c][bk] && n < MAXC) out[n++] = ACGT[bk];
        if (span[c] > 0 && 2 * insyes[c] >= span[c] && n < MAXC) {
            int ik = 0;
            for (int k = 1; k < 4; k++) if (insv[c][k] > insv[c][ik]) ik = k;
            out[n++] = ACGT[ik];
        }
    }
    return n;
}

// ---------------------------------------------------------------- statistics of the final alignment
// Cells are numbered copy * p + column; pos[] holds a base's position, -1 for a deletion, -2 outside the
// alignment. Does the alignment go on past cell k (so that the slot after it, where insertions sit, lies
// inside the alignment)?
static long *pos;
static inline int spanned(long k, long n) { return k + 1 < n && pos[k + 1] != -2; }
static void stats(int p, Rec *r) {
    static int *copy; static long capcopy;
    static long *ins0; static int *insn; static long capcells;
    if (npath > capcopy) { capcopy = npath * 2; copy = xrealloc(copy, capcopy * sizeof(int)); }
    int t = 0, prev = -1;
    long used = 0;
    for (long i = 0; i < npath; i++) {
        if (path[i].op != 'I') { if (prev >= 0 && path[i].c <= prev) t++; prev = path[i].c; used++; }
        copy[i] = t;
    }
    int ncopies = t + 1;
    r->copies = (double)used / p;
    // per (copy, column): the base's position (M), -1 (D) or -2 (not in the alignment); insertions after it
    long n = (long)ncopies * p;
    if (n > capcells) {
        capcells = n * 2;
        pos = xrealloc(pos, sizeof(long) * capcells); ins0 = xrealloc(ins0, sizeof(long) * capcells);
        insn = xrealloc(insn, sizeof(int) * capcells);
    }
    for (long k = 0; k < n; k++) { pos[k] = -2; insn[k] = 0; ins0[k] = -1; }
    for (long i = 0; i < npath; i++) {
        long k = (long)copy[i] * p + path[i].c;
        if (path[i].op == 'M') pos[k] = path[i].x;
        else if (path[i].op == 'D') pos[k] = -1;
        else { if (insn[k]++ == 0) ins0[k] = path[i].x; }
    }
    int dcount[4 * MAXC + 2]; memset(dcount, 0, sizeof dcount);
    long mt = 0, mm = 0, ind = 0;
    for (int tt = 0; tt + 1 < ncopies; tt++)
        for (int c = 0; c < p; c++) {
            long k1 = (long)tt * p + c, k2 = k1 + p, u = pos[k1], v = pos[k2];
            if (u >= 0 && v >= 0 && seq[u] == seq[v] && v - u < 4 * MAXC + 2) dcount[v - u]++;
            if (u == -2 || v == -2) continue;
            if (u == -1 && v == -1) {}
            else if (u == -1 || v == -1) ind++;
            else if (seq[u] == seq[v]) mt++;
            else mm++;
            int i1 = insn[k1], i2 = insn[k2];
            if (!spanned(k1, n) || !spanned(k2, n)) continue;     // an insertion slot at an end
            for (int z = 0; z < (i1 < i2 ? i1 : i2); z++) {
                if (seq[ins0[k1] + z] == seq[ins0[k2] + z]) mt++; else mm++;
            }
            ind += i1 > i2 ? i1 - i2 : i2 - i1;
        }
    int bd = p, bv = 0;
    for (int dd = 1; dd < 4 * MAXC + 2; dd++) if (dcount[dd] > bv) { bv = dcount[dd]; bd = dd; }
    r->per = bd;
    long tot = mt + mm + ind;
    r->pm = tot ? (int)(100 * mt / tot) : 0;
    r->pi = tot ? (int)(100 * ind / tot) : 0;
    long L = r->e - r->s + 1, comp[5] = {0, 0, 0, 0, 0};
    for (long x = r->s; x <= r->e; x++) comp[CODE[seq[x]]]++;
    double ent = 0;
    for (int k = 0; k < 4; k++) {
        r->acgt[k] = (int)(100 * comp[k] / L);
        if (comp[k]) { double f = (double)comp[k] / L; ent -= f * log2(f); }
    }
    r->ent = ent < 0.0005 ? 0.0 : ent;
}

// ---------------------------------------------------------------- analysis of one hit
// Returns the end of the alignment made at d, the realignment if there is one, else the first; later
// hits at d up to there are not analysed.
static long analyze(long i, int d) {
    unsigned char pat[MAXC], cons[MAXC];
    for (int c = 0; c < d; c++) pat[c] = seq[i - d + 1 + c];
    long a, b, a2, b2;
    align_hit(i, pat, d, &a, &b);
    long cover = b > i ? b : i;
    if (npath == 0 || path_copies(d) < MINCOPIES || !multiples_ok(a, b, d)) return cover;
    int p = consensus(d, cons);
    if (p == 0) return cover;
    int sc = align_hit(i, cons, p, &a2, &b2);
    if (npath) cover = b2 > i ? b2 : i;
    if (npath == 0 || path_copies(p) < MINCOPIES || sc < MINSCORE) return cover;
    Rec r; memset(&r, 0, sizeof r);
    r.s = a2; r.e = b2; r.score = sc; r.csz = p; r.order = nrec;
    for (int k = 0; k < p; k++) r.cons[k] = (char)cons[(path[0].c + k) % p];   // from the first column aligned
    stats(p, &r);
    if (r.per > MAXP) return cover;
    if (nrec == caprec) { caprec = caprec * 2 + 256; recs = xrealloc(recs, caprec * sizeof(Rec)); }
    recs[nrec++] = r;
    return cover;
}

// ---------------------------------------------------------------- detection
static void detect(void) {
    enum { Q = 64 };                                   // > any window's heads or tuple ends
    long runstart[64], nexttrig[64], hq[64][Q], nh[64], tq[64][Q], nt[64];
    memset(runstart, 0, sizeof runstart); memset(nexttrig, 0, sizeof nexttrig);
    memset(nh, 0, sizeof nh); memset(nt, 0, sizeof nt);
    unsigned long long prev = 0;
    for (long i = 0; i < slen; i++) {
        unsigned long long m = 0;
        unsigned char ch = seq[i];
        if (CODE[ch] < 4)
            for (int d = 1; d <= DMAX && d <= i; d++) m |= (unsigned long long)(seq[i - d] == ch) << d;
        for (unsigned long long st = m & ~prev; st; st &= st - 1) runstart[__builtin_ctzll(st)] = i;
        prev = m;
        for (unsigned long long mm = m; mm; mm &= mm - 1) {
            int d = __builtin_ctzll(mm);
            long run = i - runstart[d] + 1;
            if (run < K) continue;
            for (long x = run == K ? i - K + 1 : i; x <= i; x++) hq[d][nh[d]++ % Q] = x;   // new heads
            tq[d][nt[d]++ % Q] = i;                                                     // a tuple ends at i
            if (i < nexttrig[d]) continue;
            long w0 = i - (WIN > d ? WIN : d) + 1;                                      // window [w0, i]
            if (w0 < 0) w0 = 0;
            if (nh[d] < HEADS || hq[d][(nh[d] - HEADS) % Q] < w0) continue;            // (a)
            long first = i;
            for (long k = nt[d] - 1; k >= 0 && k >= nt[d] - Q && tq[d][k % Q] >= w0; k--) first = tq[d][k % Q];
            if (first - w0 >= WT) continue;                                             // (b)
            nexttrig[d] = analyze(i, d) + 1;
        }
    }
}

// ---------------------------------------------------------------- redundancy elimination
static int cmp_priority(const void *x, const void *y) {     // the order in which repeats are decided
    const Rec *a = x, *b = y;
    if (a->per != b->per) return a->per - b->per;
    if (a->score != b->score) return b->score - a->score;
    if (a->e - a->s != b->e - b->s) return (b->e - b->s) > (a->e - a->s) ? 1 : -1;
    if (a->s != b->s) return a->s < b->s ? -1 : 1;
    return a->order > b->order ? -1 : a->order < b->order;   // of identical repeats, the last found is kept
}
static int cmp_start(const void *x, const void *y) {
    const Rec *a = *(Rec *const *)x, *b = *(Rec *const *)y;
    return a->s < b->s ? -1 : a->s > b->s;
}
static int cmp_order(const void *x, const void *y) {
    const Rec *a = x, *b = y;
    return a->order < b->order ? -1 : a->order > b->order;
}
// A repeat B is dropped when a repeat A decided before it, and kept, has a period dividing B's, holds
// >= 90% of B and scores well enough. Only A that hold >= 90% of B can matter, so candidates are found by
// a sweep over the repeats sorted by start, cut off by the running maximum of their ends.
static void redundancy(void) {
    if (NOREDUND || nrec < 2) return;
    qsort(recs, nrec, sizeof(Rec), cmp_priority);
    static Rec **bys; static long *maxend; static long capr;
    if (nrec > capr) {
        capr = nrec * 2;
        bys = xrealloc(bys, capr * sizeof(Rec *)); maxend = xrealloc(maxend, capr * sizeof(long));
    }
    for (long j = 0; j < nrec; j++) { recs[j].dead = 0; bys[j] = &recs[j]; }
    qsort(bys, nrec, sizeof(Rec *), cmp_start);
    for (long j = 0; j < nrec; j++) maxend[j] = j && maxend[j - 1] > bys[j]->e ? maxend[j - 1] : bys[j]->e;
    for (long j = 0; j < nrec; j++) {
        Rec *B = &recs[j];
        long lb = B->e - B->s + 1, need = (9 * lb + 9) / 10;    // overlap >= 0.9 lb
        // A must start by B->e - need + 1 and end at or after B->s + need - 1
        long lo = 0, hi = nrec;                                  // first index with start > B->e - need + 1
        while (lo < hi) { long mid = (lo + hi) / 2; if (bys[mid]->s <= B->e - need + 1) lo = mid + 1; else hi = mid; }
        for (long k = lo - 1; k >= 0 && maxend[k] >= B->s + need - 1; k--) {
            Rec *A = bys[k];
            if (A >= B || A->dead || B->per % A->per || A->e < B->s) continue;   // A >= B: decided after B
            long o = (A->e < B->e ? A->e : B->e) - (A->s > B->s ? A->s : B->s) + 1;
            if (o * 10 >= lb * 9 && (A->per == B->per ? B->score <= A->score : B->score * 10 <= A->score * 11)) {
                B->dead = 1;
                break;
            }
        }
    }
    long k = 0;
    for (long j = 0; j < nrec; j++) if (!recs[j].dead) recs[k++] = recs[j];
    nrec = k;
}

// ---------------------------------------------------------------- output
static void put_range(long a, long b) {                     // seq[a, b), or '.' when empty
    if (a >= b) putchar('.');
    else fwrite(seq + a, 1, (size_t)(b - a), stdout);
}
static void report(void) {
    if (!nrec) return;
    qsort(recs, nrec, sizeof(Rec), cmp_order);
    printf("@%.*s\n", namelen, name);
    for (long j = 0; j < nrec; j++) {
        Rec *r = &recs[j];
        printf("%ld %ld %d %.1f %d %d %d %d %d %d %d %d %.2f %s ", r->s + 1, r->e + 1, r->per, r->copies,
               r->csz, r->pm, r->pi, r->score, r->acgt[0], r->acgt[1], r->acgt[2], r->acgt[3], r->ent, r->cons);
        put_range(r->s, r->e + 1);
        putchar(' ');
        put_range(r->s - FLANK < 0 ? 0 : r->s - FLANK, r->s);
        putchar(' ');
        put_range(r->e + 1, r->e + 1 + FLANK > slen ? slen : r->e + 1 + FLANK);
        putchar('\n');
    }
}

static void process(void) {
    nrec = 0;
    if (slen > 0) { detect(); redundancy(); report(); }
}

// ---------------------------------------------------------------- waiting-time threshold
// For an iid Bernoulli(pm) sequence of length n: the distribution of where the first run of k heads
// ends, given that heads in runs of >= k number at least t; returns 1 + its q-quantile. Exact DP.
static int waiting_time(int n, int k, double pm, int t, double q) {
    static double P[2][MAXC + 2][MAXC + 2];    // [run length, capped at k][heads, capped at t]
    double *pend = calloc(n + 1, sizeof(double)), total = 0;
    if (!pend) { fprintf(stderr, "strf: out of memory\n"); exit(2); }
    for (int e = k - 1; e < n; e++) {          // P(the first run of k ends at e, and >= t heads)
        memset(P, 0, sizeof P);
        P[0][0][0] = 1.0;
        int cur = 0;
        for (int x = 0; x < n; x++, cur ^= 1) {
            memset(P[cur ^ 1], 0, sizeof P[0]);
            for (int r = 0; r <= k; r++)
                for (int h = 0; h <= t; h++) {
                    double v = P[cur][r][h];
                    if (v == 0) continue;
                    int completes = r + 1 == k && h == 0;           // this head would end the first run
                    if (x != e) P[cur ^ 1][0][h] += v * (1 - pm);  // tail
                    if ((x < e && completes) || (x == e && !completes)) continue;
                    int h2 = r + 1 == k ? h + k : r == k ? h + 1 : h;
                    P[cur ^ 1][r < k ? r + 1 : k][h2 > t ? t : h2] += v * pm;
                }
        }
        for (int r = 0; r <= k; r++) pend[e] += P[cur][r][t];
        total += pend[e];
    }
    double acc = 0; int w = n;
    for (int e = 0; e < n; e++) if ((acc += pend[e]) >= q * total) { w = e + 1; break; }
    free(pend);
    return w;
}

// ---------------------------------------------------------------- input
// Reads FASTA as TRF does (see the header comment), processing each record as it ends. Returns 0, or
// 1 when the input does not start with '>'.
static int read_fasta(FILE *in) {
    static unsigned char buf[1 << 16];
    enum { START, HEADER, SEQ } state = START;
    size_t n;
    while ((n = fread(buf, 1, sizeof buf, in)) > 0) {
        for (size_t j = 0; j < n; j++) {
            unsigned char ch = buf[j];
            if (state == SEQ) {
                if (isalpha(ch)) {
                    if (slen + 1 >= scap) { scap = scap * 2 + (1 << 20); seq = xrealloc(seq, scap); }
                    seq[slen++] = (unsigned char)toupper(ch);
                } else if (ch == '>') {
                    process();
                    state = HEADER; namelen = 0; slen = 0;
                }
            } else if (state == HEADER) {
                if (ch == '\n' || ch == '\r') state = SEQ;
                else if (namelen < NAMEMAX) name[namelen++] = (char)ch;
            } else {
                if (ch != '>') return 1;
                state = HEADER; namelen = 0; slen = 0;
            }
        }
    }
    if (state == START) return 1;
    if (state == SEQ) process();       // a header with no line end after it is an empty record
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 9) {
        fprintf(stderr, "usage: strf File Match Mismatch Delta PM PI Minscore MaxPeriod [-r] [-l N] "
                        "[-h] [-ngs] [-m] [-f] [-d] [-noprune]\n");
        return 1;
    }
    MA = atoi(argv[2]); MI = atoi(argv[3]); DE = atoi(argv[4]); PM = atoi(argv[5]);
    MINSCORE = atoi(argv[7]); MAXP = atoi(argv[8]);
    if (PM != 80 && PM != 75) { fprintf(stderr, "strf: PM must be 80 or 75\n"); return 1; }
    if (MAXP < 1 || MAXP > 20) { fprintf(stderr, "strf: MaxPeriod must be 1..20\n"); return 1; }
    for (int j = 9; j < argc; j++) {
        const char *o = argv[j];
        if (!strcmp(o, "-r")) NOREDUND = 1;
        else if (!strcmp(o, "-l") || !strncmp(o, "-l=", 3)) {
            const char *v = o[2] == '=' ? o + 3 : (j + 1 < argc ? argv[++j] : "");
            long mb = atol(v);
            if (mb < 1) { fprintf(stderr, "strf: -l needs a positive number of millions\n"); return 1; }
            MAXLEN = mb * 1000000L;
        }
        else if (!strcmp(o, "-h") || !strcmp(o, "-ngs") || !strcmp(o, "-m") || !strcmp(o, "-f") ||
                 !strcmp(o, "-d") || !strcmp(o, "-noprune")) {}
        else { fprintf(stderr, "strf: unknown option %s\n", o); return 1; }
    }
    for (int c = 0; c < 256; c++) CODE[c] = 4;
    CODE['A'] = 0; CODE['C'] = 1; CODE['G'] = 2; CODE['T'] = 3;
    K = PM == 80 ? 4 : 3;
    DMAX = MAXP + 3 < 20 ? MAXP + 3 : 20;
    WT = waiting_time(WIN, K, PM / 100.0, HEADS, 0.99);
    FILE *in = strcmp(argv[1], "-") ? fopen(argv[1], "rb") : stdin;
    if (!in) { fprintf(stderr, "Error: Error opening %s\n", argv[1]); perror(argv[1]); return 1; }
    static char obuf[1 << 20];
    setvbuf(stdout, obuf, _IOFBF, sizeof obuf);
    int bad = read_fasta(in);
    if (in != stdin) fclose(in);
    if (bad) {
        fprintf(stderr, "Error: Could not load sequence. Empty file or bad format.\n");
        printf("Error processing input: Bad format.\n");
        fflush(stdout);
        return 1;
    }
    if (fflush(stdout) != 0 || ferror(stdout)) { perror("strf: writing output"); return 1; }
    return 0;
}
