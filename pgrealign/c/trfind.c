// trfind: seed-based tandem-repeat finder (periods MINP..MAXP), truth-independent.
// For every exact K-mer, the distance d to its previous occurrence is a candidate period.
// Hits with the same (drifting, ~3%-binned) period are chained while they lie within a gap
// of max(40, min(d/2, 300)) bp. A chain qualifies as a tandem repeat when its seed span
// covers >= 0.75 of a period (so the region holds >= ~2 copies), seed density >= 0.15, and the
// region [first_seed - period, last_seed + K) is >= 40 bp.
// Output (0-based, half-open): chrom start end period refined_period nseeds seedspan id_p gc
// refined_period = smallest divisor q of period whose Hamming self-identity at offset q is
// within 0.03 of the identity at offset period.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define K 10
#define NB 400

static char *seq; static long slen;
static int minp, maxp;
static const char *chrom;

typedef struct { long start, last; long n; double dsum; int active; } Run;
static Run runs[NB];

static int binof(long d) {
    if (d <= 40) return (int)d;
    int b = 40 + (int)floor(log((double)d / 40.0) / log(1.03));
    return b < NB ? b : NB - 1;
}
static long gapmax(long d) { long g = d / 2; if (g > 300) g = 300; if (g < 40) g = 40; return g; }

static double ident(long a, long b, long p) {
    long m = 0, t = 0;
    for (long i = a; i + p < b; i++) { t++; if (seq[i] == seq[i + p] && seq[i] != 'N') m++; }
    return t ? (double)m / t : 0.0;
}

static void finalize(Run *r) {
    if (!r->active) return;
    r->active = 0;
    if (r->n < 3) return;
    long p = (long)floor(r->dsum / r->n + 0.5);
    long span = r->last - r->start + K;
    if (span < 0.75 * p) return;
    if ((double)r->n / span < 0.15) return;
    long a = r->start - p; if (a < 0) a = 0;
    long b = r->last + K; if (b > slen) b = slen;
    if (b - a < 40) return;
    double idp = ident(a, b, p);
    long q = p;
    for (long dv = 2; dv <= p; dv++) {
        if (p % dv) continue;
        long qq = p / dv;
        if (qq < 1) break;
        double idq = ident(a, b, qq);
        if (idq >= idp - 0.03) q = qq;  // keep the smallest qualifying divisor
    }
    long gc = 0, acgt = 0;
    for (long i = a; i < b; i++) { char c = seq[i]; if (c == 'G' || c == 'C') gc++; if (c != 'N') acgt++; }
    printf("%s\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%.3f\t%.3f\n", chrom, a, b, p, q, r->n, span, idp,
           acgt ? (double)gc / acgt : 0.0);
}

static void scan(void) {
    static int *last = NULL;
    if (!last) last = malloc(sizeof(int) * (1 << (2 * K)));
    for (long i = 0; i < (1 << (2 * K)); i++) last[i] = -1;
    memset(runs, 0, sizeof(runs));
    unsigned code = 0; int valid = 0; unsigned mask = (1u << (2 * K)) - 1;
    for (long i = 0; i < slen; i++) {
        int v;
        switch (seq[i]) { case 'A': v = 0; break; case 'C': v = 1; break; case 'G': v = 2; break; case 'T': v = 3; break; default: v = -1; }
        if (v < 0) { valid = 0; code = 0; continue; }
        code = ((code << 2) | v) & mask; valid++;
        if (valid < K) continue;
        long pos = i - K + 1;
        int j = last[code]; last[code] = (int)pos;
        if (j < 0) continue;
        long d = pos - j;
        if (d < minp || d > maxp) continue;
        int b = binof(d);
        int chosen = -1;
        int cand[3] = { b, b - 1, b + 1 };
        for (int c = 0; c < 3; c++) {
            int bb = cand[c];
            if (bb < 0 || bb >= NB) continue;
            if (runs[bb].active && pos - runs[bb].last <= gapmax(d)) { chosen = bb; break; }
        }
        if (chosen < 0) {
            finalize(&runs[b]);
            runs[b].active = 1; runs[b].start = pos; runs[b].last = pos; runs[b].n = 1; runs[b].dsum = d;
        } else {
            Run *r = &runs[chosen];
            r->last = pos; r->n++; r->dsum += d;
        }
    }
    for (int b = 0; b < NB; b++) finalize(&runs[b]);
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: trfind in.fa minp maxp\n"); return 1; }
    minp = atoi(argv[2]); maxp = atoi(argv[3]);
    FILE *f = fopen(argv[1], "r");
    if (!f) { perror("open"); return 1; }
    fseek(f, 0, SEEK_END); long fsz = ftell(f); fseek(f, 0, SEEK_SET);
    char *buf = malloc(fsz + 1); fread(buf, 1, fsz, f); buf[fsz] = 0; fclose(f);
    seq = malloc(fsz + 1);
    char *p = buf; char namebuf[256];
    while (*p) {
        if (*p != '>') { p++; continue; }
        char *nl = strchr(p, '\n'); if (!nl) break;
        int nlen = 0; char *q = p + 1;
        while (q < nl && *q != ' ' && *q != '\t' && nlen < 255) namebuf[nlen++] = *q++;
        namebuf[nlen] = 0; chrom = namebuf;
        p = nl + 1; slen = 0;
        while (*p && *p != '>') {
            char c = *p++;
            if (c == '\n' || c == '\r') continue;
            if (c >= 'a' && c <= 'z') c -= 32;
            seq[slen++] = c;
        }
        scan();
        fflush(stdout);
    }
    return 0;
}
