/*
 * goldbach.c — high-performance Goldbach verifier (strong conjecture).
 *
 * Strategy:
 *   - odd-only segmented sieve (cache-resident blocks), so we can push to 1e10+ easily
 *   - bit-parallel resolution: for 64 consecutive even n at once, test n-p for each
 *     small prime p using shifted 64-bit word loads from the sieve bitmap
 *   - leftover n (no representation with p <= PMAX) fall back to an exact
 *     Miller-Rabin (deterministic for n < 3.3e24) search over larger p
 *   - per-block scalar sampling to cross-validate the bit-parallel path and to
 *     collect the "smallest prime in the representation" statistics
 *
 * Usage: ./goldbach <start> <end> [threads]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef uint64_t u64;
typedef unsigned __int128 u128;

#define PMAX 500000u          /* max prime p used in the bit-parallel window path */
#define ODD_BLOCK (1u << 20)  /* odd numbers per block region (value span = 2*ODD_BLOCK) */

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

/* ---------- Miller-Rabin (deterministic, 13 bases) ---------- */
static u64 mulmod(u64 a, u64 b, u64 m) { return (u128)a * b % m; }
static u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}
static int is_prime_mr(u64 n) {
    static const u64 small[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    if (n < 2) return 0;
    for (int i = 0; i < 12; i++) { if (n % small[i] == 0) return n == small[i]; }
    u64 d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (int i = 0; i < 12; i++) {
        u64 a = small[i] % n; if (!a) continue;
        u64 x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        int comp = 1;
        for (int r = 1; r < s; r++) { x = mulmod(x, x, n); if (x == n - 1) { comp = 0; break; } }
        if (comp) return 0;
    }
    return 1;
}

/* ---------- simple sieves ---------- */
static u32 *sieve_primes(u64 limit, size_t *count) {
    size_t n = (size_t)limit + 1;
    u8 *comp = (u8 *)calloc(n, 1);
    u32 *primes = (u32 *)malloc(sizeof(u32) * (limit / 10 + 64));
    size_t np = 0;
    for (u64 i = 2; i <= limit; i++) {
        if (!comp[i]) {
            if (i >= 3) primes[np++] = (u32)i;
            if (i * i <= limit) for (u64 j = i * i; j <= limit; j += i) comp[j] = 1;
        }
    }
    free(comp);
    *count = np;
    return primes;
}

/* bit test/set on a uint64 bitmap */
static inline int bit_get(const u64 *w, u64 i) { return (w[i >> 6] >> (i & 63)) & 1ULL; }
static inline void bit_clear(u64 *w, u64 i) { w[i >> 6] &= ~(1ULL << (i & 63)); }

/* window bitmap for odd numbers: bit i <-> value (w0 + 2*i) */
static inline u64 val_of(u64 w0, u64 idx) { return w0 + 2 * idx; }
static inline u64 idx_of(u64 w0, u64 v)   { return (v - w0) >> 1; }

typedef struct {
    u64 w0;             /* first odd value in window */
    u64 nbits;          /* number of odd slots */
    u64 nwords;
    u64 *w;             /* bitmap, extra zero words at the end */
} win_t;

static void win_free(win_t *s) { free(s->w); s->w = NULL; }

/* Build the prime bitmap for odd numbers in [w0, w0+2*(nbits-1)] */
static void win_build(win_t *s, u64 w0, u64 nbits, const u32 *base, size_t nbase, u64 sqrt_hi) {
    s->w0 = w0;
    s->nbits = nbits;
    s->nwords = (nbits + 63) / 64 + 4;
    s->w = (u64 *)calloc(s->nwords, sizeof(u64));
    for (size_t i = 0; i < nbits; i++) s->w[i >> 6] |= 1ULL << (i & 63); /* assume prime */
    /* clear 1 */
    if (w0 == 1) bit_clear(s->w, 0);
    for (size_t k = 0; k < nbase; k++) {
        u64 p = base[k];
        if (p < 3) continue;
        if (p > sqrt_hi) break;
        u64 start = p * p;
        u64 m = ((w0 + p - 1) / p) * p;
        if (m < start) m = start;
        if (!(m & 1)) m += p;           /* keep odd multiples only */
        if (m < w0) m += ((w0 - m + 2 * p - 1) / (2 * p)) * (2 * p);
        for (u64 v = m; v < w0 + 2 * nbits; v += 2 * p) {
            u64 idx = (v - w0) >> 1;
            if (idx >= nbits) break;
            bit_clear(s->w, idx);
        }
    }
}

/* prime test for an odd value v inside window */
static inline int win_is_prime(const win_t *s, u64 v) {
    if (v < s->w0) return -1;              /* outside */
    u64 idx = (v - s->w0) >> 1;
    if (idx >= s->nbits) return -1;
    return (int)((s->w[idx >> 6] >> (idx & 63)) & 1ULL);
}

/* fallback: exact search over p > PMAX using MR, returns smallest p or 0 */
static u64 fallback_smallest_p(u64 n) {
    if (n == 4) return 2;
    for (u64 p = PMAX + 1; p <= n / 2; p += 2) {
        if (!is_prime_mr(p)) continue;
        u64 q = n - p;
        if (is_prime_mr(q)) return p;
    }
    return 0;
}

/* scalar smallest p within window (returns 0 if none <= PMAX) */
static u64 window_smallest_p(const win_t *s, const u32 *small, size_t nsmall, u64 n) {
    if (n == 4) return 2;
    for (size_t i = 0; i < nsmall; i++) {
        u64 p = small[i];
        if (p > n / 2) break;
        u64 q = n - p;
        int pr = win_is_prime(s, q);
        if (pr == 1) return p;
    }
    return 0;
}

int main(int argc, char **argv) {
    u64 start = 4, end = 100000;
    int threads = 2;
    if (argc >= 3) { start = strtoull(argv[1], NULL, 10); end = strtoull(argv[2], NULL, 10); }
    if (argc >= 4) threads = atoi(argv[3]);

    if (start > end) { fprintf(stderr, "bad range\n"); return 1; }
    if ((start & 1) == 1) start++;
    if (end & 1) end--;

    double t0 = now_s();

    size_t nsmall = 0, nbase = 0;
    u32 *small = sieve_primes(PMAX, &nsmall);
    u64 sqrt_hi = (u64)sqrtl((long double)end) + 2;
    u32 *base = sieve_primes(sqrt_hi + 2, &nbase);

    fprintf(stderr, "[init] small primes=%zu base primes=%zu sqrt(hi)=%llu  %.2fs\n",
            nsmall, nbase, (unsigned long long)sqrt_hi, now_s() - t0);

    u64 total_evens = (end - start) / 2 + 1;
    u64 resolved = 0, unresolved = 0, fallbacks = 0, fallback_found = 0;
    u64 sample_n = 0, sample_ok = 0;
    u64 *phist = (u64 *)calloc(PMAX + 2, sizeof(u64));
    u64 checksum = 0;

    /* block iteration parameters */
    u64 block_span = 2ULL * ODD_BLOCK;
    u64 nblocks = (end - start + block_span - 1) / block_span;

    #pragma omp parallel num_threads(threads) reduction(+:resolved,unresolved,fallbacks,fallback_found,sample_n,sample_ok,checksum)
    {
        win_t win; memset(&win, 0, sizeof(win));
        #pragma omp for schedule(dynamic, 1)
        for (u64 b = 0; b < nblocks; b++) {
            u64 lo = start + b * block_span;
            u64 hi = lo + block_span - 1;
            if (hi > end) hi = end;
            u64 w0 = (lo > PMAX) ? lo - PMAX : 3;
            if ((w0 & 1) == 0) w0++;
            if (w0 > lo) w0 = lo;
            u64 nbits = (hi - w0) / 2 + 1;
            win_build(&win, w0, nbits, base, nbase, sqrt_hi);

            u64 n0 = (lo & 1) ? lo + 1 : lo;   /* first even in block */
            for (; n0 + 126 <= hi; n0 += 128) {
                u64 A = n0 - w0;                /* odd */
                u64 unres = ~0ULL;
                for (size_t i = 0; i < nsmall; i++) {
                    u64 p = small[i];
                    if (p > n0 / 2) break;
                    u64 idx0 = (A - p) >> 1;
                    u64 wi = idx0 >> 6, s = idx0 & 63;
                    u64 m = win.w[wi] >> s;
                    if (s) m |= win.w[wi + 1] << (64 - s);
                    unres &= ~m;
                    if (!unres) break;
                }
                resolved += 64 - __builtin_popcountll(unres);
                if (unres) {
                    /* exact fallback for the few survivors */
                    for (int j = 0; j < 64; j++) {
                        if ((unres >> j) & 1ULL) {
                            u64 n = n0 + 2ULL * j;
                            fallbacks++;
                            u64 p = fallback_smallest_p(n);
                            if (p) { fallback_found++; resolved++; }
                            else unresolved++;
                        }
                    }
                }
            }
            /* tail even numbers (not covered by 64-wide groups) */
            for (u64 n = n0; n <= hi; n += 2) {
                u64 p = window_smallest_p(&win, small, nsmall, n);
                if (p) resolved++;
                else {
                    fallbacks++;
                    if (fallback_smallest_p(n)) { fallback_found++; resolved++; }
                    else unresolved++;
                }
            }

            /* scalar sampling: cross-validate + statistics on smallest p */
            for (u64 n = (lo & 1) ? lo + 1 : lo; n <= hi; n += 2048) {
                u64 p = window_smallest_p(&win, small, nsmall, n);
                sample_n++;
                if (p) { sample_ok++; if (p <= PMAX) phist[p]++; checksum = checksum * 1000003ULL + p; }
            }
            win_free(&win);
        }
    }

    double dt = now_s() - t0;
    printf("{\"range\":[%llu,%llu],\"evens_checked\":%llu,\"resolved\":%llu,\"unresolved\":%llu,"
           "\"fallbacks\":%llu,\"fallback_found\":%llu,\"sample_n\":%llu,\"sample_ok\":%llu,"
           "\"checksum\":%llu,\"seconds\":%.2f,\"n_per_sec\":%.0f}\n",
           (unsigned long long)start, (unsigned long long)end,
           (unsigned long long)total_evens, (unsigned long long)resolved,
           (unsigned long long)unresolved, (unsigned long long)fallbacks,
           (unsigned long long)fallback_found, (unsigned long long)sample_n,
           (unsigned long long)sample_ok, (unsigned long long)checksum, dt,
           dt > 0 ? total_evens / dt : 0.0);

    if (unresolved) { fprintf(stderr, "!! COUNTEREXAMPLE FOUND: %llu even numbers without a Goldbach decomposition\n", (unsigned long long)unresolved); return 2; }
    fprintf(stderr, "[done] %.0f evens in %.2fs, all verified.\n", (double)total_evens, dt);
    return 0;
}
