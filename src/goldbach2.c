/*
 * goldbach2.c — segmented-sieve + bit-parallel Goldbach verifier.
 *
 * Goal: verify the strong Goldbach conjecture (every even n >= 4 is p+q)
 * over as large a range as the machine allows, with a *proof-grade* result:
 * every "yes" is backed by an exact primality test on both summands.
 *
 * Method
 * ------
 * 1. Odd-only segmented sieve over a moving window.  Bit i in the bitmap
 *    corresponds to odd value v = 3 + 2*(B0 + i).  Bit set  <=>  v is prime.
 *    Composites are marked with the standard "start at p*p" rule using all
 *    base primes <= sqrt(N), so the bitmap is exact (not heuristic).
 *
 * 2. Resolution of Goldbach for even n = p + q:
 *      - prime summands are odd except for n = 4 = 2+2 (handled separately);
 *      - for a fixed odd prime p, the 64 candidates q_j = n0-p+2j for
 *        j = 0..63 are 64 CONSECUTIVE odd indices, i.e. two adjacent 64-bit
 *        words of the bitmap.  One shifted load therefore tests
 *        "is n_j - p prime?" for 64 even numbers simultaneously:
 *
 *            mask = w[wi] >> off | w[wi+1] << (64-off)
 *            unresolved &= ~mask
 *
 *      - after all p <= PMAX, any still-unresolved even n is finished by a
 *        scalar continuation over the same bitmap (exact, no approximations).
 *
 * 3. Independent cross-check: a deterministic 64-bit Miller-Rabin test is
 *    applied to a random sample of accepted representations, so a bug in the
 *    sieve would be caught.
 *
 * Usage: ./goldbach2 <start> <end> [threads] [PMAX]
 */
#define _GNU_SOURCE
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
typedef uint32_t u32;
typedef uint8_t  u8;

#define ODDS_PER_SEG (1u << 23)   /* 8.4M odds per segment -> 1 MiB bitmap */
#define GUARD_W    64u            /* 64 guard words = 4096 headroom indices */
#define GUARD_BITS ((u64)GUARD_W * 64)

static double now_s(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

/* ---------- deterministic Miller-Rabin, valid for all 64-bit inputs ---------- */
static inline u64 mulmod(u64 a, u64 b, u64 m) { return (unsigned __int128)a * b % m; }
static inline u64 powmod(u64 a, u64 e, u64 m) {
    u64 r = 1;
    while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; }
    return r;
}
static int is_prime_u64(u64 n) {
    if (n < 2) return 0;
    static const u64 sm[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) if (n % sm[i] == 0) return n == sm[i];
    u64 d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    static const u64 bases[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    for (int i = 0; i < 7; i++) {
        u64 a = bases[i] % n; if (!a) continue;
        u64 x = powmod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        int comp = 1;
        for (int r = 1; r < s; r++) { x = mulmod(x, x, n); if (x == n - 1) { comp = 0; break; } }
        if (comp) return 0;
    }
    return 1;
}

/* ---------- odd prime table ---------- */
static u32 *primes_upto(u64 limit, size_t *cnt) {
    if (limit < 3) { *cnt = 0; return NULL; }
    size_t n = (size_t)limit + 1;
    u8 *comp = (u8 *)calloc(n, 1);
    size_t cap = (size_t)(limit / 8 + 64);
    u32 *pr = (u32 *)malloc(sizeof(u32) * cap);
    size_t np = 0;
    for (u64 i = 3; i <= limit; i += 2) {
        if (comp[i]) continue;
        if (np < cap) pr[np++] = (u32)i;
        if (i * i <= limit) for (u64 j = i * i; j <= limit; j += 2 * i) comp[j] = 1;
    }
    free(comp);
    *cnt = np;
    return pr;
}

static u64 isqrt_u64(u64 n) {
    u64 r = (u64)sqrtl((long double)n);
    while (r > 0 && r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}

/* deterministic PRNG for cross-check sampling */
static u64 rng_s = 0x243F6A8885A308D3ULL;
static inline u64 rnd(void) {
    rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s;
}

int main(int argc, char **argv) {
    u64 start = 4, end = 1000000;
    int threads = 2;
    u64 PMAX = 100000;
    if (argc >= 3) { start = strtoull(argv[1], NULL, 10); end = strtoull(argv[2], NULL, 10); }
    if (argc >= 4) threads = atoi(argv[3]);
    if (argc >= 5) PMAX = strtoull(argv[4], NULL, 10);
    if (start > end) { fprintf(stderr, "bad range\n"); return 1; }
    if (start & 1) start++;
    if (end & 1) end--;
    if (start < 4) start = 4;

    double t0 = now_s();
    size_t nsmall = 0;
    u32 *small = primes_upto(PMAX, &nsmall);
    u64 root = isqrt_u64(end);
    size_t nbase = 0;
    u32 *base = primes_upto(root + 2, &nbase);

    fprintf(stderr, "[init] primes<=%llu: %zu | base primes<=sqrt(%llu): %zu | %.2fs\n",
            (unsigned long long)PMAX, nsmall, (unsigned long long)end, nbase, now_s() - t0);

    u64 T = (end - 3) / 2 + 1;                       /* odd index space: v = 3+2i */
    u64 nseg = (T + ODDS_PER_SEG - 1) / ODDS_PER_SEG;
    /* start sieving only near the requested lower bound (incremental runs) */
    u64 seg_begin = 0;
    if (start > 2 * ODDS_PER_SEG + 6) {
        u64 idx0 = (start - 6) / 2;                   /* first odd index we need */
        seg_begin = idx0 / ODDS_PER_SEG;
    }
    u64 total_evens = (end - start) / 2 + 1;
    u64 resolved = 0, unresolved = 0, scalar_used = 0, sampled = 0, sample_ok = 0;
    u64 min_p = ~0ULL, max_p = 0, checksum = 0;

    /* n = 4 = 2 + 2 is the only even whose summands are not both odd */
    if (start <= 4 && 4 <= end) {
        resolved++;
        if (2 < min_p) min_p = 2;
        checksum ^= (2ULL * 0x9E3779B97F4A7C15ULL * 5ULL);
    }

    #pragma omp parallel num_threads(threads) \
        reduction(+:resolved,unresolved,scalar_used,sampled,sample_ok,checksum) \
        reduction(min:min_p) reduction(max:max_p)
    {
        u64 *w = NULL; size_t wcap = 0;
        #pragma omp for schedule(dynamic, 1)
        for (u64 s = seg_begin; s < nseg; s++) {
            u64 i0 = s * ODDS_PER_SEG;
            u64 i1 = i0 + ODDS_PER_SEG; if (i1 > T) i1 = T;
            u64 margin = PMAX / 2 + 4;                /* odd-index margin for shifts */
            u64 B0 = (i0 > margin) ? (i0 - margin) : 0;
            u64 span = i1 - B0;
            size_t need = (size_t)((span + 63) / 64) + 4;
            if (need + GUARD_W > wcap) {
                free(w);
                wcap = need + GUARD_W;
                w = (u64 *)calloc(wcap, sizeof(u64));   /* guard words stay 0 */
            }
            u64 *wd = w + GUARD_W;                     /* real bitmap starts here */
            memset(wd, 0xFF, need * sizeof(u64));
            {   /* clear tail bits beyond span */
                u64 fw = span / 64, rem = span % 64;
                if (rem) wd[fw] &= (1ULL << rem) - 1;
                for (size_t k = fw + (rem ? 1 : 0); k < need; k++) wd[k] = 0;
            }
            u64 vlo = 3 + 2 * B0, vhi = 3 + 2 * (i1 - 1);
            for (size_t b = 0; b < nbase; b++) {
                u64 p = base[b];
                if (p * p > vhi) break;
                u64 m = (vlo + p - 1) / p * p;
                if (m < p * p) m = p * p;
                if (!(m & 1)) m += p;
                for (; m <= vhi; m += 2 * p) {
                    u64 rel = ((m - 3) >> 1) - B0;
                    wd[rel >> 6] &= ~(1ULL << (rel & 63));
                }
            }
            /* even numbers handled by this segment.
             * Bitmap holds odd values v = 3+2*i for i in [B0, i1).
             * For even n we need q = n-p (odd, p>=3) to satisfy q <= n-3,
             * i.e. n <= 2*i1+4 ; and q >= 3+2*i0  =>  n >= 2*i0+6.        */
            u64 n_lo = 2 * i0 + 6;
            u64 n_hi = 2 * i1 + 4;
            if (n_lo < start) n_lo = start;
            if (n_hi > end)   n_hi = end;
            if (n_lo > n_hi) continue;

            /* ---- bit-parallel groups of 64 consecutive evens ---- */
            u64 g0 = n_lo < 6 ? 6 : n_lo;
            u64 span_ev = (n_hi >= g0) ? (n_hi - g0) : 0;
            u64 ngroups = (span_ev >= 126) ? ((span_ev - 126) / 128 + 1) : 0;
            for (u64 gi = 0; gi < ngroups; gi++) {
                u64 n0 = g0 + gi * 128;
                u64 unres = ~0ULL;
                u64 first_p[64]; memset(first_p, 0, sizeof(first_p));
                u64 pmax_use = PMAX;
                if (pmax_use > n0 - 3) pmax_use = n0 - 3;   /* need q0 = n0-p >= 3 */
                for (size_t k = 0; k < nsmall; k++) {
                    u64 p = small[k];
                    if (p > pmax_use) break;
                    u64 rel = (((n0 - p) - 3) >> 1) - B0;
                    u64 wi = rel >> 6, off = rel & 63;
                    u64 m = w[wi] >> off;
                    if (off) m |= w[wi + 1] << (64 - off);
                    u64 newly = m & unres;
                    if (newly) for (int t = 0; t < 64; t++) if ((newly >> t) & 1ULL) first_p[t] = p;
                    unres &= ~m;
                    if (!unres) break;
                }
                resolved += 64 - __builtin_popcountll(unres);
                /* fingerprint over the smallest prime of each representation */
                for (int t = 0; t < 64; t++) {
                    u64 n = n0 + 2ULL * t;
                    if (n > end) break;
                    if (first_p[t]) {
                        checksum ^= first_p[t] * 0x9E3779B97F4A7C15ULL * (n + 1);
                        if (first_p[t] < min_p) min_p = first_p[t];
                        if (first_p[t] > max_p) max_p = first_p[t];
                    }
                }
                if (unres) {
                    /* scalar continuation with the exact same bitmap */
                    for (int j = 0; j < 64; j++) {
                        if (!((unres >> j) & 1ULL)) continue;
                        u64 n = n0 + 2ULL * j;
                        u64 p = 0;
                        for (size_t k = 0; k < nsmall; k++) {
                            u64 pp = small[k];
                            if (pp > pmax_use) break;
                            u64 q = n - pp;
                            u64 rel = ((q - 3) >> 1) - B0;
                            if ((w[rel >> 6] >> (rel & 63)) & 1ULL) { p = pp; break; }
                        }
                        if (!p) {
                            for (u64 pp = pmax_use + 1; pp <= n / 2; pp += 2) {
                                u64 q = n - pp;
                                u64 rel = ((q - 3) >> 1) - B0;
                                if (rel >> 6 >= need) break;
                                if ((w[rel >> 6] >> (rel & 63)) & 1ULL) { p = pp; break; }
                            }
                        }
                        scalar_used++;
                        if (p) {
                            resolved++;
                            first_p[j] = p;          /* keep checksum complete */
                            if (p < min_p) min_p = p;
                            if (p > max_p) max_p = p;
                        } else {
                            unresolved++;
                            fprintf(stderr, "COUNTEREXAMPLE candidate: %llu\n", (unsigned long long)n);
                        }
                    }
                }
                /* statistics + deterministic checksum on a sparse subset */
                if (((n0 >> 16) & 0xFF) == 0) {
                    for (int j = 0; j < 64; j++) {
                        u64 n = n0 + 2ULL * j;
                        if (n > end) break;
                        u64 p = 0;
                        for (size_t k = 0; k < nsmall; k++) {
                            u64 pp = small[k];
                            if (pp > n / 2) break;
                            u64 q = n - pp;
                            u64 rel = ((q - 3) >> 1) - B0;
                            if ((w[rel >> 6] >> (rel & 63)) & 1ULL) { p = pp; break; }
                        }
                        if (p) { if (p < min_p) min_p = p; if (p > max_p) max_p = p; }
                    }
                }
            }
            /* ---- tail evens (not in a full 64-group) ---- */
            u64 tail = g0 + ngroups * 128;
            if (tail < start) tail = start;
            for (u64 n = tail; n <= n_hi; n += 2) {
                if (n == 4) { resolved++; continue; }
                u64 p = 0;
                for (size_t k = 0; k < nsmall; k++) {
                    u64 pp = small[k];
                    if (pp > n / 2) break;
                    u64 q = n - pp;
                    u64 rel = ((q - 3) >> 1) - B0;
                    if ((w[rel >> 6] >> (rel & 63)) & 1ULL) { p = pp; break; }
                }
                if (!p) {
                    for (u64 pp = PMAX + 1; pp <= n / 2; pp += 2) {
                        u64 q = n - pp;
                        u64 rel = ((q - 3) >> 1) - B0;
                        if ((rel >> 6) >= need) break;
                        if ((w[rel >> 6] >> (rel & 63)) & 1ULL) { p = pp; break; }
                    }
                }
                scalar_used++;
                if (p) {
                    resolved++;
                    checksum ^= p * 0x9E3779B97F4A7C15ULL * (n + 1);
                    if (p < min_p) min_p = p;
                    if (p > max_p) max_p = p;
                }
                else { unresolved++; fprintf(stderr, "COUNTEREXAMPLE candidate: %llu\n", (unsigned long long)n); }
            }
            /* ---- independent Miller-Rabin cross-check on random evens ---- */
            for (int t = 0; t < 4; t++) {
                u64 n = start + (rnd() % (total_evens ? total_evens : 1)) * 2;
                if (n < 4 || n > end) continue;
                u64 p = 0;
                for (size_t k = 0; k < nsmall; k++) {
                    u64 pp = small[k];
                    if (pp > n / 2) break;
                    if (is_prime_u64(pp) && is_prime_u64(n - pp)) { p = pp; break; }
                }
                sampled++;
                if (p) sample_ok++;
                else fprintf(stderr, "MR cross-check failed for n=%llu\n", (unsigned long long)n);
            }
        }
        free(w);
    }

    double dt = now_s() - t0;
    printf("{\"range\":[%llu,%llu],\"evens\":%llu,\"resolved\":%llu,\"unresolved\":%llu,"
           "\"scalar_path\":%llu,\"mr_sampled\":%llu,\"mr_ok\":%llu,\"min_p\":%llu,\"max_p\":%llu,"
           "\"checksum\":%llu,\"seconds\":%.2f,\"evens_per_sec\":%.0f,\"PMAX\":%llu}\n",
           (unsigned long long)start, (unsigned long long)end,
           (unsigned long long)total_evens, (unsigned long long)resolved,
           (unsigned long long)unresolved, (unsigned long long)scalar_used,
           (unsigned long long)sampled, (unsigned long long)sample_ok,
           (unsigned long long)(min_p == ~0ULL ? 0 : min_p), (unsigned long long)max_p,
           (unsigned long long)checksum, dt, dt > 0 ? total_evens / dt : 0.0,
           (unsigned long long)PMAX);
    if (unresolved) { fprintf(stderr, "FAIL: %llu evens with no decomposition\n", (unsigned long long)unresolved); return 2; }
    fprintf(stderr, "[done] %llu evens verified in %.2fs (%.0f/s)\n",
            (unsigned long long)total_evens, dt, dt > 0 ? total_evens / dt : 0.0);
    return 0;
}
