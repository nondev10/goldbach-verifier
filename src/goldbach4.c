/*
 * goldbach4.c — corrected + optimized Goldbach verifier.
 *
 * Fix relative to goldbach2.c:
 *   The bit-parallel group path previously restricted the prime p to
 *   p <= n0-3 (so that ALL 64 elements of the group had a valid complement).
 *   That is correct for "does a representation exist", but it makes the
 *   recorded SMALLEST prime wrong for elements further along the group
 *   (e.g. n=12 was reported as unresolved instead of 12 = 5 + 7).
 *   Here the complement q_j = n_j - p is masked per bit, so every prime
 *   p <= PMAX is tested for every element and the recorded smallest prime
 *   is the true one.
 *
 * Method:
 *   - odd-only segmented sieve  ->  exact prime bitmap (no heuristics)
 *   - bit-parallel resolution  ->  64 even numbers decided per shifted load
 *   - exact fallback with deterministic Miller-Rabin for p > PMAX
 *   - independent MR cross-check on random evens
 *
 * Usage: ./goldbach4 <start> <end> [threads] [PMAX]
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
typedef int64_t  i64;
typedef uint32_t u32;
typedef uint8_t  u8;

#define ODDS_PER_SEG (1u << 23)
#define GUARD_W    64u
#define GOLDEN     0x9E3779B97F4A7C15ULL

static double now_s(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

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

static u64 rng_s = 0x243F6A8885A308D3ULL;
static inline u64 rnd(void) { rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s; }

/* 64 bits of the bitmap starting at bit offset `off` (off may be negative:
   out-of-window positions read as 0 = "not prime", which is exactly right
   because those complements are < 3). */
static inline u64 load64(const u64 *wd, i64 off) {
    if (off < 0) {
        i64 k = -off;
        if (k >= 64) return 0;
        /* bits at absolute positions -k..-1 are invalid (complement < 3) -> 0 */
        return load64(wd, 0) >> k;
    }
    u64 wi = (u64)off >> 6, sh = (u64)off & 63;
    u64 m = wd[wi] >> sh;
    if (sh) m |= wd[wi + 1] << (64 - sh);
    return m;
}

int main(int argc, char **argv) {
    u64 start = 4, end = 1000000;
    int threads = 2;
    u64 PMAX = 2000;
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

    u64 T = (end - 3) / 2 + 1;
    u64 nseg = (T + ODDS_PER_SEG - 1) / ODDS_PER_SEG;
    u64 seg_begin = 0;
    if (start > 2 * ODDS_PER_SEG + 6) {
        u64 idx0 = (start - 6) / 2;
        seg_begin = idx0 / ODDS_PER_SEG;
    }
    u64 total_evens = (end - start) / 2 + 1;
    u64 resolved = 0, unresolved = 0, fallbacks = 0, sampled = 0, sample_ok = 0;
    u64 min_p = ~0ULL, max_p = 0, checksum = 0;
    u64 *thread_chk = (u64 *)calloc(threads > 0 ? threads : 1, sizeof(u64));

    if (start <= 4 && 4 <= end) {
        resolved++;
        if (2 < min_p) min_p = 2;
        checksum ^= 2ULL * GOLDEN * 5ULL;   /* n=4 -> p_min = 2 */
    }

    #pragma omp parallel num_threads(threads) \
        reduction(+:resolved,unresolved,fallbacks,sampled,sample_ok) \
        reduction(min:min_p) reduction(max:max_p)
    {
        u64 my_chk = 0;                 /* per-thread XOR accumulator */
        int tid = omp_get_thread_num();
        u64 *w = NULL; size_t wcap = 0;
        #pragma omp for schedule(dynamic, 1)
        for (u64 s = seg_begin; s < nseg; s++) {
            u64 i0 = s * ODDS_PER_SEG;
            u64 i1 = i0 + ODDS_PER_SEG; if (i1 > T) i1 = T;
            u64 margin = PMAX / 2 + 8;
            u64 B0 = (i0 > margin) ? (i0 - margin) : 0;
            u64 span = i1 - B0;
            size_t need = (size_t)((span + 63) / 64) + 4;
            if (need + GUARD_W + 2 > wcap) {
                free(w);
                wcap = need + GUARD_W + 2;
                w = (u64 *)calloc(wcap, sizeof(u64));
            }
            u64 *wd = w + GUARD_W;
            memset(wd, 0xFF, need * sizeof(u64));
            { u64 fw = span / 64, rem = span % 64;
              if (rem) wd[fw] &= (1ULL << rem) - 1;
              for (size_t k = fw + (rem ? 1 : 0); k < need; k++) wd[k] = 0; }

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

            u64 n_lo = 2 * i0 + 6;
            u64 n_hi = 2 * i1 + 4;
            if (n_lo < start) n_lo = start;
            if (n_hi > end)   n_hi = end;
            if (n_lo > n_hi) continue;

            u64 g0 = n_lo < 6 ? 6 : n_lo;
            u64 span_ev = (n_hi >= g0) ? (n_hi - g0) : 0;
            u64 ngroups = (span_ev >= 126) ? ((span_ev - 126) / 128 + 1) : 0;

            /* ---- bit-parallel groups of 64 consecutive evens ---- */
            for (u64 gi = 0; gi < ngroups; gi++) {
                u64 n0 = g0 + gi * 128;
                u64 unres = ~0ULL;                 /* bit j: even n0+2j still uncovered */
                u64 first_p[64]; memset(first_p, 0, sizeof(first_p));
                for (size_t k = 0; k < nsmall; k++) {
                    u64 p = small[k];              /* odd primes, ascending */
                    /* complement of n_j is q_j = n0 + 2j - p; valid iff q_j >= 3.
                       Work in signed arithmetic to avoid unsigned underflow. */
                    i64 q0 = (i64)n0 - (i64)p;
                    u64 jmin = 0;
                    if (q0 < 3) jmin = (u64)((3 - q0 + 1) >> 1);   /* ceil((3-q0)/2) */
                    if (jmin > 63) continue;      /* no element of this group is valid */
                    /* bitmap index of q_{jmin} is (q_{jmin}-3)/2 */
                    i64 qidx = (i64)(((i64)n0 + 2 * (i64)jmin - (i64)p - 3) >> 1);
                    i64 off = qidx - (i64)B0;
                    if (off < 0) {
                        /* complement below the bitmap window: those elements have
                           q < 3, i.e. are already excluded by jmin, so this can
                           only happen for jmin == 0 => skip the prime entirely */
                        continue;
                    }
                    u64 m = load64(wd, off);
                    if (!m) continue;
                    u64 unres_shr = unres >> jmin;
                    u64 newly = m & unres_shr;
                    if (newly) {
                        while (newly) {
                            int t = __builtin_ctzll(newly);
                            newly &= newly - 1;
                            int j = jmin + t;
                            if (j < 64 && !first_p[j]) first_p[j] = p;
                        }
                    }
                    unres &= ~(m << jmin);
                    if (!unres) break;
                }
                resolved += 64 - __builtin_popcountll(unres);
                for (int t = 0; t < 64; t++) {
                    u64 n = n0 + 2ULL * t;
                    if (n > end) break;
                    if (!first_p[t]) {
                        /* exact fallback: scan p > PMAX with Miller-Rabin */
                        fallbacks++;
                        u64 p = 0;
                        for (u64 pp = PMAX + 1; pp <= n / 2; pp += 2) {
                            if (is_prime_u64(pp) && is_prime_u64(n - pp)) { p = pp; break; }
                        }
                        if (!p) {
                            unresolved++;
                            fprintf(stderr, "COUNTEREXAMPLE candidate: %llu\n", (unsigned long long)n);
                            continue;
                        }
                        first_p[t] = p;
                    }
                    my_chk ^= first_p[t] * GOLDEN * (n + 1);
                    if (first_p[t] < min_p) min_p = first_p[t];
                    if (first_p[t] > max_p) max_p = first_p[t];
                }
            }

            /* ---- tail evens (not in a full 64-group) ---- */
            u64 tail = g0 + ngroups * 128;
            if (tail < start) tail = start;
            for (u64 n = tail; n <= n_hi; n += 2) {
                if (n == 4) {
                    resolved++;
                    my_chk ^= 2ULL * GOLDEN * 5ULL;
                    if (2 < min_p) min_p = 2;
                    continue;
                }
                u64 p = 0;
                for (size_t k = 0; k < nsmall; k++) {
                    u64 pp = small[k];
                    u64 q = n - pp;
                    if (q < 3) break;              /* n-pp decreasing in pp */
                    i64 off = (i64)((q - 3) >> 1) - (i64)B0;
                    if (off < 0) continue;
                    if (load64(wd, off) & 1ULL) { p = pp; break; }
                }
                if (!p) {
                    fallbacks++;
                    for (u64 pp = PMAX + 1; pp <= n / 2; pp += 2) {
                        if (is_prime_u64(pp) && is_prime_u64(n - pp)) { p = pp; break; }
                    }
                }
                if (!p) { unresolved++; fprintf(stderr, "COUNTEREXAMPLE candidate: %llu\n", (unsigned long long)n); continue; }
                resolved++;
                my_chk ^= p * GOLDEN * (n + 1);
                if (p < min_p) min_p = p;
                if (p > max_p) max_p = p;
            }

            /* ---- independent Miller-Rabin cross-check ---- */
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
        thread_chk[tid] = my_chk;
    }
    /* combine per-thread XOR accumulators (XOR is NOT an OpenMP + reduction) */
    for (int t = 0; t < threads; t++) checksum ^= thread_chk[t];
    free(thread_chk);

    double dt = now_s() - t0;
    printf("{\"range\":[%llu,%llu],\"evens\":%llu,\"resolved\":%llu,\"unresolved\":%llu,"
           "\"fallbacks\":%llu,\"mr_sampled\":%llu,\"mr_ok\":%llu,\"min_p\":%llu,\"max_p\":%llu,"
           "\"checksum\":%llu,\"seconds\":%.2f,\"evens_per_sec\":%.0f,\"PMAX\":%llu}\n",
           (unsigned long long)start, (unsigned long long)end,
           (unsigned long long)total_evens, (unsigned long long)resolved,
           (unsigned long long)unresolved, (unsigned long long)fallbacks,
           (unsigned long long)sampled, (unsigned long long)sample_ok,
           (unsigned long long)(min_p == ~0ULL ? 0 : min_p), (unsigned long long)max_p,
           (unsigned long long)checksum, dt, dt > 0 ? total_evens / dt : 0.0,
           (unsigned long long)PMAX);
    if (unresolved) { fprintf(stderr, "FAIL: %llu evens with no decomposition\n", (unsigned long long)unresolved); return 2; }
    fprintf(stderr, "[done] %llu evens verified in %.2fs (%.0f/s)\n",
            (unsigned long long)total_evens, dt, dt > 0 ? total_evens / dt : 0.0);
    return 0;
}
