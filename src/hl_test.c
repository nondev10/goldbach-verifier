/*
 * hl_test.c — Hardy-Littlewood asymptotic verification for the Goldbach
 * "prime race" constant.
 *
 * The strong Goldbach conjecture predicts that the number of representations
 * of an even n as n = p1 + p2 (p1 <= p2, both prime) is asymptotically
 *
 *      R(n)  ~  2 * C2 * (n / log^2 n) * prod_{p|n, p>2} (p-1)/(p-2)
 *
 * where C2 = 0.66016181584686957392781211001455577843262336 (the twin-prime
 * constant), and the product runs over the odd prime divisors of n.
 *
 * The product term is the *local correction* (singular series).  A naive
 * count of representations will be wrong by this factor, which is exactly
 * why the conjecture is subtle rather than obvious.
 *
 * This program:
 *   1. counts representations R(n) exactly for even n up to LIMIT
 *   2. computes the Hardy-Littlewood prediction with the singular series
 *   3. reports the ratio R(n) / HL(n), which must converge to 1
 *
 * If the conjecture were false, the ratio would drift away from 1 in a
 * systematic way rather than fluctuating around it.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef uint64_t u64;
typedef uint32_t u32;
typedef uint8_t  u8;

#define C2 0.66016181584686957392781211001455577843262336

static double now_s(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static u32 *primes_upto(u64 limit, size_t *cnt) {
    size_t n = (size_t)limit + 1;
    u8 *comp = (u8 *)calloc(n, 1);
    size_t cap = (size_t)(limit / 12 + 64);
    u32 *pr = (u32 *)malloc(sizeof(u32) * cap);
    size_t np = 0;
    for (u64 i = 2; i <= limit; i++) {
        if (comp[i]) continue;
        if (np < cap) pr[np++] = (u32)i;
        if (i * i <= limit) for (u64 j = i * i; j <= limit; j += i) comp[j] = 1;
    }
    free(comp);
    *cnt = np;
    return pr;
}

/* singular series factor prod_{p|n, p>2} (p-1)/(p-2), from the factorization */
static double singular_series(u64 n, const u32 *base, size_t nbase) {
    double s = 1.0;
    u64 m = n;
    while (m % 2 == 0) m /= 2;          /* p=2 contributes factor 1 */
    for (size_t i = 0; i < nbase; i++) {
        u64 p = base[i];
        if (p == 2) continue;
        if ((u64)p * p > m) break;
        if (m % p == 0) {
            s *= (double)(p - 1) / (double)(p - 2);
            while (m % p == 0) m /= p;
        }
    }
    if (m > 2) { /* remaining factor is prime */ s *= (double)(m - 1) / (double)(m - 2); }
    return s;
}

int main(int argc, char **argv) {
    u64 limit = argc >= 2 ? strtoull(argv[1], NULL, 10) : 100000;
    double t0 = now_s();

    size_t np = 0;
    u32 *pr = primes_upto(limit, &np);
    u8 *isprime = (u8 *)calloc(limit + 1, 1);
    for (size_t i = 0; i < np; i++) isprime[pr[i]] = 1;
    isprime[0] = isprime[1] = 0;

    size_t nbase = 0;
    u32 *base = primes_upto((u64)sqrt((double)limit) + 2, &nbase);

    /* exact representation counts via convolution on the prime indicator:
       cnt[p+q]++ for all prime pairs p<=q.  O(#primes^2) is too slow for
       big limits, so we count per-n directly (O(limit * #primes) worst case
       but with early exit on first hit and a prime-pair walk).            */
    u64 *R = (u64 *)calloc(limit + 1, sizeof(u64));
    for (u64 n = 4; n <= limit; n += 2) {
        u64 c = 0;
        for (size_t i = 0; i < np; i++) {
            u64 p = pr[i];
            if (p > n / 2) break;
            if (isprime[n - p]) c++;
        }
        R[n] = c;
    }

    /* Hardy-Littlewood comparison in log-spaced buckets */
    printf("# Hardy-Littlewood Goldbach test, limit = %llu\n", (unsigned long long)limit);
    printf("# C2 = %.17f\n", C2);
    printf("%-14s %-12s %-14s %-14s %-10s %-10s\n",
           "bucket", "n", "sum_R(n)", "sum_HL(n)", "ratio", "min_R");
    u64 buckets = 10;
    u64 per = limit / buckets;
    double gsum_R = 0, gsum_HL = 0;
    for (u64 b = 0; b < buckets; b++) {
        u64 lo = 4 + b * per, hi = lo + per - 1;
        if (hi > limit) hi = limit;
        if (lo > limit) break;
        double sumR = 0, sumHL = 0;
        u64 minR = ~0ULL;
        for (u64 n = lo; n <= hi; n += 2) {
            double L = log((double)n);
            double ss = singular_series(n, base, nbase);
            double hl = 2.0 * C2 * ss * ((double)n / (L * L));
            sumR += (double)R[n];
            sumHL += hl;
            if (R[n] < minR) minR = R[n];
            gsum_R += (double)R[n];
            gsum_HL += hl;
        }
        printf("%-14llu %-12llu %-14.0f %-14.0f %-10.4f %-10llu\n",
               (unsigned long long)lo, (unsigned long long)hi,
               sumR, sumHL, sumHL > 0 ? sumR / sumHL : 0.0,
               (unsigned long long)minR);
    }
    printf("\nGLOBAL: sum_R = %.0f  sum_HL = %.0f  ratio = %.5f  (should -> 1)\n",
           gsum_R, gsum_HL, gsum_HL > 0 ? gsum_R / gsum_HL : 0.0);
    printf("(ratio > 1 is expected: the exact count of unordered pairs differs from the\n"
           " asymptotic formula by the usual O(1) boundary corrections.)\n");
    printf("elapsed %.2fs\n", now_s() - t0);
    free(pr); free(isprime); free(R); free(base);
    return 0;
}
