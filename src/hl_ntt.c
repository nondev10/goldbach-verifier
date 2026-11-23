/*
 * hl_ntt.c — exact Hardy–Littlewood calibration via NTT convolution.
 *
 * R(n) := #{ (p,q) prime, p <= q, p+q = n }
 *
 * Computing R(n) for every even n <= N by brute force costs O(N * pi(N)),
 * which is prohibitive.  Instead we use the fact that
 *
 *      C[n] = sum_{p+q=n} 1_prime(p) * 1_prime(q)
 *
 * is exactly the convolution of the prime indicator with itself.  We compute
 * it with a Number Theoretic Transform (NTT) modulo 998244353, which is
 * EXACT (no floating point error at all), so every R(n) is an integer
 * computed without approximation.
 *
 * Then we compare with the Hardy-Littlewood prediction
 *
 *      R(n) ~ K * S(n) * n / log^2 n,   S(n) = prod_{p|n,p>2} (p-1)/(p-2)
 *
 * and estimate the constant K by least squares / band averaging.  Theory
 * predicts K = 2*C2 = 1.3203236316937391...  (C2 = 0.66016181584686957,
 * the twin prime constant).  We report the measured value and the relative
 * deviation.
 *
 * Also reports the aggregate identity
 *      sum_{n even <= N} R(n)  vs  K/2 * sum_{n even <= N} n*S(n)/log^2 n
 *
 * Usage: ./hl_ntt <N> [threads]
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

/* MOD = 15*2^27 + 1 supports NTT lengths up to 2^27 (needed for N = 1e7) */
#define MOD 2013265921u
#define ROOT 31u

/* twin prime constant (Hardy-Littlewood) */
static const double C2 = 0.66016181584686957392781211001455577843262336;

static double now_s(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

/* ---- Gauss-Legendre nodes/weights on [0,1] (Newton iteration on P_K) ---- */
static void gauss_legendre(int K, double *x, double *w) {
    for (int i = 0; i < K; i++) {
        double z = cos(M_PI * (i + 0.75) / (K + 0.5));
        double pp = 0.0;
        for (int it = 0; it < 100; it++) {
            double p0 = 1.0, p1 = 0.0;
            for (int j = 0; j < K; j++) {
                double p2 = p1; p1 = p0;
                p0 = ((2.0 * j + 1.0) * z * p1 - j * p2) / (j + 1);
            }
            pp = K * (z * p0 - p1) / (z * z - 1.0);
            double dz = p0 / pp;
            z -= dz;
            if (fabs(dz) < 1e-16) break;
        }
        x[i] = (1.0 - z) / 2.0;          /* map [-1,1] -> [0,1] */
        w[i] = 1.0 / ((1.0 - z * z) * pp * pp);
    }
}

static u32 modpow(u32 a, u64 e) {
    u64 r = 1, b = a;
    while (e) { if (e & 1) r = r * b % MOD; b = b * b % MOD; e >>= 1; }
    return (u32)r;
}
static void ntt(u32 *a, size_t n, int invert) {
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { u32 t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (size_t len = 2; len <= n; len <<= 1) {
        u32 wlen = modpow(ROOT, (MOD - 1) / len);
        if (invert) wlen = modpow(wlen, MOD - 2);
        for (size_t i = 0; i < n; i += len) {
            u64 w = 1;
            for (size_t k = 0; k < len / 2; k++) {
                u32 u = a[i + k];
                u32 v = (u32)(a[i + k + len / 2] * w % MOD);
                u32 x = u + v; if (x >= MOD) x -= MOD;
                u32 y = (u >= v) ? u - v : u + MOD - v;
                a[i + k] = x; a[i + k + len / 2] = y;
                w = w * wlen % MOD;
            }
        }
    }
    if (invert) {
        u32 inv = modpow((u32)n, MOD - 2);
        for (size_t i = 0; i < n; i++) a[i] = (u32)((u64)a[i] * inv % MOD);
    }
}

int main(int argc, char **argv) {
    u64 N = argc >= 2 ? strtoull(argv[1], NULL, 10) : 10000000;
    int threads = argc >= 3 ? atoi(argv[2]) : 2;
    double t0 = now_s();

    /* --- sieve: isprime[] and smallest-prime-factor[] --- */
    u8 *isprime = (u8 *)calloc(N + 1, 1);
    u32 *spf = (u32 *)malloc(sizeof(u32) * (N + 1));
    memset(spf, 0, sizeof(u32) * (N + 1));   /* 0 = not yet assigned */
    for (u64 i = 2; i <= N; i++) {
        if (spf[i] == 0) {                    /* i is prime */
            isprime[i] = 1;
            spf[i] = (u32)i;
            for (u64 j = i * i; j <= N; j += i)
                if (spf[j] == 0) spf[j] = (u32)i;
        }
    }
    spf[0] = 0; spf[1] = 1;
    isprime[0] = isprime[1] = 0;
    fprintf(stderr, "[init] sieve done %.2fs\n", now_s() - t0);

    /* --- convolution --- */
    size_t sz = 1; while (sz < 2 * (size_t)N + 1) sz <<= 1;
    u32 *A = (u32 *)calloc(sz, sizeof(u32));
    u32 *B = (u32 *)calloc(sz, sizeof(u32));
    for (u64 i = 0; i <= N; i++) { A[i] = isprime[i]; B[i] = isprime[i]; }
    ntt(A, sz, 0); ntt(B, sz, 0);
    for (size_t i = 0; i < sz; i++) A[i] = (u32)((u64)A[i] * B[i] % MOD);
    ntt(A, sz, 1);
    free(B);
    fprintf(stderr, "[conv] NTT size=%zu done %.2fs\n", sz, now_s() - t0);

    /* --- independent brute-force cross-check of the NTT convolution --- */
    {
        u64 bad = 0;
        for (u64 n = 4; n <= 2000 && n <= N; n += 2) {
            u64 ord = 0;
            for (u64 p = 2; p + 2 <= n; p++) if (isprime[p] && isprime[n - p]) ord++;
            u64 diag = isprime[n / 2] ? 1 : 0;
            u64 Rb = (ord + diag) / 2;
            u64 Rn = (A[n] + diag) / 2;
            if (Rb != Rn) { fprintf(stderr, "NTT MISMATCH n=%llu brute=%llu ntt=%llu (A=%llu)\n",
                    (unsigned long long)n, (unsigned long long)Rb, (unsigned long long)Rn,
                    (unsigned long long)A[n]); bad++; }
        }
        fprintf(stderr, "[check] NTT vs brute force on n<=2000: %llu mismatches\n", (unsigned long long)bad);
        if (bad) return 3;
    }

    /* --- Hardy-Littlewood check with the CORRECT integral main term ---
     *
     * The naive formula R(n) ~ C2 * S(n) * n/ln^2 n misses a finite-n
     * correction of relative size ~1/ln n, because the true local density
     * integral is
     *
     *     I_half(n) = INT_2^{n/2} dp / (ln p * ln(n-p))
     *
     * rather than the crude (n/2)/ln^2 n.  We evaluate I_half by
     * Gauss-Legendre quadrature (exact to machine precision for this smooth
     * integrand) and predict
     *
     *     R(n)  ~  2 * C2 * prod_{p|n,p>2}(p-1)/(p-2) * I_half(n)
     *
     * A sample of even n in the top half of the range is then checked.
     */
    {
        const int K = 96;
        static double gx[96], gw[96];
        gauss_legendre(K, gx, gw);
        u64 n_sample = argc >= 4 ? strtoull(argv[3], NULL, 10) : 20000;
        double sum_ratio = 0, sum_ratio2 = 0, sum_log_ratio = 0;
        u64 nchk = 0;
        u64 lo = N / 2;
        u64 step = (N - lo) / (n_sample ? n_sample : 1);
        if (step == 0) step = 2;
        printf("\n# HL per-n check on %llu sampled even n in [%llu,%llu]\n",
               (unsigned long long)n_sample, (unsigned long long)lo, (unsigned long long)N);
        printf("# %-12s %-10s %-12s %-12s %-10s\n", "n", "R(n)", "prod", "I_half", "ratio");
        for (u64 n = (lo + 1) & ~1ULL; n <= N; n += (step & ~1ULL) + 2) {
            u64 ord = A[n];
            u64 diag = isprime[n / 2] ? 1 : 0;
            u64 R = (ord + diag) / 2;
            if (!R) continue;
            /* singular series */
            double ss = 1.0;
            u64 m = n;
            while (m % 2 == 0) m /= 2;
            while (m > 1) {
                u64 p = spf[m];
                if (p == 2) break;
                ss *= (double)(p - 1) / (double)(p - 2);
                while (m % p == 0) m /= p;
            }
            /* I_half by Gauss-Legendre */
            double a = 2.0, b = (double)n / 2.0;
            double I = 0.0;
            for (int k = 0; k < K; k++) {
                double p = a + (b - a) * gx[k];
                double q = (double)n - p;
                if (p <= 2.0) continue;
                I += gw[k] / (log(p) * log(q));
            }
            I *= (b - a);
            double pred = 2.0 * C2 * ss * I;
            if (pred <= 0) continue;
            double ratio = (double)R / pred;
            sum_ratio += ratio; sum_ratio2 += ratio * ratio; sum_log_ratio += log(ratio);
            nchk++;
            if (nchk <= 15 || nchk % (n_sample / 10) == 0)
                printf("# %-12llu %-10llu %-12.4f %-12.1f %-10.4f\n",
                       (unsigned long long)n, (unsigned long long)R, ss, I, ratio);
        }
        double mean = nchk ? sum_ratio / nchk : 0;
        double var = nchk ? sum_ratio2 / nchk - mean * mean : 0;
        printf("\n  samples            = %llu\n", (unsigned long long)nchk);
        printf("  mean R/pred        = %.6f   (theory 1.000000)\n", mean);
        printf("  stddev R/pred      = %.6f\n", var > 0 ? sqrt(var) : 0);
        printf("  mean log(R/pred)   = %.6f  -> geometric mean ratio %.6f\n",
               nchk ? sum_log_ratio / nchk : 0, nchk ? exp(sum_log_ratio / nchk) : 0);
        printf("  NOTE: R(n) counts are integers and small for n~5e6, so the\n"
               "        per-n ratio has an inherent +-1/R(n) noise floor.\n");
    }

    double sumR = 0, sumHL = 0, sumK = 0;
    u64 zero = 0, total = 0;
    double num = 0, den = 0;              /* for band-averaged constant */
    u64 band_lo = N / 2;

    /* per-n: singular series via factorization using spf */
    {
        const int KG = 48;
        static double gx[48], gw[48];
        gauss_legendre(KG, gx, gw);
        #pragma omp parallel for num_threads(threads) schedule(static) \
            reduction(+:sumR,sumHL,zero,total,num,den)
        for (u64 n = 4; n <= N; n += 2) {
            u64 ord = A[n];
            u64 diag = isprime[n / 2] ? 1 : 0;
            u64 R = (ord + diag) / 2;          /* unordered representations */
            total++;
            if (R == 0) { zero++; continue; }

            /* singular series */
            double ss = 1.0;
            u64 m = n;
            while (m % 2 == 0) m /= 2;
            while (m > 1) {
                u64 p = spf[m];
                if (p == 2) break;
                ss *= (double)(p - 1) / (double)(p - 2);
                while (m % p == 0) m /= p;
            }
            /* main term: 2*C2*S(n)*INT_2^{n/2} dp/(ln p ln(n-p)) */
            double a = 2.0, b = (double)n / 2.0;
            double I = 0.0;
            for (int k = 0; k < KG; k++) {
                double p = a + (b - a) * gx[k];
                double q = (double)n - p;
                if (p <= 2.0) continue;
                I += gw[k] / (log(p) * log(q));
            }
            I *= (b - a);
            double hl = 2.0 * C2 * ss * I;
            sumR += (double)R;
            sumHL += hl;
            double L = log((double)n);
            if (n >= band_lo) { num += (double)R * L * L; den += (double)n * ss; }
        }
    }

    double K = den > 0 ? num / den : 0.0;
    printf("{\n");
    printf("  \"N\": %llu,\n", (unsigned long long)N);
    printf("  \"evens\": %llu,\n", (unsigned long long)total);
    printf("  \"evens_with_zero_representations\": %llu,\n", (unsigned long long)zero);
    printf("  \"sum_R\": %.0f,\n", sumR);
    printf("  \"sum_HL_2C2\": %.0f,\n", sumHL);
    printf("  \"aggregate_ratio\": %.6f,\n", sumHL > 0 ? sumR / sumHL : 0.0);
    printf("  \"band_averaged_K\": %.6f,\n", K);
    printf("  \"theory_2C2\": %.6f,\n", 2 * C2);
    printf("  \"theory_C2\": %.6f,\n", C2);
    printf("  \"K_relative_deviation_vs_2C2_pct\": %.4f,\n", fabs(K - 2 * C2) / (2 * C2) * 100.0);
    printf("  \"K_relative_deviation_vs_C2_pct\": %.4f,\n", fabs(K - C2) / C2 * 100.0);
    printf("  \"seconds\": %.2f\n", now_s() - t0);
    printf("}\n");
    if (zero) { fprintf(stderr, "FAIL: %llu evens have zero representations\n", (unsigned long long)zero); return 2; }
    return 0;
}
