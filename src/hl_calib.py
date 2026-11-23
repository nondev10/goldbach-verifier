#!/usr/bin/env python3
"""
hl_calib.py — rigorous Hardy–Littlewood calibration for the Goldbach
representation count R(n).

The strong Goldbach conjecture predicts

    R(n) := #{ (p,q) prime, p <= q, p + q = n }
      ~  2 * C2 * (n / log^2 n) * prod_{p | n, p > 2} (p-1)/(p-2)

with C2 = 0.66016181584686957392781211001455577843262336 (twin prime
constant).  The product is the singular series.

This script:
  1. computes R(n) EXACTLY for every even n <= N (via a prime-pair count),
  2. computes the HL prediction with the singular series,
  3. reports the empirical per-n constant K = R(n) * log^2(n) / (n * S(n)),
     which must converge to 2*C2 = 1.32032363169373914785562422002911155686524672,
  4. reports the aggregate identity  sum_{n even <= N} R(n) ~ C2 * N^2 / (2 log^2 N),
     which must converge to ratio 1.

If Goldbach were false for large n, K would collapse towards 0.

Complexity note: R(n) is computed by walking primes p <= n/2 with a lookup,
so the whole sweep is O(N * pi(N)) which is fine to N ~ 1e6..1e7.
"""
import math
import sys
from array import array

C2 = 0.66016181584686957392781211001455577843262336


def sieve(limit):
    """bytearray sieve; byte 1 <=> prime."""
    s = bytearray(b"\x01") * (limit + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(math.isqrt(limit)) + 1):
        if s[i]:
            s[i * i::i] = b"\x00" * ((limit - i * i) // i + 1)
    return s


def primes_from(s):
    return [i for i in range(2, len(s)) if s[i]]


def singular_series(n, small_primes):
    """prod_{p | n, p > 2} (p-1)/(p-2)  — exact rational factors in float."""
    ss = 1.0
    m = n
    while m % 2 == 0:
        m //= 2
    for p in small_primes:
        if p == 2:
            continue
        if p * p > m:
            break
        if m % p == 0:
            ss *= (p - 1) / (p - 2)
            while m % p == 0:
                m //= p
    if m > 2:
        ss *= (m - 1) / (m - 2)
    return ss


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1_000_000
    mode = sys.argv[2] if len(sys.argv) > 2 else "both"

    s = sieve(N + 1)
    primes = primes_from(s)
    small = primes_from(sieve(int(math.isqrt(N)) + 2))
    print(f"# N = {N}   pi(N) = {len(primes)}   C2 = {C2:.17f}", flush=True)

    # ---- exact representation counts -------------------------------------
    R = array("Q", bytes(8 * (N + 1)))
    zero, counter = 0, 0
    for n in range(4, N + 1, 2):
        c = 0
        for p in primes:
            if p > n // 2:
                break
            if s[n - p]:
                c += 1
        R[n] = c
        if c == 0:
            zero += 1
            counter += 1
    print(f"# evens with ZERO representations (would be counterexamples): {zero}", flush=True)

    # ---- per-n empirical constant ----------------------------------------
    if mode in ("both", "pern"):
        print(f"\n{'n':>12} {'R(n)':>10} {'S(n)':>10} {'HL pred':>14} {'R/HL':>8}")
        probes = [N // 8, N // 4, N // 2, N]
        for n in probes:
            if n < 4:
                continue
            L = math.log(n)
            ss = singular_series(n, small)
            hl = 2 * C2 * ss * n / (L * L)
            print(f"{n:>12} {R[n]:>10} {ss:>10.4f} {hl:>14.1f} {R[n]/hl:>8.4f}")

        # band-averaged constant: K -> 2*C2
        num = den = 0.0
        for n in range(N // 2, N + 1, 2):
            L = math.log(n)
            ss = singular_series(n, small)
            num += R[n] * L * L
            den += n * ss
        K = num / den
        print(f"\nband-averaged per-n constant K = {K:.6f}")
        print(f"theoretical 2*C2                = {2*C2:.6f}")
        print(f"relative deviation              = {abs(K - 2*C2)/(2*C2)*100:.3f}%")

    # ---- aggregate identity ---------------------------------------------
    if mode in ("both", "agg"):
        tot = sum(R[n] for n in range(4, N + 1, 2))
        L = math.log(N)
        pred = C2 * N * N / (2 * L * L)
        print(f"\nsum_{{n even <= N}} R(n) = {tot}")
        print(f"HL aggregate C2*N^2/(2 log^2 N) = {pred:.0f}")
        print(f"ratio = {tot/pred:.5f}   (should -> 1)")


if __name__ == "__main__":
    main()
