#!/usr/bin/env python3
"""
independent_check.py — a deliberately NAIVE, independent implementation of the
strong Goldbach check, used to cross-validate the fast C code.

Design goal: maximum independence.  It shares no code, no data structures and
no algorithmic idea with goldbach2.c / goldbach3.c:
  - plain bytearray sieve (not a segmented bitmap)
  - straightforward nested loop over primes (no bit parallelism)
  - for every even n it finds the SMALLEST prime p with n-p prime, using a
    plain early-exit scan over the prime list

The checksum it produces is the same functional one the C code computes:
    chk = XOR over even n of  p_min(n) * 0x9E3779B97F4A7C15 * (n+1)  (mod 2^64)

so the two programs' checksums can be compared directly.

A second mode recomputes the FULL representation count R(n) for a sample of
even n and compares against the NTT-based C program.

Usage:  python3 independent_check.py <limit> [out.json] [--count-every K]
"""
import json
import math
import sys
import time

GOLDEN = 0x9E3779B97F4A7C15
MASK = (1 << 64) - 1


def build_sieve(limit):
    s = bytearray(b"\x01") * (limit + 1)
    s[0] = s[1] = 0
    for i in range(2, int(math.isqrt(limit)) + 1):
        if s[i]:
            start = i * i
            s[start::i] = b"\x00" * (((limit - start) // i) + 1)
    primes = [i for i in range(2, limit + 1) if s[i]]
    return s, primes


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    count_every = 1
    for a in sys.argv[1:]:
        if a.startswith("--count-every"):
            count_every = int(a.split("=")[1])
    limit = int(args[0]) if args else 200000
    outfile = args[1] if len(args) > 1 else None
    lo = 4
    for a in sys.argv[1:]:
        if a.startswith("--lo"):
            lo = int(a.split("=")[1])

    t0 = time.time()
    isp, primes = build_sieve(limit)
    print(f"[init] sieve up to {limit} in {time.time()-t0:.2f}s, pi={len(primes)}", file=sys.stderr, flush=True)

    # primes up to limit/2 are the only ones we ever need
    pmax = 2000
    for a in sys.argv[1:]:
        if a.startswith("--pmax"):
            pmax = int(a.split("=")[1])
    half = [p for p in primes if p <= limit // 2]

    smallest = {}
    counterexamples = []
    chk = 0
    total_repr = 0
    sampled_counts = {}

    for n in range(lo, limit + 1, 2):
        sp = 0
        for p in half:
            if p > n // 2:
                break
            if isp[n - p]:
                sp = p
                break
        if sp == 0:
            counterexamples.append(n)
            smallest[n] = 0
        else:
            smallest[n] = sp
            chk ^= (sp * GOLDEN * (n + 1)) & MASK

        if count_every > 1 and (n // 2) % count_every == 0:
            c = 0
            for p in half:
                if p > n // 2:
                    break
                if isp[n - p]:
                    c += 1
            sampled_counts[n] = c
            total_repr += c

    out = {
        "limit": limit,
        "lo": lo,
        "evens": (limit - lo) // 2 + 1,
        "counterexamples": len(counterexamples),
        "counterexample_list": counterexamples[:20],
        "checksum_smallest_p": chk,
        "min_smallest_p": min(v for v in smallest.values() if v),
        "max_smallest_p": max(v for v in smallest.values() if v),
        "sampled_full_counts": {str(k): v for k, v in list(sampled_counts.items())[:500]},
        "sampled_total_representations": total_repr,
        "seconds": round(time.time() - t0, 3),
    }
    print(json.dumps({k: v for k, v in out.items() if k != "sampled_full_counts"}, indent=2))
    if outfile:
        with open(outfile, "w") as f:
            json.dump(out, f, indent=2)
    return 0 if not counterexamples else 2


if __name__ == "__main__":
    sys.exit(main())
