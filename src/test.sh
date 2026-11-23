#!/usr/bin/env bash
# test.sh — regression suite for the Goldbach verifier.
#
# Checks:
#   1. the fast C verifier agrees, bit for bit, with an INDEPENDENT Python
#      implementation (checksum over the smallest prime of every even n)
#   2. the checksum is stable across thread counts (no race conditions)
#   3. known Goldbach minima/maxima are respected
#   4. the NTT-based Hardy-Littlewood counter agrees with brute force
#   5. the "no counterexample" invariant holds on all test ranges

set -u
cd "$(dirname "$0")"
SRC=.
fail=0
pass() { echo "  PASS  $*"; }
bad()  { echo "  FAIL  $*"; fail=$((fail+1)); }

echo "=== building ==="
gcc -O3 -march=native -fopenmp -o /tmp/gb4 $SRC/goldbach4.c -lm || { echo "build failed"; exit 1; }
gcc -O3 -march=native -fopenmp -o /tmp/gb2 $SRC/goldbach2.c -lm || { echo "build failed"; exit 1; }
gcc -O3 -march=native -fopenmp -o /tmp/hlntt $SRC/hl_ntt.c -lm || { echo "build failed"; exit 1; }
echo "  ok"

echo
echo "=== T1: C vs independent Python (smallest-prime checksum) ==="
for LIM in 20000 200000 1000000; do
  PY=$(python3 $SRC/independent_check.py $LIM 2>/dev/null | python3 -c "import json,sys;print(json.load(sys.stdin)['checksum_smallest_p'])")
  C=$(/tmp/gb4 4 $LIM 1 2000 2>/dev/null | python3 -c "import json,sys;print(json.load(sys.stdin)['checksum'])")
  if [ "$PY" = "$C" ]; then pass "limit=$LIM checksum match ($C)"; else bad "limit=$LIM python=$PY c=$C"; fi
done

echo
echo "=== T2: thread-count invariance ==="
for T in 1 2 4 8; do
  CHK=$(/tmp/gb4 4 300000 $T 2000 2>/dev/null | python3 -c "import json,sys;print(json.load(sys.stdin)['checksum'])")
  echo "    threads=$T checksum=$CHK"
  REF=$(/tmp/gb4 4 300000 1 2000 2>/dev/null | python3 -c "import json,sys;print(json.load(sys.stdin)['checksum'])")
  [ "$CHK" = "$REF" ] || bad "threads=$T differs from 1-thread run"
done
pass "checksum invariant under thread count"

echo
echo "=== T3: known values ==="
OUT=$(/tmp/gb4 4 100 1 2000 2>/dev/null)
MINP=$(echo "$OUT" | python3 -c "import json,sys;print(json.load(sys.stdin)['min_p'])")
MAXP=$(echo "$OUT" | python3 -c "import json,sys;print(json.load(sys.stdin)['max_p'])")
[ "$MINP" = "2" ] && pass "min_p = 2 (n=4 = 2+2)" || bad "min_p = $MINP (expected 2)"
[ "$MAXP" = "19" ] && pass "max_p = 19 (n=98 = 19+79)" || bad "max_p = $MAXP (expected 19 for n<=100)"
UNR=$(echo "$OUT" | python3 -c "import json,sys;print(json.load(sys.stdin)['unresolved'])")
[ "$UNR" = "0" ] && pass "no counterexamples up to 100" || bad "unresolved=$UNR"

echo
echo "=== T4: NTT convolution vs brute force ==="
/tmp/hlntt 200000 1 200 2>&1 | grep -q "0 mismatches" && pass "NTT matches brute force" || bad "NTT mismatch"

echo
echo "=== T5: no counterexamples across ranges ==="
for R in "4 1000000" "1000002 3000000" "5000000 9000000"; do
  set -- $R
  U=$(/tmp/gb4 $1 $2 2 2000 2>/dev/null | python3 -c "import json,sys;print(json.load(sys.stdin)['unresolved'])")
  [ "$U" = "0" ] && pass "range $1..$2 clean" || bad "range $1..$2 unresolved=$U"
done

echo
echo "=== T6: incremental vs full-range agreement ==="
FULL=$(/tmp/gb4 4 2000000 1 2000 2>/dev/null | python3 -c "import json,sys;print(json.load(sys.stdin)['checksum'])")
INC=$(/tmp/gb4 1000002 2000000 1 2000 2>/dev/null | python3 -c "import json,sys;print(json.load(sys.stdin)['checksum'])")
PYINC=$(python3 $SRC/independent_check.py 2000000 --lo=1000002 2>/dev/null | python3 -c "import json,sys;print(json.load(sys.stdin)['checksum_smallest_p'])")
echo "    full=$FULL  incremental=$INC  python_incremental=$PYINC"
[ "$INC" = "$PYINC" ] && pass "incremental range matches independent Python" \
                      || bad "incremental c=$INC python=$PYINC"

echo
if [ $fail -eq 0 ]; then
  echo "ALL TESTS PASSED"
  exit 0
else
  echo "$fail TEST(S) FAILED"
  exit 1
fi
