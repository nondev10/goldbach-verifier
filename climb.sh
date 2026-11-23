#!/usr/bin/env bash
# climb.sh — continuously extend the verified Goldbach range, one power of
# ten at a time, forever.
#
# Each stage verifies the NEW range only (10^(k-1) < n <= 10^k); the union of
# all stages is the verified interval.  Every stage writes:
#   logs/run_1eNN.json   machine-readable result (range + counts + checksum)
#   logs/run_1eNN.err   progress log
#   logs/climb.log      stage index
#
# Env:
#   THREADS (default 2)
#   PMAX    (default 2000)
#   FROM    (default 11)     exponent to start at
#   BIN     (default ./src/gb4)

set -u
cd /root/goldbach
BIN=${BIN:-./src/gb4}
THREADS=${THREADS:-2}
PMAX=${PMAX:-2000}
FROM=${FROM:-11}

if [ ! -x "$BIN" ]; then
  echo "missing $BIN" >&2
  exit 1
fi
mkdir -p logs

echo "[climb] start $(date -Is) bin=$BIN threads=$THREADS pmax=$PMAX from=1e$FROM" | tee -a logs/climb.log

exp=$FROM
while true; do
  N=$(python3 -c "print(10**$exp)")
  LO=$(python3 -c "print(10**($exp-1)+2)")
  tag="1e$exp"
  if [ -s "logs/run_${tag}.json" ] && grep -q '"unresolved":0' "logs/run_${tag}.json" 2>/dev/null; then
    echo "[climb] $tag already verified, skipping" | tee -a logs/climb.log
    exp=$((exp+1))
    continue
  fi
  echo "[climb] === $tag : $LO .. $N  start $(date -Is) ===" | tee -a logs/climb.log
  t0=$(date +%s)
  "$BIN" "$LO" "$N" "$THREADS" "$PMAX" > "logs/run_${tag}.json" 2> "logs/run_${tag}.err"
  rc=$?
  t1=$(date +%s)
  dur=$((t1 - t0))
  if [ $rc -eq 0 ] && grep -q '"unresolved":0' "logs/run_${tag}.json"; then
    echo "[climb] $tag VERIFIED in ${dur}s  $(head -c 400 logs/run_${tag}.json)" | tee -a logs/climb.log
  else
    echo "[climb] $tag FAILED rc=$rc after ${dur}s" | tee -a logs/climb.log
    tail -5 "logs/run_${tag}.err" | tee -a logs/climb.log
  fi
  exp=$((exp+1))
done
