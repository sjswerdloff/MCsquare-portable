#!/bin/bash
# Step-size follow-up 2 (see PREDICTION_followup2.md). Sequential, nice 19, same pinned archive as the first test.
#   run_followup.sh [threads=12]
set -uo pipefail
THREADS=${1:-12}
TOP="/Volumes/T7 Shield/SMBWritable/topas/stepsize/6c4aaff81994"
T="$TOP/validation/topas"
[ -f "$T/make_run.sh" ] || { echo "pinned archive missing at $T"; exit 3; }
ROOT="$TOP/followup2"
mkdir "$ROOT" || { echo "refusing: $ROOT exists"; exit 3; }
N=1000000
# energy:arm:maxstep_mm(or -):seed
SPECS="100:ms0p05:0.05:913131 100:ms0p05:0.05:913132
150:ms0p05:0.05:913231 150:ms0p05:0.05:913232"
for spec in $SPECS; do
  IFS=: read -r E ARM MS SEED <<< "$spec"
  RT=$(bash "$T/make_run.sh" "$E" opt0 "$SEED" "$N" "$THREADS" "$ROOT/E${E}_$ARM") || { echo "make_run failed $spec"; exit 4; }
  D=$(dirname "$RT")
  if [ "$MS" != "-" ]; then echo "d:Ge/Phantom/MaxStepSize = $MS mm" >> "$D/run.txt" || exit 4; fi
  start=$(date +%s)
  nice -n 19 bash "$T/run_topas.sh" "$D"; rc=$?
  echo "$(date '+%F %T') E=$E arm=$ARM maxstep=$MS seed=$SEED rc=$rc wall=$(( $(date +%s) - start ))s"
  [ "$rc" = 0 ] || { echo "run_topas.sh returned $rc for $spec: stopping"; exit 5; }
done
echo "ALL DONE $(date '+%F %T')"
