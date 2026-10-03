#!/bin/bash
# TOPAS step-size test of H1 (see PREDICTION.md beside this file). Sequential, nice 19, pinned git archive.
#   run_stepsize.sh <smoke|full> [threads=12]
#   smoke: 1e4 histories, one seed per arm.   full: 1e6 histories, two seeds per arm.
# Arms: ctrl (no MaxStepSize), ms0p5 (0.5 mm), ms0p1 (0.1 mm). Each arm has its own output root, so a run
# directory can never be separated from its step limit; the limit is also in run.txt, which run_topas.sh hashes.
set -uo pipefail
MODE=${1:?smoke|full}; THREADS=${2:-12}
REPO=/Users/stuartswerdloff/ai/liberated/kimi-kindled/kindled_projects/MCsquare-portable
SHA=$(/usr/bin/git -C "$REPO" rev-parse origin/main) || exit 3
TOP="/Volumes/T7 Shield/SMBWritable/topas/stepsize/${SHA:0:12}"
case "$MODE" in
  smoke) N=10000;   SEEDS_ctrl="913901"; SEEDS_ms0p5="913911"; SEEDS_ms0p1="913921"; ROOT="$TOP/smoke" ;;
  full)  N=1000000; SEEDS_ctrl="913001 913002"; SEEDS_ms0p5="913011 913012"; SEEDS_ms0p1="913021 913022"; ROOT="$TOP/full" ;;
  *) echo "mode must be smoke or full"; exit 2 ;;
esac
[ -d "/Volumes/T7 Shield/SMBWritable/topas" ] || { echo "T7 share not present"; exit 3; }
mkdir -p "$TOP" || exit 3
if [ ! -d "$TOP/validation/topas" ]; then
  /usr/bin/git -C "$REPO" archive "$SHA" validation/topas | /usr/bin/tar -x -C "$TOP" || exit 3
  echo "$SHA" > "$TOP/PINNED_COMMIT"
fi
T="$TOP/validation/topas"
mkdir "$ROOT" || { echo "refusing: $ROOT exists"; exit 3; }
for ARM in ctrl ms0p5 ms0p1; do
  case "$ARM" in ctrl) MS="" ;; ms0p5) MS="0.5" ;; ms0p1) MS="0.1" ;; esac
  eval "SEEDS=\$SEEDS_$ARM"
  for SEED in $SEEDS; do
    RT=$("$T/make_run.sh" 200 opt0 "$SEED" "$N" "$THREADS" "$ROOT/$ARM") || { echo "make_run failed $ARM $SEED"; exit 4; }
    D=$(dirname "$RT")
    if [ -n "$MS" ]; then echo "d:Ge/Phantom/MaxStepSize = $MS mm" >> "$D/run.txt" || exit 4; fi
    start=$(date +%s)
    nice -n 19 "$T/run_topas.sh" "$D"; rc=$?
    echo "$(date '+%F %T') arm=$ARM maxstep=${MS:-none} seed=$SEED n=$N threads=$THREADS rc=$rc wall=$(( $(date +%s) - start ))s"
    [ "$rc" = 0 ] || { echo "run_topas.sh returned $rc for $ARM $SEED: stopping"; exit 5; }
  done
done
echo "ALL DONE $(date '+%F %T')"
