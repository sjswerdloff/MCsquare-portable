#!/bin/bash
# Nuclear elastic off (see PREDICTION_elastic.md). Sequential, nice 19, same pinned archive as the step-size tests.
#   run_elastic.sh [threads=12]
set -uo pipefail
THREADS=${1:-12}
TOP="/Volumes/T7 Shield/SMBWritable/topas/stepsize/6c4aaff81994"
T="$TOP/validation/topas"
[ -f "$T/make_run.sh" ] || { echo "pinned archive missing at $T"; exit 3; }
ROOT="$TOP/elastic_off"
mkdir "$ROOT" || { echo "refusing: $ROOT exists"; exit 3; }
N=10000000
FROZEN='sv:Ph/Default/Modules = 6 "g4em-standard_opt0" "g4h-phy_QGSP_BIC_HP" "g4decay" "g4ion-binarycascade" "g4h-elastic_HP" "g4stopping"'
NOEL='sv:Ph/Default/Modules = 5 "g4em-standard_opt0" "g4h-phy_QGSP_BIC_HP" "g4decay" "g4ion-binarycascade" "g4stopping"'
for SEED in 914201 914202; do
  RT=$(bash "$T/make_run.sh" 200 opt0 "$SEED" "$N" "$THREADS" "$ROOT/E200_noelastic") || { echo "make_run failed $SEED"; exit 4; }
  D=$(dirname "$RT")
  [ "$(grep -cxF "$FROZEN" "$D/run.txt")" = 1 ] || { echo "frozen module line not found once in $D/run.txt"; exit 4; }
  python3 - "$D/run.txt" "$FROZEN" "$NOEL" <<'PY' || { echo "module edit failed"; exit 4; }
import sys
p, a, b = sys.argv[1:]
s = open(p).read()
assert s.count(a + "\n") == 1
open(p, "w").write(s.replace(a + "\n", b + "\n"))
PY
  [ "$(grep -cxF "$NOEL" "$D/run.txt")" = 1 ] && [ "$(grep -c 'g4h-elastic' "$D/run.txt")" = 0 ] || { echo "edit not verified"; exit 4; }
  start=$(date +%s)
  nice -n 19 bash "$T/run_topas.sh" "$D"; rc=$?
  echo "$(date '+%F %T') E=200 arm=noelastic seed=$SEED rc=$rc wall=$(( $(date +%s) - start ))s"
  [ "$rc" = 0 ] || { echo "run_topas.sh returned $rc for $SEED: stopping"; exit 5; }
done
echo "ALL DONE $(date '+%F %T')"
