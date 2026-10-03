#!/bin/bash
# Halo check without new runs (Fable test 1): ring fractions from the TOPAS halo-addendum runs, Dose vs DoseAll.
# Dose excludes neutrons, gammas and their descendants; DoseAll includes them. Reads only; writes JSON lines here.
set -uo pipefail
RUNS="/Volumes/T7 Shield/SMBWritable/topas/halo_addendum/5aedf3ac49d7/runs"
EP=/Users/stuartswerdloff/ai/liberated/kimi-kindled/kindled_projects/MCsquare-portable/validation/pencil_endpoints.py
OUT=/Users/stuartswerdloff/ai/ClaudeInstanceHomeOffices/clement-7074f29f/mc_validation/stepsize/doseall_rings.jsonl
: > "$OUT"
for spec in "100:40,60:912001 912002 912003 912004" "150:80,125:912011 912012 912013 912014" "200:100,200:912021 912022 912023 912024"; do
  E=${spec%%:*}; rest=${spec#*:}; DEPTHS=${rest%%:*}; SEEDS=${rest#*:}
  for SEED in $SEEDS; do
    D="$RUNS/E${E}_opt0_seed${SEED}_n10000000_th12_a1"
    for F in dose dose_all; do
      L=$(/opt/homebrew/bin/uv run --no-project --with numpy --with scipy python "$EP" topas "$D/$F.bin" --depths "$DEPTHS" --label "E${E}_${SEED}_${F}" 2>&1 | /usr/bin/grep '^ENDPOINTS ' | cut -c11-)
      [ -n "$L" ] || { echo "no endpoints for $D/$F.bin"; exit 4; }
      echo "$L" >> "$OUT"
      echo "$(date '+%T') E=$E seed=$SEED $F ok"
    done
  done
done
echo "ALL DONE $(date '+%T')"
