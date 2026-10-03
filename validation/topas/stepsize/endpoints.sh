#!/bin/bash
# Endpoints for every dose.bin under a run root, with the slab depths of each energy (as in the halo addendum).
#   endpoints.sh <run_root> <out.jsonl>     (runs pencil_endpoints.py from the MCsquare-portable worktree)
set -uo pipefail
ROOT="$1"; OUT="$2"
V=/Users/stuartswerdloff/ai/ClaudeInstanceHomeOffices/clement-7074f29f/worktrees/mcsq-stepsize/validation
: > "$OUT"
find "$ROOT" -name dose.bin | sort | while IFS= read -r f; do
  case "$f" in
    */E100_*) dp=3,40,60 ;;
    */E150_*) dp=3,80,125 ;;
    */E200_*) dp=3,100,200 ;;
    *) echo "no energy in $f" >&2; exit 4 ;;
  esac
  (cd "$V" && uv run --quiet --no-project --with numpy --with scipy python pencil_endpoints.py topas "$f" \
     --label "$(basename "$(dirname "$f")")" --depths "$dp") >> "$OUT" || { echo "failed: $f" >&2; exit 5; }
done
