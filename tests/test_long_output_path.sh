#!/bin/bash
# Regression test for issue #5: output paths were built in 200-byte buffers, so a long
# Output_Directory overflowed them (a fortify trap on macOS, silent stack corruption without
# fortify). Paths now use PATH_SIZE buffers, and the config refuses a path longer than
# CONFIG_PATH_MAX instead of overflowing.
#
# Usage: tests/test_long_output_path.sh <MCsquare binary>   (run from the repository root)
set -u
BIN=${1:?usage: $0 <MCsquare binary>}
WORK=$(mktemp -d "${TMPDIR:-/tmp}/mcsq_longpath.XXXXXX")
trap 'rm -rf "$WORK"' EXIT
fail=0

# A directory path of about 900 characters: far past the old 200, inside the OS path limit
# (1024 on macOS) with room for the file names MCsquare appends.
long="$WORK"
while [ ${#long} -lt 880 ]; do long="$long/segment_of_a_long_output_directory"; done
mkdir -p "$long"

config() {  # config <output dir> <file>
  sed -e 's/^Num_Primaries .*/Num_Primaries 1000/' \
      -e "s#^Output_Directory .*#Output_Directory $1#" \
      Sample_input_data/smoke_test_config.txt > "$2"
}

# 1. A long but valid output directory must run to completion and write the dose.
config "$long" "$WORK/long.txt"
"$BIN" "$WORK/long.txt" > "$WORK/long.log" 2>&1
rc=$?
if [ "$rc" -ne 0 ] || [ ! -f "$long/Dose.mhd" ]; then
  echo "FAIL: ${#long}-character Output_Directory: rc=$rc, Dose.mhd $( [ -f "$long/Dose.mhd" ] && echo present || echo missing )"
  tail -5 "$WORK/long.log"
  fail=1
else
  echo "ok: ${#long}-character Output_Directory ran and wrote Dose.mhd"
fi

# 2. An output directory longer than CONFIG_PATH_MAX (1024) must be refused cleanly, not crash.
over=$(printf 'x%.0s' $(seq 1 1025))
config "$over" "$WORK/over.txt"
"$BIN" "$WORK/over.txt" > "$WORK/over.log" 2>&1
rc=$?
if [ "$rc" -ge 128 ] || ! grep -q 'Output_Directory value must contain between' "$WORK/over.log"; then
  echo "FAIL: 1025-character Output_Directory: rc=$rc, expected a clean refusal"
  tail -5 "$WORK/over.log"
  fail=1
else
  echo "ok: 1025-character Output_Directory refused (rc=$rc)"
fi

exit $fail
