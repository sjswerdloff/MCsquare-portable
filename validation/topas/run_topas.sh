#!/bin/bash
# Run ONE stage-1 TOPAS run directory (made by make_run.sh) and record provenance (issue #32).
#
#   run_topas.sh <run_dir>
#
# Launch heavy runs detached:  nohup run_topas.sh <dir> > /dev/null 2>&1 < /dev/null &
# The environment is set inside this script because macOS strips DYLD_* from anything launched through a
# SIP-protected binary (/usr/bin/time, env -i); wrapping THIS script is therefore safe.
#
# Fails closed. Exit status:
#   0  COMPLETE: TOPAS exited 0, all four outputs present with the size the grid implies, and the logged EM
#      physics matches the arm requested in run.txt
#   2  usage / not a run directory
#   3  refused: the directory was already claimed, or holds previous attempt or output state (never clobbered)
#   4  refused: free-space guard (unreadable free space, malformed MIN_FREE_GB, or below it)
#   5  INCOMPLETE: TOPAS exited 0 but an output is missing or has the wrong size
#   6  UNVERIFIED: the logged EM physics is missing, or does not match the requested arm
#   7  TOPAS itself exited nonzero (its status is recorded in provenance.txt as topas_exit)
# A single run is attempted once. A retry is a NEW directory (make_run.sh ... a2), so the records of a failed
# attempt are kept.
#
# TOPAS_BIN overrides the executable (tests use a fake). Whatever runs, its path and sha256 are recorded.
set -uo pipefail

die() { echo "run_topas.sh: $1" >&2; exit "$2"; }

[ $# -eq 1 ] || die "usage: $0 <run_dir>" 2
DIR=$1
[ -d "$DIR" ] && [ -f "$DIR/run.txt" ] && [ -f "$DIR/stage1_base.txt" ] || die "not a run directory: $DIR" 2
cd "$DIR" || die "cannot cd into $DIR" 2

# --- Claim: mkdir is atomic, so of two concurrent invocations exactly one gets past this line.
mkdir .claimed 2>/dev/null || die "refusing: $DIR was already claimed by an earlier or concurrent attempt" 3
for f in provenance.txt topas.log dose.bin dose.binheader dose_all.bin dose_all.binheader; do
    [ -e "$f" ] && die "refusing: $DIR already holds $f from a previous attempt" 3
done

MC=/Users/stuartswerdloff/MonteCarlo
TOPAS_BIN=${TOPAS_BIN:-$MC/opentopas-4.3.0/bin/topas}
export TOPAS_G4_DATA_DIR=$MC/geant4-data-11.3.2
export DYLD_LIBRARY_PATH=$MC/opentopas-4.3.0/lib:$MC/geant4-11.3.2/lib:$MC/gdcm-2.6.8/lib

sha() { /usr/bin/shasum -a 256 "$1" | /usr/bin/cut -d' ' -f1; }
# noclobber: provenance.txt is created here or not at all.
set -o noclobber
{
    echo "host: $(/bin/hostname)"
    echo "start: $(/bin/date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "runner: $0"
    echo "runner sha256: $(sha "$0")"
    echo "topas_bin: $TOPAS_BIN"
    echo "topas_bin sha256: $( [ -f "$TOPAS_BIN" ] && sha "$TOPAS_BIN" || echo MISSING)"
    echo "geant4_data: $TOPAS_G4_DATA_DIR"
    echo "run.txt sha256: $(sha run.txt)"
    echo "stage1_base.txt sha256: $(sha stage1_base.txt)"
} > provenance.txt || die "cannot create provenance.txt" 3
set +o noclobber

# --- Free-space guard, fails closed.
MIN_FREE_GB=${MIN_FREE_GB-100}   # default only when UNSET: an explicit empty value is malformed and refused
if ! [[ "$MIN_FREE_GB" =~ ^[0-9]{1,6}$ ]]; then
    echo "verdict: REFUSED malformed MIN_FREE_GB='$MIN_FREE_GB'" >> provenance.txt
    die "refusing: MIN_FREE_GB must be a whole number of GB, got '$MIN_FREE_GB'" 4
fi
FREE_GB=$(/bin/df -g . | /usr/bin/awk 'NR==2 {print $4}')
if ! [[ "$FREE_GB" =~ ^[0-9]+$ ]]; then
    echo "verdict: REFUSED unreadable free space '$FREE_GB'" >> provenance.txt
    die "refusing: could not read free space for $DIR (got '$FREE_GB')" 4
fi
if (( FREE_GB < MIN_FREE_GB )); then
    echo "verdict: REFUSED ${FREE_GB} GB free < MIN_FREE_GB=${MIN_FREE_GB}" >> provenance.txt
    die "refusing: ${FREE_GB} GB free on the volume holding $DIR, below MIN_FREE_GB=${MIN_FREE_GB}" 4
fi

# --- What this run must produce, derived from ITS OWN frozen inputs rather than assumed.
bins() { /usr/bin/awk -F= -v k="i:Ge/Phantom/$1" '{g=$1; gsub(/[ \t]/,"",g)} g==k {v=$2; gsub(/[ \t]/,"",v); print v}' stage1_base.txt; }
NX=$(bins XBins); NY=$(bins YBins); NZ=$(bins ZBins)
for v in "$NX" "$NY" "$NZ"; do [[ "$v" =~ ^[1-9][0-9]*$ ]] || die "cannot read the scoring grid from stage1_base.txt" 2; done
EXPECT_BYTES=$(( NX * NY * NZ * 8 ))   # one double per voxel; the base reports "Sum" only
if /usr/bin/grep -q '"g4em-standard_opt0"' run.txt; then EXPECT_ICRU90=0
elif /usr/bin/grep -q '"g4em-standard_opt4"' run.txt; then EXPECT_ICRU90=1
else die "run.txt requests neither opt0 nor opt4" 2
fi

START=$(/bin/date +%s)
"$TOPAS_BIN" run.txt > topas.log 2>&1
TOPAS_EXIT=$?
END=$(/bin/date +%s)

VERDICT=COMPLETE; STATUS=0
{
    echo "topas_exit: $TOPAS_EXIT"
    echo "wall_seconds: $((END - START))"
    echo "end: $(/bin/date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "expected_bin_bytes: $EXPECT_BYTES"
} >> provenance.txt

if (( TOPAS_EXIT != 0 )); then
    VERDICT="FAILED topas exited $TOPAS_EXIT"; STATUS=7
fi

for f in dose.bin dose_all.bin dose.binheader dose_all.binheader; do
    if [ ! -s "$f" ]; then
        echo "$f: MISSING_OR_EMPTY" >> provenance.txt
        (( STATUS == 0 )) && { VERDICT="INCOMPLETE $f missing or empty"; STATUS=5; }
        continue
    fi
    SIZE=$(/usr/bin/stat -f %z "$f")
    echo "$f bytes: $SIZE sha256: $(sha "$f")" >> provenance.txt
    if [[ "$f" == *.bin ]] && (( SIZE != EXPECT_BYTES )); then
        (( STATUS == 0 )) && { VERDICT="INCOMPLETE $f is $SIZE bytes, expected $EXPECT_BYTES"; STATUS=5; }
    fi
done

EM_LINE=$(/usr/bin/grep -m1 'Use ICRU90 data' topas.log)
EM_VAL=$(echo "$EM_LINE" | /usr/bin/awk '{print $NF}')
echo "em_icru90_logged: ${EM_VAL:-NOT_FOUND} expected: $EXPECT_ICRU90" >> provenance.txt
if [ "$EM_VAL" != "$EXPECT_ICRU90" ]; then
    (( STATUS == 0 )) && { VERDICT="UNVERIFIED EM physics: logged '${EM_VAL:-NOT_FOUND}', expected $EXPECT_ICRU90"; STATUS=6; }
fi

echo "verdict: $VERDICT" >> provenance.txt
(( STATUS == 0 )) || echo "run_topas.sh: $DIR: $VERDICT" >&2
exit $STATUS
