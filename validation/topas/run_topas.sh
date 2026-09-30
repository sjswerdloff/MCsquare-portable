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
#   8  PROVENANCE FAILURE: a mandatory hash or provenance write failed (nothing is reported as COMPLETE
#      without every hash present)
# A single run is attempted once. A retry is a NEW directory (make_run.sh ... a2), so the records of a failed
# attempt are kept.
#
# TOPAS_BIN overrides the executable (tests use a fake). Whatever runs, its path and sha256 are recorded.
set -uo pipefail

die() { echo "run_topas.sh: $1" >&2; exit "$2"; }

[ $# -eq 1 ] || die "usage: $0 <run_dir>" 2
DIR=$1
[ -d "$DIR" ] && [ -f "$DIR/run.txt" ] && [ -f "$DIR/stage1_base.txt" ] || die "not a run directory: $DIR" 2

# Resolve THIS script's real path BEFORE changing directory: a relative "$0" means something else after cd.
RUNNER="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
[ -f "$RUNNER" ] || die "cannot resolve the runner's own path from '$0'" 2

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

# sha <file>: prints a 64-hex digest or NOTHING. Callers capture it into a variable and check it, because a
# failure inside $( ) cannot stop this script and would otherwise be recorded as an empty hash.
sha() { local h; h=$(/usr/bin/shasum -a 256 "$1" 2>/dev/null | /usr/bin/cut -d' ' -f1); [[ "$h" =~ ^[0-9a-f]{64}$ ]] && echo "$h"; }
need_sha() { local h; h=$(sha "$2"); [ -n "$h" ] || die "provenance failure: cannot hash $1 ($2)" 8; echo "$h"; }
# append <line>: every provenance write is checked.
append() { echo "$1" >> provenance.txt || die "provenance failure: cannot append to $DIR/provenance.txt" 8; }

# Mandatory input hashes, computed and checked BEFORE anything is written.
[ -f "$TOPAS_BIN" ] || die "topas executable not found: $TOPAS_BIN" 2
H_RUNNER=$(need_sha runner "$RUNNER") || exit 8
H_TOPAS=$(need_sha topas_bin "$TOPAS_BIN") || exit 8
H_RUNTXT=$(need_sha run.txt run.txt) || exit 8
H_BASE=$(need_sha stage1_base.txt stage1_base.txt) || exit 8

# noclobber: provenance.txt is created here or not at all.
set -o noclobber
{
    echo "host: $(/bin/hostname)"
    echo "start: $(/bin/date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "runner: $RUNNER"
    echo "runner sha256: $H_RUNNER"
    echo "topas_bin: $TOPAS_BIN"
    echo "topas_bin sha256: $H_TOPAS"
    echo "geant4_data: $TOPAS_G4_DATA_DIR"
    echo "run.txt sha256: $H_RUNTXT"
    echo "stage1_base.txt sha256: $H_BASE"
} > provenance.txt || die "cannot create provenance.txt" 3
set +o noclobber

# --- Free-space guard, fails closed.
MIN_FREE_GB=${MIN_FREE_GB-100}   # default only when UNSET: an explicit empty value is malformed and refused
# Decimal only: no sign, no leading zeros (bash reads 010 as octal 8, and 008 as an error), at most 6 digits.
if ! [[ "$MIN_FREE_GB" =~ ^(0|[1-9][0-9]{0,5})$ ]]; then
    append "verdict: REFUSED malformed MIN_FREE_GB='$MIN_FREE_GB'"
    die "refusing: MIN_FREE_GB must be a whole number of GB, got '$MIN_FREE_GB'" 4
fi
FREE_GB=$(/bin/df -g . | /usr/bin/awk 'NR==2 {print $4}')
if ! [[ "$FREE_GB" =~ ^[0-9]+$ ]]; then
    append "verdict: REFUSED unreadable free space '$FREE_GB'"
    die "refusing: could not read free space for $DIR (got '$FREE_GB')" 4
fi
if (( 10#$FREE_GB < 10#$MIN_FREE_GB )); then
    append "verdict: REFUSED ${FREE_GB} GB free < MIN_FREE_GB=${MIN_FREE_GB}"
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
} >> provenance.txt || die "provenance failure: cannot append to $DIR/provenance.txt" 8

if (( TOPAS_EXIT != 0 )); then
    VERDICT="FAILED topas exited $TOPAS_EXIT"; STATUS=7
fi

for f in dose.bin dose_all.bin dose.binheader dose_all.binheader; do
    if [ ! -s "$f" ]; then
        append "$f: MISSING_OR_EMPTY"
        (( STATUS == 0 )) && { VERDICT="INCOMPLETE $f missing or empty"; STATUS=5; }
        continue
    fi
    SIZE=$(/usr/bin/stat -f %z "$f")
    if ! [[ "$SIZE" =~ ^[0-9]+$ ]]; then
        append "$f: SIZE_UNREADABLE '$SIZE'"; VERDICT="PROVENANCE FAILURE cannot read the size of $f"; STATUS=8; continue
    fi
    H_OUT=$(sha "$f")
    append "$f bytes: $SIZE sha256: ${H_OUT:-HASH_FAILED}"
    if [ -z "$H_OUT" ]; then VERDICT="PROVENANCE FAILURE cannot hash $f"; STATUS=8; fi
    if [[ "$f" == *.bin ]] && (( SIZE != EXPECT_BYTES )); then
        (( STATUS == 0 )) && { VERDICT="INCOMPLETE $f is $SIZE bytes, expected $EXPECT_BYTES"; STATUS=5; }
    fi
done

EM_LINE=$(/usr/bin/grep -m1 'Use ICRU90 data' topas.log)
EM_VAL=$(echo "$EM_LINE" | /usr/bin/awk '{print $NF}')
append "em_icru90_logged: ${EM_VAL:-NOT_FOUND} expected: $EXPECT_ICRU90"
if [ "$EM_VAL" != "$EXPECT_ICRU90" ]; then
    (( STATUS == 0 )) && { VERDICT="UNVERIFIED EM physics: logged '${EM_VAL:-NOT_FOUND}', expected $EXPECT_ICRU90"; STATUS=6; }
fi

append "verdict: $VERDICT"
(( STATUS == 0 )) || echo "run_topas.sh: $DIR: $VERDICT" >&2
exit $STATUS
