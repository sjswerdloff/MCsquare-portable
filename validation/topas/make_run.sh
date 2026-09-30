#!/bin/bash
# Write ONE stage-1 TOPAS run directory (issue #32) and print the path of its run.txt.
#
#   make_run.sh <energy_MeV> <em: opt0|opt4> <seed> <histories> <threads> <out_root> [attempt]
#
#   energy_MeV   integer, 10..400    (the MCsquare/Geant4 stopping-power tables stop at 400 MeV)
#   seed         integer, 1..2147483647   (TOPAS i:Ts/Seed is a signed 32-bit int)
#   histories    integer, 1..2147483647   (i:So/Beam/NumberOfHistoriesInRun, same type)
#   threads      integer, 1..32      (Mac Studio M3 Ultra: 24 performance + 8 efficiency cores)
#   attempt      a1, a2, ...  default a1. A retry of the same configuration gets a new attempt tag;
#                a directory is never reused, so a failed attempt's records are kept.
#
# The directory name carries every per-run choice, so an output can never be separated from the parameters
# that produced it. The leaf directory is claimed with a plain mkdir (no -p), which is atomic: of two
# concurrent calls for the same configuration exactly one succeeds, and the other refuses before writing.
set -uo pipefail

die() { echo "make_run.sh: $1" >&2; exit "${2:-2}"; }

[ $# -eq 6 ] || [ $# -eq 7 ] || die "usage: $0 <energy_MeV> <opt0|opt4> <seed> <histories> <threads> <out_root> [attempt]"
E=$1; EM=$2; SEED=$3; N=$4; THREADS=$5; ROOT=$6; ATTEMPT=${7:-a1}
HERE="$(cd "$(dirname "$0")" && pwd)"

# int_in_range <name> <value> <min> <max>: decimal integer, no sign, no leading zeros, within [min, max].
int_in_range() {
    [[ "$2" =~ ^[1-9][0-9]{0,9}$ ]] || die "$1 must be a positive decimal integer, got '$2'"
    (( 10#$2 >= $3 && 10#$2 <= $4 )) || die "$1 must be in $3..$4, got $2"
}
int_in_range energy_MeV "$E" 10 400
int_in_range seed "$SEED" 1 2147483647
int_in_range histories "$N" 1 2147483647
int_in_range threads "$THREADS" 1 32
[[ "$ATTEMPT" =~ ^a[1-9][0-9]{0,2}$ ]] || die "attempt must be a1..a999, got '$ATTEMPT'"

case "$EM" in
    opt0) EMMOD="g4em-standard_opt0" ;;   # primary arm: matches MCsquare's stopping-power table
    opt4) EMMOD="g4em-standard_opt4" ;;   # secondary arm: TOPAS default / clinical list
    *) die "em must be opt0 or opt4, got '$EM'" ;;
esac

[ -f "$HERE/stage1_base.txt" ] || die "missing $HERE/stage1_base.txt"
mkdir -p "$ROOT" || die "cannot create output root $ROOT"
DIR="$ROOT/E${E}_${EM}_seed${SEED}_n${N}_th${THREADS}_${ATTEMPT}"
mkdir "$DIR" 2>/dev/null || die "refusing: $DIR already exists (or cannot be created)" 3

cp "$HERE/stage1_base.txt" "$DIR/stage1_base.txt" || die "cannot copy base into $DIR"
cat > "$DIR/run.txt" <<EOF || die "cannot write $DIR/run.txt"
includeFile = stage1_base.txt
d:So/Beam/BeamEnergy = ${E} MeV
sv:Ph/Default/Modules = 6 "${EMMOD}" "g4h-phy_QGSP_BIC_HP" "g4decay" "g4ion-binarycascade" "g4h-elastic_HP" "g4stopping"
i:Ts/Seed = ${SEED}
i:So/Beam/NumberOfHistoriesInRun = ${N}
i:Ts/NumberOfThreads = ${THREADS}
s:Sc/Dose/OutputFile = "dose"
s:Sc/DoseAll/OutputFile = "dose_all"
EOF
echo "$DIR/run.txt"
