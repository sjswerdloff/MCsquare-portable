#!/bin/bash
# Write ONE stage-1 TOPAS run file (MCsquare-portable #32) into its own directory and print the path.
#
#   make_run.sh <energy_MeV> <em: opt0|opt4> <seed> <histories> <threads> <out_root>
#
# The run directory name carries every per-run choice, so an output file can never be separated from
# the parameters that produced it. Refuses to reuse a directory: a pre-registered run is written once.
set -euo pipefail

if [ $# -ne 6 ]; then
    echo "usage: $0 <energy_MeV> <opt0|opt4> <seed> <histories> <threads> <out_root>" >&2
    exit 2
fi
E=$1; EM=$2; SEED=$3; N=$4; THREADS=$5; ROOT=$6
HERE="$(cd "$(dirname "$0")" && pwd)"

case "$EM" in
    opt0) EMMOD="g4em-standard_opt0" ;;   # primary arm: matches MCsquare's stopping-power table
    opt4) EMMOD="g4em-standard_opt4" ;;   # secondary arm: TOPAS default / clinical list
    *) echo "em must be opt0 or opt4, got '$EM'" >&2; exit 2 ;;
esac
for v in "$E" "$SEED" "$N" "$THREADS"; do
    [[ "$v" =~ ^[0-9]+$ ]] || { echo "not a positive integer: '$v'" >&2; exit 2; }
done

DIR="$ROOT/E${E}_${EM}_seed${SEED}_n${N}"
if [ -e "$DIR" ]; then
    echo "refusing: $DIR already exists" >&2
    exit 3
fi
mkdir -p "$DIR"
cp "$HERE/stage1_base.txt" "$DIR/stage1_base.txt"   # frozen copy: the run records the base it used

cat > "$DIR/run.txt" <<EOF
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
