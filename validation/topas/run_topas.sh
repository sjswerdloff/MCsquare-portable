#!/bin/bash
# Run one stage-1 TOPAS run directory (made by make_run.sh) on Mac-Studio and record provenance.
#
#   run_topas.sh <run_dir>
#
# Launch heavy runs detached:  nohup run_topas.sh <dir> > /dev/null 2>&1 < /dev/null &
# The environment lives in this script because macOS strips DYLD_* from anything launched through a
# SIP-protected binary (/usr/bin/time, env -i): topas then dies with "Library not loaded" (exit 134).
set -uo pipefail

DIR=${1:?usage: run_topas.sh <run_dir>}
MC=/Users/stuartswerdloff/MonteCarlo
export TOPAS_G4_DATA_DIR=$MC/geant4-data-11.3.2
export DYLD_LIBRARY_PATH=$MC/opentopas-4.3.0/lib:$MC/geant4-11.3.2/lib:$MC/gdcm-2.6.8/lib
cd "$DIR" || exit 2

# Free-space guard. Two 400x400x350 double scorers write ~0.86 GB per run, and the obvious output volume
# (/Volumes/T7 Shield) also holds mlxcel-cold-storage, the family cold store. Refuse to start rather than
# let a batch of runs fill a disk someone's memory lives on. Override with MIN_FREE_GB.
MIN_FREE_GB=${MIN_FREE_GB:-100}
FREE_GB=$(/bin/df -g . | /usr/bin/awk 'NR==2 {print $4}')
if ! [[ "$FREE_GB" =~ ^[0-9]+$ ]]; then
    echo "refusing: could not read free space for $DIR (got '$FREE_GB')" >&2
    exit 4
fi
if [ "$FREE_GB" -lt "$MIN_FREE_GB" ]; then
    echo "refusing: ${FREE_GB} GB free on the volume holding $DIR, below MIN_FREE_GB=${MIN_FREE_GB}" >&2
    exit 4
fi

{
    echo "host: $(/bin/hostname)"
    echo "start: $(/bin/date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "topas: $MC/opentopas-4.3.0/bin/topas (OpenTOPAS 4.3.0, commit 134451a)"
    echo "geant4: 11.3.2, data $TOPAS_G4_DATA_DIR"
    echo "run.txt sha256: $(/usr/bin/shasum -a 256 run.txt | cut -d' ' -f1)"
    echo "stage1_base.txt sha256: $(/usr/bin/shasum -a 256 stage1_base.txt | cut -d' ' -f1)"
} > provenance.txt

START=$(/bin/date +%s)
"$MC/opentopas-4.3.0/bin/topas" run.txt > topas.log 2>&1
RC=$?
END=$(/bin/date +%s)
{
    echo "exit: $RC"
    echo "wall_seconds: $((END - START))"
    echo "end: $(/bin/date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "em_icru90: $(/usr/bin/grep -m1 'Use ICRU90 data' topas.log || echo 'NOT FOUND - EM physics unverified')"
    for f in dose.bin dose.binheader dose_all.bin dose_all.binheader; do
        if [ -f "$f" ]; then echo "$f sha256: $(/usr/bin/shasum -a 256 "$f" | cut -d' ' -f1)"; else echo "$f: MISSING"; fi
    done
} >> provenance.txt
exit $RC
