#!/usr/bin/env bats
# Contract tests for make_run.sh and run_topas.sh (issue #32). No transport: TOPAS_BIN is a fake, and the
# scoring grid is shrunk to 2x3x4 in a private copy of the base so outputs are a few hundred bytes.

setup() {
    SRC="$(cd "$BATS_TEST_DIRNAME/.." && pwd)"
    W="$BATS_TEST_TMPDIR/w"
    mkdir -p "$W/bin"
    cp "$SRC/make_run.sh" "$SRC/run_topas.sh" "$W/bin/"
    sed -e 's/^i:Ge\/Phantom\/XBins = .*/i:Ge\/Phantom\/XBins = 2/' \
        -e 's/^i:Ge\/Phantom\/YBins = .*/i:Ge\/Phantom\/YBins = 3/' \
        -e 's/^i:Ge\/Phantom\/ZBins = .*/i:Ge\/Phantom\/ZBins = 4/' \
        "$SRC/stage1_base.txt" > "$W/bin/stage1_base.txt"
    grep -q '^i:Ge/Phantom/ZBins = 4$' "$W/bin/stage1_base.txt"   # the shrink must have applied
    export TOPAS_BIN="$BATS_TEST_DIRNAME/fake_topas.sh"
    export MIN_FREE_GB=0
    OUT="$W/out"
}

mk() { local a=("$@"); "$W/bin/make_run.sh" "${a[@]:0:5}" "$OUT" "${a[@]:5}"; }
rundir() { dirname "$(mk "$@")"; }
prov() { grep "^$2" "$1/provenance.txt"; }

# ---------- make_run.sh
@test "make_run: a valid configuration writes run.txt, a frozen base, and threads in the name" {
    run mk 150 opt0 1 100 2
    [ "$status" -eq 0 ]
    [[ "$output" == */E150_opt0_seed1_n100_th2_a1/run.txt ]]
    grep -q '"g4em-standard_opt0"' "$output"
    cmp -s "$W/bin/stage1_base.txt" "$(dirname "$output")/stage1_base.txt"
}

@test "make_run: zero or out-of-range energy, histories, threads, seed are refused and write nothing" {
    for args in "0 opt0 1 100 2" "9 opt0 1 100 2" "401 opt0 1 100 2" "150 opt0 1 0 2" "150 opt0 1 100 0" \
                "150 opt0 1 100 33" "150 opt0 0 100 2" "150 opt0 2147483648 100 2" "150 opt0 1 -5 2" \
                "150 opt0 1 007 2" "150 opt0 1 1e6 2" "150 opt2 1 100 2"; do
        # shellcheck disable=SC2086
        run mk $args
        [ "$status" -eq 2 ] || { echo "accepted: $args"; false; }
    done
    [ ! -d "$OUT" ] || [ -z "$(ls -A "$OUT")" ]
}

@test "make_run: the same configuration twice is refused; a new attempt tag is allowed" {
    run mk 150 opt0 1 100 2;      [ "$status" -eq 0 ]
    run mk 150 opt0 1 100 2;      [ "$status" -eq 3 ]
    run mk 150 opt0 1 100 2 a2;   [ "$status" -eq 0 ]
    run mk 150 opt0 1 100 2 b1;   [ "$status" -eq 2 ]
}

@test "make_run: two concurrent claims of one configuration, exactly one wins" {
    ( rc=0; mk 150 opt0 7 100 2 > "$W/c1" 2>&1 || rc=$?; echo $rc > "$W/r1" ) &
    ( rc=0; mk 150 opt0 7 100 2 > "$W/c2" 2>&1 || rc=$?; echo $rc > "$W/r2" ) &
    wait
    r="$(cat "$W/r1") $(cat "$W/r2")"
    [ "$r" = "0 3" ] || [ "$r" = "3 0" ]
}

# ---------- run_topas.sh
@test "run_topas: complete output and matching EM physics -> 0, verdict COMPLETE" {
    d=$(rundir 150 opt0 1 100 2)
    run "$W/bin/run_topas.sh" "$d"
    [ "$status" -eq 0 ]
    prov "$d" "verdict: COMPLETE"
    prov "$d" "topas_exit: 0"
    prov "$d" "em_icru90_logged: 0 expected: 0"
    prov "$d" "topas_bin sha256: "
}

@test "run_topas: opt4 arm expects ICRU90 = 1" {
    d=$(rundir 150 opt4 1 100 2)
    run "$W/bin/run_topas.sh" "$d"
    [ "$status" -eq 0 ]
    prov "$d" "em_icru90_logged: 1 expected: 1"
}

@test "run_topas: TOPAS nonzero exit -> 7, its status preserved" {
    d=$(rundir 150 opt0 1 100 2)
    FAKE_MODE=exit3 run "$W/bin/run_topas.sh" "$d"
    [ "$status" -eq 7 ]
    prov "$d" "topas_exit: 3"
    prov "$d" "verdict: FAILED"
}

@test "run_topas: exit 0 with a missing, empty or wrong-size output -> 5" {
    for m in missing empty wrongsize; do
        d=$(rundir 150 opt0 1 100 2 "a$(( ${#m} ))")
        FAKE_MODE=$m run "$W/bin/run_topas.sh" "$d"
        [ "$status" -eq 5 ] || { echo "mode $m gave $status"; false; }
        prov "$d" "topas_exit: 0"
        prov "$d" "verdict: INCOMPLETE"
    done
}

@test "run_topas: exit 0 with EM line missing or mismatched -> 6" {
    for m in noem wrongem; do
        d=$(rundir 150 opt0 1 100 2 "a$(( ${#m} ))")
        FAKE_MODE=$m run "$W/bin/run_topas.sh" "$d"
        [ "$status" -eq 6 ] || { echo "mode $m gave $status"; false; }
        prov "$d" "verdict: UNVERIFIED"
    done
}

@test "run_topas: a repeat invocation is refused and leaves the first attempt's records byte-identical" {
    d=$(rundir 150 opt0 1 100 2)
    run "$W/bin/run_topas.sh" "$d";  [ "$status" -eq 0 ]
    before=$(shasum -a 256 "$d"/provenance.txt "$d"/topas.log "$d"/dose.bin)
    FAKE_MODE=exit3 run "$W/bin/run_topas.sh" "$d"
    [ "$status" -eq 3 ]
    [ "$(shasum -a 256 "$d"/provenance.txt "$d"/topas.log "$d"/dose.bin)" = "$before" ]
}

@test "run_topas: a directory holding stale outputs but no claim is refused without writing" {
    d=$(rundir 150 opt0 1 100 2)
    echo old > "$d/topas.log"
    run "$W/bin/run_topas.sh" "$d"
    [ "$status" -eq 3 ]
    [ "$(cat "$d/topas.log")" = old ]
    [ ! -e "$d/provenance.txt" ]
}

@test "run_topas: two concurrent invocations on one directory, exactly one runs" {
    d=$(rundir 150 opt0 1 100 2)
    ( rc=0; FAKE_MODE=slow "$W/bin/run_topas.sh" "$d" > /dev/null 2>&1 || rc=$?; echo $rc > "$W/q1" ) &
    ( rc=0; FAKE_MODE=slow "$W/bin/run_topas.sh" "$d" > /dev/null 2>&1 || rc=$?; echo $rc > "$W/q2" ) &
    wait
    r="$(cat "$W/q1") $(cat "$W/q2")"
    [ "$r" = "0 3" ] || [ "$r" = "3 0" ]
    prov "$d" "verdict: COMPLETE"
}

@test "run_topas: malformed or excessive MIN_FREE_GB fails closed without running TOPAS" {
    i=10
    for v in abc -1 "" 1e3 10000000 008 010 00 " 5"; do
        i=$(( i + 1 ))
        d=$(rundir 150 opt0 1 100 2 "a$i")
        MIN_FREE_GB=$v run "$W/bin/run_topas.sh" "$d"
        [ "$status" -eq 4 ] || { echo "MIN_FREE_GB='$v' gave $status"; false; }
        [ ! -e "$d/topas.log" ]
        prov "$d" "verdict: REFUSED"
    done
}

@test "run_topas: not a run directory -> 2" {
    run "$W/bin/run_topas.sh" "$W"
    [ "$status" -eq 2 ]
}

@test "run_topas: MIN_FREE_GB=0 is a valid decimal and runs" {
    d=$(rundir 150 opt0 1 100 2)
    MIN_FREE_GB=0 run "$W/bin/run_topas.sh" "$d"
    [ "$status" -eq 0 ]
}

@test "run_topas: relative invocation records the runner's real digest" {
    d=$(rundir 150 opt0 1 100 2)
    expected=$(shasum -a 256 "$W/bin/run_topas.sh" | cut -d' ' -f1)
    [[ "$expected" =~ ^[0-9a-f]{64}$ ]]
    cd "$W"
    run bin/run_topas.sh "$d"
    [ "$status" -eq 0 ]
    prov "$d" "runner sha256: $expected\$"
    prov "$d" "runner: $W/bin/run_topas.sh\$"
}

@test "run_topas: every recorded digest is a full 64-hex hash" {
    d=$(rundir 150 opt0 1 100 2)
    run "$W/bin/run_topas.sh" "$d"
    [ "$status" -eq 0 ]
    n=$(grep -c 'sha256: ' "$d/provenance.txt")
    [ "$n" -eq 8 ]   # runner, topas_bin, run.txt, base, and 4 outputs
    [ "$(grep -cE 'sha256: [0-9a-f]{64}$' "$d/provenance.txt")" -eq 8 ]
}
