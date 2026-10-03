"""Tests for field_followup_analyse.py (part E) on SYNTHETIC inputs only: no run output from any machine is read.

The run trees are written the way the two workflows write them: the Lenovo's files with CRLF line ends, compact JSON
and a UTF-16 log; the HP's with LF, spaced JSON and a UTF-8 log. Doses are full size (150 voxels per side), because
the endpoints need the whole grid: the 80-run tree is about 1.1 GB and is built once per test session.

Run: uv run --no-project --with numpy --with scipy --with pytest pytest validation/tests/test_field_followup_analyse.py
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import random
import shutil
import statistics as st
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import t as student_t

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import apples_analyse as aa
import field_edge_make_cases as fm
import field_followup_analyse as fa
import field_followup_runs as fr
import platform_study_analyse as psa
import platform_study_metrics as psm
from field_edge_analyse import SP, C, N

COMMIT = "ab" * 20
WINDOWS = {"A-up", "A-port"}  # the arms whose files the Lenovo writes
PEAK_ROW = {100: 111, 150: 71}  # row of the synthetic Bragg peak; both slab depths of each energy lie short of it
REPO = fa.REPO


# ------------------------------------------------------------------------------------------------ synthetic trees
def working_tree_blobs(_commit: str, _repo: Path) -> dict[str, str]:
    """Stand-in for `git ls-tree` at the acquisition commit: the blob ids of the working tree's files."""
    files = [p for top in fa.INPUT_PATHS for p in ([REPO / top] if (REPO / top).is_file() else sorted((REPO / top).rglob("*")))
             if p.is_file()] + [REPO / rel for rel in (*fa.RUN_LIST_FILES, *fa.ANALYSIS_FILES) if (REPO / rel).is_file()]
    return {p.relative_to(REPO).as_posix(): fa.blob_id(p.read_bytes()) for p in files}


@pytest.fixture(autouse=True)
def _blobs_from_the_working_tree(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fa, "committed_blobs", working_tree_blobs)


def dose(energy: int, rng: random.Random, *, far_zero: bool = False, tail: float = 1.0, scale: float = 1.0) -> np.ndarray:
    """d[z, y, x]: a 150 mm field with an exponential tail outside, times a depth curve with a peak at PEAK_ROW."""
    c = np.abs((np.arange(N) - C + 0.5) * SP)
    lat = np.where(c < 75, 1.0, tail * (1 + rng.gauss(0, 0.002)) * 0.5 * np.exp(-(c - 75) / 10.0))
    if far_zero:
        lat[c > 140] = 0.0
    iy = np.arange(N)
    depth = scale * (1 + rng.gauss(0, 0.001)) * (1.0 + 3.0 * np.exp(-(((iy - PEAK_ROW[energy] - rng.gauss(0, 0.02)) / 4.0) ** 2)))
    depth[:PEAK_ROW[energy] - 8] = 0.0  # beyond the distal edge (smaller iy is deeper)
    return (lat[:, None, None] * depth[None, :, None] * lat[None, None, :]).astype("<f4")


class Template:
    """The files every run directory shares, made once: the generated case and the committed inputs."""

    def __init__(self, home: Path) -> None:
        self.home = home
        self.case = home / "case"
        self.case.mkdir(parents=True)
        cwd = Path.cwd()
        os.chdir(self.case)
        try:
            fm.main([float(e) for e in fr.ENERGIES], [fr.SIDE_MM])
        finally:
            os.chdir(cwd)
        for top in ("Materials", "Scanners/default"):
            shutil.copytree(REPO / top, home / "inputs" / top)
        (home / "inputs" / "BDL").mkdir()
        shutil.copy2(REPO / fa.INPUT_PATHS[2], home / "inputs" / fa.INPUT_PATHS[2])

    def link(self, src: Path, dst: Path) -> None:
        dst.parent.mkdir(parents=True, exist_ok=True)
        try:
            os.link(src, dst)
        except OSError:
            shutil.copy2(src, dst)


def text(arm: str, lines: list[str]) -> bytes:
    return ("\r\n" if arm in WINDOWS else "\n").join([*lines, ""]).encode("ascii")


def write_run(template: Template, cell: Path, arm: str, energy: int, seed: int, mode: str, d: np.ndarray,
              override: dict[str, object] | None = None) -> Path:
    """One run directory as the arm's workflow writes it. `override` replaces run.json fields."""
    run = cell / f"s{seed}"
    (run / "out_seed").mkdir(parents=True)
    win = arm in WINDOWS
    requested = fr.HISTORIES[mode]
    overshoot = 0 if fr.ARMS[arm].portable else 3
    rec = {"seed": seed, "requested": requested, "simulated": requested + overshoot, "overshoot": overshoot,
           "threads": fr.ARMS[arm].threads, "energy_mev": energy, "case": "E", "arm": arm, "mode": mode,
           "host": fr.ARMS[arm].host, "cpu": "synthetic", "binary_sha256": fr.ARMS[arm].binary_sha256, "commit": COMMIT,
           "start_utc": "2026-10-03T00:00:00Z", "end_utc": "2026-10-03T00:06:00Z", "wall_s": 360, "transport_status": "ok", **(override or {})}
    (run / "run.json").write_bytes(text(arm, [json.dumps(rec, separators=(",", ":")) if win else json.dumps(rec)]))
    (run / "cfg.txt").write_bytes(text(arm, fa.expected_cfg(arm, energy, seed, mode)))
    log = f"\nRange shifter initialized for beam 0:\n\nNbr primaries simulated: {rec['simulated']} \n"
    (run / "log.txt").write_bytes(log.replace("\n", "\r\n").encode("utf-16") if win else log.encode())
    for name in ("cube.mhd", f"{fr.plan_name(energy)}.txt"):
        data = (template.case / name).read_bytes()
        (run / name).write_bytes(data.replace(b"\n", b"\r\n") if win else data)
    template.link(template.case / "cube.raw", run / "cube.raw")
    for src in (template.home / "inputs").rglob("*"):
        if src.is_file():
            template.link(src, run / src.relative_to(template.home / "inputs"))
    (run / "out_seed" / "Dose.raw").write_bytes(d.tobytes())
    (run / "out_seed" / "Dose.mhd").write_bytes(text(arm, [
        "ObjectType = Image", "NDims = 3", f"DimSize = {N} {N} {N}", "ElementSpacing = 2.000000 2.000000 2.000000",
        "Offset = 0.000000 0.000000 0.000000", "ElementType = MET_FLOAT", "ElementByteOrderMSB = False", "ElementDataFile = Dose.raw"]))
    record_hashes(run, arm)
    return run


def record_hashes(run: Path, arm: str) -> None:
    out = run / "out_seed"
    (out / "sha256.txt").write_bytes(text(arm, [f"{hashlib.sha256((out / n).read_bytes()).hexdigest()}  {n}" for n in ("Dose.raw", "Dose.mhd")]))


def write_cell(template: Template, root: Path, arm: str, energy: int, mode: str, **dose_args: object) -> Path:
    cell = root / fr.ARMS[arm].directory / f"e{energy}"
    for rel in fa.RUN_LIST_FILES:
        (cell / "snapshot" / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(REPO / rel, cell / "snapshot" / rel)
    for seed in fr.seeds(arm, energy, mode):
        write_run(template, cell, arm, energy, seed, mode, dose(energy, random.Random(seed), **dose_args))  # type: ignore[arg-type]
    return cell


@pytest.fixture(scope="session")
def template(tmp_path_factory: pytest.TempPathFactory) -> Template:
    return Template(tmp_path_factory.mktemp("field_e_template"))


@pytest.fixture(scope="session")
def trees(template: Template, tmp_path_factory: pytest.TempPathFactory) -> list[Path]:
    """The complete population as it reaches the share: one root per host. At 100 MeV nothing is scored 70 mm out."""
    home = tmp_path_factory.mktemp("field_e_full")
    roots = [home / "lenovo", home / "hp"]
    for arm in fr.ARMS:
        for energy in fr.ENERGIES:
            write_cell(template, roots[0 if arm in WINDOWS else 1], arm, energy, "full", far_zero=energy == 100)
    return roots


@pytest.fixture(scope="session")
def document(trees: list[Path]) -> dict:
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(fa, "committed_blobs", working_tree_blobs)
        return fa.run(trees, COMMIT)


@pytest.fixture
def one(template: Template, tmp_path: Path, request: pytest.FixtureRequest) -> tuple[fa.Run, fa.Expected, aa.Ledger]:
    """One run directory (zero dose) of the arm in request.param, with what it is held to."""
    arm = getattr(request, "param", "A-port")
    seed = fr.seeds(arm, 100, "full")[0]
    path = write_run(template, tmp_path / fr.ARMS[arm].directory / "e100", arm, 100, seed, "full", np.zeros((N, N, N), dtype="<f4"))
    return fa.Run(arm, 100, seed, path), fa.expected(COMMIT, "full"), aa.Ledger(tmp_path)


def problems(one: tuple[fa.Run, fa.Expected, aa.Ledger]) -> str:
    run, want, ledger = one
    fa.verify_run(run, want, ledger)
    return " | ".join(run.problems)


def rewrite(path: Path, old: bytes, new: bytes) -> None:
    """Replace the file (never in place: input files are hard links shared with other runs)."""
    data = path.read_bytes()
    assert data.count(old) == 1, (path.name, data.count(old))
    path.unlink()
    path.write_bytes(data.replace(old, new))


# ------------------------------------------------------------------------------------------------- frozen design
def test_the_family_is_the_plans_34_endpoints() -> None:
    specs = fa.specs()
    assert len(specs) == 34 and [s.energy for s in specs] == [100] * 17 + [150] * 17
    assert [s.key for s in specs[:17]] == [
        *(f"lateral_{d}_{o}" for d in (39, 61) for o in (5, 10, 20, 30, 50, 70)), "cax_39", "cax_61", "cax_max", "r80_mm", "r20_mm"]
    assert [s.key for s in specs[17:29]] == [f"lateral_{d}_{o}" for d in (79, 125) for o in (5, 10, 20, 30, 50, 70)]
    assert (fa.RATIO_OFFSETS_MM, fa.ALPHA, fa.GRID_N) == ((50, 70), 0.05, 150)


def test_margins_are_those_of_the_200_mev_field() -> None:
    for new, old in (("lateral_39_5", "lateral_127_5"), ("lateral_61_10", "lateral_127_10"), ("lateral_125_30", "lateral_201_30"),
                     ("cax_39", "cax_127"), ("cax_max", "cax_i23"), ("r80_mm", "r80_mm"), ("r20_mm", "r20_mm")):
        assert fa.margin(new)[:2] == psa.margin(old)
    assert fa.margin("lateral_39_50") == fa.margin("lateral_125_70") == (-0.20, 0.20, False)  # the two new offsets
    assert fa.margin("lateral_125_5") == (-0.50, 0.50, False)
    assert fa.margin("cax_79") == (math.log(0.995), math.log(1.005), True)
    with pytest.raises(ValueError, match="no margin"):
        fa.margin("sigma_40")


# ---------------------------------------------------------------------------------------- verification of one run
@pytest.mark.parametrize("one", ["A-port", "B-pgcc"], indirect=True)
def test_a_run_as_either_host_writes_it_verifies(one: tuple[fa.Run, fa.Expected, aa.Ledger]) -> None:
    assert problems(one) == ""
    run, _want, ledger = one
    assert run.dose_sha256 == hashlib.sha256((run.path / "out_seed/Dose.raw").read_bytes()).hexdigest()
    assert len(ledger.files) == 9 + 291  # every file read is in the fingerprint: nine of the run, 291 inputs


@pytest.mark.parametrize(("field", "value", "message"), [
    ("seed", 961042, "run.json seed is 961042, expected 961041"),
    ("requested", 6e7, "run.json requested is 60000000.0, expected 60000000"),
    ("requested", 30_000_000, "run.json requested is 30000000"),
    ("threads", 3, "run.json threads is 3, expected 4"),
    ("energy_mev", 150, "run.json energy_mev is 150, expected 100"),
    ("arm", "A-up", "run.json arm is 'A-up'"),
    ("mode", "smoke", "run.json mode is 'smoke', expected 'full'"),
    ("host", "DESKTOP-5H86O9N", "run.json host is 'DESKTOP-5H86O9N'"),
    ("binary_sha256", "f3a28398a399224c1e20d6e54e945571f6ee17693f4380ee7ad9c0e3ef84b966", "run.json binary_sha256 is"),
    ("commit", "cd" * 20, "run.json commit is"),
    ("transport_status", "dose_missing", "run.json transport_status is 'dose_missing', expected 'ok'"),
    ("simulated", 59_999_999, "run.json simulated is 59999999, expected an integer of at least 60000000"),
    ("simulated", None, "run.json simulated is None"),
    ("overshoot", 5, "run.json overshoot is 5, expected 0"),
    ("case", "F", "run.json case is 'F', expected 'E'"),
])
def test_run_json_must_be_the_frozen_run(template: Template, tmp_path: Path, field: str, value: object, message: str) -> None:
    path = write_run(template, tmp_path / "apt" / "e100", "A-port", 100, 961041, "full", np.zeros((N, N, N), dtype="<f4"), {field: value})
    run = fa.Run("A-port", 100, 961041, path)
    fa.verify_run(run, fa.expected(COMMIT, "full"), aa.Ledger(tmp_path))
    assert message in " | ".join(run.problems)


MUTATIONS = [
    ("run.json", b'{"seed"', b'["seed"', "run.json is not a JSON object"),
    ("run.json", b'"transport_status":"ok"', b'"transport_status":"failed","transport_status":"ok"', "run.json is not a JSON object with distinct keys"),
    ("cfg.txt", b"Num_Primaries 60000000", b"Num_Primaries 6e7", "cfg.txt is not the configuration"),
    ("cfg.txt", b"Dose_MHD_Output True", b"Dose_MHD_Output True\r\nSimulate_Nuclear_Interactions False", "cfg.txt is not the configuration"),
    ("cfg.txt", b"Scanners/default/HU_Density", b"Scanners/other/HU_Density", "cfg.txt is not the configuration"),
    ("E100_S150.txt", b"100.000000", b"101.000000", "E100_S150.txt is not what field_edge_make_cases.py writes"),
    ("cube.mhd", b"DimSize = 150 150 150", b"DimSize = 150 150 151", "cube.mhd is not what field_edge_make_cases.py writes"),
    ("BDL/BDL_default_DN_RangeShifter.txt", b"\n", b"\n\n", "1 input file(s) differ from the committed bytes: ['BDL/BDL_default_DN_RangeShifter.txt']"),
    ("Scanners/default/HU_Density_Conversion.txt", b"\n", b" \n", "1 input file(s) differ from the committed bytes"),
    ("out_seed/Dose.mhd", b"MET_FLOAT", b"MET_DOUBLE", "Dose.mhd is not the file hashed on the run host"),
]


@pytest.mark.parametrize(("name", "old", "new", "message"), MUTATIONS)
def test_a_changed_file_is_a_problem(one: tuple[fa.Run, fa.Expected, aa.Ledger], name: str, old: bytes, new: bytes, message: str) -> None:
    data = (one[0].path / name).read_bytes()
    first = data.find(old)
    assert first >= 0
    (one[0].path / name).unlink()
    (one[0].path / name).write_bytes(data[:first] + new + data[first + len(old):])
    assert message in problems(one)


def test_the_utf16_log_of_the_lenovo_is_read(one: tuple[fa.Run, fa.Expected, aa.Ledger]) -> None:
    log = one[0].path / "log.txt"
    assert log.read_bytes().startswith(b"\xff\xfe")
    rewrite(log, "60000000".encode("utf-16-le"), "60000001".encode("utf-16-le"))
    assert "log.txt reports 60000001 primaries simulated, run.json 60000000" in problems(one)


@pytest.mark.parametrize("one", ["B-up"], indirect=True)
def test_the_log_must_carry_the_simulated_count_and_no_unknown_tag(one: tuple[fa.Run, fa.Expected, aa.Ledger]) -> None:
    log = one[0].path / "log.txt"
    log.write_bytes(b"Unknown tag: Num_Primarys\nno count here\n")
    msg = problems(one)
    assert "log.txt has 0 'Nbr primaries simulated' lines, expected exactly 1" in msg
    assert "a configuration tag MCsquare did not recognise" in msg


@pytest.mark.parametrize("one", ["A-port", "B-up"], indirect=True)  # the Lenovo's UTF-16 CRLF log and the HP's LF log
@pytest.mark.parametrize(("line", "message"), [
    ("Nbr primaries simulated: {n}garbage ", "a malformed 'Nbr primaries simulated' line"),
    ("Nbr primaries simulated: {n} \nNbr primaries simulated: {m} ", "2 'Nbr primaries simulated' lines, expected exactly 1"),
    ("Nbr primaries simulated: {n} \nNbr primaries simulated: {n} ", "2 'Nbr primaries simulated' lines, expected exactly 1"),
    ("Nbr primaries simulated: {n} Nbr primaries simulated: {m} ", "a malformed 'Nbr primaries simulated' line"),
    ("xNbr primaries simulated: {n} ", "a malformed 'Nbr primaries simulated' line"),
    ("Nbr primaries simulated: {n} (7 generated outside the geometry) extra", "a malformed 'Nbr primaries simulated' line"),
    ("Nbr primaries simulated:  {n} ", "a malformed 'Nbr primaries simulated' line"),
])
def test_the_completion_count_must_be_one_whole_unambiguous_line(
        one: tuple[fa.Run, fa.Expected, aa.Ledger], line: str, message: str) -> None:
    """#59 review 7116: a numeric prefix ('60000000garbage') and the first of two conflicting counts both passed."""
    run = one[0]
    n = json.loads(fa._text((run.path / "run.json").read_bytes()))["simulated"]
    win = run.arm in WINDOWS
    log = "\nRange shifter initialized for beam 0:\n\n" + line.format(n=n, m=n + 1) + "\n"
    (run.path / "log.txt").unlink()
    (run.path / "log.txt").write_bytes(log.replace("\n", "\r\n").encode("utf-16") if win else log.encode())
    assert message in problems(one)


@pytest.mark.parametrize("one", ["A-port", "B-up"], indirect=True)
@pytest.mark.parametrize("line", ["Nbr primaries simulated: {n} ", "Nbr primaries simulated: {n} (19 generated outside the geometry) ",
                                  "Nbr primaries simulated: {n}"])
def test_both_forms_of_the_completion_line_that_mcsquare_prints_verify(one: tuple[fa.Run, fa.Expected, aa.Ledger], line: str) -> None:
    run = one[0]
    n = json.loads(fa._text((run.path / "run.json").read_bytes()))["simulated"]
    log = "\nRange shifter initialized for beam 0:\n\n" + line.format(n=n) + "\n"
    (run.path / "log.txt").unlink()
    (run.path / "log.txt").write_bytes(log.replace("\n", "\r\n").encode("utf-16") if run.arm in WINDOWS else log.encode())
    assert problems(one) == ""


def test_cube_raw_must_be_the_generated_ct(one: tuple[fa.Run, fa.Expected, aa.Ledger]) -> None:
    cube = one[0].path / "cube.raw"
    data = bytearray(cube.read_bytes())
    data[100] = 1
    cube.unlink()
    cube.write_bytes(bytes(data))
    assert "cube.raw is not what field_edge_make_cases.py writes" in problems(one)


def test_an_input_file_added_or_missing_is_a_problem(one: tuple[fa.Run, fa.Expected, aa.Ledger]) -> None:
    run = one[0].path
    (run / "Materials" / "extra.txt").write_text("x")
    next(p for p in sorted((run / "Materials").rglob("*")) if p.is_file() and p.name != "extra.txt").unlink()
    assert "inputs differ from the committed set: 1 added ['Materials/extra.txt'], 1 missing" in problems(one)


def test_a_missing_file_is_named(one: tuple[fa.Run, fa.Expected, aa.Ledger]) -> None:
    (one[0].path / "out_seed" / "Dose.raw").unlink()
    (one[0].path / "cfg.txt").unlink()
    assert problems(one) == "missing: cfg.txt, out_seed/Dose.raw"


def test_a_dose_changed_after_it_was_hashed_is_a_problem(one: tuple[fa.Run, fa.Expected, aa.Ledger]) -> None:
    raw = one[0].path / "out_seed" / "Dose.raw"
    data = bytearray(raw.read_bytes())
    data[0] = 1
    raw.write_bytes(bytes(data))
    assert problems(one) == "Dose.raw is not the file hashed on the run host"


@pytest.mark.parametrize(("old", "new", "message"), [
    (b"DimSize = 150 150 150", b"DimSize = 150 150 149", "Dose.mhd DimSize is '150 150 149', expected '150 150 150'"),
    (b"2.000000 2.000000 2.000000", b"2.000000 1.000000 2.000000", "Dose.mhd ElementSpacing is '2.000000 1.000000 2.000000'"),
    (b"MET_FLOAT", b"MET_DOUBLE", "Dose.mhd ElementType is 'MET_DOUBLE'"),
    (b"ElementByteOrderMSB = False", b"ElementByteOrderMSB = True", "Dose.mhd ElementByteOrderMSB is 'True'"),
    (b"NDims = 3", b"NDims = 3\r\nDimSize = 150 150 150", "Dose.mhd sets DimSize more than once"),
])
def test_a_dose_header_that_is_not_the_assumed_grid_is_a_problem_though_correctly_hashed(
    one: tuple[fa.Run, fa.Expected, aa.Ledger], old: bytes, new: bytes, message: str
) -> None:
    rewrite(one[0].path / "out_seed" / "Dose.mhd", old, new)
    record_hashes(one[0].path, "A-port")  # identity holds: the header is the one its record names
    msg = problems(one)
    assert "is not the file hashed" not in msg and message in msg


def test_a_dose_of_another_size_is_a_problem_though_correctly_hashed(one: tuple[fa.Run, fa.Expected, aa.Ledger]) -> None:
    raw = one[0].path / "out_seed" / "Dose.raw"
    raw.write_bytes(raw.read_bytes()[:-4])
    record_hashes(one[0].path, "A-port")
    assert problems(one) == f"Dose.raw is {4 * N**3 - 4} bytes, expected {4 * N**3}"


@pytest.mark.parametrize(("keep", "extra", "listed"), [
    (1, b"", "['Dose.raw']"),
    (2, b"ab" * 32 + b"  Dose.raw\r\n", "['Dose.mhd', 'Dose.raw', 'Dose.raw']"),  # a second hash for one file
    (2, b"ab" * 32 + b"  Energy.raw\r\n", "['Dose.mhd', 'Dose.raw', 'Energy.raw']"),
])
def test_the_hash_record_must_name_both_dose_files_once_each(
    one: tuple[fa.Run, fa.Expected, aa.Ledger], keep: int, extra: bytes, listed: str
) -> None:
    rec = one[0].path / "out_seed" / "sha256.txt"
    rec.write_bytes(b"".join(rec.read_bytes().splitlines(keepends=True)[:keep]) + extra)
    assert f"out_seed/sha256.txt lists {listed}, expected Dose.mhd and Dose.raw once each" in problems(one)


def test_measure_refuses_a_dose_that_changed_after_verification(one: tuple[fa.Run, fa.Expected, aa.Ledger], monkeypatch: pytest.MonkeyPatch) -> None:
    assert problems(one) == ""
    raw = one[0].path / "out_seed" / "Dose.raw"
    raw.write_bytes(b"\x01" + raw.read_bytes()[1:])
    monkeypatch.setattr(psm, "field_endpoints", lambda *_a: pytest.fail("endpoints were computed from an unverified dose"))
    fa.measure(one[0])
    assert one[0].problems == ["Dose.raw changed between verification and reading"]


# ------------------------------------------------------------------------------------------------ what is expected
def test_the_run_lists_this_analysis_imports_must_be_the_acquisition_commits(monkeypatch: pytest.MonkeyPatch) -> None:
    for rel in fa.RUN_LIST_FILES:
        monkeypatch.setattr(fa, "committed_blobs", lambda c, r, rel=rel: {**working_tree_blobs(c, r), rel: "0" * 40})
        with pytest.raises(fa.InputError, match=f"{rel} as this analysis imports it is not the file of the acquisition commit"):
            fa.expected(COMMIT, "full")


@pytest.mark.parametrize(("commit", "mode", "message"), [("ab" * 6, "full", "full 40-hex sha"), (COMMIT, "pilot", "mode must be one of")])
def test_expected_refuses_an_abbreviated_commit_and_an_unknown_mode(commit: str, mode: str, message: str) -> None:
    with pytest.raises(fa.InputError, match=message):
        fa.expected(commit, mode)


def test_committed_blobs_reads_the_tree_of_a_commit(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.undo()  # the real function
    import subprocess

    head = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    blobs = fa.committed_blobs(head, REPO)
    assert sum(p.startswith("Materials/") for p in blobs) == 288 and fa.INPUT_PATHS[2] in blobs
    want = subprocess.run(["git", "-C", str(REPO), "rev-parse", f"{head}:{fa.INPUT_PATHS[2]}"], capture_output=True, text=True, check=True)
    assert blobs[fa.INPUT_PATHS[2]] == want.stdout.strip()
    assert fa.blob_id(b"") == "e69de29bb2d1d6434b8b29ae775ad8c2e48c5391"  # git's empty blob
    with pytest.raises(fa.InputError, match="cannot list the tree"):
        fa.committed_blobs("0" * 40, REPO)


# ------------------------------------------------------------------------------------------- endpoint validity
def test_endpoint_validity_rules() -> None:
    ep = {k: 1.0 for k in fa.endpoint_keys(100)} | {"r80_crossings": 1, "r20_crossings": 1}
    assert all(v == 1.0 for v in fa.endpoint_values(ep, 100).values())
    assert fa.endpoint_values(ep | {"lateral_39_70": 0.0}, 100)["lateral_39_70"] == 0.0  # zero is a value
    assert fa.endpoint_values(ep | {"lateral_39_70": float("nan")}, 100)["lateral_39_70"] is None
    assert fa.endpoint_values(ep | {"lateral_39_5": float("inf")}, 100)["lateral_39_5"] is None
    assert fa.endpoint_values(ep | {"cax_39": 0.0}, 100)["cax_39"] is None
    assert fa.endpoint_values(ep | {"cax_max_invalid": True}, 100)["cax_max"] is None
    assert fa.endpoint_values(ep | {"r80_crossings": 2}, 100)["r80_mm"] is None
    assert fa.endpoint_values(ep | {"r20_crossings": 0, "r20_mm": None}, 100)["r20_mm"] is None
    assert fa.endpoint_values(ep | {"r80_mm": True}, 100)["r80_mm"] is None
    assert len(fa.endpoint_values(ep, 150)) == 17 and fa.endpoint_values(ep, 150)["cax_79"] is None  # absent is invalid


def test_peak_three_row_mean() -> None:
    d = dose(100, random.Random(1))
    cax = d[C - 5:C + 5, :, C - 5:C + 5].mean(axis=(0, 2)).astype(np.float64)
    i = int(np.argmax(cax))
    assert psm.peak_three_row_mean(d) == pytest.approx(float(cax[i - 1:i + 2].mean()), rel=1e-12)
    assert cax[i] == pytest.approx(psm.field_endpoints(d, 150.0, (39,), (5,))["cax_max"], rel=1e-12)
    edge = np.zeros((N, N, N), dtype=np.float32)
    edge[:, 0, :] = 1.0
    assert psm.peak_three_row_mean(edge) is None  # maximum on the first row
    assert psm.peak_three_row_mean(np.zeros((N, N, N), dtype=np.float32)) is None
    with pytest.raises(ValueError, match="expected"):
        psm.peak_three_row_mean(np.zeros((N, N, N - 1), dtype=np.float32))


# ------------------------------------------------------------------------------------- statistics, on values alone
def cells_of(value: object, *, issue: dict[tuple[str, int], list[str]] | None = None) -> tuple[fa.Cells, fa.Issues]:
    """In-memory population: value(arm, energy, i, key) gives each run's endpoint value."""
    cells: fa.Cells = {}
    for arm in fr.ARMS:
        for energy in fr.ENERGIES:
            runs = []
            for i, seed in enumerate(fr.seeds(arm, energy, "full")):
                run = fa.Run(arm, energy, seed, Path("unused"), record={"simulated": 60_000_000})
                run.values = {k: value(arm, energy, i, k) for k in fa.endpoint_keys(energy)}  # type: ignore[operator]
                run.peak3 = 4.0 * (1 + 0.0002 * (i - 3.5))
                runs.append(run)
            cells[(arm, energy)] = runs
    return cells, {k: [] for k in cells} | (issue or {})


def noisy(arm: str, energy: int, i: int, key: str) -> float:
    """Every arm the same but for a deterministic spread: all 34 equivalent."""
    wobble = (i - 3.5) / 3.5  # -1 .. 1 over the 8 runs
    if key.startswith("lateral_"):
        return 5.0 + 0.02 * wobble
    if key.startswith("cax_"):
        return 2.0 * (1 + 0.0005 * wobble)
    return 77.0 + 0.01 * wobble


def test_all_equivalent_gives_both_joint_claims_and_holm_over_34() -> None:
    cells, issues = cells_of(noisy)
    c = fa.analyse_contrast(cells, issues, "A-port", "A-up", confirmatory=True)
    assert len(c["rows"]) == 34 and {r["outcome"] for r in c["rows"]} == {"equivalent"}
    assert c["claims"] == {
        "joint_claim_equivalent_on_all_endpoints": True, "secondary_joint_claim_holm_excluding_all_zero_rows": True,
        "all_zero_rows_excluded_from_the_secondary_claim": [], "n_equivalent_unadjusted": 34, "n_equivalent_holm": 34, "n_not_established": 0}
    smallest = min(c["rows"], key=lambda r: r["p_tost"])
    assert smallest["holm_p"] == pytest.approx(min(1.0, 34 * smallest["p_tost"]))


def test_a_difference_row_is_hand_checkable() -> None:
    def value(arm: str, energy: int, i: int, key: str) -> float:
        if key == "lateral_39_10":
            return [1.9, 2.0, 2.1, 2.0, 1.9, 2.1, 2.0, 2.0][i] if arm == "A-port" else 1.5
        return noisy(arm, energy, i, key)

    cells, issues = cells_of(value)
    row = next(r for r in fa.analyse_contrast(cells, issues, "A-port", "A-up", confirmatory=True)["rows"] if r["endpoint"] == "E100/lateral_39_10")
    a = [1.9, 2.0, 2.1, 2.0, 1.9, 2.1, 2.0, 2.0]
    se = math.sqrt(st.variance(a) / 8)  # the reference has zero variance, so df = 7
    assert (row["arm_mean"], row["reference_mean"]) == (pytest.approx(2.0), 1.5)
    assert row["estimate"] == pytest.approx(0.5) and row["se"] == pytest.approx(se) and row["df"] == pytest.approx(7.0)
    half = float(student_t.ppf(0.95, 7)) * se
    assert row["ci90"] == pytest.approx([0.5 - half, 0.5 + half])
    assert row["outcome"] == "not_equivalent" and row["margin"] == [-0.2, 0.2]  # the whole interval is above +0.20
    assert row["p_tost"] == pytest.approx(float(student_t.sf((0.2 - 0.5) / se, 7)))


def test_a_central_axis_row_is_a_geometric_mean_ratio() -> None:
    def value(arm: str, energy: int, i: int, key: str) -> float:
        if key == "cax_max":
            return (2.0 if i % 2 else 2.002) * (1.0049 if arm == "B-pgcc" else 1.0)
        return noisy(arm, energy, i, key)

    cells, issues = cells_of(value)
    row = next(r for r in fa.analyse_contrast(cells, issues, "B-pgcc", "B-up", confirmatory=True)["rows"] if r["endpoint"] == "E150/cax_max")
    assert row["scale"] == "ratio" and row["estimate"] == pytest.approx(1.0049)
    assert row["margin"] == pytest.approx([0.995, 1.005])
    assert row["ci90"][0] < 1.005 < row["ci90"][1] and row["outcome"] == "inconclusive"  # the interval straddles the margin


def test_all_runs_zero_in_both_arms_is_not_established_and_only_the_secondary_claim_can_hold() -> None:
    def value(arm: str, energy: int, i: int, key: str) -> float:
        return 0.0 if key.endswith("_70") and energy == 100 else noisy(arm, energy, i, key)

    cells, issues = cells_of(value)
    c = fa.analyse_contrast(cells, issues, "B-picc", "B-up", confirmatory=True)
    zero = [r for r in c["rows"] if r["all_runs_zero"]]
    assert [r["endpoint"] for r in zero] == ["E100/lateral_39_70", "E100/lateral_61_70"]
    for r in zero:
        assert r["outcome"] == "not_established" and r["p_tost"] is None and r["holm_p"] is None and r["holm_decision"] == "not_established"
        assert "ci90" not in r and r["estimate"] == 0.0
        assert r["reason"] == "all runs zero in both arms (simulated histories: arm 480000000, reference 480000000)"
    assert c["claims"]["joint_claim_equivalent_on_all_endpoints"] is False
    assert c["claims"]["secondary_joint_claim_holm_excluding_all_zero_rows"] is True
    assert c["claims"]["all_zero_rows_excluded_from_the_secondary_claim"] == ["E100/lateral_39_70", "E100/lateral_61_70"]
    assert (c["claims"]["n_equivalent_unadjusted"], c["claims"]["n_not_established"]) == (32, 2)
    # the two rows stay in the Holm family: the smallest p is still multiplied by 34
    smallest = min((r for r in c["rows"] if r["p_tost"] is not None), key=lambda r: r["p_tost"])
    assert smallest["holm_p"] == pytest.approx(min(1.0, 34 * smallest["p_tost"]))


def test_the_secondary_claim_needs_every_retained_row_to_survive_holm_over_all_34() -> None:
    """#59 review 7116: rows are set aside by the data, so the unadjusted tests give the secondary claim no error
    control. One retained row equivalent at an unadjusted p near 0.03 does not survive Holm, and the claim fails."""
    def value(arm: str, energy: int, i: int, key: str) -> float:
        if key.endswith("_70") and energy == 100:
            return 0.0
        if key == "lateral_79_30":
            return 5.0 + 0.2815 * (i - 3.5) / 3.5  # SE about 0.0985 against a margin of 0.2
        return noisy(arm, energy, i, key)

    cells, issues = cells_of(value)
    c = fa.analyse_contrast(cells, issues, "B-picc", "B-up", confirmatory=True)
    weak = next(r for r in c["rows"] if r["endpoint"] == "E150/lateral_79_30")
    assert 0.05 / 34 < weak["p_tost"] < 0.05 and weak["outcome"] == "equivalent"
    assert weak["holm_p"] > 0.05 and weak["holm_decision"] == "not shown"
    retained = [r for r in c["rows"] if not r["all_runs_zero"]]
    assert len(retained) == 32 and all(r["outcome"] == "equivalent" for r in retained)  # unadjusted, every retained row passes
    assert c["claims"]["secondary_joint_claim_holm_excluding_all_zero_rows"] is False
    assert c["claims"]["joint_claim_equivalent_on_all_endpoints"] is False
    assert (c["claims"]["n_equivalent_unadjusted"], c["claims"]["n_equivalent_holm"]) == (32, 31)
    md = "\n".join(fa.markdown({"acquisition_commit": COMMIT, "label": fa.LABEL, "contrasts": [c],
                                "collection": {"population_issues_before_any_dose_was_read": 0, "declared_closed": None}}))
    assert "is equivalent after Holm over all 34): **NOT ESTABLISHED**" in md


def test_identical_non_zero_values_in_both_arms_are_not_established_and_not_excluded() -> None:
    def value(arm: str, energy: int, i: int, key: str) -> float:
        return 0.25 if key == "lateral_79_50" else noisy(arm, energy, i, key)

    cells, issues = cells_of(value)
    c = fa.analyse_contrast(cells, issues, "A-port", "A-up", confirmatory=True)
    row = next(r for r in c["rows"] if r["endpoint"] == "E150/lateral_79_50")
    assert row["outcome"] == "not_established" and row["all_runs_zero"] is False
    assert row["reason"].startswith("zero sample variance in both arms (simulated histories: arm 480000000")
    assert c["claims"]["joint_claim_equivalent_on_all_endpoints"] is False
    assert c["claims"]["secondary_joint_claim_holm_excluding_all_zero_rows"] is False  # only all-zero rows are set aside


def test_zero_variance_in_one_arm_only_is_evaluated() -> None:
    def value(arm: str, energy: int, i: int, key: str) -> float:
        return 5.0 if key == "lateral_39_5" and arm == "A-up" else noisy(arm, energy, i, key)

    cells, issues = cells_of(value)
    row = next(r for r in fa.analyse_contrast(cells, issues, "A-port", "A-up", confirmatory=True)["rows"] if r["endpoint"] == "E100/lateral_39_5")
    assert row["outcome"] == "equivalent" and row["df"] == pytest.approx(7.0)


def test_an_invalid_value_in_one_run_makes_the_endpoint_not_established() -> None:
    def value(arm: str, energy: int, i: int, key: str) -> float | None:
        return None if (arm, energy, i, key) == ("A-up", 150, 2, "r20_mm") else noisy(arm, energy, i, key)

    cells, issues = cells_of(value)
    c = fa.analyse_contrast(cells, issues, "A-port", "A-up", confirmatory=True)
    row = next(r for r in c["rows"] if r["endpoint"] == "E150/r20_mm")
    assert row["outcome"] == "not_established" and row["reason"] == "reference: invalid or unusable in 1 of 8 runs: A-up 150 MeV s960053"
    assert c["claims"]["joint_claim_equivalent_on_all_endpoints"] is False and c["claims"]["n_not_established"] == 1
    assert c["claims"]["secondary_joint_claim_holm_excluding_all_zero_rows"] is False


def test_a_population_issue_withholds_the_claims_of_the_contrasts_it_touches_only() -> None:
    cells, issues = cells_of(noisy, issue={("B-up", 150): ["B-up 150 MeV: run s962058 is missing"]})
    cells[("B-up", 150)].pop()
    for arm, ref in (("B-pgcc", "B-up"), ("B-picc", "B-up")):
        c = fa.analyse_contrast(cells, issues, arm, ref, confirmatory=True)
        assert c["claims_withheld"] is True and c["claims"] == {} and c["partial_reasons"] == ["B-up 150 MeV: run s962058 is missing"]
        assert {r["outcome"] for r in c["rows"]} == {"descriptive"} and all("holm_p" not in r and r["p_tost"] is None for r in c["rows"])
    c = fa.analyse_contrast(cells, issues, "A-port", "A-up", confirmatory=True)
    assert c["claims_withheld"] is False and c["claims"]["joint_claim_equivalent_on_all_endpoints"] is True


def test_the_descriptive_contrast_makes_no_claim() -> None:
    cells, issues = cells_of(noisy)
    c = fa.analyse_contrast(cells, issues, "B-pgcc", "B-picc", confirmatory=False)
    assert c["claims"] == {} and c["claims_withheld"] is False and {r["outcome"] for r in c["rows"]} == {"descriptive"}


def test_ratio_of_arm_means_by_the_delta_method_hand_checked() -> None:
    a, b = [0.2, 0.4], [0.1, 0.1, 0.4]  # means 0.3 and 0.2
    out = fa.ratio_of_means(a, b)
    ra = (st.variance(a) / 2) / 0.3**2  # squared SE of ln(mean a): 0.01 / 0.09
    rb = (st.variance(b) / 3) / 0.2**2  # 0.01 / 0.04
    assert (ra, rb) == (pytest.approx(1 / 9), pytest.approx(1 / 4))
    df = (ra + rb) ** 2 / (ra**2 / 1 + rb**2 / 2)
    h = float(student_t.ppf(0.975, df)) * math.sqrt(ra + rb)
    assert out["status"] == "ratio" and out["estimate"] == pytest.approx(1.5) and out["df"] == pytest.approx(df)
    assert out["ci95"] == pytest.approx([1.5 * math.exp(-h), 1.5 * math.exp(h)])
    assert fa.ratio_of_means([0.0, 0.0], [0.1, 0.2])["status"] == "not given"
    assert fa.ratio_of_means([0.1, 0.1], [0.2, 0.2]) == {
        "status": "ratio without interval", "estimate": 0.5, "arm_mean": 0.1, "reference_mean": 0.2, "reason": "zero sample variance in both arms"}


def test_the_ratio_rows_are_the_50_and_70_mm_offsets_and_enter_no_claim() -> None:
    cells, issues = cells_of(noisy)
    c = fa.analyse_contrast(cells, issues, "B-pgcc", "B-up", confirmatory=True)
    assert [r["endpoint"] for r in c["ratio_of_arm_means_at_50_and_70_mm"]] == [
        f"E{e}/lateral_{d}_{o}" for e, depths in ((100, (39, 61)), (150, (79, 125))) for d in depths for o in (50, 70)]
    assert all(r["status"] == "ratio" and r["estimate"] == pytest.approx(1.0) for r in c["ratio_of_arm_means_at_50_and_70_mm"])


def test_the_three_row_ratio_sits_beside_the_maximum_and_the_r80_difference() -> None:
    cells, issues = cells_of(noisy)
    for r in cells[("A-port", 100)]:
        r.peak3 *= 1.002  # type: ignore[operator]
    c = fa.analyse_contrast(cells, issues, "A-port", "A-up", confirmatory=True)
    m100, m150 = c["maximum"]
    assert (m100["energy"], m150["energy"]) == (100, 150) and m100["maximum_outcome"] == "equivalent"
    assert m100["three_row_mean_ratio"]["estimate"] == pytest.approx(1.002) and m100["three_row_mean_ratio"]["interval_within_margin"] is True
    assert m150["three_row_mean_ratio"]["estimate"] == pytest.approx(1.0)
    assert m100["r80_difference_mm"] == pytest.approx(0.0, abs=1e-12) and len(m100["r80_ci90"]) == 2
    for r in cells[("A-port", 150)]:
        r.peak3 *= 1.0049  # type: ignore[operator]  # the 95% interval now reaches beyond 1.005
    wide = fa.analyse_contrast(cells, issues, "A-port", "A-up", confirmatory=True)["maximum"][1]["three_row_mean_ratio"]
    assert wide["ci95"][0] < 1.005 < wide["ci95"][1] and wide["interval_within_margin"] is False
    cells[("A-up", 100)][0].peak3 = None
    assert fa.analyse_contrast(cells, issues, "A-port", "A-up", confirmatory=True)["maximum"][0]["three_row_mean_ratio"]["status"] == "not given"


def test_the_reading_of_a_maximum_is_printed_only_when_it_falls_short_and_the_three_row_interval_is_inside() -> None:
    def value(arm: str, energy: int, i: int, key: str) -> float:
        if key == "cax_max" and energy == 100:
            return (2.0 if i % 2 else 2.002) * (1.0049 if arm == "A-port" else 1.0)
        return noisy(arm, energy, i, key)

    cells, issues = cells_of(value)
    doc = {"acquisition_commit": COMMIT, "label": fa.LABEL, "contrasts": [fa.analyse_contrast(cells, issues, "A-port", "A-up", confirmatory=True)],
           "collection": {"population_issues_before_any_dose_was_read": 0, "declared_closed": None}}
    md = "\n".join(fa.markdown(doc))
    assert "At 100 MeV the maximum is inconclusive while the three-row ratio's 95% interval lies within [0.995, 1.005]" in md
    assert "At 150 MeV the maximum" not in md and "No cause is attributed." in md
    assert "Joint claim (all 34 endpoints equivalent): **NOT ESTABLISHED**" in md


@pytest.mark.parametrize(("x", "digits", "text"), [
    (0.0, 4, "0"), (-0.0, 4, "0"), (0, 4, "0"),
    (3.2e-7, 4, "3.2e-07"), (-4.56789e-5, 4, "-4.57e-05"), (0.000123456, 4, "0.000123"),
    (0.00999, 4, "0.00999"), (0.01, 4, "0.0100"), (-0.01, 4, "-0.0100"), (1.23456, 4, "1.2346"), (5, 4, "5.0000"),
    (0.000999, 5, "0.000999"), (0.001, 5, "0.00100"), (1.000049, 5, "1.00005"),
    (None, 4, ""), (True, 4, ""), ("1.0", 4, ""),
])
def test_a_table_number_keeps_three_significant_figures_when_small_and_only_exact_zero_prints_as_zero(
        x: object, digits: int, text: str) -> None:
    assert fa._num(x, digits) == text


def test_a_table_interval_is_formatted_as_its_numbers_are() -> None:
    assert fa._pair([-3.2e-7, 4.1e-6]) == "[-3.2e-07, 4.1e-06]"
    assert fa._pair([-0.2, 0.2]) == "[-0.2000, 0.2000]" and fa._pair([0.995, 1.005], 5) == "[0.99500, 1.00500]"
    assert fa._pair(None) == ""


def test_very_small_lateral_means_are_legible_in_the_table() -> None:
    """#60: at 100 MeV, 70 mm, the dose can be a few millionths of a point; the table must show how small."""
    def value(arm: str, energy: int, i: int, key: str) -> float:
        if key.endswith("_70") and energy == 100:
            return (3.0e-6 if arm == "B-picc" else 2.0e-6) * (1 + 0.1 * (i - 3.5) / 3.5)
        return noisy(arm, energy, i, key)

    cells, issues = cells_of(value)
    doc = {"acquisition_commit": COMMIT, "label": fa.LABEL, "contrasts": [fa.analyse_contrast(cells, issues, "B-picc", "B-up", confirmatory=True)],
           "collection": {"population_issues_before_any_dose_was_read": 0, "declared_closed": None}}
    lines = fa.markdown(doc)
    row = next(line for line in lines if line.startswith("| E100/lateral_39_70 |")).split(" | ")
    assert row[1:4] == ["3e-06", "2e-06", "1e-06"] and row[6] == "equivalent"
    assert row[4].startswith("[") and "e-0" in row[4] and row[5] == "[-0.2000, 0.2000] (difference)"
    ratio = [line for line in lines if line.startswith("| E100/lateral_39_70 |")][1].split(" | ")
    assert ratio[1:4] == ["3e-06", "2e-06", "1.5000"]
    assert not any("0.0000" in line for line in lines if line.startswith("| E100/lateral_"))


# ------------------------------------------------------------------------------------------ the whole population
def test_status_verifies_every_run_without_reading_a_dose_value(trees: list[Path], monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fa, "measure", lambda _run: pytest.fail("a dose value was read in --status"))
    ok, lines = fa.status(trees, COMMIT, "full")
    assert ok and lines == [f"{a} {e} MeV: 8 of 8 verified" for a in fr.ARMS for e in fr.ENERGIES]


def test_end_to_end_population_fingerprint_and_claims(document: dict) -> None:
    assert document["acquisition_commit"] == COMMIT and document["histories_per_run"] == 60_000_000
    assert len(document["runs"]) == 80 and all(r["usable"] for r in document["runs"])
    assert all(p == {"expected": 8, "verified": 8, "issues": []} for p in document["population"].values())
    assert document["dataset_fingerprint"]["n_files"] == 80 * (9 + 291) + 10 * 2  # every file read, the doses included
    assert [(c["arm"], c["reference"], c["confirmatory"]) for c in document["contrasts"]] == list(fr.CONTRASTS)
    for c in document["contrasts"][:3]:
        # the synthetic 100 MeV field scores nothing 70 mm out: the foreseeable case of the plan
        assert c["claims"]["all_zero_rows_excluded_from_the_secondary_claim"] == ["E100/lateral_39_70", "E100/lateral_61_70"]
        assert c["claims"]["joint_claim_equivalent_on_all_endpoints"] is False
        assert c["claims"]["secondary_joint_claim_holm_excluding_all_zero_rows"] is True
        assert c["claims"]["n_equivalent_unadjusted"] == 32
        notes = {r["endpoint"]: r for r in c["ratio_of_arm_means_at_50_and_70_mm"]}
        assert notes["E100/lateral_39_70"]["status"] == "not given" and notes["E150/lateral_79_70"]["status"] == "ratio"
    assert document["contrasts"][3]["claims"] == {}
    assert document["analysis_files"]["validation/field_followup_analyse.py"]["unchanged"] is True  # blobs are the working tree's here


def test_end_to_end_rows_match_an_independent_computation(document: dict, trees: list[Path]) -> None:
    def metrics(arm: str, energy: int) -> list[dict]:
        root = trees[0 if arm in WINDOWS else 1] / fr.ARMS[arm].directory / f"e{energy}"
        return [psm.field_endpoints(np.fromfile(root / f"s{s}" / "out_seed" / "Dose.raw", dtype=np.float32).reshape(N, N, N),
                                    150.0, fa.DEPTHS_MM[energy], fa.OFFSETS_MM) for s in fr.seeds(arm, energy, "full")]

    rows = {r["endpoint"]: r for r in document["contrasts"][1]["rows"]}  # B-pgcc against B-up
    a, b = metrics("B-pgcc", 150), metrics("B-up", 150)
    for key in ("lateral_79_5", "lateral_125_50", "r80_mm"):
        xa, xb = [m[key] for m in a], [m[key] for m in b]
        se = math.sqrt(st.variance(xa) / 8 + st.variance(xb) / 8)
        df = se**4 / ((st.variance(xa) / 8) ** 2 / 7 + (st.variance(xb) / 8) ** 2 / 7)
        row = rows[f"E150/{key}"]
        assert row["estimate"] == pytest.approx(st.mean(xa) - st.mean(xb), abs=1e-12)
        half = float(student_t.ppf(0.95, df)) * se
        assert row["ci90"] == pytest.approx([row["estimate"] - half, row["estimate"] + half])
        assert row["arm_mean"] == pytest.approx(st.mean(xa)) and row["outcome"] == "equivalent"
    la, lb = [math.log(m["cax_125"]) for m in a], [math.log(m["cax_125"]) for m in b]
    assert rows["E150/cax_125"]["estimate"] == pytest.approx(math.exp(st.mean(la) - st.mean(lb)))
    assert 75.0 < a[0]["r80_mm"] < 170.0 and a[0]["lateral_79_5"] > 10  # the synthetic field is where the estimators look


def test_an_incomplete_population_is_refused_before_any_dose_value_is_read(trees: list[Path], monkeypatch: pytest.MonkeyPatch) -> None:
    """#59 review 7116: with one frozen seed missing and both hosts' trees present, 79 runs were measured."""
    monkeypatch.setattr(fa, "measure", lambda _run: pytest.fail("a dose value was read although the population is incomplete"))
    cell = trees[0] / "apt" / "e150"
    (cell / "s961058").rename(cell.parent / "held_s961058")  # one frozen run is absent; nothing else is wrong
    try:
        with pytest.raises(fa.InputError, match=r"not complete and verified \(\d+ issue\(s\)\); no dose value was read") as e:
            fa.run(trees, COMMIT)
        with pytest.raises(fa.InputError, match="needs a statement"):
            fa.run(trees, COMMIT, closed="  ")
    finally:
        (cell.parent / "held_s961058").rename(cell / "s961058")
    assert "A-port 150 MeV: run s961058 is missing" in str(e.value)


def test_the_whole_population_is_verified_before_the_first_dose_value_is_read(trees: list[Path], monkeypatch: pytest.MonkeyPatch) -> None:
    order: list[str] = []
    verify, measure = fa.verify_run, fa.measure
    monkeypatch.setattr(fa, "verify_run", lambda run, *a: (order.append("verify"), verify(run, *a))[1])
    monkeypatch.setattr(fa, "measure", lambda run: (order.append("measure"), measure(run))[1])
    fa.run(trees, COMMIT)
    assert order == ["verify"] * 80 + ["measure"] * 80


def test_a_missing_and_an_unexpected_run_make_the_lenovo_contrast_partial_only(trees: list[Path]) -> None:
    cell = trees[0] / "apt" / "e150"
    (cell / "s961058").rename(cell / "s961059")  # the frozen seed is gone and an unlisted one is there
    statement = "connor-227743e6: both workflows finished (runs 1 and 2); archives verified on the share"
    try:
        doc = fa.run(trees, COMMIT, closed=statement)
    finally:
        (cell / "s961059").rename(cell / "s961058")
    assert doc["collection"] == {"state": "closed by declaration, with population issues",
                                 "population_issues_before_any_dose_was_read": 2, "declared_closed": statement}
    assert f"The collection was declared closed: {statement}**" in "\n".join(fa.markdown(doc))
    a, pg, pi, _d = doc["contrasts"]
    assert a["claims_withheld"] is True and a["claims"] == {}
    assert a["partial_reasons"] == ["A-port 150 MeV: run s961058 is missing", "A-port 150 MeV: run s961059 is not in the frozen list"]
    assert {r["outcome"] for r in a["rows"]} == {"descriptive", "not_established"}
    assert pg["claims_withheld"] is False and pi["claims"]["n_equivalent_unadjusted"] == 32
    assert doc["population"]["A-port 150 MeV"]["verified"] == 7


def test_a_run_hosts_snapshot_of_the_run_lists_must_be_the_acquisition_commits(trees: list[Path], monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fa, "verify_run", lambda *_a: None)  # the cell-level check alone
    snap = trees[1] / "bpi" / "e100" / "snapshot" / fa.RUN_LIST_FILES[0]
    original = snap.read_bytes()
    snap.write_bytes(original + b"# changed\n")
    try:
        _cells, issues, _ledger, _notes = fa.load(trees, fa.expected(COMMIT, "full"))
    finally:
        snap.write_bytes(original)
    assert issues[("B-picc", 100)] == [f"B-picc 100 MeV: the run host's snapshot of {fa.RUN_LIST_FILES[0]} is not the acquisition commit's"]
    assert not any(v for k, v in issues.items() if k != ("B-picc", 100))


@pytest.mark.parametrize(("where", "name", "holds_run", "cells_touched"), [
    ("cell", "s9610417", True, [("A-port", 150)]),      # a seven-digit seed directory holding run.json (review 7116)
    ("cell", "s96105", False, [("A-port", 150)]),       # named like a run, empty
    ("cell", "old", True, [("A-port", 150)]),           # not named like a run, but it holds one
    ("cell", "snapshot/copy", True, [("A-port", 150)]), # run data hidden inside a known entry
    ("arm", "e150_first_try", True, [("A-port", 100), ("A-port", 150)]),
    ("root", "apt_old", True, [("A-port", 100), ("A-port", 150), ("A-up", 100), ("A-up", 150)]),
])
def test_run_data_that_is_not_a_frozen_run_is_a_population_issue_wherever_it_sits(
        trees: list[Path], monkeypatch: pytest.MonkeyPatch, where: str, name: str, holds_run: bool, cells_touched: list[tuple[str, int]]) -> None:
    monkeypatch.setattr(fa, "verify_run", lambda *_a: None)  # the layout check alone
    base = {"cell": trees[0] / "apt" / "e150", "arm": trees[0] / "apt", "root": trees[0]}[where]
    stray = base / name
    (stray / "deeper").mkdir(parents=True)
    if holds_run:
        (stray / "deeper" / "run.json").write_text("{}")
    try:
        _cells, issues, _ledger, notes = fa.load(trees, fa.expected(COMMIT, "full"))
    finally:
        shutil.rmtree(base / name.split("/")[0] if "/" not in name else stray)
    top = name.split("/")[0]
    for key, why in issues.items():
        assert (len(why) == 1 and f"{top} looks like run data and is not a run of the frozen lists" in why[0]) if key in cells_touched else why == []
    assert not any(top in n for n in notes)


def test_an_entry_that_is_not_run_data_is_noted_and_is_no_issue(trees: list[Path], monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fa, "verify_run", lambda *_a: None)
    extra = [trees[0] / "apt" / "e150" / "notes.txt", trees[0] / "apt" / "README", trees[0] / "listing.txt"]
    for f in extra:
        f.write_text("not a run\n")
    try:
        _cells, issues, _ledger, notes = fa.load(trees, fa.expected(COMMIT, "full"))
    finally:
        for f in extra:
            f.unlink()
    assert not any(issues.values())
    assert sum("is ignored (not run data)" in n for n in notes) == 3


def test_the_real_layout_of_both_hosts_has_no_stray_entry(trees: list[Path], template: Template, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The Lenovo's cells also hold inputs.tar; every cell holds the case, the snapshot and their hash files."""
    monkeypatch.setattr(fa, "verify_run", lambda *_a: None)
    made = []
    for root in trees:
        for cell in (c for arm in root.iterdir() for c in arm.iterdir()):
            for name in ("fcase_sha256.txt", "snapshot_sha256.txt", "run_root.txt", "inputs.tar"):
                (cell / name).write_text("x\n")
                made.append(cell / name)
            (cell / "fcase").mkdir()
            made.append(cell / "fcase")
    try:
        _cells, issues, _ledger, notes = fa.load(trees, fa.expected(COMMIT, "full"))
    finally:
        for m in made:
            m.rmdir() if m.is_dir() else m.unlink()
    assert not any(issues.values()) and notes == []


def test_two_runs_with_the_same_dose_are_a_duplicated_run(trees: list[Path], monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fa, "measure", lambda _run: pytest.fail("load() read a dose value"))
    a, b = trees[1] / "bpg" / "e100" / "s963041" / "out_seed", trees[1] / "bup" / "e150" / "s962052" / "out_seed"
    kept = {p: p.read_bytes() for p in (b / "Dose.raw", b / "sha256.txt")}
    (b / "Dose.raw").write_bytes((a / "Dose.raw").read_bytes())
    record_hashes(b.parent, "B-up")  # each run agrees with its own record; only the pair gives it away
    try:
        _cells, issues, _ledger, _notes = fa.load(trees, fa.expected(COMMIT, "full"))
    finally:
        for p, data in kept.items():
            p.write_bytes(data)
    why = "B-pgcc 100 MeV s963041 and B-up 150 MeV s962052 have the same Dose.raw"
    assert issues[("B-up", 150)] == [why] and issues[("B-pgcc", 100)] == [why]
    assert not any(v for k, v in issues.items() if k not in {("B-up", 150), ("B-pgcc", 100)})


def test_one_hosts_tree_alone_is_refused_before_any_run_is_read(trees: list[Path], monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fa, "verify_run", lambda *_a: pytest.fail("a run was read although a host's tree is absent"))
    with pytest.raises(fa.InputError, match=r"arm B-up \(directory bup\) is in 0 of the roots"):
        fa.run(trees[:1], COMMIT)
    with pytest.raises(fa.InputError, match=r"arm A-up \(directory aup\) is in 2 of the roots"):
        fa.run([trees[0], trees[0], trees[1]], COMMIT)
    with pytest.raises(fa.InputError, match="is not a directory"):
        fa.run([trees[0] / "nowhere", trees[1]], COMMIT)


def test_smoke_trees_are_verified_for_completion_only(template: Template, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    roots = [tmp_path / "lenovo" / "smoke", tmp_path / "hp" / "smoke"]
    for arm in fr.ARMS:
        for energy in fr.ENERGIES:
            cell = roots[0 if arm in WINDOWS else 1] / fr.ARMS[arm].directory / f"e{energy}"
            for rel in fa.RUN_LIST_FILES:
                (cell / "snapshot" / rel).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(REPO / rel, cell / "snapshot" / rel)
            seed = fr.seeds(arm, energy, "smoke")[0]
            write_run(template, cell, arm, energy, seed, "smoke", np.full((N, N, N), seed, dtype="<f4"))
    args = ["--root", str(roots[0]), "--root", str(roots[1]), "--commit", COMMIT]
    assert fa.main([*args, "--mode", "smoke", "--status"]) == 0
    out = capsys.readouterr().out
    assert "A-up 100 MeV: 1 of 1 verified" in out and "VERIFIED: mode smoke; no dose value was read" in out
    assert fa.main([*args, "--mode", "full", "--status"]) == 1  # a smoke tree is not the full population
    assert "run.json mode is 'smoke', expected 'full'" not in capsys.readouterr().out  # the smoke seeds are simply not the full list
    with pytest.raises(SystemExit):
        fa.main([*args, "--mode", "smoke", "--json", str(tmp_path / "o.json"), "--md", str(tmp_path / "o.md")])
    assert not (tmp_path / "o.json").exists()


def test_main_writes_the_document_and_the_table(trees: list[Path], tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    out, md = tmp_path / "o.json", tmp_path / "o.md"
    args = ["--root", str(trees[0]), "--root", str(trees[1]), "--commit", COMMIT]
    assert fa.main([*args, "--json", str(out), "--md", str(md)]) == 0
    doc = json.loads(out.read_text())
    assert len(doc["contrasts"]) == 4 and len(doc["contrasts"][0]["rows"]) == 34
    text_md = md.read_text()
    assert "## A-port against A-up (confirmatory)" in text_md and "## B-pgcc against B-picc (descriptive)" in text_md
    assert "rows excluded: E100/lateral_39_70, E100/lateral_61_70" in text_md
    assert "not_established (all runs zero in both arms (simulated histories: arm 480000000, reference 480000024))" in text_md
    assert fa.main(["--root", str(trees[0]), "--commit", COMMIT, "--status"]) == 2
    assert "REFUSED" in capsys.readouterr().err


# ------------------------------------------------------------------------------------------------ acquisition amendment 1
RERUN = "cd" * 20


@pytest.fixture(scope="session")
def rerun_trees(template: Template, tmp_path_factory: pytest.TempPathFactory) -> list[Path]:
    """The rerun of the two cells cut by the timeout, as the two workflows write it at the rerun commit."""
    home = tmp_path_factory.mktemp("field_e_rerun")
    roots = [home / "lenovo", home / "hp"]
    for arm, energy in fa.RERUN_CELLS:
        cell = roots[0 if arm in WINDOWS else 1] / fr.ARMS[arm].directory / f"e{energy}"
        for rel in fa.RUN_LIST_FILES:
            (cell / "snapshot" / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(REPO / rel, cell / "snapshot" / rel)
        for seed in fr.seeds(arm, energy, "full"):
            write_run(template, cell, arm, energy, seed, "full", dose(energy, random.Random(seed + 7)), {"commit": RERUN})
    return roots


class Interrupted:
    """The original trees as the timeout left them: in each rerun cell the 8th run absent and the 7th without run.json
    or dose. What is taken out is parked beside the roots, outside both; restored on exit."""

    def __init__(self, trees: list[Path]) -> None:
        self.moves: list[tuple[Path, Path]] = []
        park = trees[0].parent / "held"
        park.mkdir(exist_ok=True)
        for arm, energy in fa.RERUN_CELLS:
            cell = trees[0 if arm in WINDOWS else 1] / fr.ARMS[arm].directory / f"e{energy}"
            *_, seventh, eighth = fr.seeds(arm, energy, "full")
            self.moves.append((cell / f"s{eighth}", park / f"s{eighth}"))
            for name in ("run.json", "out_seed"):
                self.moves.append((cell / f"s{seventh}" / name, park / f"{seventh}_{name}"))

    def __enter__(self) -> None:
        for src, dst in self.moves:
            src.rename(dst)

    def __exit__(self, *_exc: object) -> None:
        for src, dst in reversed(self.moves):
            dst.rename(src)


@pytest.fixture
def rerun_ok(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fa, "rerun_commit_problems", lambda *_a: [])


def test_without_the_rerun_the_interrupted_population_is_not_verified(trees: list[Path]) -> None:
    with Interrupted(trees):
        ok, lines = fa.status(trees, COMMIT, "full")
    assert not ok
    assert "  A-port 150 MeV: run s961058 is missing" in lines and "B-pgcc 150 MeV: 6 of 8 verified" in lines


@pytest.mark.usefixtures("rerun_ok")
def test_the_rerun_cells_come_from_the_rerun_and_only_from_it(trees: list[Path], rerun_trees: list[Path]) -> None:
    with Interrupted(trees):
        ok, lines = fa.status(trees, COMMIT, "full", rerun_roots=rerun_trees, rerun_commit=RERUN)
    assert ok, lines
    assert "A-port 150 MeV: 8 of 8 verified" in lines and "B-pgcc 150 MeV: 8 of 8 verified" in lines
    assert any(x.startswith("note: A-port 150 MeV: read from the rerun at cdcdcdcdcdcd") for x in lines)


@pytest.mark.usefixtures("rerun_ok")
def test_a_rerun_run_at_the_acquisition_commit_and_an_original_run_at_the_rerun_commit_are_problems(
        trees: list[Path], rerun_trees: list[Path]) -> None:
    # The rerun's runs are checked against the rerun commit: the original roots, read as the rerun, fail the commit check.
    ok, lines = fa.status(trees, COMMIT, "full", rerun_roots=trees, rerun_commit=RERUN)
    assert not ok
    assert any("A-port 150 MeV s961051: run.json commit is" in x and "expected 'cdcd" in x for x in lines)
    # and the rerun's runs, read under the acquisition commit, fail it
    ok, lines = fa.status(trees, COMMIT, "full", rerun_roots=rerun_trees, rerun_commit=COMMIT[:-1] + "c")
    assert not ok and any("B-pgcc 150 MeV s963051: run.json commit is 'cdcd" in x for x in lines)


@pytest.mark.usefixtures("rerun_ok")
@pytest.mark.parametrize("extra", ["apt/e100", "aup/e150", "bpg/e100/s963001"])
def test_anything_but_the_rerun_cells_in_a_rerun_root_is_a_population_issue(
        trees: list[Path], rerun_trees: list[Path], tmp_path: Path, extra: str) -> None:
    roots = [tmp_path / "lenovo", tmp_path / "hp"]
    for src, dst in zip(rerun_trees, roots, strict=True):
        shutil.copytree(src, dst, copy_function=os.link)
    host = roots[0] if extra.startswith("a") else roots[1]
    (host / extra).mkdir(parents=True)
    with Interrupted(trees):
        ok, lines = fa.status(trees, COMMIT, "full", rerun_roots=roots, rerun_commit=RERUN)
    assert not ok
    assert any(("is not a rerun" in x or "run data" in x) for x in lines), lines


def test_rerun_root_and_rerun_commit_go_together(trees: list[Path]) -> None:
    with pytest.raises(fa.InputError, match="go together"):
        fa.status(trees, COMMIT, "full", rerun_roots=trees)
    with pytest.raises(fa.InputError, match="go together"):
        fa.status(trees, COMMIT, "full", rerun_commit=RERUN)


def _third_root(tmp_path: Path, content: str) -> Path:
    """A further root holding `content` and no arm of its own: known-arm run data, or an unrecognised run directory."""
    third = tmp_path / "third"
    if content == "known-arm run":
        (third / "aup" / "e150" / "s960051").mkdir(parents=True)
        (third / "aup" / "e150" / "s960051" / "run.json").write_text("{}")
    else:
        (third / "stuff" / "x").mkdir(parents=True)
        for name in ("run.json", "Dose.raw"):
            (third / "stuff" / "x" / name).write_bytes(b"0")
    return third


@pytest.mark.usefixtures("rerun_ok")
@pytest.mark.parametrize("content", ["known-arm run", "unrecognised run"])
def test_a_rerun_root_that_holds_no_rerun_arm_is_refused(
        trees: list[Path], rerun_trees: list[Path], tmp_path: Path, content: str) -> None:
    """alden-ec2221c7, review of #66: its contents were attributed to no arm and so reported by nobody."""
    third = _third_root(tmp_path, content)
    with Interrupted(trees), pytest.raises(fa.InputError, match="rerun root .*third holds none of the arms"):
        fa.status(trees, COMMIT, "full", rerun_roots=[*rerun_trees, third], rerun_commit=RERUN)


@pytest.mark.parametrize("content", ["known-arm run", "unrecognised run"])
def test_an_original_root_that_holds_no_arm_is_refused(trees: list[Path], tmp_path: Path, content: str) -> None:
    third = _third_root(tmp_path, content)
    # Known-arm data puts that arm in two roots, which the one-home check already refuses; anything else is armless.
    message = "arm A-up .* is in 2 of the roots" if content == "known-arm run" else "root .*third holds none of the arms"
    with pytest.raises(fa.InputError, match=message):
        fa.status([*trees, third], COMMIT, "full")


@pytest.mark.usefixtures("rerun_ok")
def test_with_the_rerun_the_analysis_measures_80_runs_none_from_an_interrupted_cell(
        trees: list[Path], rerun_trees: list[Path], monkeypatch: pytest.MonkeyPatch) -> None:
    measured: list[Path] = []
    measure = fa.measure
    monkeypatch.setattr(fa, "measure", lambda run: (measured.append(run.path), measure(run))[1])
    with Interrupted(trees):
        doc = fa.run(trees, COMMIT, rerun_roots=rerun_trees, rerun_commit=RERUN)
    assert len(measured) == 80
    rerun_base = {Path(r).resolve() for r in rerun_trees}
    from_rerun = [p for p in measured if any(b in p.resolve().parents for b in rerun_base)]
    assert len(from_rerun) == 16
    assert all(p.resolve().parent.parent.name in {"apt", "bpg"} and p.resolve().parent.name == "e150" for p in from_rerun)
    assert doc["collection"]["state"] == "complete and verified"
    assert doc["rerun"] == {"amendment": "acquisition amendment 1", "commit": RERUN, "cells": ["A-port 150 MeV", "B-pgcc 150 MeV"]}


def test_rerun_commit_problems_reads_the_repository(tmp_path: Path) -> None:
    """A real repository: descent, the files changed, the commit form."""
    def git(*args: str) -> str:
        return subprocess.run(["git", "-C", str(tmp_path), *args], capture_output=True, text=True, check=True).stdout.strip()

    git("init", "-q"); git("config", "user.email", "t@t"); git("config", "user.name", "t")
    for rel in (*fa.RERUN_MAY_CHANGE, "validation/field_followup_runs.py"):
        (tmp_path / rel).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / rel).write_text("v1\n")
    git("add", "-A"); git("commit", "-q", "-m", "acq"); acq = git("rev-parse", "HEAD")
    (tmp_path / fa.RERUN_MAY_CHANGE[0]).write_text("v2\n")
    git("commit", "-q", "-am", "request"); good = git("rev-parse", "HEAD")
    assert fa.rerun_commit_problems(acq, good, tmp_path) == []
    (tmp_path / "validation/field_followup_runs.py").write_text("v2\n")
    git("commit", "-q", "-am", "seeds"); bad = git("rev-parse", "HEAD")
    assert "other than the run requests" in fa.rerun_commit_problems(acq, bad, tmp_path)[0]
    assert "does not descend" in fa.rerun_commit_problems(good, acq, tmp_path)[0]
    assert "is the acquisition commit" in fa.rerun_commit_problems(acq, acq, tmp_path)[0]
    assert "full 40-hex" in fa.rerun_commit_problems(acq, good[:12], tmp_path)[0]


def test_a_rerun_commit_with_problems_is_refused_before_any_run_is_read(trees: list[Path], rerun_trees: list[Path],
                                                                         monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(fa, "rerun_commit_problems", lambda *_a: ["the rerun commit changes files other than the run requests"])
    monkeypatch.setattr(fa, "verify_run", lambda *_a: pytest.fail("a run was read although the rerun commit is refused"))
    with pytest.raises(fa.InputError, match="other than the run requests"):
        fa.status(trees, COMMIT, "full", rerun_roots=rerun_trees, rerun_commit=RERUN)
