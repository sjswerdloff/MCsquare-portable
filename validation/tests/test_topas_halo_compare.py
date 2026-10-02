"""Tests for topas_halo_compare on SYNTHETIC inputs only (no run output from any machine is read).

The TOPAS run directories are made by the REAL make_run.sh and run_topas.sh with the fake TOPAS of the wrapper tests
and a shrunk grid, so the provenance and run.txt this analysis parses are in the format the wrappers write. The fake
is followed by a step that writes dose.binheader as OpenTOPAS 4.3.0 writes it (REAL_HEADER is one of its own), so
run_topas.sh records that header's hash itself.

Run: uv run --no-project --with numpy --with scipy --with pytest pytest validation/tests/test_topas_halo_compare.py
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import random
import re
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
import pencil_endpoints as pe
import topas_halo_compare as thc

# isort: split
from apples_fixtures import (  # noqa: F401 - the autouse fixture must be in this module's namespace
    ALL_ARMS,
    _frozen_from_the_working_tree,
    attest,
    build,
)

pytestmark = pytest.mark.skipif(sys.platform != "darwin", reason="run_topas.sh uses macOS stat and shasum paths")

TOPAS_DIR = Path(__file__).resolve().parents[1] / "topas"
PEAK = {100: 75, 150: 155, 200: 204}  # depth bin of the synthetic Bragg peak; each energy's slabs lie short of it
# The shrunk phantom: 1 mm voxels, X and Y unequal so that an exchange of the two is a change.
BINS = {"X": 24, "Y": 26, "Z": 210}
HALF_CM = {"X": "1.2", "Y": "1.3", "Z": "10.5"}
HEADER = (
    "# TOPAS Version: 4.3\n"
    "# Parameter File: run.txt\n"
    "# Results for scorer: Dose\n"
    '# Filtered by: OnlyIncludeIfParticleOrAncestorNotNamed = 2 "neutron" "gamma"\n'
    "# Scored in component: Phantom\n"
    "# X in {X} bins of 0.1 cm\n"
    "# Y in {Y} bins of 0.1 cm\n"
    "# Z in {Z} bins of 0.1 cm\n"
    "# DoseToMedium ( Gy ) : Sum   \n"
    "# Binary file: dose.bin\n"
)
# dose.binheader of a run OpenTOPAS 4.3.0 made from the committed base (planning run E100 seed 900011, 2026-09-30),
# and the size and hash run_topas.sh recorded for it.
REAL_HEADER = HEADER.format(X=400, Y=400, Z=350)
REAL_HEADER_RECORD = (315, "e0acf0430f0b4c08fecc743a286dc8b92415a2877d433cbe6e42d36165aed334")


# ------------------------------------------------------------------------------------------------ frozen design
def test_seed_blocks_are_the_frozen_ones_and_outside_tds_confirmatory_blocks() -> None:
    assert thc.TOPAS_SEEDS == {
        100: tuple(range(912001, 912009)), 150: tuple(range(912011, 912019)), 200: tuple(range(912021, 912029))
    }
    seeds = [s for e in thc.ENERGIES for s in thc.TOPAS_SEEDS[e]]
    assert len(seeds) == len(set(seeds)) == 24
    assert all(912000 < s < 913000 for s in seeds)  # TD: 910001+ opt0 and 911001+ opt4 are the confirmatory arms


def test_design_constants() -> None:
    assert thc.PORTABLE_ARMS == ("A-port", "B-pgcc", "B-picc")
    assert (thc.TOPAS_HISTORIES, thc.TOPAS_EM, thc.LEVELS) == (10_000_000, "opt0", (0.90, 0.95))
    assert aa.SLAB_DEPTHS == {100: (40, 60), 150: (80, 125), 200: (100, 200)}
    assert aa.RINGS == ((5, 10), (10, 20), (20, 40), (40, 80), (80, 200))


# --------------------------------------------------------------------------------------------------- estimates
def test_difference_row_hand_checked() -> None:
    # MCsquare 1, 2, 3 (mean 2, variance 1); TOPAS 1, 1, 1. d = 1, SE = sqrt(1/3), Welch df = 2.
    row = thc.difference_row([1.0, 2.0, 3.0], [1.0, 1.0, 1.0])
    assert row["status"] == "difference"
    assert row["estimate"] == pytest.approx(1.0)
    assert row["se"] == pytest.approx(math.sqrt(1 / 3))
    assert row["df"] == pytest.approx(2.0)
    assert row["ci95"] == pytest.approx([1 - 4.302653 * math.sqrt(1 / 3), 1 + 4.302653 * math.sqrt(1 / 3)], abs=1e-5)
    assert row["ci90"] == pytest.approx([1 - 2.919986 * math.sqrt(1 / 3), 1 + 2.919986 * math.sqrt(1 / 3)], abs=1e-5)
    assert row["uncertainty"] == "welch"  # zero sample variance in ONE code is an ordinary Welch row


def test_difference_is_mcsquare_minus_topas() -> None:
    assert thc.difference_row([5.0, 5.2], [4.0, 4.2])["estimate"] == pytest.approx(1.0)


def test_ring_ratio_is_the_geometric_mean_ratio_not_the_ratio_of_arithmetic_means() -> None:
    # geometric means 0.02 and 0.02: ratio 1. Arithmetic means 0.025 and 0.02 would give 1.25.
    row = thc.ring_row([0.01, 0.04], [0.02, 0.02])
    assert row["status"] == "ratio"
    assert row["estimate"] == pytest.approx(1.0)
    se = math.log(4) / 2  # sd of the logs is ln(4)/sqrt(2); /sqrt(n=2); TOPAS has zero variance, so df = 1
    assert row["se_log"] == pytest.approx(se)
    assert row["df"] == pytest.approx(1.0)
    half = float(student_t.ppf(0.975, 1)) * se
    assert row["ci95"] == pytest.approx([math.exp(-half), math.exp(half)])


def test_ring_ratio_is_mcsquare_over_topas() -> None:
    assert thc.ring_row([0.02, 0.08], [0.01, 0.04])["estimate"] == pytest.approx(2.0)


@pytest.mark.parametrize(("mc", "tp"), [([0.0, 0.0], [0.01, 0.02]), ([0.01, 0.02], [0.0, 0.03]), ([0.0, 0.0], [0.0, 0.0])])
def test_a_zero_run_in_either_code_gives_counts_and_means_and_no_ratio(mc: list[float], tp: list[float]) -> None:
    row = thc.ring_row(mc, tp)
    assert row["status"] == "no ratio (zero runs)"
    assert "estimate" not in row and "ci95" not in row
    assert row["mcsquare"] == {"n_nonzero": sum(v > 0 for v in mc), "n": 2, "mean_fraction": st.mean(mc)}
    assert row["topas"] == {"n_nonzero": sum(v > 0 for v in tp), "n": 2, "mean_fraction": st.mean(tp)}


def test_constant_values_in_both_codes_give_the_difference_and_no_interval() -> None:
    """Amendment 2 (#54 review 7114): zero sample variance in both codes is not zero variance."""
    row = thc.difference_row([1.0] * 8, [1.5] * 8)
    assert row == {"status": "difference", "estimate": -0.5, "n_mcsquare": 8, "n_topas": 8,
                   "uncertainty": "not estimated: zero sample variance in both codes"}


def test_constant_positive_values_in_both_codes_give_the_ratio_and_no_interval() -> None:
    row = thc.ring_row([0.02] * 8, [0.01] * 8)
    assert row["estimate"] == pytest.approx(2.0)
    assert {k: v for k, v in row.items() if k != "estimate"} == {
        "status": "ratio", "n_mcsquare": 8, "n_topas": 8, "uncertainty": "not estimated: zero sample variance in both codes"}


def test_the_ordinary_and_the_zero_run_rows_are_unchanged_beside_the_constant_ones() -> None:
    ordinary = thc.ring_row([0.02, 0.021] * 4, [0.01] * 8)  # variable in one code: Welch, with intervals
    assert ordinary["uncertainty"] == "welch" and ordinary["ci90"][0] < ordinary["estimate"] < ordinary["ci90"][1]
    assert ordinary["ci95"][0] < ordinary["ci90"][0] and ordinary["se_log"] > 0 and ordinary["df"] == pytest.approx(7.0)
    zero = thc.ring_row([0.0] * 8, [0.0] * 8)  # constant AND zero: the zero-run rule, not a ratio
    assert zero["status"] == "no ratio (zero runs)" and "uncertainty" not in zero and "estimate" not in zero


def test_a_row_without_an_interval_says_so_in_the_table() -> None:
    rows = [thc._row(100, "A-port", "R80", thc.difference_row([1.0] * 8, [1.5] * 8), None),
            thc._row(100, "A-port", "ring_40_5_10", thc.ring_row([0.02] * 8, [0.01] * 8), None),
            thc._row(100, "A-port", "sigma_40", thc.difference_row([1.0, 2.0, 3.0], [1.0, 1.0, 1.0]), None)]
    md = thc.markdown({"label": "label", "rows": rows})
    assert "| R80 (difference) | A-port | -0.5000 | not estimated: zero sample variance in both codes | not estimated |  |" in md
    assert "| ring_40_5_10 (ratio) | A-port | 2.000 | not estimated: zero sample variance in both codes | not estimated |  |" in md
    assert sum("[" in line for line in md if line.startswith(("| R80", "| ring_"))) == 0
    assert next(line for line in md if line.startswith("| sigma_40")).count("[") == 2


def test_an_invalid_run_keeps_the_endpoint_visible_as_not_computed() -> None:
    for row in (thc.difference_row([1.0, None, 3.0], [1.0, 2.0]), thc.ring_row([0.1, 0.2], [None, 0.2])):
        assert row["status"] == "not computed"
        assert "estimate" not in row
    assert thc.difference_row([1.0, None, 3.0], [1.0, 2.0])["reason"] == "invalid in 1 of 3 MCsquare runs and 0 of 2 TOPAS runs"


def test_value_validity_rules() -> None:
    assert thc.r80_value({"R80": 77.4, "R80_multiple_crossings": False}) == 77.4
    assert thc.r80_value({"R80": 77.4, "R80_multiple_crossings": True}) is None
    assert thc.r80_value({"R80": 77.4}) is None  # the flag must be present and exactly False
    assert thc.sigma_value({"sigma_40": 3.1}, 40) == 3.1
    assert thc.sigma_value({"sigma_40": 0.5}, 40) is None
    assert thc.sigma_value({"sigma_40": None}, 40) is None
    assert thc.ring_value({"ring_40_80_200": 0.0}, 40, 80, 200) == 0.0  # zero is a value, not an invalid run
    assert thc.ring_value({"ring_40_80_200": float("nan")}, 40, 80, 200) is None
    assert thc.ring_value({"ring_40_80_200": 1.5}, 40, 80, 200) is None
    assert thc.ring_value({"ring_40_80_200": True}, 40, 80, 200) is None


def test_parse_provenance_two_verdicts_is_not_complete() -> None:
    assert thc.parse_provenance("host: h\nverdict: COMPLETE\n")["verdict"] == "COMPLETE"
    assert thc.parse_provenance("verdict: FAILED topas exited 3\nverdict: COMPLETE\n")["verdict"] != "COMPLETE"


# ------------------------------------------------------------------- TOPAS run directories from the real wrappers
class Wrappers:
    """The real make_run.sh and run_topas.sh, a shrunk base and the fake TOPAS."""

    def __init__(self, home: Path) -> None:
        self.bin = home / "bin"
        self.bin.mkdir(parents=True)
        for name in ("make_run.sh", "run_topas.sh"):
            shutil.copy2(TOPAS_DIR / name, self.bin / name)
        base, n = re.subn(
            r"(?m)^(i:Ge/Phantom/([XYZ])Bins) = .*$", lambda m: f"{m[1]} = {BINS[m[2]]}",
            (TOPAS_DIR / "stage1_base.txt").read_text(encoding="utf-8"),
        )
        assert n == 3  # the shrink must have applied to all three axes
        base, n = re.subn(r"(?m)^(d:Ge/Phantom/HL([XYZ])) = .*$", lambda m: f"{m[1]} = {HALF_CM[m[2]]} cm", base)
        assert n == 3
        self.base = self.bin / "stage1_base.txt"
        self.base.write_text(base, encoding="utf-8")
        header = home / "dose.binheader"
        header.write_text(HEADER.format(**BINS), encoding="utf-8")
        fake = home / "fake_topas_then_header.sh"
        fake.write_text(
            "#!/bin/bash\n"
            f"'{TOPAS_DIR / 'tests' / 'fake_topas.sh'}' \"$@\" || exit $?\n"
            f"[ -e dose.binheader ] && /bin/cat '{header}' > dose.binheader\n"
            "exit 0\n",
            encoding="utf-8",
        )
        fake.chmod(0o755)
        self.env = {**os.environ, "TOPAS_BIN": str(fake), "MIN_FREE_GB": "0"}

    def run(self, out: Path, energy: int, seed: int, *, em: str = "opt0", histories: int = thc.TOPAS_HISTORIES,
            threads: int = 12, attempt: str = "a1", mode: str = "ok") -> Path:
        made = subprocess.run(
            [str(self.bin / "make_run.sh"), str(energy), em, str(seed), str(histories), str(threads), str(out), attempt],
            env=self.env, capture_output=True, text=True, check=True,
        )
        run_dir = Path(made.stdout.strip()).parent
        done = subprocess.run([str(self.bin / "run_topas.sh"), str(run_dir)], env={**self.env, "FAKE_MODE": mode},
                              capture_output=True, text=True, check=False)
        assert (done.returncode == 0) == (mode == "ok"), done.stderr
        return run_dir


@pytest.fixture(scope="module")
def wrappers(tmp_path_factory: pytest.TempPathFactory) -> Wrappers:
    return Wrappers(tmp_path_factory.mktemp("wrappers"))


@pytest.fixture(scope="module")
def complete_set(wrappers: Wrappers, tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("complete") / "halo_addendum"
    for energy in thc.ENERGIES:
        for i, seed in enumerate(thc.TOPAS_SEEDS[energy]):
            wrappers.run(out, energy, seed, threads=12 if i < 4 else 16)
    return out


@pytest.fixture
def topas(complete_set: Path, wrappers: Wrappers, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A private copy of the complete set of 24, with the analysis pointed at the shrunk base the runs used."""
    monkeypatch.setattr(thc, "BASE", wrappers.base)
    root = tmp_path / "halo_addendum"
    shutil.copytree(complete_set, root)
    return root


def run_dir(root: Path, energy: int, seed: int) -> Path:
    found = sorted(root.glob(f"E{energy}_opt0_seed{seed}_n{thc.TOPAS_HISTORIES}_th*_a1"))
    assert len(found) == 1
    return found[0]


def refusal(root: Path) -> str:
    with pytest.raises(thc.InputError) as e:
        thc.discover(root)
    return str(e.value)


def rerecord(run: Path, name: str) -> None:
    """Make provenance.txt record `name` as it is now, the way run_topas.sh records it: a changed file that still
    agrees with its own record, so that only a check of CONFORMITY can refuse it."""
    f = run / name
    digest = hashlib.sha256(f.read_bytes()).hexdigest()
    line = f"{name} sha256: {digest}" if name.endswith(".txt") else f"{name} bytes: {f.stat().st_size} sha256: {digest}"
    prov = run / "provenance.txt"
    text, n = re.subn(rf"(?m)^{re.escape(line.split(': ')[0])}: .*$", line, prov.read_text())
    assert n == 1
    prov.write_text(text)


def edit(path: Path, old: str, new: str) -> None:
    text = path.read_text()
    assert text.count(old) == 1
    path.write_text(text.replace(old, new))


def test_the_wrappers_format_is_the_one_parsed(topas: Path) -> None:
    prov = thc.parse_provenance((run_dir(topas, 100, 912001) / "provenance.txt").read_text())
    assert prov["verdict"] == "COMPLETE"
    assert thc._DOSE_LINE.match(prov["dose.bin bytes"])
    assert re.fullmatch(r"[0-9a-f]{64}", prov["run.txt sha256"])
    assert re.fullmatch(r"[0-9a-f]{64}", prov["topas_bin sha256"])


def test_the_complete_set_is_accepted_with_thread_counts_recorded(topas: Path) -> None:
    runs, notes = thc.discover(topas)
    assert sorted(runs) == [(e, s) for e in thc.ENERGIES for s in thc.TOPAS_SEEDS[e]]
    assert sorted({r.threads for r in runs.values()}) == [12, 16]
    assert notes == {"other_attempts": [], "ignored": []}


def test_a_missing_run_refuses_and_names_it(topas: Path) -> None:
    shutil.rmtree(run_dir(topas, 150, 912013))
    msg = refusal(topas)
    assert "E150 seed 912013: 0 COMPLETE run directories" in msg
    assert "no dose file was opened" in msg


def test_refusal_happens_before_any_dose_or_mcsquare_file_is_read(topas: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    shutil.rmtree(run_dir(topas, 200, 912028))

    def forbidden(*_a: object, **_k: object) -> None:
        raise AssertionError("an input was opened although the TOPAS set is incomplete")

    monkeypatch.setattr(pe, "read_topas_bin", forbidden)
    monkeypatch.setattr(thc, "load_mcsquare", forbidden)
    with pytest.raises(thc.InputError):
        thc.run(topas, topas / "no-such-mcsquare-root")


def test_a_run_that_did_not_complete_refuses_and_shows_its_verdict(topas: Path, wrappers: Wrappers) -> None:
    shutil.rmtree(run_dir(topas, 100, 912004))
    wrappers.run(topas, 100, 912004, attempt="a2", mode="exit3")
    msg = refusal(topas)
    assert "E100 seed 912004: 0 COMPLETE run directories" in msg
    assert "_a2: FAILED topas exited 3" in msg


def test_a_run_still_in_progress_is_not_complete(topas: Path) -> None:
    (run_dir(topas, 100, 912004) / "provenance.txt").write_text("host: h\nstart: now\n")
    assert "E100 seed 912004: 0 COMPLETE run directories" in refusal(topas)


def test_a_failed_attempt_beside_a_complete_one_is_noted_not_fatal(topas: Path, wrappers: Wrappers) -> None:
    wrappers.run(topas, 100, 912004, attempt="a2", mode="exit3")
    runs, notes = thc.discover(topas)
    assert len(runs) == 24
    assert len(notes["other_attempts"]) == 1 and "FAILED topas exited 3" in notes["other_attempts"][0]


def test_two_complete_attempts_of_one_seed_refuse(topas: Path, wrappers: Wrappers) -> None:
    wrappers.run(topas, 150, 912011, attempt="a2")
    assert "E150 seed 912011: 2 COMPLETE run directories" in refusal(topas)


@pytest.mark.parametrize("name", [
    "E100_opt0_seed912009_n10000000_th12_a1",  # a seed outside the frozen block
    "E100_opt4_seed912001_n10000000_th12_a1",  # the other EM option
    "E100_opt0_seed912001_n1000000_th12_a2",  # other histories
    "E150_opt0_seed912001_n10000000_th12_a2",  # a frozen seed at another energy
    "E70_opt0_seed912001_n10000000_th12_a1",  # an energy outside the design
])
def test_a_run_directory_outside_the_frozen_list_refuses(topas: Path, name: str) -> None:
    (topas / name).mkdir()
    assert f"{name}: not a run of the frozen list" in refusal(topas)


def test_entries_that_are_not_run_directories_are_listed_as_ignored(topas: Path) -> None:
    (topas / "launch.log").write_text("x")
    (topas / "notes").mkdir()
    _runs, notes = thc.discover(topas)
    assert notes["ignored"] == ["launch.log", "notes"]


def test_run_txt_changed_after_the_run_refuses(topas: Path) -> None:
    with (run_dir(topas, 200, 912021) / "run.txt").open("a") as f:
        f.write("# edited\n")
    assert "run.txt is not the file recorded in provenance.txt" in refusal(topas)


def test_a_base_other_than_this_commits_refuses(topas: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(thc, "BASE", TOPAS_DIR / "stage1_base.txt")  # the runs used the shrunk base
    msg = refusal(topas)
    assert msg.count("stage1_base.txt differs from this commit's validation/topas/stage1_base.txt") == 24


def test_run_txt_must_request_what_the_directory_name_says(topas: Path) -> None:
    victim = run_dir(topas, 100, 912002)
    name = victim.name
    shutil.rmtree(victim)
    run_dir(topas, 100, 912001).rename(topas / name)  # seed 912001's files under seed 912002's name
    msg = refusal(topas)
    assert (f"{name}: run.txt is not what make_run.sh writes for this directory: line 4 is 'i:Ts/Seed = 912001', "
            "expected 'i:Ts/Seed = 912002'") in msg
    assert "E100 seed 912001: 0 COMPLETE run directories" in msg


def test_missing_dose_file_refuses(topas: Path) -> None:
    (run_dir(topas, 100, 912001) / "dose.bin").unlink()
    assert "dose.bin is missing" in refusal(topas)


def test_runs_from_two_topas_executables_refuse(topas: Path) -> None:
    prov = run_dir(topas, 100, 912001) / "provenance.txt"
    text, n = re.subn(r"(?m)^topas_bin sha256: .*$", "topas_bin sha256: " + "ab" * 32, prov.read_text())
    assert n == 1
    prov.write_text(text)
    assert "2 different TOPAS executables" in refusal(topas)


# ------------------------------------------------------------ amendment 1 (review 7101): conformity, not only identity
def test_expected_run_txt_is_what_make_run_sh_writes(topas: Path) -> None:
    for (energy, seed), r in thc.discover(topas)[0].items():
        assert (r.path / "run.txt").read_text() == thc.expected_run_txt(energy, seed, r.threads)
    assert '6 "g4em-standard_opt0" "g4h-phy_QGSP_BIC_HP" "g4decay" "g4ion-binarycascade" "g4h-elastic_HP" "g4stopping"' in (
        thc.expected_run_txt(100, 912001, 12))


def test_a_set_with_no_topas_executable_hash_refuses(topas: Path) -> None:
    for prov in topas.glob("E*/provenance.txt"):
        text, n = re.subn(r"(?m)^topas_bin sha256: .*\n", "", prov.read_text())
        assert n == 1
        prov.write_text(text)
    assert refusal(topas).count("provenance.txt has no well-formed 'topas_bin sha256' line") == 24


@pytest.mark.parametrize("value", ["", "absent", "ab" * 31, "AB" * 32, "ab" * 32 + " x"])
def test_a_malformed_topas_executable_hash_refuses(topas: Path, value: str) -> None:
    prov = run_dir(topas, 150, 912015) / "provenance.txt"
    text, n = re.subn(r"(?m)^topas_bin sha256: .*$", f"topas_bin sha256: {value}", prov.read_text())
    assert n == 1
    prov.write_text(text)
    msg = refusal(topas)
    assert "th16_a1: provenance.txt has no well-formed 'topas_bin sha256' line" in msg
    assert msg.count("no well-formed") == 1


@pytest.mark.parametrize("key", ["topas_bin sha256", "run.txt sha256", "dose.bin bytes", "dose.binheader bytes", "host"])
def test_a_provenance_key_written_twice_refuses_even_with_equal_values(topas: Path, key: str) -> None:
    prov = run_dir(topas, 100, 912003) / "provenance.txt"
    line = next(ln for ln in prov.read_text().splitlines() if ln.startswith(key + ": "))
    with prov.open("a") as f:
        f.write(line + "\n")
    assert f"provenance.txt has more than one {key!r} line" in refusal(topas)


def test_a_run_not_written_by_this_commits_runner_refuses(topas: Path) -> None:
    edit(run_dir(topas, 100, 912001) / "provenance.txt", f"runner sha256: {thc.sha256_file(thc.RUNNER)}",
         "runner sha256: " + "cd" * 32)
    msg = refusal(topas)
    assert msg.count("provenance.txt was not written by this commit's validation/topas/run_topas.sh") == 1


@pytest.mark.parametrize(("old", "new", "line"), [
    # review 7101: the EM option alone instead of the six frozen modules
    ('6 "g4em-standard_opt0" "g4h-phy_QGSP_BIC_HP" "g4decay" "g4ion-binarycascade" "g4h-elastic_HP" "g4stopping"',
     '1 "g4em-standard_opt0"', 3),
    ('"g4h-phy_QGSP_BIC_HP"', '"g4h-phy_QGSP_BERT_HP"', 3),
    ("includeFile = stage1_base.txt", "includeFile = other_base.txt", 1),
    ('s:Sc/Dose/OutputFile = "dose"', 's:Sc/Dose/OutputFile = "dose"\nd:Ph/Default/CutForAllParticles = 1 mm', 8),
    ('s:Sc/DoseAll/OutputFile = "dose_all"\n', 's:Sc/DoseAll/OutputFile = "dose_all"\ni:Ts/Seed = 5\n', None),
    ("i:Ts/Seed = 912021\n", "i:Ts/Seed = 912021\ni:Ts/Seed = 912021\n", 5),
    ('s:Sc/DoseAll/OutputFile = "dose_all"\n', 's:Sc/DoseAll/OutputFile = "dose_all"', None),
])
def test_a_run_txt_that_is_not_the_frozen_configuration_refuses_though_correctly_recorded(
    topas: Path, old: str, new: str, line: int | None
) -> None:
    run = run_dir(topas, 200, 912021)
    edit(run / "run.txt", old, new)
    rerecord(run, "run.txt")  # identity holds: the file is the one its provenance records
    msg = refusal(topas)
    assert "run.txt is not the file recorded" not in msg
    assert msg.count("run.txt is not what make_run.sh writes for this directory") == 1
    if line is not None:
        assert f"for this directory: line {line} is " in msg


def test_a_header_whose_recorded_byte_count_is_not_its_size_refuses(topas: Path) -> None:
    run = run_dir(topas, 100, 912001)
    n = (run / "dose.binheader").stat().st_size
    edit(run / "provenance.txt", f"dose.binheader bytes: {n} sha256:", f"dose.binheader bytes: {n + 1} sha256:")
    msg = refusal(topas)
    assert f"dose.binheader is {n} bytes, provenance.txt records {n + 1}" in msg
    assert "no dose file was opened" in msg


def test_a_header_changed_after_the_run_refuses(topas: Path) -> None:
    run = run_dir(topas, 100, 912001)
    edit(run / "dose.binheader", f"# X in {BINS['X']} bins", f"# X in {BINS['Y']} bins")
    edit(run / "dose.binheader", f"# Y in {BINS['Y']} bins", f"# Y in {BINS['X']} bins")
    msg = refusal(topas)
    assert "dose.binheader is not the file recorded in provenance.txt" in msg
    assert "no dose file was opened" in msg


@pytest.mark.parametrize(("old", "new", "what"), [
    ("# X in 24 bins of 0.1 cm\n# Y in 26 bins", "# X in 26 bins of 0.1 cm\n# Y in 24 bins", "line 6 is ('X', 26, 1.0)"),
    ("# Z in 210 bins of 0.1 cm", "# Z in 105 bins of 0.2 cm", "line 8 is ('Z', 105, 2.0)"),
    ("# Z in 210 bins of 0.1 cm", "# Z in 210 bins of 0.2 cm", "line 8 is ('Z', 210, 2.0)"),
    ("# X in 24 bins of 0.1 cm\n# Y in 26 bins of 0.1 cm\n", "# Y in 26 bins of 0.1 cm\n# X in 24 bins of 0.1 cm\n",
     "line 6 is"),
    ("# Results for scorer: Dose", "# Results for scorer: DoseAll", "line 3 is"),
    ('# Filtered by: OnlyIncludeIfParticleOrAncestorNotNamed = 2 "neutron" "gamma"\n', "", "has 9 lines"),
    ("# Scored in component: Phantom", "# Scored in component: World", "line 5 is"),
    ("( Gy ) : Sum", "( Gy ) : Sum Mean", "line 9 is"),
    ("# DoseToMedium ( Gy )", "# DoseToWater ( Gy )", "line 9 is"),
    ("# Binary file: dose.bin", "# Binary file: dose_all.bin", "line 10 is"),
    ("# TOPAS Version: 4.3", "# TOPAS Version: 4.2", "line 1 is"),
    ("# Binary file: dose.bin\n", "# Binary file: dose.bin\n# X in 24 bins of 0.1 cm\n", "has 11 lines"),
])
def test_a_header_that_is_not_what_the_base_requests_refuses_though_correctly_recorded(
    topas: Path, monkeypatch: pytest.MonkeyPatch, old: str, new: str, what: str
) -> None:
    run = run_dir(topas, 150, 912012)
    edit(run / "dose.binheader", old, new)
    rerecord(run, "dose.binheader")
    monkeypatch.setattr(pe, "read_topas_bin", lambda *_a: pytest.fail("a dose file was opened"))
    msg = refusal(topas)
    assert "dose.binheader is not the file recorded" not in msg
    assert f"th12_a1: dose.binheader {what}" in msg
    with pytest.raises(thc.InputError, match="dose.binheader"):  # and the reader of one run refuses it by itself
        thc.topas_record(thc.TopasRun(150, 912012, run, 12, thc.parse_provenance((run / "provenance.txt").read_text())))


def test_the_header_template_is_a_real_opentopas_header_and_matches_the_committed_base() -> None:
    assert len(REAL_HEADER.encode()) == REAL_HEADER_RECORD[0]
    assert hashlib.sha256(REAL_HEADER.encode()).hexdigest() == REAL_HEADER_RECORD[1]
    want = thc.expected_header(TOPAS_DIR / "stage1_base.txt")
    assert want[5:8] == [("X", 400, 1.0), ("Y", 400, 1.0), ("Z", 350, 1.0)]
    assert [w[0] for w in want if len(w) == 1] == [ln.rstrip() for ln in REAL_HEADER.splitlines() if " bins of " not in ln]
    assert thc.base_grid(TOPAS_DIR / "stage1_base.txt") == ([400, 400, 350], [1.0, 1.0, 1.0])


def test_an_axis_width_in_mm_is_the_same_statement_as_in_cm(topas: Path) -> None:
    run = run_dir(topas, 100, 912002)
    edit(run / "dose.binheader", "# Z in 210 bins of 0.1 cm", "# Z in 210 bins of 1 mm")
    assert thc.header_problems(run / "dose.binheader", thc.expected_header(thc.BASE)) == []


def test_a_dose_file_of_another_size_than_the_grid_refuses_though_correctly_recorded(topas: Path) -> None:
    run = run_dir(topas, 200, 912024)
    dose = run / "dose.bin"
    dose.write_bytes(dose.read_bytes()[:-8])
    rerecord(run, "dose.bin")
    n = 8 * BINS["X"] * BINS["Y"] * BINS["Z"]
    assert f"dose.bin is not the {n} bytes the base's grid implies (recorded {n - 8}, on disk {n - 8})" in refusal(topas)


def test_a_base_that_sets_a_parameter_twice_is_refused(tmp_path: Path) -> None:
    base = tmp_path / "base.txt"
    base.write_text((TOPAS_DIR / "stage1_base.txt").read_text() + "i:Ge/Phantom/XBins = 200\n")
    with pytest.raises(thc.InputError, match="Ge/Phantom/XBins is set more than once"):
        thc.base_grid(base)


@pytest.mark.parametrize(("name", "kind"), [
    ("E100_opt0_seed912001_n10000000_th12_a1.bak", "dir"),
    ("E100_opt0_seed0912001_n10000000_th12_a1", "dir"),  # a leading zero: make_run.sh refuses to write it
    ("E100_opt0_seed912001_n10000000_th12", "dir"),  # no attempt tag
    ("E100_opt0_seed912001_n10000000_th012_a2", "dir"),
    ("E100_opt0_seed912001_n10000000_th12_a2", "file"),  # a well-formed name that is not a directory
    ("e100-copy", "dir with run files"),
])
def test_a_run_like_entry_that_is_not_named_as_a_run_refuses(topas: Path, name: str, kind: str) -> None:
    if kind == "file":
        (topas / name).write_text("x")
    else:
        (topas / name).mkdir()
    if kind == "dir with run files":
        shutil.copy2(run_dir(topas, 100, 912001) / "provenance.txt", topas / name / "provenance.txt")
    assert f"{name}: looks like a run but is not a directory named as make_run.sh names one" in refusal(topas)


def test_status_reports_completeness_without_opening_a_dose_file(
    topas: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(pe, "read_topas_bin", lambda *_a: pytest.fail("--status opened a dose file"))
    assert thc.main(["--topas-root", str(topas), "--status"]) == 0
    assert "all 24 frozen TOPAS runs are present" in capsys.readouterr().out
    shutil.rmtree(run_dir(topas, 100, 912001))
    assert thc.main(["--topas-root", str(topas), "--status"]) == 2
    assert "REFUSED" in capsys.readouterr().err


# ------------------------------------------------------------------------------------------------ end to end
def synthetic_dose(energy: int, rng: random.Random) -> np.ndarray:
    """v[ix, iy, kz] as TOPAS writes it: a Gaussian spot (sigma about 3 mm) with a Bragg-like peak at PEAK[energy]."""
    cx, cy = (np.arange(BINS[a]) + 0.5 - BINS[a] / 2 for a in "XY")
    sigma = 3.0 * (1 + rng.gauss(0, 0.01))
    lateral = np.exp(-(cx[:, None] ** 2 + cy[None, :] ** 2) / (2 * sigma**2))
    k = np.arange(BINS["Z"])
    idd = 1 + 3 * np.exp(-(((k - PEAK[energy] - rng.gauss(0, 0.2)) / 5) ** 2))
    idd[PEAK[energy] + 6:] = 0.0  # nothing beyond the distal edge
    canonical = idd[:, None, None] * lateral[None, :, :]  # D[k, x, y]
    return np.ascontiguousarray(canonical.transpose(1, 2, 0)[:, :, ::-1])


def put_dose(run: Path, v: np.ndarray) -> None:
    """Replace a wrapper-made run's dose with `v` on the run's own grid, recorded as run_topas.sh records it.
    The header stays the one the run wrote."""
    assert v.shape == (BINS["X"], BINS["Y"], BINS["Z"])
    (run / "dose.bin").write_bytes(np.asfortranarray(v, dtype="<f8").tobytes(order="F"))
    rerecord(run, "dose.bin")


@pytest.fixture
def both(topas: Path, tmp_path: Path) -> tuple[Path, Path, dict[int, list[dict]]]:
    """(TOPAS root with synthetic doses, attested synthetic MCsquare tree, the TOPAS records computed directly)."""
    direct: dict[int, list[dict]] = {e: [] for e in thc.ENERGIES}
    for energy in thc.ENERGIES:
        for seed in thc.TOPAS_SEEDS[energy]:
            v = synthetic_dose(energy, random.Random(seed))
            put_dose(run_dir(topas, energy, seed), v)
            direct[energy].append(pe.endpoints(pe.canonical_from_topas(v), aa.SLAB_DEPTHS[energy]))

    def one_mcsquare_zero(arm: str, case: str, energy: int | None, i: int, rec: dict) -> None:
        if (arm, case, energy) == ("B-picc", "P", 150):
            rec["ring_125_80_200"] = 0.0 if i < 3 else rec["ring_125_80_200"]

    mc = build(tmp_path / "mc", ALL_ARMS, mutate=one_mcsquare_zero)
    attest(mc, ("A", "B"))
    return topas, mc, direct


def mc_records(mc: Path, arm: str, energy: int) -> list[dict]:
    seeds = sorted(s for s, e in aa.expected_population(aa.ACQUISITION_COMMIT).seeds[(arm, "P")].items() if e == energy)
    return [json.loads((mc / arm / "P" / f"s{s}" / "endpoints.json").read_text()[len("ENDPOINTS "):]) for s in seeds]


def row_of(doc: dict, energy: int, arm: str, endpoint: str) -> dict:
    found = [r for r in doc["rows"] if (r["energy"], r["arm"], r["endpoint"]) == (energy, arm, endpoint)]
    assert len(found) == 1
    return found[0]


def test_end_to_end_rows_match_an_independent_computation(both: tuple[Path, Path, dict[int, list[dict]]]) -> None:
    topas, mc, direct = both
    doc = thc.run(topas, mc, expect_fingerprint=None)
    assert len(doc["rows"]) == 3 * 3 * 13
    assert doc["label"] == thc.LABEL and "read before this design was fixed" in thc.LABEL
    assert len(doc["topas"]["runs"]) == 24
    for energy in thc.ENERGIES:
        d0 = aa.SLAB_DEPTHS[energy][0]
        tp = direct[energy]
        for arm in thc.PORTABLE_ARMS:
            recs = mc_records(mc, arm, energy)
            assert len(recs) == 8
            r80 = row_of(doc, energy, arm, "R80")
            assert r80["estimate"] == pytest.approx(st.mean(r["R80"] for r in recs) - st.mean(r["R80"] for r in tp))
            key = f"ring_{d0}_5_10"
            ratio = row_of(doc, energy, arm, key)
            want = math.exp(st.mean(math.log(r[key]) for r in recs) - st.mean(math.log(r[key]) for r in tp))
            assert ratio["status"] == "ratio"
            assert ratio["estimate"] == pytest.approx(want)
            assert ratio["ci90"][0] > ratio["ci95"][0] and ratio["ci90"][1] < ratio["ci95"][1]
            # the synthetic TOPAS grid is 24 x 26 mm, so it scores nothing beyond 20 mm
            far = row_of(doc, energy, arm, f"ring_{d0}_40_80")
            assert far["status"] == "no ratio (zero runs)"
            assert far["topas"]["n_nonzero"] == 0 and far["mcsquare"]["n_nonzero"] == 8
    zero = row_of(doc, 150, "B-picc", "ring_125_80_200")
    assert zero["mcsquare"]["n_nonzero"] == 5 and zero["mcsquare"]["n"] == 8


def test_td_bands_are_shown_at_200_mev_only(both: tuple[Path, Path, dict[int, list[dict]]]) -> None:
    topas, mc, _direct = both
    doc = thc.run(topas, mc, expect_fingerprint=None)
    banded = {r["endpoint"]: r["td_reference_band"] for r in doc["rows"] if r["td_reference_band"] is not None}
    assert all(r["energy"] == 200 for r in doc["rows"] if r["td_reference_band"] is not None)
    assert banded == {
        "R80": [-0.3, 0.3], "sigma_100": [-0.1, 0.1], "sigma_200": [-0.15, 0.15],
        "ring_100_20_40": [0.90, 1.10], "ring_100_40_80": [0.90, 1.10], "ring_100_80_200": [0.75, 1.25],
        "ring_200_20_40": [0.90, 1.10], "ring_200_40_80": [0.90, 1.10], "ring_200_80_200": [0.75, 1.25],
    }
    md = "\n".join(thc.markdown(doc))
    assert "Rings with a zero run in either code (no ratio):" in md
    assert "| R80 (difference) | A-port |" in md


def test_a_dose_file_changed_after_the_run_refuses(both: tuple[Path, Path, dict[int, list[dict]]]) -> None:
    topas, mc, _direct = both
    dose = run_dir(topas, 100, 912001) / "dose.bin"
    data = bytearray(dose.read_bytes())
    data[0] ^= 1
    dose.write_bytes(bytes(data))
    with pytest.raises(thc.InputError, match="dose.bin is not the file recorded in provenance.txt"):
        thc.run(topas, mc, expect_fingerprint=None)


def test_a_header_changed_between_the_gate_and_the_read_refuses_before_the_dose_is_interpreted(
    both: tuple[Path, Path, dict[int, list[dict]]], monkeypatch: pytest.MonkeyPatch
) -> None:
    topas, _mc, _direct = both
    runs, _notes = thc.discover(topas)
    victim = runs[(100, 912001)]
    edit(victim.path / "dose.binheader", f"# X in {BINS['X']} bins", f"# X in {BINS['Y']} bins")
    edit(victim.path / "dose.binheader", f"# Y in {BINS['Y']} bins", f"# Y in {BINS['X']} bins")
    monkeypatch.setattr(pe, "read_topas_bin", lambda *_a: pytest.fail("the dose was interpreted with a changed header"))
    with pytest.raises(thc.InputError, match="dose.binheader is not the file recorded in provenance.txt"):
        thc.topas_record(victim)


def test_main_refuses_any_mcsquare_dataset_but_the_reused_one(
    both: tuple[Path, Path, dict[int, list[dict]]], tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    topas, mc, _direct = both
    out = tmp_path / "out.json"
    assert thc.main(["--topas-root", str(topas), "--mcsquare-root", str(mc), "--json", str(out), "--md", str(tmp_path / "o.md")]) == 2
    assert f"is not the reused dataset {thc.MCSQUARE_FINGERPRINT}" in capsys.readouterr().err
    assert not out.exists()


def test_mcsquare_tree_without_a_manifest_refuses(tmp_path: Path) -> None:
    mc = build(tmp_path / "mc", ALL_ARMS)
    with pytest.raises(thc.InputError, match="no collection manifest"):
        thc.load_mcsquare(mc, expect_fingerprint=None)


def test_an_unusable_portable_run_refuses(tmp_path: Path) -> None:
    def failed(arm: str, case: str, energy: int | None, i: int, run: dict) -> None:
        if (arm, case, energy, i) == ("B-pgcc", "P", 100, 0):
            run["transport_status"] = "failed"

    mc = build(tmp_path / "mc", ALL_ARMS, run_mutate=failed)
    attest(mc, ("A", "B"))
    with pytest.raises(thc.InputError, match="B-pgcc at 100 MeV: 8 runs, 1 unusable"):
        thc.load_mcsquare(mc, expect_fingerprint=None)
