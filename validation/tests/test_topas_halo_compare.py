"""Tests for topas_halo_compare on SYNTHETIC inputs only (no run output from any machine is read).

The TOPAS run directories are made by the REAL make_run.sh and run_topas.sh with the fake TOPAS of the wrapper tests
and a shrunk grid, so the provenance and run.txt this analysis parses are in the format the wrappers write.

Run: uv run --no-project --with numpy --with scipy --with pytest pytest validation/tests/test_topas_halo_compare.py
"""

from __future__ import annotations

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
SHRUNK = {"X": 2, "Y": 3, "Z": 4}
PEAK = {100: 75, 150: 155, 200: 204}  # depth bin of the synthetic Bragg peak; each energy's slabs lie short of it
N_LATERAL = 24


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
            r"(?m)^(i:Ge/Phantom/([XYZ])Bins) = .*$", lambda m: f"{m[1]} = {SHRUNK[m[2]]}",
            (TOPAS_DIR / "stage1_base.txt").read_text(encoding="utf-8"),
        )
        assert n == 3  # the shrink must have applied to all three axes
        self.base = self.bin / "stage1_base.txt"
        self.base.write_text(base, encoding="utf-8")
        self.env = {**os.environ, "TOPAS_BIN": str(TOPAS_DIR / "tests" / "fake_topas.sh"), "MIN_FREE_GB": "0"}

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
    assert f"{name}: run.txt lacks the line 'i:Ts/Seed = 912002'" in msg
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
    c = np.arange(N_LATERAL) + 0.5 - N_LATERAL / 2
    sigma = 3.0 * (1 + rng.gauss(0, 0.01))
    lateral = np.exp(-(c[:, None] ** 2 + c[None, :] ** 2) / (2 * sigma**2))
    k = np.arange(PEAK[energy] + 6)
    idd = 1 + 3 * np.exp(-(((k - PEAK[energy] - rng.gauss(0, 0.2)) / 5) ** 2))
    canonical = idd[:, None, None] * lateral[None, :, :]  # D[k, x, y]
    return np.ascontiguousarray(canonical.transpose(1, 2, 0)[:, :, ::-1])


def put_dose(run: Path, v: np.ndarray) -> None:
    """Replace a wrapper-made run's dose with `v`, recorded in its provenance as run_topas.sh records it."""
    nx, ny, nz = v.shape
    payload = np.asfortranarray(v, dtype="<f8").tobytes(order="F")
    (run / "dose.bin").write_bytes(payload)
    (run / "dose.binheader").write_text(
        "# TOPAS Version: 4.3\n# Results for scorer: Dose\n"
        f"# X in {nx} bins of 0.1 cm\n# Y in {ny} bins of 0.1 cm\n# Z in {nz} bins of 0.1 cm\n"
        "# DoseToMedium ( Gy ) : Sum   \n# Binary file: dose.bin\n"
    )
    prov = run / "provenance.txt"
    text, n = re.subn(r"(?m)^dose\.bin bytes: .*$", f"dose.bin bytes: {len(payload)} sha256: {thc.sha256_file(run / 'dose.bin')}",
                      prov.read_text())
    assert n == 1
    prov.write_text(text)


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
            # the synthetic TOPAS grid is 24 mm wide, so it scores nothing beyond 20 mm
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
