"""Tests for apples_collect on SYNTHETIC acquisition-shaped trees (layout, file names and hash-list formats follow the two
workflow files; no run output from any machine is read).

Run: uv run --no-project --with numpy --with scipy --with pytest pytest validation/tests/test_apples_collect.py
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import apples_analyse as aa
import apples_collect as ac

# isort: split
from apples_fixtures import (  # noqa: F401 - the autouse fixture must be in this module's namespace
    _REAL_EXPECTED,
    ACQ,
    NATIVE_ARM,
    SOURCES,
    _frozen_from_the_working_tree,
    edit_run_json,
    interrupted,
    make_native,
    native_seed_dir,
    sha,
    stop_job_after,
    transport_failure,
)

NATIVE_P = "apt/p"
NATIVE_F = "apt/f"


@pytest.fixture(scope="module")
def native_template(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """One full acquisition-shaped tree, built once; each test copies it (the copy is the system under test)."""
    root = tmp_path_factory.mktemp("native")
    with pytest.MonkeyPatch.context() as mp:  # the autouse fixture is function-scoped: patch the same way here
        mp.setattr(aa, "expected_population", lambda commit=ACQ, repo=aa.REPO, sources=None: _REAL_EXPECTED(commit, repo, SOURCES))
        make_native(root)
    return root


@pytest.fixture
def native(native_template: Path, tmp_path: Path) -> dict[str, Path]:
    import shutil

    shutil.copytree(native_template, tmp_path / "n")
    return {"win": tmp_path / "n" / "win" / ACQ[:12], "lin": tmp_path / "n" / "lin" / ACQ[:12]}


def run_collect(native: dict[str, Path], out: Path, **kw: object) -> dict[str, object]:
    return ac.collect([native["win"], native["lin"]], out, **kw)  # type: ignore[arg-type]


def refuses(native: dict[str, Path], tmp_path: Path, *needles: str, **kw: object) -> str:
    out = tmp_path / "out"
    with pytest.raises(ac.CollectionError) as e:
        run_collect(native, out, **kw)
    assert not out.exists(), "a refused collection must write nothing"
    for needle in needles:
        assert needle in str(e.value), (needle, str(e.value))
    return str(e.value)


def tree_state(root: Path) -> dict[str, tuple[str, float]]:
    return {p.relative_to(root).as_posix(): (sha(p.read_bytes()), p.stat().st_mtime) for p in sorted(root.rglob("*")) if p.is_file()}


# ----------------------------------------------------------------------------------------- the happy path


def test_acquisition_shaped_trees_collect_and_the_analysis_reads_the_result(native: dict[str, Path], tmp_path: Path) -> None:
    out = tmp_path / "out"
    assert ac.main(["--source", str(native["win"]), "--source", str(native["lin"]), "--out", str(out)]) == 0
    assert aa.main([str(out)]) == 0
    _, doc = aa.run_analysis(out, ("A", "B"))
    assert doc["collection_manifest"]["present"] is True
    assert all(not c["partial"] for c in doc["contrasts"])
    assert all(c["claims"]["joint_claim_equivalent_on_all_endpoints"] is True for c in doc["contrasts"] if c["kind"] == "confirmatory")
    assert doc["runs"]["A-port/F"]["legacy_bindings_used"] == ["pe1-at-acquisition-2f9dab40"]


def test_output_layout_is_arm_case_seed_with_provenance_beside_the_runs(native: dict[str, Path], tmp_path: Path) -> None:
    out = tmp_path / "out"
    run_collect(native, out)
    assert sorted(p.name for p in out.iterdir()) == sorted([*NATIVE_ARM, aa.MANIFEST_NAME])
    seed_dir = out / "A-port" / "P" / "s961011"  # P's e150/ layer is flattened; the energy lives in run.json
    assert sorted(p.name for p in seed_dir.iterdir()) == ["config.txt", "endpoints.json", "log.txt", "out_sha256.txt", "run.json"]
    assert json.loads((seed_dir / "run.json").read_text())["energy_mev"] == 150
    prov = sorted(p.name for p in (out / "A-port" / "P" / aa.PROVENANCE_DIR).iterdir())
    assert prov == ["e100_inputs_sha256.txt", "e150_inputs_sha256.txt", "e200_inputs_sha256.txt", "run_root.txt",
                    "snapshot_binary_sha256.txt", "snapshot_build.txt", "snapshot_sha256.txt"]
    assert sorted(p.name for p in (out / "B-up" / "F").iterdir() if p.name.startswith("s")) == [f"s96203{n}" for n in range(1, 9)]
    assert "fcase_sha256.txt" in {p.name for p in (out / "B-up" / "F" / aa.PROVENANCE_DIR).iterdir()}


def test_copies_are_byte_identical_and_the_manifest_records_source_and_sha256(native: dict[str, Path], tmp_path: Path) -> None:
    out = tmp_path / "out"
    manifest = run_collect(native, out)
    on_disk = json.loads((out / aa.MANIFEST_NAME).read_text())
    assert on_disk["files"] == manifest["files"] and on_disk["acquisition_commit"] == ACQ
    assert len(on_disk["files"]) == 160 * 5 + 10 * 4 + 5 * 3 + 5 * 1  # runs; 4 per arm/case; 3 inputs lists per P; 1 fcase list per F
    for entry in on_disk["files"]:
        src, dest = Path(entry["source"]), out / entry["path"]
        assert src.read_bytes() == dest.read_bytes() and sha(dest.read_bytes()) == entry["sha256"]
        assert any(str(src).startswith(str(r)) for r in native.values())


def test_sources_are_untouched(native: dict[str, Path], tmp_path: Path) -> None:
    before = {k: tree_state(v) for k, v in native.items()}
    run_collect(native, tmp_path / "out")
    assert {k: tree_state(v) for k, v in native.items()} == before


def test_a_pycache_in_the_snapshot_is_tolerated_and_recorded(native: dict[str, Path], tmp_path: Path) -> None:
    manifest = run_collect(native, tmp_path / "out")
    assert any("__pycache__/x.pyc" in t for t in manifest["tolerated"])  # type: ignore[attr-defined]


def test_one_part_only(native: dict[str, Path], tmp_path: Path) -> None:
    import shutil

    shutil.rmtree(native["lin"])
    out = tmp_path / "out"
    ac.collect([native["win"]], out, ("A",))
    assert sorted(p.name for p in out.iterdir()) == ["A-port", "A-up", aa.MANIFEST_NAME]
    assert aa.main([str(out), "--parts", "A"]) == 0


# ------------------------------------------------------------------------------- refusals: layout and population


def test_existing_output_directory_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    (tmp_path / "out").mkdir()
    with pytest.raises(ac.CollectionError, match="must not exist"):
        run_collect(native, tmp_path / "out")


def test_run_root_txt_beside_the_cases_is_part_of_the_layout_not_an_error(native: dict[str, Path], tmp_path: Path) -> None:
    """The reviewer's trigger: the analysis loader rejected a native root holding run_root.txt. The collector reads it."""
    assert (native["win"] / NATIVE_P / "run_root.txt").is_file()
    run_collect(native, tmp_path / "out")


def test_missing_dose_file_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    (native["lin"] / "bpg" / "f" / "s963031" / "out_seed" / "Dose.raw").unlink()
    refuses(native, tmp_path, "Dose.raw")


def test_missing_expected_seed_directory_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    import shutil

    shutil.rmtree(native["win"] / NATIVE_P / "e200" / "s961021")
    refuses(native, tmp_path, "expected seed 961021 has no directory")


def test_missing_arm_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    import shutil

    shutil.rmtree(native["lin"] / "bpi")
    refuses(native, tmp_path, "arm B-picc is not present")


def test_smoke_directory_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    (native["win"] / "smoke").mkdir()
    refuses(native, tmp_path, "unexpected entry 'smoke'")


def test_stray_file_at_each_level_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    for where in (native["win"], native["win"] / "apt", native["win"] / NATIVE_P, native["win"] / NATIVE_P / "e100",
                  native["win"] / NATIVE_P / "e100" / "s961001", native["lin"] / "bup" / "f" / "s962031"):
        stray = where / "notes.txt"
        stray.write_text("x")
        refuses(native, tmp_path, "notes.txt")
        stray.unlink()


def test_unexpected_seed_directory_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    import shutil

    shutil.copytree(native["win"] / NATIVE_P / "e100" / "s961001", native["win"] / NATIVE_P / "e100" / "s999999")
    refuses(native, tmp_path, "s999999 is not in the frozen list")


def test_seed_directory_under_the_wrong_energy_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    import shutil

    shutil.move(str(native["win"] / NATIVE_P / "e100" / "s961001"), str(native["win"] / NATIVE_P / "e150" / "s961001"))
    refuses(native, tmp_path, "seed 961001 is frozen at 100 MeV but sits under e150")


def test_leading_zero_alias_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    (native["lin"] / "bup" / "f" / "s962031").rename(native["lin"] / "bup" / "f" / "s0962031")
    refuses(native, tmp_path, "s0962031", "962031")


def test_source_directory_must_be_named_for_the_commit(native: dict[str, Path], tmp_path: Path) -> None:
    moved = native["win"].rename(native["win"].parent / "123456789abc")
    refuses({"win": moved, "lin": native["lin"]}, tmp_path, "is named '123456789abc'")


def test_arm_in_two_sources_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    import shutil

    shutil.copytree(native["win"] / "apt", native["lin"] / "apt")
    refuses(native, tmp_path, "arm A-port is also present")


def test_arm_outside_the_requested_parts_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    refuses(native, tmp_path, "not in parts A", parts=("A",))


# ---------------------------------------------------------------------------- refusals: provenance and bindings


def test_run_root_must_name_the_arm_case_mode_and_commit(native: dict[str, Path], tmp_path: Path) -> None:
    f = native["win"] / NATIVE_P / "run_root.txt"
    original = f.read_text()
    for bad in (original.replace("mode full", "mode smoke"), original.replace(ACQ, "0" * 40),
                original.replace("A-port", "A-up"), "garbage\n"):
        f.write_text(bad)
        refuses(native, tmp_path, "run_root.txt")
    f.write_text(original)


def test_snapshot_file_changed_after_the_hash_list_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    (native["win"] / NATIVE_P / "snapshot" / "validation" / "pencil_endpoints.py").write_text("# changed\n")
    refuses(native, tmp_path, "pencil_endpoints.py", "hashes differently")


def test_snapshot_extra_file_outside_pycache_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    (native["lin"] / "bup" / "p" / "snapshot" / "extra.txt").write_text("x")
    refuses(native, tmp_path, "extra.txt", "not in the list")


def test_binary_hash_in_run_json_must_equal_the_snapshots(native: dict[str, Path], tmp_path: Path) -> None:
    f = native["win"] / "apt" / "p" / "e100" / "s961001" / "run.json"
    run = json.loads(f.read_text())
    run["binary_sha256"] = "0" * 64
    f.write_text(json.dumps(run))
    refuses(native, tmp_path, "binary_sha256 differs from snapshot/binary_sha256.txt")


def test_case_input_file_changed_after_its_hash_list_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    (native["win"] / NATIVE_P / "e200" / "plan_E200.txt").write_text("another plan")
    refuses(native, tmp_path, "plan_E200.txt")


def test_materials_copy_that_differs_from_the_snapshot_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    (native["lin"] / "bpg" / "f" / "s963031" / "Materials" / "Water.txt").write_text("changed")
    refuses(native, tmp_path, "Water.txt")


def test_fcase_file_changed_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    (native["lin"] / "bpi" / "f" / "fcase" / "cube.raw").write_text("another cube")
    refuses(native, tmp_path, "cube.raw")


def test_dose_file_that_does_not_match_its_sha256_txt_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    (native["win"] / NATIVE_P / "e100" / "s961001" / "out" / "Dose.raw").write_text("tampered")
    refuses(native, tmp_path, "Dose.raw", "hashes differently")


def test_f_record_hashes_must_match_the_files(native: dict[str, Path], tmp_path: Path) -> None:
    f = native["win"] / NATIVE_F / "s961031" / "record.json"
    record = json.loads(f.read_text())
    original = json.dumps(record)
    for edit in (lambda r: r["sha256"].update({"Dose.raw": "0" * 64}), lambda r: r.update(cfg_sha256="0" * 64),
                 lambda r: r["sha256"].update({"cube.raw": "0" * 64}), lambda r: r["sha256"].update({"binary": "0" * 64})):
        r = json.loads(original)
        edit(r)
        f.write_text(json.dumps(r))
        refuses(native, tmp_path, "record.json")
    f.write_text(original)


@pytest.mark.parametrize(
    ("key", "value"),
    [("arm", "A-up"), ("case", "F"), ("seed", 961002), ("energy_mev", 150)],
)
def test_run_json_must_agree_with_its_directory(native: dict[str, Path], tmp_path: Path, key: str, value: object) -> None:
    f = native["win"] / NATIVE_P / "e100" / "s961001" / "run.json"
    run = json.loads(f.read_text())
    run[key] = value
    f.write_text(json.dumps(run))
    refuses(native, tmp_path, "run.json", key)


def test_config_must_agree_with_run_json(native: dict[str, Path], tmp_path: Path) -> None:
    f = native["lin"] / "bup" / "p" / "e100" / "s962001" / "config.txt"
    for old, new, needle in (("Num_Threads 3", "Num_Threads 4", "Num_Threads"), ("RNG_Seed 962001", "RNG_Seed 962002", "RNG_Seed"),
                             ("Num_Primaries 1e7", "Num_Primaries 1e5", "Num_Primaries")):
        original = f.read_text()
        f.write_text(original.replace(old, new))
        refuses(native, tmp_path, needle)
        f.write_text(original)


def test_log_primaries_must_agree_with_run_json(native: dict[str, Path], tmp_path: Path) -> None:
    (native["lin"] / "bpi" / "p" / "e200" / "s964021" / "log.txt").write_text("Nbr primaries simulated: 5\n")
    refuses(native, tmp_path, "primaries line")


def test_malformed_hash_list_is_refused_not_crashed(native: dict[str, Path], tmp_path: Path) -> None:
    (native["win"] / NATIVE_P / "snapshot_sha256.txt").write_text("not a hash list\n")
    refuses(native, tmp_path, "unusable hash list")


def test_all_problems_are_reported_together(native: dict[str, Path], tmp_path: Path) -> None:
    (native["win"] / NATIVE_P / "e150" / "s961011" / "notes.txt").write_text("x")
    (native["lin"] / "bpg" / "f" / "s963031" / "out_seed" / "Dose.raw").unlink()
    (native["win"] / "smoke").mkdir()
    message = refuses(native, tmp_path, "problem(s)")
    assert "refusing, 3 problem(s)" not in message  # a missing product is reported more than once; all three appear
    assert "notes.txt" in message and "Dose.raw" in message and "smoke" in message


# ------------------------------------------------------------------------------------- the manifest and the CLI


def test_analysis_refuses_a_collected_file_changed_afterwards(native: dict[str, Path], tmp_path: Path) -> None:
    out = tmp_path / "out"
    run_collect(native, out)
    f = out / "B-pgcc" / "P" / "s963001" / "endpoints.json"
    f.write_text(f.read_text() + " ")
    assert aa.main([str(out)]) == 2


def test_collected_trees_with_a_modified_native_run_are_not_confirmatory(native: dict[str, Path], tmp_path: Path) -> None:
    """End to end through the seam: a run whose run.json says smoke is copied faithfully and the analysis refuses its claim."""
    f = native["win"] / NATIVE_P / "e100" / "s961001" / "run.json"
    run = json.loads(f.read_text())
    run.update(mode="smoke")
    f.write_text(json.dumps(run))
    out = tmp_path / "out"
    run_collect(native, out)
    _, doc = aa.run_analysis(out, ("A", "B"))
    contrast = next(c for c in doc["contrasts"] if c["name"] == "A-port vs A-up")
    assert contrast["claims_withheld"] is True and contrast["claims"] == {}


def test_cli_exit_codes(native: dict[str, Path], tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    out = tmp_path / "out"
    argv = ["--source", str(native["win"]), "--source", str(native["lin"]), "--out", str(out)]
    assert ac.main([*argv, "--parts", "Z"]) == 2
    assert ac.main(argv) == 0
    assert "collected" in capsys.readouterr().out
    assert ac.main(argv) == 2  # the output now exists
    assert "must not exist" in capsys.readouterr().err
    assert ac.main(["--source", str(tmp_path / "nowhere"), "--out", str(tmp_path / "o2")]) == 2


def test_hash_list_parsers_read_both_hosts_forms() -> None:
    h = hashlib.sha256(b"x").hexdigest()
    assert ac.parse_hash_list(f"{h}  .\\a\\b.txt\r\n{h}  ./c/d.txt\n{h}  e\n") == {"a/b.txt": h, "c/d.txt": h, "e": h}
    for bad in ("", "zz  a", f"{h} a b".replace(" a b", ""), f"{h}  a\n{h}  a\n"):
        with pytest.raises(ValueError):
            ac.parse_hash_list(bad)


# ------------------------------------------------------------------- failed and absent runs: collected, not refused


def ne_entries(manifest: dict[str, object]) -> dict[tuple[str, str, int], dict]:
    return {(e["arm"], e["case"], e["seed"]): e for e in manifest["not_established"]}  # type: ignore[attr-defined]


def contrast(doc: dict[str, object], name: str) -> dict:
    return next(c for c in doc["contrasts"] if c["name"] == name)  # type: ignore[attr-defined]


def test_transport_failure_collects_with_the_later_seeds_absent_and_the_contrast_partial(
    native: dict[str, Path], tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    transport_failure(native, "A-port", "P", 961003, "mcsquare_exit_139")
    absent = stop_job_after(native, "A-port", "P", 961003)
    assert len(absent) == 21  # 961004..961008 at 100 MeV, then the 150 and 200 MeV cells, never reached
    out = tmp_path / "out"
    argv = ["--source", str(native["win"]), "--source", str(native["lin"]), "--out", str(out)]
    assert ac.main(argv) == 0
    assert "collected 138, failed 1, absent 21" in capsys.readouterr().out
    manifest = json.loads((out / aa.MANIFEST_NAME).read_text())
    assert manifest["outcomes"] == {"collected": 138, "failed": 1, "absent": 21}
    ne = ne_entries(manifest)
    assert set(ne) == {("A-port", "P", s) for s in [961003, *absent]}
    failed = ne[("A-port", "P", 961003)]
    assert failed["outcome"] == "failed" and failed["energy"] == 100
    assert failed["reason"] == "transport failure: transport_status 'mcsquare_exit_139'"
    assert sorted(Path(f["path"]).name for f in failed["files"]) == ["config.txt", "log.txt", "run.json"]
    for f in failed["files"]:
        assert f["path"].startswith(".not_established/A-port/P/s961003/")
        assert sha((out / f["path"]).read_bytes()) == f["sha256"] == sha(Path(f["source"]).read_bytes())
    assert not (out / "A-port" / "P" / "s961003").exists()  # never where the analysis reads a run
    for s in absent:
        entry = ne[("A-port", "P", s)]
        assert entry["outcome"] == "absent" and entry["files"] == []
        assert entry["energy"] == aa.expected_population(ACQ).seeds[("A-port", "P")][s]
        assert entry["reason"].startswith("no directory: the job's last run directory is s961003 (failed: transport")
    assert "e150_inputs_sha256.txt" not in {p.name for p in (out / "A-port" / "P" / aa.PROVENANCE_DIR).iterdir()}

    _, doc = aa.run_analysis(out, ("A", "B"))
    c = contrast(doc, "A-port vs A-up")
    assert c["partial"] is True and c["claims_withheld"] is True and c["claims"] == {}
    assert {e["outcome"] for e in c["endpoints"]} <= {"descriptive", "not_established"}
    assert all(e["holm_decision"] == "" and e["p_tost"] is None for e in c["endpoints"])
    assert "A-port/P/s961003 (100 MeV): not established, failed: transport failure: transport_status 'mcsquare_exit_139'" in c[
        "partial_reasons"]
    assert all(not contrast(doc, n)["partial"] for n in ("B-pgcc vs B-up", "B-picc vs B-up"))
    listed = {(n["arm"], n["case"], n["seed"]): n for n in doc["not_established"]}
    assert set(listed) == set(ne) and listed[("A-port", "P", 961003)]["outcome"] == "failed"
    assert listed[("A-port", "P", 961028)]["energy_mev"] == 200
    md, _ = aa.run_analysis(out, ("A", "B"))
    text = "\n".join(md)
    assert "Runs not established (listed by the collector, never read as runs; their contrasts are PARTIAL): 22" in text
    assert "- A-port/P/s961003 (100 MeV): failed: transport failure: transport_status 'mcsquare_exit_139'" in text
    # the dataset fingerprint covers the .not_established copies: 138 runs x (run.json + endpoint record), the manifest
    # and the failed run's three diagnostic files
    assert doc["dataset_fingerprint"]["n_files"] == 138 * 2 + 1 + 3


@pytest.mark.parametrize("status", ["simulated_count_below_requested", "config_tag_not_recognised", "primaries_outside_geometry"])
def test_transport_failure_after_the_primaries_line_was_read_still_checks_the_log(
    native: dict[str, Path], tmp_path: Path, status: str
) -> None:
    d = transport_failure(native, "B-up", "P", 962028, status, simulated=9_999_000, keep_dose=True)
    (d / "log.txt").write_text("Nbr primaries simulated: 9999000\n")
    manifest = run_collect(native, tmp_path / "first")
    assert ne_entries(manifest)[("B-up", "P", 962028)]["reason"] == f"transport failure: transport_status {status!r}"
    assert any("not verified, no sha256.txt" in t for t in manifest["tolerated"])  # type: ignore[attr-defined]
    (d / "log.txt").write_text("Nbr primaries simulated: 5\n")  # the log must still give run.json simulated
    refuses(native, tmp_path, "primaries line")


def test_interrupted_first_run_of_a_job_is_failed_and_the_rest_absent(native: dict[str, Path], tmp_path: Path) -> None:
    interrupted(native, "A-port", "F", 961031)
    absent = stop_job_after(native, "A-port", "F", 961031)
    out = tmp_path / "out"
    manifest = run_collect(native, out)
    ne = ne_entries(manifest)
    assert ne[("A-port", "F", 961031)]["outcome"] == "failed"
    assert ne[("A-port", "F", 961031)]["reason"] == "interrupted: the run directory has no run.json"
    assert sorted(Path(f["path"]).name for f in ne[("A-port", "F", 961031)]["files"]) == ["cfg.txt", "log.txt"]
    assert all(ne[("A-port", "F", s)]["outcome"] == "absent" for s in absent) and len(absent) == 7
    assert sorted(p.name for p in (out / "A-port" / "F").iterdir()) == [aa.PROVENANCE_DIR]  # no run at all
    _, doc = aa.run_analysis(out, ("A",))
    c = contrast(doc, "A-port vs A-up")
    assert c["partial"] is True and c["claims"] == {}
    assert all(e["outcome"] == "not_established" for e in c["endpoints"] if e["endpoint"].startswith("F/"))


def test_job_that_stopped_before_its_first_run_lists_every_seed_absent(native: dict[str, Path], tmp_path: Path) -> None:
    stop_job_after(native, "B-picc", "P", None)
    manifest = run_collect(native, tmp_path / "out")
    ne = ne_entries(manifest)
    assert len(ne) == 24 and {e["outcome"] for e in ne.values()} == {"absent"}
    assert {e["reason"] for e in ne.values()} == {"no directory: the job stopped before its first run"}


@pytest.mark.parametrize(("arm", "seed", "content", "why"), [
    ("B-up", 962028, "", "endpoints.json is empty"),  # Linux: the redirect created it, pipefail killed the job
    ("B-up", 962028, "Traceback (most recent call last):\n", "endpoints.json does not hold exactly one ENDPOINTS line"),
    ("A-up", 960028, None, "endpoints.json is missing"),  # Windows: $ep is written only when the script exited 0
])
def test_p_endpoint_failure_after_transport_ok_is_failed(
    native: dict[str, Path], tmp_path: Path, arm: str, seed: int, content: str | None, why: str
) -> None:
    f = native_seed_dir(native, arm, "P", seed) / "endpoints.json"
    if content is None:
        f.unlink()
    else:
        f.write_text(content)
    manifest = run_collect(native, tmp_path / "out")
    entry = ne_entries(manifest)[(arm, "P", seed)]
    assert entry["outcome"] == "failed" and entry["reason"].startswith(f"endpoint failure: {why}")
    names = sorted(Path(x["path"]).name for x in entry["files"])
    assert names == sorted(["config.txt", "log.txt", "run.json", "out_sha256.txt", *([] if content is None else ["endpoints.json"])])


def test_f_endpoint_status_error_is_failed_and_the_job_continues(native: dict[str, Path], tmp_path: Path) -> None:
    f = native_seed_dir(native, "B-pgcc", "F", 963032) / "record.json"
    record = json.loads(f.read_text())
    record.update(metrics=None, metrics_sha256=None, endpoint_status="error", endpoint_error="ValueError: no dose")
    f.write_text(json.dumps(record))
    manifest = run_collect(native, tmp_path / "out")
    ne = ne_entries(manifest)
    assert set(ne) == {("B-pgcc", "F", 963032)}  # platform_study_record exits 0 on an endpoint error: no absent seeds
    assert ne[("B-pgcc", "F", 963032)]["reason"] == "endpoint failure: record.json endpoint_status 'error' (ValueError: no dose)"
    _, doc = aa.run_analysis(tmp_path / "out", ("B",))
    assert contrast(doc, "B-pgcc vs B-up")["partial"] is True and contrast(doc, "B-picc vs B-up")["partial"] is False


def test_f_record_missing_after_transport_ok_is_failed(native: dict[str, Path], tmp_path: Path) -> None:
    (native_seed_dir(native, "B-picc", "F", 964038) / "record.json").unlink()
    manifest = run_collect(native, tmp_path / "out")
    assert ne_entries(manifest)[("B-picc", "F", 964038)]["reason"] == "endpoint failure: record.json is missing"


# ---------------------------------------------------------------- failed runs do not open a route around the checks


def test_transport_ok_without_dose_or_its_hash_list_still_refuses(native: dict[str, Path], tmp_path: Path) -> None:
    out = native_seed_dir(native, "A-port", "P", 961005) / "out"
    for name in ("sha256.txt", "Dose.mhd"):
        data = (out / name).read_bytes()
        (out / name).unlink()
        refuses(native, tmp_path, f"expected product '{name}' is missing (transport_status is ok)")
        (out / name).write_bytes(data)


def test_failed_run_whose_run_json_names_another_seed_refuses(native: dict[str, Path], tmp_path: Path) -> None:
    d = transport_failure(native, "A-port", "P", 961003, "dose_missing")
    stop_job_after(native, "A-port", "P", 961003)
    edit_run_json(d, seed=961004)
    refuses(native, tmp_path, "seed is 961004, the directory says 961003")


def test_failed_run_whose_config_names_another_seed_refuses(native: dict[str, Path], tmp_path: Path) -> None:
    d = interrupted(native, "B-up", "F", 962038)
    (d / "cfg.txt").write_text((d / "cfg.txt").read_text().replace("RNG_Seed 962038", "RNG_Seed 962037"))
    refuses(native, tmp_path, "RNG_Seed")


def test_failed_run_with_an_unexpected_entry_refuses(native: dict[str, Path], tmp_path: Path) -> None:
    d = transport_failure(native, "B-up", "P", 962028, "no_primaries_line", keep_dose=True)
    (d / "notes.txt").write_text("x")
    refuses(native, tmp_path, "unexpected entry 'notes.txt'")


@pytest.mark.parametrize("left", ["out/sha256.txt", "endpoints.json"])
def test_non_ok_run_holding_what_only_ok_writes_refuses(native: dict[str, Path], tmp_path: Path, left: str) -> None:
    d = native_seed_dir(native, "B-up", "P", 962028)
    data = (d / left).read_bytes()
    transport_failure(native, "B-up", "P", 962028, "dose_missing")
    (d / left).write_bytes(data)
    refuses(native, tmp_path, "the workflows write it only after ok")


def test_endpoint_record_without_run_json_refuses(native: dict[str, Path], tmp_path: Path) -> None:
    (native_seed_dir(native, "A-up", "F", 960038) / "run.json").unlink()
    refuses(native, tmp_path, "record.json exists without run.json")


def test_missing_transport_status_refuses(native: dict[str, Path], tmp_path: Path) -> None:
    edit_run_json(native_seed_dir(native, "A-up", "P", 960028), transport_status=None)
    refuses(native, tmp_path, "transport_status None is not a status the workflows write")


def test_energy_directory_and_its_inputs_list_go_together(native: dict[str, Path], tmp_path: Path) -> None:
    (native["lin"] / "bpg" / "p" / "e150_inputs_sha256.txt").unlink()
    refuses(native, tmp_path, "e150 and e150_inputs_sha256.txt must exist together")


# ----------------------------------------------------- the analysis re-verifies the not_established entries


@pytest.fixture
def failed_collection(native: dict[str, Path], tmp_path: Path) -> Path:
    transport_failure(native, "A-port", "P", 961003, "mcsquare_exit_139")
    stop_job_after(native, "A-port", "P", 961003)
    out = tmp_path / "out"
    run_collect(native, out)
    return out


def edit_manifest(out: Path, edit: object) -> None:
    path = out / aa.MANIFEST_NAME
    manifest = json.loads(path.read_text())
    edit(manifest)  # type: ignore[operator]
    path.write_text(json.dumps(manifest))


def analysis_refuses(out: Path, needle: str, capsys: pytest.CaptureFixture[str]) -> None:
    assert aa.main([str(out)]) == 2
    err = capsys.readouterr().err
    assert needle in err, err


def test_failed_collection_analyses(failed_collection: Path) -> None:
    assert aa.main([str(failed_collection)]) == 0


def test_analysis_refuses_an_altered_not_established_copy(failed_collection: Path, capsys: pytest.CaptureFixture[str]) -> None:
    f = failed_collection / ".not_established" / "A-port" / "P" / "s961003" / "log.txt"
    f.write_text(f.read_text() + "edited\n")
    analysis_refuses(failed_collection, ".not_established/A-port/P/s961003/log.txt is changed since collection", capsys)


def test_analysis_refuses_a_seed_both_collected_and_listed(failed_collection: Path, capsys: pytest.CaptureFixture[str]) -> None:
    edit_manifest(failed_collection, lambda m: m["not_established"].append(
        {"arm": "A-port", "case": "P", "seed": 961001, "energy": 100, "outcome": "absent", "reason": "x", "files": []}))
    analysis_refuses(failed_collection, "961001 both collected as runs and listed as not established", capsys)


def test_analysis_refuses_a_frozen_seed_neither_collected_nor_listed(
    failed_collection: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    edit_manifest(failed_collection, lambda m: m.update(not_established=[e for e in m["not_established"] if e["seed"] != 961028]))
    analysis_refuses(failed_collection, "frozen seed(s) 961028 neither collected nor listed", capsys)


def test_analysis_refuses_a_seed_listed_twice(failed_collection: Path, capsys: pytest.CaptureFixture[str]) -> None:
    edit_manifest(failed_collection, lambda m: m["not_established"].append(dict(m["not_established"][-1])))
    analysis_refuses(failed_collection, "is listed as not established twice", capsys)


@pytest.mark.parametrize(("change", "needle"), [
    ({"seed": 999999}, "is not in the frozen population"),
    ({"arm": "A-up", "seed": 961028}, "is not in the frozen population"),
    ({"energy": 150}, "energy 150, frozen 200"),
    ({"outcome": "collected"}, "needs an outcome"),
    ({"reason": ""}, "needs an outcome"),
])
def test_analysis_refuses_a_malformed_not_established_entry(
    failed_collection: Path, capsys: pytest.CaptureFixture[str], change: dict, needle: str
) -> None:
    def edit(m: dict) -> None:
        last = next(e for e in m["not_established"] if e["seed"] == 961028)
        last.update(change)

    edit_manifest(failed_collection, edit)
    analysis_refuses(failed_collection, needle, capsys)


def test_analysis_refuses_a_failed_runs_file_outside_its_not_established_directory(
    failed_collection: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """A FAILED entry (files allowed) naming a real, correctly hashed file of a collected run: only the path rule
    can refuse it, so the hash check and the absent-has-no-files rule cannot stand in for it."""
    real = "A-port/P/s961001/run.json"
    good = sha((failed_collection / real).read_bytes())

    def edit(m: dict) -> None:
        entry = next(e for e in m["not_established"] if e["seed"] == 961003)
        entry["files"].append({"path": real, "source": "x", "sha256": good})

    edit_manifest(failed_collection, edit)
    analysis_refuses(failed_collection, "lists file(s) outside .not_established/A-port/P/s961003/", capsys)


def test_analysis_refuses_files_listed_for_an_absent_run(failed_collection: Path, capsys: pytest.CaptureFixture[str]) -> None:
    def edit(m: dict) -> None:
        entry = next(e for e in m["not_established"] if e["seed"] == 961028)
        failed = next(e for e in m["not_established"] if e["seed"] == 961003)
        entry["files"] = [{**failed["files"][0], "path": ".not_established/A-port/P/s961028/run.json"}]

    edit_manifest(failed_collection, edit)
    analysis_refuses(failed_collection, "or for an absent run", capsys)


# ----------------------------------------------------------------------------------------- part A's inputs.tar


def _rewrite_tar(path: Path, files: dict[str, bytes]) -> None:
    from apples_fixtures import write_git_archive

    path.unlink()
    write_git_archive(path, files, ACQ)


def _tar_files(path: Path) -> dict[str, bytes]:
    import tarfile

    with tarfile.open(path) as tf:
        return {m.name: tf.extractfile(m).read() for m in tf.getmembers() if m.isfile()}  # type: ignore[union-attr]


def test_part_a_inputs_tar_matching_the_snapshot_is_verified_and_collects(native: dict[str, Path], tmp_path: Path) -> None:
    manifest = run_collect(native, tmp_path / "out")
    assert manifest["files"]
    assert not any("inputs.tar" in t for t in manifest["tolerated"])  # type: ignore[union-attr]


def test_part_a_without_inputs_tar_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    (native["win"] / NATIVE_P / "inputs.tar").unlink()
    refuses(native, tmp_path, "expected product 'inputs.tar' is missing")


def test_inputs_tar_on_a_part_b_arm_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    import shutil

    shutil.copy2(native["win"] / NATIVE_P / "inputs.tar", native["lin"] / "bup" / "p" / "inputs.tar")
    refuses(native, tmp_path, "unexpected entry 'inputs.tar'")


def test_inputs_tar_member_differing_from_the_snapshot_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    tar = native["win"] / NATIVE_F / "inputs.tar"
    files = _tar_files(tar)
    files["Materials/Water.txt"] = b"not water\n"
    _rewrite_tar(tar, files)
    refuses(native, tmp_path, "Materials/Water.txt differs between the archive and the snapshot")


def test_inputs_tar_with_an_extra_or_missing_member_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    tar = native["win"] / NATIVE_P / "inputs.tar"
    files = _tar_files(tar)
    files["validation/extra.py"] = b"# extra\n"
    del files["BDL/BDL_default_DN_RangeShifter.txt"]
    _rewrite_tar(tar, files)
    refuses(native, tmp_path, "validation/extra.py is in the archive but not in the snapshot",
            "BDL/BDL_default_DN_RangeShifter.txt is in the snapshot but not in the archive")


def test_unreadable_inputs_tar_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    (native["win"] / NATIVE_P / "inputs.tar").write_bytes(b"not a tar")
    refuses(native, tmp_path, "unreadable archive")
