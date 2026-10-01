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
    make_native,
    sha,
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


def test_missing_expected_product_is_refused(native: dict[str, Path], tmp_path: Path) -> None:
    (native["win"] / NATIVE_P / "e150" / "s961011" / "endpoints.json").unlink()
    refuses(native, tmp_path, "endpoints.json", "missing")


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
    (native["win"] / NATIVE_P / "e150" / "s961011" / "endpoints.json").unlink()
    (native["lin"] / "bpg" / "f" / "s963031" / "out_seed" / "Dose.raw").unlink()
    (native["win"] / "smoke").mkdir()
    message = refuses(native, tmp_path, "problem(s)")
    assert "refusing, 3 problem(s)" not in message  # a missing product is reported more than once; all three appear
    assert "endpoints.json" in message and "Dose.raw" in message and "smoke" in message


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
