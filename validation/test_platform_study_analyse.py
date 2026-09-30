"""Behavioural tests for platform_study_analyse.py on SYNTHETIC records only (no study data is read).

Run: uv run --no-project --with numpy --with scipy --with pytest pytest validation/test_platform_study_analyse.py
"""

import copy
import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import platform_study_analyse as A  # noqa: E402

PLATFORMS = ("windows", "linux", "macos")
BIN = {"windows": "a" * 64, "linux": "b" * 64, "macos": "c" * 64}


def syn_manifest() -> dict:
    """A SYNTHETIC frozen manifest: sentinel strings stand in for hashes (tests never read the real manifest)."""
    return {
        "study_commit": A.STUDY_COMMIT,
        "text": {n: {"LF": f"lf_{n}", "CRLF": f"crlf_{n}"} for n in A.TEXT_INPUTS},
        "cube.raw": "raw",
        "materials": {"LF": {"posix": "M_lf", "windows": "M_lf"}, "CRLF": {"posix": "M_crlf", "windows": "M_crlf"}},
        "config": {p: {str(s): f"cfg_{p}_{s}" for w in A.WAVE_SEEDS.values() for s in w[p]} for p in PLATFORMS},
    }


def rec(plat: str, seed: int, shift: float = 0.0, commit: str | None = None, eol: str = "LF") -> dict:
    """One manifest-consistent synthetic record; values are platform-independent unless shifted, plus small noise."""
    noise = ((seed * 7919) % 97 - 48) / 48 * 0.01
    metrics = {ep: 3.0 + shift + noise for ep in A.LATERAL}
    metrics |= {ep: 1.0e-3 * (1 + noise / 10) for ep in A.CAX}
    metrics |= {"r80_mm": 250.0 + noise, "r20_mm": 255.0 + noise, "r80_crossings": 1, "r20_crossings": 1}
    if commit is None:
        commit = A.RERUN_COMMIT if (plat, seed) in A.RERUN_SEEDS else A.STUDY_COMMIT
    sha = {n: f"{eol.lower()}_{n}" for n in A.TEXT_INPUTS}
    sha |= {"cube.raw": "raw", "config": f"cfg_{plat}_{seed}", "binary": BIN[plat]}
    return {"study": "pe1", "platform": plat, "seed": seed, "host": "h", "commit": commit, "sha256": sha,
            "materials": {"combined_sha256": f"M_{eol.lower()}"}, "metrics": metrics}


def wave(w: int, plats=PLATFORMS, commit: str | None = None, eol: str = "LF") -> list[dict]:
    return [rec(p, s, commit=commit, eol=eol) for p in plats for s in A.WAVE_SEEDS[w][p]]


def write_manifest(tmp_path: Path, manifest: dict | None = None) -> Path:
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(syn_manifest() if manifest is None else manifest, indent=1, sort_keys=True))
    return path


def write_wave2(tmp_path: Path, commit: str | None) -> Path:
    path = tmp_path / "wave2.json"
    path.write_text(json.dumps({"wave2_commit": commit}))
    return path


def run(tmp_path: Path, recs: list[dict], look: int, verdicts: Path | None = None, manifest: Path | None = None,
        analyzer: Path | None = None, wave2: Path | None = None) -> subprocess.CompletedProcess:
    data = tmp_path / f"records_{look}.jsonl"
    data.write_text("".join(json.dumps(r) + "\n" for r in recs))
    v = verdicts or tmp_path / "look1.json"
    m = manifest or (tmp_path / "manifest.json")
    if not m.exists():
        write_manifest(tmp_path)
    w2 = wave2 or (tmp_path / "wave2.json")
    if not w2.exists():
        write_wave2(tmp_path, None)
    return subprocess.run([sys.executable, str(analyzer or HERE / "platform_study_analyse.py"), str(data), str(look), str(v),
                           "--manifest", str(m), "--wave2", str(w2)], capture_output=True, text=True)


def refused(p: subprocess.CompletedProcess, why: str) -> bool:
    return p.returncode != 0 and "REFUSING TO ANALYSE" in p.stderr and why in p.stderr


def look1_state(manifest: Path, recs: list[dict], verdicts: dict[str, str]) -> dict:
    """A look-1 state file built with the analyzer's own binding functions (in-process tests)."""
    return {"verdicts": verdicts, "analyzer_sha256": A.sha_of(A.__file__), "manifest_sha256": A.sha_of(manifest),
            "wave1_fingerprint": A.fingerprint(A.wave1_records(recs))}


def test_complete_wave1_runs_and_writes_verdicts(tmp_path):
    p = run(tmp_path, wave(1), 1)
    assert p.returncode == 0, p.stderr
    state = json.loads((tmp_path / "look1.json").read_text())
    assert set(state) == {"verdicts", "analyzer_sha256", "manifest_sha256", "wave1_fingerprint"}
    assert set(state["verdicts"]) == {"windows", "macos"}
    assert state["analyzer_sha256"] == A.sha_of(HERE / "platform_study_analyse.py")
    assert state["manifest_sha256"] == A.sha_of(tmp_path / "manifest.json")
    assert state["wave1_fingerprint"] == A.fingerprint(wave(1))
    assert "VERDICT windows vs linux (look 1)" in p.stdout
    assert "distinct binary sha256 per platform (not source attestation; build time is embedded): linux 1, macos 1, windows 1" in p.stdout


def test_look1_refuses_to_overwrite_its_verdicts(tmp_path):
    assert run(tmp_path, wave(1), 1).returncode == 0
    assert refused(run(tmp_path, wave(1), 1), "already analysed")


@pytest.mark.parametrize("look", [0, 3, -1])
def test_invalid_look(tmp_path, look):
    assert refused(run(tmp_path, wave(1), look), "look must be 1 or 2")


def test_no_data(tmp_path):
    assert refused(run(tmp_path, [], 1), "no records")


def test_missing_seed(tmp_path):
    assert refused(run(tmp_path, wave(1)[1:], 1), "missing")


def test_unexpected_seed(tmp_path):
    assert refused(run(tmp_path, wave(1) + [rec("linux", 2999)], 1), "unexpected")


def test_duplicate_seed(tmp_path):
    recs = wave(1)
    assert refused(run(tmp_path, recs + [copy.deepcopy(recs[0])], 1), "duplicate")


def test_missing_platform(tmp_path):
    assert refused(run(tmp_path, wave(1, ("windows", "linux")), 1), "missing")


def test_foreign_study(tmp_path):
    recs = wave(1)
    recs[5]["study"] = "other"
    assert refused(run(tmp_path, recs, 1), "study")


def test_foreign_commit(tmp_path):
    recs = wave(1)
    recs[3]["commit"] = "0" * 40
    assert refused(run(tmp_path, recs, 1), "commit")


def test_rerun_commit_allowed_only_for_the_recorded_seeds(tmp_path):
    recs = wave(1)
    assert any(r["commit"] == A.RERUN_COMMIT for r in recs)  # control: the two re-runs carry it and pass
    assert run(tmp_path, recs, 1).returncode == 0
    recs2 = wave(1)
    next(r for r in recs2 if (r["platform"], r["seed"]) == ("macos", 3001))["commit"] = A.RERUN_COMMIT
    other = tmp_path / "b"
    other.mkdir()
    assert refused(run(other, recs2, 1), "commit")


def test_input_differs_within_a_platform(tmp_path):
    recs = wave(1)
    victim = next(r for r in recs if r["platform"] == "linux")
    victim["sha256"]["BDL"] = "other"
    p = run(tmp_path, recs, 1)
    assert refused(p, "sha256['BDL']") and f"({victim['platform']}, {victim['seed']})" in p.stderr


def test_cube_raw_differs_across_platforms(tmp_path):
    recs = wave(1)
    for r in recs:
        if r["platform"] == "macos":
            r["sha256"]["cube.raw"] = "other"
    assert refused(run(tmp_path, recs, 1), "sha256['cube.raw']")


def test_cube_raw_has_no_crlf_allowance_even_on_windows(tmp_path):
    recs = wave(1)
    next(r for r in recs if r["platform"] == "windows")["sha256"]["cube.raw"] = "crlf_raw"
    assert refused(run(tmp_path, recs, 1), "sha256['cube.raw']")


def test_old_differs_but_content_fixed_by_commit_note_is_gone(tmp_path):
    recs = wave(1, eol="CRLF")  # every platform CRLF: linux/macos must refuse, not note
    p = run(tmp_path, recs, 1)
    assert p.returncode != 0 and "note:" not in p.stdout


@pytest.mark.parametrize("field", A.TEXT_INPUTS)
def test_arbitrary_windows_text_input_change_refuses(tmp_path, field):
    recs = wave(1)
    next(r for r in recs if r["platform"] == "windows" and r["seed"] == 1005)["sha256"][field] = "edited"
    p = run(tmp_path, recs, 1)
    assert refused(p, f"sha256[{field!r}]") and "(windows, 1005)" in p.stderr


def test_arbitrary_windows_bdl_change_refuses_even_if_consistent_across_windows_seeds(tmp_path):
    recs = wave(1)
    for r in recs:
        if r["platform"] == "windows":
            r["sha256"]["BDL"] = "edited"  # the old within-platform-consistency check would have passed this
    assert refused(run(tmp_path, recs, 1), "sha256['BDL']")


def test_arbitrary_windows_materials_change_refuses(tmp_path):
    recs = wave(1)
    for r in recs:
        if r["platform"] == "windows":
            r["materials"]["combined_sha256"] = "edited"
    p = run(tmp_path, recs, 1)
    assert refused(p, "materials.combined_sha256") and "(windows," in p.stderr


def test_wrong_config_refuses(tmp_path):
    recs = wave(1)
    victim = next(r for r in recs if r["platform"] == "macos" and r["seed"] == 3004)
    victim["sha256"]["config"] = "cfg_with_other_primaries"
    p = run(tmp_path, recs, 1)
    assert refused(p, "sha256['config']") and "(macos, 3004)" in p.stderr


def test_config_of_another_seed_refuses(tmp_path):
    recs = wave(1)
    next(r for r in recs if r["platform"] == "linux" and r["seed"] == 2004)["sha256"]["config"] = "cfg_linux_2005"
    assert refused(run(tmp_path, recs, 1), "sha256['config']")


def test_crlf_only_variant_on_windows_passes(tmp_path):
    recs = wave(1, ("linux", "macos")) + wave(1, ("windows",), eol="CRLF")
    p = run(tmp_path, recs, 1)
    assert p.returncode == 0, p.stderr


def test_crlf_and_lf_windows_records_are_each_accepted(tmp_path):
    recs = wave(1, ("linux", "macos")) + wave(1, ("windows",))
    recs[-1] = rec("windows", 1016, eol="CRLF")
    assert run(tmp_path, recs, 1).returncode == 0


@pytest.mark.parametrize("plat", ["linux", "macos"])
def test_crlf_text_input_on_posix_platform_refuses(tmp_path, plat):
    recs = wave(1)
    next(r for r in recs if r["platform"] == plat and r["seed"] == A.WAVE_SEEDS[1][plat][3])["sha256"]["BDL"] = "crlf_BDL"
    p = run(tmp_path, recs, 1)
    assert refused(p, "sha256['BDL']") and f"({plat}," in p.stderr


@pytest.mark.parametrize("plat", ["linux", "macos"])
def test_crlf_materials_on_posix_platform_refuses(tmp_path, plat):
    recs = wave(1)
    next(r for r in recs if r["platform"] == plat)["materials"]["combined_sha256"] = "M_crlf"
    assert refused(run(tmp_path, recs, 1), "materials.combined_sha256")


def test_mixed_binary_hashes_within_a_platform_are_reported_not_refused(tmp_path):
    recs = wave(1)
    next(r for r in recs if r["platform"] == "windows" and r["seed"] == 1012)["sha256"]["binary"] = "d" * 64
    p = run(tmp_path, recs, 1)
    assert p.returncode == 0 and "windows 2" in p.stdout


def test_different_binaries_across_platforms_are_fine_and_printed(tmp_path):
    p = run(tmp_path, wave(1), 1)
    assert p.returncode == 0 and "distinct binary sha256 per platform" in p.stdout


def test_missing_binary_hash_refuses(tmp_path):
    recs = wave(1)
    for r in recs:
        if r["platform"] == "macos":
            del r["sha256"]["binary"]
    assert refused(run(tmp_path, recs, 1), "binary sha256 missing or not 64 lowercase hex")


def test_manifest_absent_refuses(tmp_path):
    p = run(tmp_path, wave(1), 1, manifest=tmp_path / "nowhere.json")
    assert refused(p, "is absent")


def test_manifest_for_another_study_commit_refuses(tmp_path):
    m = syn_manifest() | {"study_commit": "0" * 40}
    assert refused(run(tmp_path, wave(1), 1, manifest=write_manifest(tmp_path, m)), "study commit")


def test_manifest_without_config_for_a_seed_refuses(tmp_path):
    m = syn_manifest()
    del m["config"]["linux"]["2007"]
    p = run(tmp_path, wave(1), 1, manifest=write_manifest(tmp_path, m))
    assert refused(p, "no config hash for (linux, 2007)")


def test_malformed_manifest_refuses(tmp_path):
    m = syn_manifest()
    del m["text"]["BDL"]
    assert refused(run(tmp_path, wave(1), 1, manifest=write_manifest(tmp_path, m)), "malformed")


def test_default_manifest_path_is_next_to_the_analyzer():
    assert A.MANIFEST_PATH == HERE / "platform_study_manifest.json"


def test_fingerprint_is_order_independent_and_sensitive_to_any_record_change():
    recs = wave(1)
    f = A.fingerprint(recs)
    assert A.fingerprint(list(reversed(recs))) == f
    changed = copy.deepcopy(recs)
    changed[7]["metrics"]["cax_127"] *= 1.0000001
    assert A.fingerprint(changed) != f
    assert A.fingerprint(recs[1:]) != f


def test_nonfinite_endpoint_cannot_pass(tmp_path):
    recs = wave(1)
    next(r for r in recs if r["platform"] == "windows")["metrics"]["lateral_127_10"] = math.nan
    p = run(tmp_path, recs, 1)
    assert p.returncode == 0
    assert "lateral_127_10" in p.stdout and "CANNOT PASS: non-finite" in p.stdout
    assert "EQUIVALENT" not in json.loads((tmp_path / "look1.json").read_text())["verdicts"]["windows"]


def test_look2_needs_look1_verdicts(tmp_path):
    assert refused(run(tmp_path, wave(1) + wave(2), 2), "needs look 1")


def test_look2_needs_frozen_wave2_commit(tmp_path):
    m = write_manifest(tmp_path)
    cont = {"windows": "CONTINUE to wave 2", "macos": "CONTINUE to wave 2"}
    (tmp_path / "look1.json").write_text(json.dumps(look1_state(m, wave(1), cont)))
    assert refused(run(tmp_path, wave(1) + wave(2), 2), "wave2_commit not frozen")


def test_stopped_comparison_is_frozen_at_look2(tmp_path, monkeypatch):
    """windows stopped at look 1; look 2 needs only linux+macos wave 2 and never recomputes windows."""
    m = write_manifest(tmp_path)
    monkeypatch.setattr(A, "MANIFEST_PATH", m)
    monkeypatch.setattr(A, "WAVE2_PATH", write_wave2(tmp_path, "f" * 40))
    recs = wave(1) + wave(2, ("linux", "macos"), commit="f" * 40)
    verdicts = {"windows": "EQUIVALENT", "macos": "CONTINUE to wave 2"}
    (tmp_path / "look1.json").write_text(json.dumps(look1_state(m, wave(1), verdicts)))
    data = tmp_path / "r.jsonl"
    data.write_text("".join(json.dumps(r) + "\n" for r in recs))
    with pytest.raises(SystemExit) as e:  # the in-process run exits 0 at the end only if nothing refused
        A.main(str(data), 2, str(tmp_path / "look1.json"))
        raise SystemExit(0)
    assert e.value.code == 0
    # windows wave-2 records present -> refused as unexpected
    recs_w = recs + wave(2, ("windows",), commit="f" * 40)
    data.write_text("".join(json.dumps(r) + "\n" for r in recs_w))
    with pytest.raises(SystemExit, match="unexpected"):
        A.main(str(data), 2, str(tmp_path / "look1.json"))


def test_no_comparison_continued_means_no_look2(tmp_path):
    m = write_manifest(tmp_path)
    verdicts = {"windows": "EQUIVALENT", "macos": "NON-EQUIVALENT on cax_127"}
    (tmp_path / "look1.json").write_text(json.dumps(look1_state(m, wave(1), verdicts)))
    assert refused(run(tmp_path, wave(1), 2), "there is no look 2")


# ---- look-1 state binding (end to end: a COPY of the analyzer runs both looks so its own hash can be changed) ----

def two_looks(tmp_path: Path) -> tuple[Path, Path, list[dict]]:
    """Run look 1 with a private analyzer copy (wave2.json still null); force both verdicts to CONTINUE, then freeze wave 2.

    Returns (analyzer copy, manifest, wave-2 records). The bindings in look1.json are the real ones the copy wrote."""
    ana = tmp_path / "ana" / "platform_study_analyse.py"
    ana.parent.mkdir()
    ana.write_text((HERE / "platform_study_analyse.py").read_text())
    m = write_manifest(tmp_path)
    assert run(tmp_path, wave(1), 1, analyzer=ana).returncode == 0
    state = json.loads((tmp_path / "look1.json").read_text())
    state["verdicts"] = {"windows": "CONTINUE to wave 2", "macos": "CONTINUE to wave 2"}
    (tmp_path / "look1.json").write_text(json.dumps(state))
    write_wave2(tmp_path, "f" * 40)  # frozen AFTER look 1: must not disturb the look-1 binding
    return ana, m, wave(2, commit="f" * 40)


def test_look2_baseline_unchanged_state_is_accepted(tmp_path):
    ana, _, w2 = two_looks(tmp_path)
    p = run(tmp_path, wave(1) + w2, 2, analyzer=ana)
    assert p.returncode == 0, p.stderr
    assert "VERDICT windows vs linux (look 2)" in p.stdout


def test_look2_with_altered_wave1_record_refuses(tmp_path):
    ana, _, w2 = two_looks(tmp_path)
    recs = wave(1)
    next(r for r in recs if r["platform"] == "linux" and r["seed"] == 2009)["metrics"]["cax_127"] *= 1.001
    assert refused(run(tmp_path, recs + w2, 2, analyzer=ana), "wave-1 records differ")


def test_look2_with_any_field_of_a_wave1_record_changed_refuses(tmp_path):
    ana, _, w2 = two_looks(tmp_path)
    recs = wave(1)
    recs[0]["host"] = "another host"  # any field of the record is bound, not only the metrics
    assert refused(run(tmp_path, recs + w2, 2, analyzer=ana), "wave-1 records differ")


def test_look2_with_changed_analyzer_refuses(tmp_path):
    ana, _, w2 = two_looks(tmp_path)
    ana.write_text(ana.read_text() + "\n# an edit after look 1\n")
    assert refused(run(tmp_path, wave(1) + w2, 2, analyzer=ana), "analyzer file changed")


def test_look2_with_changed_manifest_refuses(tmp_path):
    ana, m, w2 = two_looks(tmp_path)
    m.write_text(json.dumps(json.loads(m.read_text()), indent=4))  # same content, different bytes
    assert refused(run(tmp_path, wave(1) + w2, 2, analyzer=ana), "manifest file changed")


def test_look2_with_state_lacking_the_bindings_refuses(tmp_path):
    ana, _, w2 = two_looks(tmp_path)
    state = json.loads((tmp_path / "look1.json").read_text())
    del state["wave1_fingerprint"]
    (tmp_path / "look1.json").write_text(json.dumps(state))
    assert refused(run(tmp_path, wave(1) + w2, 2, analyzer=ana), "lacks 'wave1_fingerprint'")


def test_editing_wave2_json_after_look1_does_not_break_the_look2_binding(tmp_path):
    ana, _, w2 = two_looks(tmp_path)  # look 1 ran with wave2_commit null; it is frozen only now
    state = json.loads((tmp_path / "look1.json").read_text())
    assert "wave2" not in json.dumps(state).lower()  # nothing about it is bound
    p = run(tmp_path, wave(1) + w2, 2, analyzer=ana)
    assert p.returncode == 0, p.stderr


def test_wave2_records_refused_while_wave2_commit_is_null(tmp_path):
    ana, _, w2 = two_looks(tmp_path)
    write_wave2(tmp_path, None)
    assert refused(run(tmp_path, wave(1) + w2, 2, analyzer=ana), "wave2_commit not frozen")


def test_wave2_records_refused_while_wave2_file_is_absent(tmp_path):
    ana, _, w2 = two_looks(tmp_path)
    p = run(tmp_path, wave(1) + w2, 2, analyzer=ana, wave2=tmp_path / "missing.json")
    assert refused(p, "wave2_commit not frozen")


def test_wave2_records_at_another_commit_refused(tmp_path):
    ana, _, _ = two_looks(tmp_path)
    assert refused(run(tmp_path, wave(1) + wave(2, commit="e" * 40), 2, analyzer=ana), "built at commit")


@pytest.mark.parametrize("content", ['{"wave2_commit": "abc"}', '{"other": 1}', "not json", '{"wave2_commit": 5}'])
def test_malformed_wave2_json_refuses(tmp_path, content):
    ana, _, w2 = two_looks(tmp_path)
    (tmp_path / "wave2.json").write_text(content)
    assert refused(run(tmp_path, wave(1) + w2, 2, analyzer=ana), "wave2")


@pytest.mark.parametrize("bad", ["", "abc", "a" * 63, "a" * 65, "g" * 64, "A" * 64, " " + "a" * 63, None, 5])
def test_malformed_binary_sha256_refuses(tmp_path, bad):
    recs = wave(1)
    next(r for r in recs if r["platform"] == "linux" and r["seed"] == 2005)["sha256"]["binary"] = bad
    assert refused(run(tmp_path, recs, 1), "binary sha256 missing or not 64 lowercase hex")


def test_distinct_valid_binary_hashes_within_a_platform_pass(tmp_path):
    recs = wave(1)
    for i, r in enumerate(x for x in recs if x["platform"] == "linux"):
        r["sha256"]["binary"] = f"{i:064x}"
    p = run(tmp_path, recs, 1)
    assert p.returncode == 0 and "linux 16" in p.stdout
