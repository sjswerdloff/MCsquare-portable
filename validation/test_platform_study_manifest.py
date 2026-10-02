"""Tests for platform_study_manifest.py on SYNTHETIC inputs only (no study data is read).

Run: uv run --no-project --with numpy --with scipy --with pytest pytest validation/test_platform_study_manifest.py
"""

import hashlib
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import platform_study_analyse as A  # noqa: E402
import platform_study_manifest as G  # noqa: E402
import platform_study_record as R  # noqa: E402


def test_generator_constants_equal_the_analyzers():
    assert (G.STUDY_COMMIT, G.RERUN_COMMIT) == (A.STUDY_COMMIT, A.RERUN_COMMIT)
    assert G.RERUN_SEEDS == A.RERUN_SEEDS
    assert {w: {p: list(r) for p, r in d.items()} for w, d in G.WAVE_SEEDS.items()} == \
           {w: {p: list(r) for p, r in d.items()} for w, d in A.WAVE_SEEDS.items()}
    assert set(G.TEXT_PATHS) | {"plan E200_S150.txt", "cube.mhd"} == set(A.TEXT_INPUTS)


@pytest.mark.parametrize(("raw", "expected"), [
    (b"a\nb\n", b"a\r\nb\r\n"),
    (b"a\r\nb\n", b"a\r\nb\r\n"),          # an LF already preceded by CR is left alone
    (b"a\rb\n", b"a\rb\r\n"),              # a lone CR is left alone
    (b"", b""),
    (b"no newline", b"no newline"),
    (b"a\0\nb\n", b"a\0\nb\n"),            # a NUL byte: treated as binary, unchanged
])
def test_to_crlf(raw, expected):
    assert G.to_crlf(raw) == expected


def test_git_binary_heuristic_agrees_with_git_on_the_basic_cases():
    assert not G.git_says_binary(b"plain text\nwith lines\r\n\ttabs\n")
    assert G.git_says_binary(b"nul\0byte")
    assert G.git_says_binary(b"lone\rcr")
    assert G.git_says_binary(bytes(range(1, 8)) * 40)


FILES = {"Materials/B/x.dat": b"1\n2\n", "Materials/A/y.dat": b"y\n", "Materials/A2/z.dat~": b"z\r\nzz\n",
         "Materials/top.txt": b"t\n"}


def record_digest(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, files: dict[str, bytes]) -> str:
    for p, b in files.items():
        (tmp_path / p).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / p).write_bytes(b)
    monkeypatch.chdir(tmp_path)
    return R.tree_digest("Materials/**/*")["combined_sha256"]


def test_tree_digest_reproduces_the_record_lf(tmp_path, monkeypatch):
    assert G.tree_digest(FILES, "posix") == record_digest(tmp_path, monkeypatch, FILES)


def test_tree_digest_reproduces_the_record_crlf(tmp_path, monkeypatch):
    crlf = {p: G.to_crlf(b) for p, b in FILES.items()}
    assert G.tree_digest(crlf, "posix") == record_digest(tmp_path, monkeypatch, crlf)
    assert G.tree_digest(crlf, "posix") != G.tree_digest(FILES, "posix")


def test_windows_order_differs_from_posix_when_a_directory_name_prefixes_another():
    # "A/y" < "A2/z" in posix order but "A\\y" > "A2\\z" in backslash order: the digests must differ
    assert G.tree_digest(FILES, "posix") != G.tree_digest(FILES, "windows")


def test_windows_order_matches_the_record_run_with_backslash_paths(monkeypatch):
    """Drive the record's own tree_digest with backslash-form paths, as a Windows glob returns them."""
    digests = {p: hashlib.sha256(b).hexdigest() for p, b in FILES.items()}
    win = {p.replace("/", "\\"): p for p in FILES}
    monkeypatch.setattr(R.glob, "glob", lambda pattern, recursive: list(win))
    monkeypatch.setattr(R.os.path, "isfile", lambda p: True)
    monkeypatch.setattr(R, "sha256", lambda p: digests[win[p]])
    monkeypatch.setattr(R.os, "sep", "\\")
    got = R.tree_digest("Materials/**/*")["combined_sha256"]
    monkeypatch.undo()
    assert got == G.tree_digest(FILES, "windows")
    assert got != G.tree_digest(FILES, "posix")


WORKFLOW = '''
env:
  PRIMARIES: "3e7"
jobs:
  linux:
    steps:
      - run: |
          printf '%s\\n' "Num_Threads 3" "Num_Primaries $PRIMARIES" "RNG_Seed $SEED" \\
            "Output_Directory out_seed" > cfg.txt
  macos:
    steps:
      - run: |
          printf '%s\\n' "Num_Threads 3" "Num_Primaries $PRIMARIES" "RNG_Seed $SEED" \\
            "Output_Directory out_seed" > cfg.txt
  windows:
    steps:
      - run: |
          @("Num_Threads 3", "Num_Primaries $env:PRIMARIES", "RNG_Seed $seed",
            "Output_Directory out_seed") |
            Set-Content -Encoding Ascii cfg.txt
'''
LINES = ["Num_Threads 3", "Num_Primaries 3e7", "RNG_Seed {SEED}", "Output_Directory out_seed"]


@pytest.mark.parametrize("plat", ["linux", "macos", "windows"])
def test_parse_cfg_template(plat):
    assert G.parse_cfg_template(WORKFLOW, plat) == LINES


def test_cfg_bytes_line_endings_per_platform():
    assert G.cfg_bytes(LINES, "linux", 2001) == b"Num_Threads 3\nNum_Primaries 3e7\nRNG_Seed 2001\nOutput_Directory out_seed\n"
    assert G.cfg_bytes(LINES, "macos", 3001).endswith(b"out_seed\n") and b"\r" not in G.cfg_bytes(LINES, "macos", 3001)
    win = G.cfg_bytes(LINES, "windows", 1001)
    assert win == b"Num_Threads 3\r\nNum_Primaries 3e7\r\nRNG_Seed 1001\r\nOutput_Directory out_seed\r\n"


def test_template_with_run_id_in_the_output_directory_is_refused_by_the_parser():
    wf = WORKFLOW.replace('"Output_Directory out_seed" > cfg.txt', '"Output_Directory out_${SEED}_${GITHUB_RUN_ID}" > cfg.txt')
    with pytest.raises(G.ManifestError, match="unresolved variable"):
        G.parse_cfg_template(wf, "linux")


def test_parser_refuses_a_workflow_without_the_expected_block():
    with pytest.raises(G.ManifestError):
        G.parse_cfg_template(WORKFLOW.replace("> cfg.txt", "> other.txt"), "linux")
    with pytest.raises(G.ManifestError):
        G.parse_cfg_template(WORKFLOW.replace("Set-Content -Encoding Ascii", "Set-Content -Encoding UTF8"), "windows")
    with pytest.raises(G.ManifestError):
        G.parse_cfg_template(WORKFLOW.replace('PRIMARIES: "3e7"', "PRIMARIES: 3e7"), "linux")
