"""Behavioural tests for mcsquare_pencil_case.py: a valid case is written; an invalid one writes nothing."""

import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "mcsquare_pencil_case.py"


def _run(tmp_path: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), *args], cwd=tmp_path, capture_output=True, text=True, check=False)


def test_valid_case_writes_all_inputs_with_lf_and_exact_bdl_header(tmp_path):
    r = _run(tmp_path, "200", "--nx", "4", "--ny", "6", "--nz", "4", "--iso-y", "6")
    assert r.returncode == 0, r.stderr
    assert sorted(p.name for p in tmp_path.iterdir()) == ["bdl_E200.txt", "plan_E200.txt", "water.mhd", "water.raw"]
    bdl = (tmp_path / "bdl_E200.txt").read_bytes()
    assert bdl.startswith(b"--UPenn beam model (double gaussian)--\n")  # data_beam_model.c strcmp, LF only
    assert b"\r" not in bdl
    assert (tmp_path / "water.raw").stat().st_size == 4 * 4 * 6 * 4
    assert b"2.000000\t 6.000000\t 2.000000" in (tmp_path / "plan_E200.txt").read_bytes()


@pytest.mark.parametrize(
    "args",
    [
        ["0", "--iso-y", "6"],  # energy below range
        ["nan", "--iso-y", "6"],
        ["inf", "--iso-y", "6"],
        ["200.5", "--iso-y", "6"],  # file names round it
        ["301", "--iso-y", "6"],
        ["200", "--iso-y", "6", "--nozzle", "0"],  # the documented broken case: source on the entrance face
        ["200", "--iso-y", "3", "--nozzle", "1"],  # source inside the phantom
        ["200", "--iso-y", "7", "--nozzle", "1"],  # isocentre outside the CT
        ["200", "--iso-y", "nan"],
        ["200", "--iso-y", "6", "--nozzle", "-1"],
        ["200", "--iso-y", "6", "--nx", "5"],  # odd lateral size
        ["200", "--iso-y", "0", "--ny", "0"],
    ],
)
def test_invalid_settings_are_refused_and_write_nothing(tmp_path, args):
    full = args + (["--nx", "4"] if "--nx" not in args else []) + (["--ny", "6"] if "--ny" not in args else []) + ["--nz", "4"]
    r = _run(tmp_path, *full)
    assert r.returncode != 0
    assert list(tmp_path.iterdir()) == []
