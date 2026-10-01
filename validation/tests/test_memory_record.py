"""Tests for memory_record.py: valid controls taken from task 10538's real output, and each refusal."""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import memory_record as mr

# The two lines as /usr/bin/time -l printed them for the pencil case in task 10538 (run 305).
TIME_OK = (
    "        0.50 real         1.20 user         0.30 sys\n"
    "          2137325568  maximum resident set size\n"
    "                   0  average shared memory size\n"
    "          2137343872  peak memory footprint\n"
)
LOG_OK = "Simulation started\nNbr primaries simulated: 100000 \nSimulation completed\n"
CT_OK = "ObjectType = Image\nNDims = 3\nDimSize = 400 350 400\nElementSpacing = 1.000000 1.000000 1.000000\n"


def _args(tmp_path, time=TIME_OK, log=LOG_OK, ct=CT_OK, requested="1e5"):
    for name, text in (("time.txt", time), ("log.txt", log), ("ct.mhd", ct)):
        (tmp_path / name).write_text(text)
    return [
        "--case", "pencil", "--time", str(tmp_path / "time.txt"), "--log", str(tmp_path / "log.txt"),
        "--ct", str(tmp_path / "ct.mhd"), "--requested", requested, "--threads", "4",
        "--binary-sha256", "71111c17", "--material", "water", "--os", "macOS 26.6.1",
    ]  # fmt: skip


def test_valid_control_gives_typed_record(tmp_path, capsys):
    assert mr.main(_args(tmp_path)) == 0
    rec = json.loads(capsys.readouterr().out)
    assert rec["macos_max_rss_bytes"] == 2137325568 and rec["macos_peak_footprint_bytes"] == 2137343872
    assert rec["simulated"] == 100000 and rec["requested"] == 100000
    assert rec["ct_dims"] == [400, 350, 400] and rec["voxels"] == 56_000_000


def test_missing_footprint_is_null_not_refused(tmp_path, capsys):
    t = TIME_OK.replace("          2137343872  peak memory footprint\n", "")
    assert mr.main(_args(tmp_path, time=t)) == 0
    assert json.loads(capsys.readouterr().out)["macos_peak_footprint_bytes"] is None


@pytest.mark.parametrize(
    "time",
    [
        TIME_OK.replace("2137325568  maximum", "not_a_number  maximum"),
        TIME_OK + "          1000  maximum resident set size\n",
        TIME_OK.replace("2137325568  maximum", "0  maximum"),
        TIME_OK.replace("          2137325568  maximum resident set size\n", ""),
        TIME_OK.replace("2137343872  peak", "-5  peak"),
        TIME_OK + "          1  peak memory footprint\n",
        "",
    ],
    ids=["rss-nonnumeric", "rss-twice", "rss-zero", "rss-missing", "footprint-negative", "footprint-twice", "empty"],
)
def test_malformed_time_output_is_refused(tmp_path, capsys, time):
    assert mr.main(_args(tmp_path, time=time)) == 1
    assert capsys.readouterr().out == ""


@pytest.mark.parametrize(
    "log",
    [
        "Nbr primaries simulated: 99999\n",
        "",
        LOG_OK + "Nbr primaries simulated: 100000\n",
        LOG_OK + "WARNING: Unknown tag Foo\n",
        LOG_OK + "12 primaries generated outside the geometry\n",
    ],
    ids=["short", "absent", "twice", "unknown-tag", "outside-geometry"],
)
def test_incomplete_or_wrong_run_is_refused(tmp_path, log):
    assert mr.main(_args(tmp_path, log=log)) == 1


@pytest.mark.parametrize(
    "ct",
    [CT_OK.replace("DimSize = 400 350 400", "DimSize = 400 0 400"), CT_OK.replace("ElementSpacing", "Spacing"), "x"],
    ids=["zero-dim", "no-spacing", "garbage"],
)
def test_invalid_ct_header_is_refused(tmp_path, ct):
    assert mr.main(_args(tmp_path, ct=ct)) == 1


@pytest.mark.parametrize("requested", ["1.5", "0", "-1"])
def test_requested_must_be_a_positive_whole_count(tmp_path, requested):
    assert mr.main(_args(tmp_path, requested=requested)) == 1
