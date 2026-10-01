"""Build one fail-closed peak-memory record for an MCsquare run timed with macOS `/usr/bin/time -l`.

Refuses (exit 1, message on stderr, nothing on stdout) unless:
  - the time output has exactly one "maximum resident set size" line, a positive integer;
  - "peak memory footprint" is either absent (recorded as null) or exactly one positive integer;
  - the MCsquare log has exactly one "Nbr primaries simulated: N" line with N >= the requested count;
  - the log has no "Unknown tag" (unrecognised config key) and no primaries "generated outside the geometry";
  - the CT header gives three positive dimensions and spacings.

The record states what was measured and no more: macOS RSS and macOS physical footprint, in bytes, of one short run,
with the conditions needed to repeat it. It is not a general peak for longer runs or other modes, and not a
measurement of any other OS.

Usage:
    python memory_record.py --case NAME --time time.txt --log log.txt --ct CT.mhd --requested 1e5 \
        --threads 4 --binary-sha256 HEX --material TEXT [--os TEXT]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


class Refused(ValueError):
    """The inputs do not establish a valid measurement."""


def _one_count(text: str, label: str, required: bool) -> int | None:
    """The single integer printed before `label` by `time -l`; None if absent and not required."""
    lines = [line for line in text.splitlines() if line.strip().endswith(label)]
    if not lines:
        if required:
            raise Refused(f"no '{label}' line")
        return None
    if len(lines) > 1:
        raise Refused(f"{len(lines)} '{label}' lines")
    value = lines[0].strip()[: -len(label)].strip()
    if not value.isdigit() or int(value) <= 0:
        raise Refused(f"'{label}' is not a positive integer: {value!r}")
    return int(value)


def parse_time(text: str) -> tuple[int, int | None]:
    """(max RSS bytes, peak footprint bytes or None) from macOS `/usr/bin/time -l` output."""
    rss = _one_count(text, "maximum resident set size", required=True)
    assert rss is not None
    return rss, _one_count(text, "peak memory footprint", required=False)


def parse_log(text: str, requested: int) -> int:
    """Simulated primaries from the MCsquare log, after checking the run did what was asked."""
    if "Unknown tag" in text:
        raise Refused("config key not recognised (Unknown tag)")
    if "generated outside the geometry" in text:
        raise Refused("primaries generated outside the geometry")
    # Every line carrying the label counts, then the whole field must be one integer: a regex on the digits alone read
    # "100000garbage" and "100000.5" as 100000 and skipped a malformed labelled line beside a valid one.
    label = "Nbr primaries simulated:"
    lines = [line for line in text.splitlines() if label in line]
    if len(lines) != 1:
        raise Refused(f"{len(lines)} '{label[:-1]}' lines")
    before, value = (part.strip() for part in lines[0].split(label, 1))
    if before or not value.isdigit():
        raise Refused(f"'{label[:-1]}' line is not one whole count: {lines[0].strip()!r}")
    simulated = int(value)
    if simulated < requested:
        raise Refused(f"simulated {simulated} < requested {requested}")
    return simulated


def parse_ct(text: str) -> tuple[list[int], list[float]]:
    """DimSize and ElementSpacing from an MHD header."""
    header = dict(line.split("=", 1) for line in text.splitlines() if "=" in line)
    header = {k.strip(): v.strip() for k, v in header.items()}
    try:
        dims = [int(x) for x in header["DimSize"].split()]
        spacing = [float(x) for x in header["ElementSpacing"].split()]
    except (KeyError, ValueError) as e:
        raise Refused(f"CT header lacks valid DimSize/ElementSpacing: {e}") from e
    if len(dims) != 3 or len(spacing) != 3 or min(dims) <= 0 or min(spacing) <= 0:
        raise Refused(f"CT geometry invalid: {dims} {spacing}")
    return dims, spacing


def record(a: argparse.Namespace) -> dict:
    """The measurement record, or Refused."""
    requested = int(float(a.requested))
    if requested <= 0 or float(a.requested) != requested:
        raise Refused(f"requested must be a positive whole count, got {a.requested!r}")
    rss, footprint = parse_time(Path(a.time).read_text())
    simulated = parse_log(Path(a.log).read_text(), requested)
    dims, spacing = parse_ct(Path(a.ct).read_text())
    return {
        "case": a.case,
        "material": a.material,
        "ct_dims": dims,
        "ct_spacing_mm": spacing,
        "voxels": dims[0] * dims[1] * dims[2],
        "threads": a.threads,
        "requested": requested,
        "simulated": simulated,
        "macos_max_rss_bytes": rss,
        "macos_peak_footprint_bytes": footprint,
        "binary_sha256": a.binary_sha256,
        "os": a.os,
        "scope": "one short run on this host and OS; not a general peak for longer runs or other modes, "
        "and not a measurement of any other OS",
    }


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name in ("--case", "--time", "--log", "--ct", "--requested", "--binary-sha256", "--material"):
        p.add_argument(name, required=True)
    p.add_argument("--threads", type=int, required=True)
    p.add_argument("--os", default="")
    a = p.parse_args(argv)
    try:
        rec = record(a)
    except (Refused, OSError) as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 1
    print(json.dumps(rec, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
