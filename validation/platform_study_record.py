"""Write one seed's provenance record for the cross-platform study and print it as a STUDY_RESULT line.

Usage: python platform_study_record.py <retain dir> <output dir> <platform> <seed> <binary> <compiler text> <commit>

Optional environment variable STUDY_ID names the study (default "pe1", the original cross-platform study, whose records
are unchanged); other studies that reuse this writer set it so their records cannot pass as pe1. record.json also holds
cfg_sha256 (cfg.txt), metrics_sha256 (the endpoints) and endpoint_status: "ok", or "error" with endpoint_error when the
endpoint computation raised. An endpoint failure still writes a record (metrics null); it is distinct from transport
status, which this script does not see.

Run from the repository root after the simulation. Computes the design's 13 endpoints from <output dir>/Dose.raw,
sha256 of every input and output the design lists, and writes <retain dir>/record.json. The retain dir must already
hold the copied outputs; this script never deletes or overwrites anything but record.json in that new directory.
"""

import glob
import hashlib
import json
import os
import platform as pyplatform
import re
import socket
import sys

from field_edge_analyse import load
from platform_study_metrics import endpoints

SIDE_MM = 150.0
DEFAULT_STUDY = "pe1"
STUDY_ID_FORM = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def study_id() -> str:
    """STUDY_ID from the environment (default pe1); refuses an empty or oddly shaped value."""
    value = os.environ.get("STUDY_ID", DEFAULT_STUDY)
    if not STUDY_ID_FORM.match(value):
        raise SystemExit(f"STUDY_ID {value!r} is not a plain identifier")
    return value


def metrics_digest(metrics: dict) -> str:
    return hashlib.sha256(json.dumps(metrics, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def compute_endpoints(outdir: str) -> dict:
    """The record's endpoint fields: metrics, metrics_sha256, endpoint_status (+ endpoint_error on failure)."""
    try:
        metrics = endpoints(load(outdir), SIDE_MM)
        digest = metrics_digest(metrics)
    except Exception as exc:  # noqa: BLE001 - any failure is recorded, never lost
        print(f"endpoint computation failed: {exc!r}", file=sys.stderr)
        return {"metrics": None, "metrics_sha256": None, "endpoint_status": "error",
                "endpoint_error": f"{type(exc).__name__}: {exc}"}
    return {"metrics": metrics, "metrics_sha256": digest, "endpoint_status": "ok"}


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def tree_digest(pattern: str) -> dict:
    files = sorted(p for p in glob.glob(pattern, recursive=True) if os.path.isfile(p))
    per = {p.replace(os.sep, "/"): sha256(p) for p in files}
    combined = hashlib.sha256("".join(f"{k} {v}\n" for k, v in per.items()).encode()).hexdigest()
    return {"files": len(per), "combined_sha256": combined}


def main(retain: str, outdir: str, plat: str, seed: str, binary: str, compiler: str, commit: str) -> None:
    rec = {
        "study": study_id(),
        "platform": plat,
        "seed": int(seed),
        "host": socket.gethostname(),
        "machine": pyplatform.machine(),
        "commit": commit,
        "compiler": compiler.strip(),
        "sha256": {
            "Dose.raw": sha256(os.path.join(outdir, "Dose.raw")),
            "Dose.mhd": sha256(os.path.join(outdir, "Dose.mhd")),
            "config": sha256("cfg.txt"),
            "plan E200_S150.txt": sha256("E200_S150.txt"),
            "cube.mhd": sha256("cube.mhd"),
            "cube.raw": sha256("cube.raw"),
            "BDL": sha256("BDL/BDL_default_DN_RangeShifter.txt"),
            "HU_Density": sha256("Scanners/default/HU_Density_Conversion.txt"),
            "HU_Material": sha256("Scanners/default/HU_Material_Conversion.txt"),
            "binary": sha256(binary),
        },
        "materials": tree_digest("Materials/**/*"),
        "cfg_sha256": sha256("cfg.txt"),
        **compute_endpoints(outdir),
    }
    path = os.path.join(retain, "record.json")
    if os.path.exists(path):
        raise SystemExit(f"{path} exists: not overwriting")
    with open(path, "w") as f:
        json.dump(rec, f, indent=1, sort_keys=True)
    print("STUDY_RESULT " + json.dumps(rec, sort_keys=True))


if __name__ == "__main__":
    main(*sys.argv[1:8])
