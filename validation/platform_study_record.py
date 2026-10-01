"""Write one seed's provenance record for the cross-platform study and print it as a STUDY_RESULT line.

Usage: python platform_study_record.py <retain dir> <output dir> <platform> <seed> <binary> <compiler text> <commit>

Run from the repository root after the simulation. Computes the design's 13 endpoints from <output dir>/Dose.raw,
sha256 of every input and output the design lists, and writes <retain dir>/record.json. The retain dir must already
hold the copied outputs; this script never deletes or overwrites anything but record.json in that new directory.
"""

import glob
import hashlib
import json
import os
import platform as pyplatform
import socket
import sys

from platform_study_metrics import endpoints
from field_edge_analyse import load

SIDE_MM = 150.0


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
        "study": "pe1",
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
        "metrics": endpoints(load(outdir), SIDE_MM),
    }
    path = os.path.join(retain, "record.json")
    if os.path.exists(path):
        raise SystemExit(f"{path} exists: not overwriting")
    with open(path, "w") as f:
        json.dump(rec, f, indent=1, sort_keys=True)
    print("STUDY_RESULT " + json.dumps(rec, sort_keys=True))


if __name__ == "__main__":
    main(*sys.argv[1:8])
