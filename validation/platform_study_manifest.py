"""Generate the FROZEN EXPECTED-INPUT MANIFEST for the cross-platform study (platform_study_manifest.json).

Usage: python platform_study_manifest.py <repo path> [<output json>]     (needs numpy and git; default output is
platform_study_manifest.json next to this file; refuses to overwrite a manifest whose content would differ)

The analyzer (platform_study_analyse.py) refuses any record whose input hashes are not one of the values computed here
from the frozen study commit. Nothing here reads a STUDY_RESULT. Everything is computed from git objects at
STUDY_COMMIT (checked identical at RERUN_COMMIT), plus the output of validation/field_edge_make_cases.py run from
that commit's own copy of the script.

WHAT EACH EXPECTATION ASSUMES (each is also written into the manifest's "assumptions" list)
  text inputs from git (BDL, HU_Density, HU_Material)
    LF variant    = the git blob bytes (Linux and macOS checkouts, no eol conversion; no .gitattributes exists at the commit).
    CRLF variant  = every LF not preceded by CR becomes CRLF, applied only to blobs without a NUL byte (the reviewer's
                    rule for "git would treat as text"). The generator also evaluates git's own binary heuristic
                    (NUL, lone CR, or more non-printables than a 128th of the printables) and REPORTS every file on which
                    it disagrees with the NUL-only rule, and every file that already holds CR; a nonzero count is an
                    ambiguity that stays in the manifest.
    Which variant a Windows runner produces (autocrlf true or false) is not knowable from here, so Windows accepts both.
  generated inputs (plan E200_S150.txt, cube.mhd, cube.raw)
    Produced by running `field_edge_make_cases.py 200 150` from the study commit, in an empty directory, with FE_N
    unset (the workflow sets none). Run twice in separate directories; the manifest is written only if both runs agree.
    LF variant is what a POSIX Python writes. CRLF variant is what Python text mode writes on Windows
    (os.linesep); the two text files are written in text mode, so on Windows CRLF is the expected variant and LF is
    accepted only because the reviewer's rule says so. cube.raw is binary (numpy tofile, little-endian float32
    zeros; cube.mhd says ElementByteOrderMSB = False) and has exactly one expected hash on every platform.
    ASSUMPTION: the numpy of each platform formats/generates identical floats (np.arange with a float step and "%f"
    formatting); the manifest cannot check this for platforms other than the one it is generated on.
  Materials tree digest (record.tree_digest over Materials/**/*)
    Same path strings, same "<path> <sha256>\\n" lines, same combined sha256. The record sorts the OS-form paths BEFORE
    replacing separators. On Windows the OS form uses backslashes, which can change the order (e.g. "A/x" sorts
    before "A2/x" but "A\\x" sorts after "A2\\x"). The generator computes BOTH orders for both eol variants and reports
    whether they differ for this tree. Linux/macOS accept the posix order; Windows accepts both orders (whether the
    MSYS2 Python's os.sep is a backslash is not verifiable from here; if the two orders coincide the question is moot).
    Every Materials file (including *.dat~ backups) is included; the tree has no dotfiles (globbing would skip them),
    no symlinks and no empty directories at the commit; the generator refuses if that stops being true.
  config (cfg.txt), per (platform, seed)
    Reconstructed from each platform's workflow template at the study commit: Linux/macOS `printf '%s\\n' ... > cfg.txt`
    (LF after every line, including the last); Windows `@(...) | Set-Content -Encoding Ascii cfg.txt` (CRLF after every
    element including the last, no BOM: true of Windows PowerShell 5.1 and PowerShell 7 with -Encoding Ascii, an
    assumption about the runner's PowerShell that cannot be checked from here). PRIMARIES is read from the workflow's
    env block. Seeds cover BOTH waves (the manifest must be frozen before look 1 and cannot change afterwards);
    ASSUMPTION: wave 2 is run from a workflow whose cfg template is identical to this one. The re-run workflow is
    parsed the same way and must produce identical bytes for its two seeds, or the generator refuses.
    NOT REPRODUCIBLE: a template whose Output_Directory depends on GITHUB_RUN_ID cannot be frozen; this generator
    refuses any template that is not the literal expected line set.
  binary
    No expectation: the analyzer requires one hash per platform and prints it.
"""

import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

STUDY_COMMIT = "7d07db0145671924428936d5a7fa896d8c3f8e22"
RERUN_COMMIT = "600c66e0fe19adaaee49e51fa023edf89c47f299"
WAVE_SEEDS = {1: {"windows": range(1001, 1017), "linux": range(2001, 2017), "macos": range(3001, 3017)},
              2: {"windows": range(1017, 1033), "linux": range(2017, 2033), "macos": range(3017, 3033)}}
RERUN_SEEDS = {("windows", 1008), ("linux", 2003)}
TEXT_PATHS = {"BDL": "BDL/BDL_default_DN_RangeShifter.txt",
              "HU_Density": "Scanners/default/HU_Density_Conversion.txt",
              "HU_Material": "Scanners/default/HU_Material_Conversion.txt"}
GENERATED = {"plan E200_S150.txt": "E200_S150.txt", "cube.mhd": "cube.mhd", "cube.raw": "cube.raw"}
MAKE_CASES = "validation/field_edge_make_cases.py"
WORKFLOW = ".gitea/workflows/platform-study.yml"
RERUN_WORKFLOW = ".gitea/workflows/platform-study-rerun.yml"
CFG_TEMPLATE = ["Num_Threads 3", "Num_Primaries 3e7", "RNG_Seed {SEED}", "CT_File cube.mhd",
                "HU_Density_Conversion_File Scanners/default/HU_Density_Conversion.txt",
                "HU_Material_Conversion_File Scanners/default/HU_Material_Conversion.txt",
                "BDL_Machine_Parameter_File BDL/BDL_default_DN_RangeShifter.txt",
                "BDL_Plan_File E200_S150.txt", "Output_Directory out_seed", "Dose_MHD_Output True"]


class ManifestError(Exception):
    """The manifest cannot be produced faithfully; nothing is written."""


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def to_crlf(data: bytes) -> bytes:
    """Every LF not preceded by CR becomes CRLF; blobs holding a NUL byte are left alone (treated as binary)."""
    if b"\0" in data:
        return data
    return re.sub(rb"(?<!\r)\n", b"\r\n", data)


def git_says_binary(data: bytes) -> bool:
    """git's own auto-detection (convert.c gather_stats/convert_is_binary): lone CR, NUL, or too many non-printables."""
    nul = lonecr = printable = nonprintable = 0
    i = 0
    while i < len(data):
        c = data[i]
        if c == 0x0D:
            if i + 1 < len(data) and data[i + 1] == 0x0A:
                printable += 1
                i += 1
            else:
                lonecr += 1
        elif c == 0x0A:
            printable += 1
        elif c == 0x08 or c == 0x09 or c == 0x1B or c == 0x0C:
            printable += 1
        elif c < 32 or c == 127:
            if c == 0:
                nul += 1
            nonprintable += 1
        else:
            printable += 1
        i += 1
    return bool(lonecr or nul or (printable >> 7) < nonprintable)


def tree_digest(files: dict[str, bytes], order: str) -> str:
    """platform_study_record.tree_digest(...)["combined_sha256"] for these files (keys are posix paths)."""
    key = (lambda p: p) if order == "posix" else (lambda p: p.replace("/", "\\"))
    per = [(p, sha256(files[p])) for p in sorted(files, key=key)]
    return sha256("".join(f"{k} {v}\n" for k, v in per).encode())


def parse_cfg_template(text: str, platform: str) -> list[str]:
    """The cfg.txt lines a workflow writes for `platform`, PRIMARIES substituted and the seed left as {SEED}."""
    prim = re.search(r'^\s+PRIMARIES:\s*"([^"]+)"\s*$', text, re.M)
    if not prim:
        raise ManifestError("no PRIMARIES in the workflow env block")
    job = re.search(rf"^  {platform}:\n(.*?)(?=^  \w+:\n|\Z)", text, re.M | re.S)
    if not job:
        raise ManifestError(f"no job {platform!r} in the workflow")
    if platform == "windows":
        m = re.search(r"@\((.*?)\)\s*\|\s*\n\s*Set-Content -Encoding Ascii cfg\.txt", job.group(1), re.S)
        if not m:
            raise ManifestError("windows cfg block not found (expected `@(...) | Set-Content -Encoding Ascii cfg.txt`)")
        toks = re.findall(r'"([^"]*)"', m.group(1))
        toks = [t.replace("$env:PRIMARIES", prim.group(1)).replace("$seed", "{SEED}") for t in toks]
    else:
        m = re.search(r"printf '%s\\n' (.*?) > cfg\.txt", job.group(1), re.S)
        if not m:
            raise ManifestError(f"{platform} cfg block not found (expected `printf '%s\\n' ... > cfg.txt`)")
        toks = shlex.split(re.sub(r"\\\n\s*", " ", m.group(1)))
        toks = [t.replace("$PRIMARIES", prim.group(1)).replace("$SEED", "{SEED}") for t in toks]
    if any("$" in t for t in toks):
        raise ManifestError(f"{platform}: unresolved variable in cfg template {toks}")
    return toks


def cfg_bytes(lines: list[str], platform: str, seed: int) -> bytes:
    nl = "\r\n" if platform == "windows" else "\n"
    return "".join(ln.replace("{SEED}", str(seed)) + nl for ln in lines).encode("ascii")


def git(repo: str, *args: str, text: bool = False) -> bytes | str:
    exe = shutil.which("git")
    if not exe:
        raise ManifestError("git not found")
    p = subprocess.run([exe, "-C", repo, *args], capture_output=True, check=False)
    if p.returncode:
        raise ManifestError(f"git {' '.join(args)} failed: {p.stderr.decode(errors='replace').strip()}")
    return p.stdout.decode() if text else p.stdout


def generate_inputs(script: bytes) -> dict[str, bytes]:
    """Run the make-cases script (bytes from the study commit) twice in empty directories; return its outputs."""
    runs = []
    for _ in range(2):
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / "make_cases.py").write_bytes(script)
            env = {k: v for k, v in os.environ.items() if k != "FE_N"}
            p = subprocess.run([sys.executable, "make_cases.py", "200", "150"], cwd=d, env=env,
                               capture_output=True, text=True, check=False)
            if p.returncode:
                raise ManifestError(f"field_edge_make_cases.py failed: {p.stderr.strip()}")
            names = sorted(x.name for x in Path(d).iterdir() if x.name != "make_cases.py")
            if names != sorted(GENERATED.values()):
                raise ManifestError(f"make_cases wrote {names}, expected {sorted(GENERATED.values())}")
            runs.append({n: (Path(d) / n).read_bytes() for n in names})
    if runs[0] != runs[1]:
        raise ManifestError("field_edge_make_cases.py output differs between two runs: not deterministic")
    return runs[0]


def build_manifest(repo: str) -> dict:
    same = ["BDL", "Scanners", "Materials", MAKE_CASES]
    for path in same:
        a, b = (git(repo, "rev-parse", f"{c}:{path}", text=True).strip() for c in (STUDY_COMMIT, RERUN_COMMIT))
        if a != b:
            raise ManifestError(f"{path} differs between the study commit and the re-run commit")
    wf_ids = [git(repo, "rev-parse", f"{c}:{WORKFLOW}", text=True) for c in (STUDY_COMMIT, RERUN_COMMIT)]
    if wf_ids[0] != wf_ids[1]:
        raise ManifestError("platform-study.yml differs between the study commit and the re-run commit")
    counts: dict = {}
    text = {}
    for name, path in TEXT_PATHS.items():
        blob = git(repo, "show", f"{STUDY_COMMIT}:{path}")
        if b"\0" in blob:
            raise ManifestError(f"{path} holds a NUL byte: the text/binary rule cannot classify it as text")
        text[name] = {"LF": sha256(blob), "CRLF": sha256(to_crlf(blob))}
        counts[name] = {"bytes": len(blob), "lf": blob.count(b"\n"), "cr": blob.count(b"\r"),
                        "git_says_binary": git_says_binary(blob)}
    listing = git(repo, "ls-tree", "-r", "-z", STUDY_COMMIT, "--", "Materials").split(b"\0")
    files: dict[str, bytes] = {}
    for entry in filter(None, listing):
        meta, path = entry.decode().split("\t", 1)
        mode, kind, _ = meta.split(" ")
        if mode != "100644" or kind != "blob":
            raise ManifestError(f"Materials entry {path} is {mode} {kind}: only plain files are modelled")
        if "/." in path:
            raise ManifestError(f"Materials entry {path} is a dotfile: glob would skip it")
        files[path] = git(repo, "show", f"{STUDY_COMMIT}:{path}")
    crlf_files = {p: to_crlf(b) for p, b in files.items()}
    mats = {eol: {order: tree_digest(fs, order) for order in ("posix", "windows")}
            for eol, fs in (("LF", files), ("CRLF", crlf_files))}
    binary = [p for p, b in files.items() if b"\0" in b]
    counts["materials"] = {
        "files": len(files), "with_NUL_left_unconverted": len(binary),
        "with_any_LF_so_CRLF_differs": sum(1 for p in files if crlf_files[p] != files[p]),
        "already_holding_CR": sum(1 for b in files.values() if b"\r" in b),
        "git_binary_heuristic_disagrees_with_NUL_rule": sum(1 for b in files.values() if git_says_binary(b) != (b"\0" in b)),
        "posix_order_equals_windows_order": mats["LF"]["posix"] == mats["LF"]["windows"],
    }
    gen = generate_inputs(git(repo, "show", f"{STUDY_COMMIT}:{MAKE_CASES}"))
    for name, fname in GENERATED.items():
        data = gen[fname]
        if name == "cube.raw":
            continue
        if b"\r" in data or b"\0" in data:
            raise ManifestError(f"generated {fname} unexpectedly holds CR or NUL")
        text[name] = {"LF": sha256(data), "CRLF": sha256(to_crlf(data))}
        counts[name] = {"bytes": len(data), "lf": data.count(b"\n")}
    counts["cube.raw"] = {"bytes": len(gen["cube.raw"])}
    templates = {}
    study_wf = git(repo, "show", f"{STUDY_COMMIT}:{WORKFLOW}", text=True)
    rerun_wf = git(repo, "show", f"{RERUN_COMMIT}:{RERUN_WORKFLOW}", text=True)
    for plat in ("linux", "macos", "windows"):
        templates[plat] = parse_cfg_template(study_wf, plat)
        if templates[plat] != CFG_TEMPLATE:
            raise ManifestError(f"{plat} cfg template is not the literal expected one: {templates[plat]}")
    for plat, seed in sorted(RERUN_SEEDS):
        if cfg_bytes(parse_cfg_template(rerun_wf, plat), plat, seed) != cfg_bytes(templates[plat], plat, seed):
            raise ManifestError(f"re-run workflow writes a different cfg.txt for {plat} seed {seed}")
    config = {plat: {str(s): sha256(cfg_bytes(templates[plat], plat, s))
                     for w in WAVE_SEEDS.values() for s in w[plat]} for plat in templates}
    assumptions = [
        "LF variant = git blob bytes; CRLF variant = LF not preceded by CR -> CRLF, only for blobs without NUL",
        "Linux/macOS accept LF only; Windows accepts LF or CRLF for text inputs and Materials (any autocrlf setting)",
        "Windows Materials digest accepted in posix and in backslash-sorted order (os.sep of MSYS2 python unverifiable here)",
        "generated plan/cube.mhd: POSIX python writes LF, Windows text mode CRLF; numpy formats identically everywhere",
        "cube.raw: one exact hash on all platforms (little-endian float32 zeros)",
        "cfg.txt: Linux/macOS LF after every line; Windows Set-Content -Encoding Ascii = CRLF after every line, no BOM",
        "cfg.txt seeds cover both waves and assume wave 2 uses this same template (Output_Directory out_seed)",
        "no expectation for the binary: one hash per platform is required and printed",
    ]
    return {"study_commit": STUDY_COMMIT, "rerun_commit": RERUN_COMMIT, "generator": "validation/platform_study_manifest.py",
            "text": text, "cube.raw": sha256(gen["cube.raw"]), "materials": mats, "config": config,
            "config_template": CFG_TEMPLATE, "counts": counts, "assumptions": assumptions,
            "not_reproducible": ["Windows checkout eol behaviour and MSYS2 python os.sep are assumed, not observed",
                                 "a cfg.txt whose Output_Directory embeds GITHUB_RUN_ID cannot be frozen"]}


def main(repo: str, out: str) -> None:
    text = json.dumps(build_manifest(repo), indent=1, sort_keys=True) + "\n"
    path = Path(out)
    if path.exists():
        if path.read_text() != text:
            raise SystemExit(f"{path} exists with different content: not overwriting (delete it deliberately to regenerate)")
        print(f"{path} already up to date")
        return
    path.write_text(text)
    print(f"wrote {path}")


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        raise SystemExit(__doc__.split("\n")[2])
    main(sys.argv[1], sys.argv[2] if len(sys.argv) == 3 else str(Path(__file__).with_name("platform_study_manifest.json")))
