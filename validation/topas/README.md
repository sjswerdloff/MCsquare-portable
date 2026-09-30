# TOPAS arm inputs (issue #32, design in `../topas_design.md`)

The TOPAS side of the MCsquare-vs-TOPAS comparison. It runs on the Mac Studio against OpenTOPAS 4.3.0 / Geant4
11.3.2, installed under `/Users/stuartswerdloff/MonteCarlo/`.

| file | what it does |
|---|---|
| `stage1_base.txt` | Shared parameters: water phantom, 1 mm grid, beam axis on the voxel corner, Gaussian source, the non-EM physics modules, and the two scorers (`Dose` excludes neutron and gamma descendants; `DoseAll` is unfiltered). |
| `make_run.sh <E_MeV> <opt0\|opt4> <seed> <histories> <threads> <out_root> [attempt]` | Writes one run directory, `E<E>_<em>_seed<s>_n<N>_th<T>_<attempt>`, containing `run.txt` and a frozen copy of the base. Validates ranges (energy 10–400 MeV; seed and histories 1–2147483647; threads 1–32; attempt a1–a999) and claims the directory atomically. Refuses to reuse one; a retry is a new attempt tag. |
| `run_topas.sh <run_dir>` | Runs one directory once, and fails closed. Claims the directory atomically and refuses any previous attempt or output state, so a finished run's records are never overwritten. Writes `provenance.txt`: host, UTC times, runner and executable paths with their sha256, sha256 of the inputs, TOPAS's own exit status, output sizes and sha256s, the logged EM physics against the requested arm, and a verdict. Exit 0 only when TOPAS exits 0, all four outputs have the size the grid implies, and the logged EM physics matches; otherwise 5 (incomplete), 6 (EM unverified) or 7 (TOPAS failed). It refuses (4) if `MIN_FREE_GB` (default 100, whole GB) is malformed or not met. |

Study outputs go to `/Volumes/T7 Shield/MonteCarlo_runs/` (the home volume is nearly full). Launch heavy runs
detached, e.g. `nohup ./run_topas.sh <dir> > /dev/null 2>&1 < /dev/null &`. Never launch the `topas` binary itself through
`/usr/bin/time` or `env -i`: macOS strips `DYLD_*` from those, and topas then fails to load its libraries.
Wrapping `run_topas.sh` is safe, because it sets the environment inside the script (that is how the timings were taken).

Output layout: the binary is x-fastest (`reshape((400, 400, 350), order="F")` gives `v[ix, iy, kz]`), and
depth = 349.5 − k mm.

Tests: `bats tests/test_wrappers.bats` (no transport; a fake TOPAS and a 2×3×4 grid) cover the wrapper contract:
valid and invalid inputs, repeat and concurrent claims, TOPAS failure, missing, empty or wrong-size outputs,
missing or mismatched EM physics, and malformed free-space overrides.
