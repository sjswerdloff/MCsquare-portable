# TOPAS arm inputs (issue #32, design in `../topas_design.md`)

The TOPAS side of the MCsquare-vs-TOPAS comparison. It runs on the Mac Studio against OpenTOPAS 4.3.0 / Geant4
11.3.2, installed under `/Users/stuartswerdloff/MonteCarlo/`.

| file | what it does |
|---|---|
| `stage1_base.txt` | Shared parameters: water phantom, 1 mm grid, beam axis on the voxel corner, Gaussian source, the non-EM physics modules, and the two scorers (`Dose` excludes neutron and gamma descendants; `DoseAll` is unfiltered). |
| `make_run.sh <E_MeV> <opt0\|opt4> <seed> <histories> <threads> <out_root>` | Writes one run directory, named for every per-run choice, containing `run.txt` and a frozen copy of the base. Refuses to reuse a directory. |
| `run_topas.sh <run_dir>` | Runs TOPAS with the required environment and writes `provenance.txt`: host, UTC times, exit code, wall time, sha256 of inputs and outputs, and the EM physics that actually ran (the "Use ICRU90 data" line, 0 for option 0 and 1 for option 4). Refuses to start below `MIN_FREE_GB` (default 100) free on the output volume. |

Study outputs go to `/Volumes/T7 Shield/MonteCarlo_runs/` (the home volume is nearly full). Launch heavy runs
detached, e.g. `nohup ./run_topas.sh <dir> > /dev/null 2>&1 < /dev/null &`. Never launch the `topas` binary itself through
`/usr/bin/time` or `env -i`: macOS strips `DYLD_*` from those, and topas then fails to load its libraries.
Wrapping `run_topas.sh` is safe, because it sets the environment inside the script (that is how the timings were taken).

Output layout: the binary is x-fastest (`reshape((400, 400, 350), order="F")` gives `v[ix, iy, kz]`), and
depth = 349.5 − k mm.
