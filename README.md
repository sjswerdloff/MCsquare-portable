# MCsquare (portable)

Fast Monte Carlo dose calculation for pencil-beam scanning proton therapy. This is a portable fork of
[OpenMCsquare](https://gitlab.com/openmcsquare/MCsquare) (from upstream commit `85bf2911`) that builds with
free compilers on Linux, macOS on Apple Silicon and Windows.

**No Intel software is needed.** Upstream requires the Intel compiler and Intel MKL. Here Cilk Plus is
replaced by OpenMP SIMD, MKL's random number generator by PCG, and the MKL vector helpers by standard C.
The physics models, data tables and beam model are upstream's.

**This is research software. It is not a medical device and is not cleared for clinical use.**

## Status

- **Compared with upstream.** On the same hosts, and with the same fix applied to both (see below),
  this fork and Intel-built upstream agree within margins fixed in advance on range, spot size, the core
  and near halo of a pencil beam at 100, 150 and 200 MeV, and every endpoint of a broad field. Of 52
  endpoints, 46 are equivalent and none is not equivalent. The other six, all in the pencil beam's far halo
  (40 mm or more from the axis) at 100 and 150 MeV, are inconclusive or not established, so equivalence on
  all 52 together is not established.
- **Across platforms.** For one broad field, Windows and macOS are each equivalent to Linux on 13 of 13
  endpoints.
- **Compared with TOPAS/Geant4.** A descriptive comparison (no margins) shows a range offset of 0.42 to
  0.54 mm, most of which depends on the TOPAS step size, and less dose than TOPAS in the far halo at 100
  and 150 MeV. By inference both are common to this fork and upstream.
- **Not compared with measurement.**
- **Speed.** Built with gcc, this fork took 2.4 to 2.7 times as long as Intel-built upstream on the two
  Intel hosts used. Built with the Intel compiler it took about as long as upstream.

The methods, results, limitations and the data behind every number are in the draft report,
[`docs/validation_report_draft.md`](docs/validation_report_draft.md). It is a draft and has not been
peer reviewed.

## Platforms

| Platform | Compiler | Status |
|---|---|---|
| Linux x86-64 | GCC with OpenMP | Built and tested in CI on every pull request (Ubuntu 24.04) |
| macOS on Apple Silicon (arm64) | Apple clang + Homebrew `libomp` | Built and smoke-tested in CI on every pull request |
| macOS on Apple Silicon (arm64) | Homebrew GCC (`gcc-15`) | Builds and runs the smoke test; not run in CI |
| Windows x86-64 | MinGW-w64 GCC (MSYS2 UCRT64) | Built and smoke-tested on the development server's CI; not in the GitHub CI |

## Build on Linux (x86-64)

```
sudo apt install build-essential    # gcc with OpenMP
make MCsquare_portable
```

Plain `make` builds upstream's Intel targets, so always name a target.

## Build on Apple Silicon

**Xcode is not needed.** The Apple Command Line Tools are, for either compiler: they supply the macOS
SDK headers and the linker, which Homebrew's GCC also uses.

```
xcode-select --install          # Command Line Tools (skip if already installed)
```

**Option 1: Apple clang + libomp** (what CI uses):

```
brew install libomp
make MCsquare_arm64
```

**Option 2: GNU GCC** (uses GCC's own OpenMP runtime, `libgomp`):

```
brew install gcc
gcc-15 src/*.c -fopenmp -lm -O3 -DVERSION='"portable, gcc-15"' -o MCsquare_arm64
```

On macOS, plain `gcc` is Apple clang under another name and does not accept `-fopenmp`, so call the
versioned Homebrew binary (`gcc-15`, or whichever version `brew` installed).

## Build on Windows (x86-64)

Install [MSYS2](https://www.msys2.org), then in its **UCRT64** shell:

```
pacman -S --needed mingw-w64-ucrt-x86_64-gcc mingw-w64-ucrt-x86_64-libgomp make
make MCsquare_win_portable
```

This produces `MCsquare_win_portable.exe`. It uses GCC's OpenMP runtime, `libgomp-1.dll`, which MSYS2
provides only as a DLL. So run it from the UCRT64 shell, add `C:\msys64\ucrt64\bin` to `PATH`, or copy
the DLLs it needs from there next to the `.exe`. Microsoft's compiler (MSVC) is not supported: its C
OpenMP support stops at version 2.0.

## Run

Run the binary directly with a configuration file. Upstream's `MCsquare` launcher script selects the
Intel builds and is not used here. The binary is `MCsquare_portable` on Linux, `MCsquare_arm64` on Apple
Silicon and `MCsquare_win_portable.exe` on Windows:

```
./MCsquare_portable Sample_input_data/config.txt    # Linux
./MCsquare_arm64 Sample_input_data/config.txt       # Apple Silicon
./MCsquare_win_portable.exe Sample_input_data/config.txt     # Windows, MSYS2 UCRT64 shell
.\MCsquare_win_portable.exe Sample_input_data\config.txt     # Windows, PowerShell
```

Materials are read from `./Materials` if present, otherwise from the directory in the environment
variable `MCsquare_Materials_Dir`. A quick check with 1000 primaries (use `./MCsquare_arm64` on Apple
Silicon):

```
mkdir -p smoke_test_output
./MCsquare_portable Sample_input_data/smoke_test_config.txt
```

The configuration file, beam model (`BDL/`), CT calibration (`Scanners/`) and plan formats are
upstream's and are unchanged.

## Tests

The C regression tests run with UndefinedBehaviorSanitizer. They build with gcc on Linux and with clang +
libomp on macOS:

```
make test_transport test_remove_tmp test_material_labels test_angle
tests/test_long_output_path.sh ./MCsquare_portable
```

The validation and analysis tools under `validation/` have their own tests, run with
[uv](https://docs.astral.sh/uv/):

```
uv run --project validation --frozen --with scipy==1.18.1 python -m pytest validation/tests
```

CI (`.github/workflows/ci.yml`) runs all of these on Ubuntu 24.04, and the build and smoke test on macOS.
The workflows under `.gitea/` belong to the development server that produced the validation runs; they
are kept as the record of how each result was made and do not run here.

## Differences from upstream

Apart from portability, this fork fixes defects found while porting. Each has a regression test.

- **Secondary emission angles** (`make test_angle`). Upstream interpolates the angular table of
  nuclear-inelastic secondaries at an energy in eV between brackets in MeV, so the sampled angles do not
  follow the ICRU data, and the angle index can run one past the end of the table. The measured effect,
  for a 200 MeV, 15 × 15 cm field, is that upstream computes 13 to 18% less dose than the corrected code
  between 5 and 30 mm outside the field edge. Reported upstream as
  [OpenMCsquare work item 42](https://gitlab.com/openmcsquare/MCsquare/-/work_items/42), with the patch
  (`validation/topas/openmcsquare/fix_secondary_angle_energy.patch`). Whether the corrected dose is closer
  to physical dose has not been shown.
- **Out-of-range emission angle.** An angle index outside the table now aborts with a message instead
  of producing an angle.
- **Path buffers** (`tests/test_long_output_path.sh`). A buffer overflow in beamlet mode, because path
  buffers were too short.
- **Rays that miss the CT** (`make test_transport`). `Transport_to_CT` read uninitialised data.
- **Temporary folders** (`make test_remove_tmp`). They were removed through `system()` with a path from
  the configuration; they are now removed without a shell.
- **Material labels** (`make test_material_labels`). An inner loop in `CT_Transport_SPR` reused the
  outer lane index. That code is not reached in default builds, because `define.h` sets
  `InterfaceCrossing` to `VoxelInterface`.

## Licence and attribution

MCsquare was developed by Kevin Souris at Université catholique de Louvain (UCLouvain, Louvain-la-Neuve,
Belgium) in a collaboration with IBA s.a., and is released under the Apache 2.0 licence; see
[`LICENSE`](LICENSE), which also states upstream's attribution requirement. This fork keeps that licence.
Please cite the original work when using it:

- K. Souris, J. A. Lee, E. Sterpin, "Fast multipurpose Monte Carlo simulation for proton therapy using
  multi- and many-core CPU architectures", *Medical Physics* 43(4), 2016.

PCG random number generation is from the [PCG](https://www.pcg-random.org) minimal C implementation
(`src/pcg_basic.c`), Apache 2.0.

Contributors to this fork and to the validation are listed in section 9 of the draft report.
