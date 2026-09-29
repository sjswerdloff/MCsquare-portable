# MCsquare (portable)

Fast Monte Carlo dose calculation for pencil-beam scanning proton therapy. This is a portable fork of
[OpenMCsquare](https://gitlab.com/openmcsquare/MCsquare) that builds with open-source compilers and
runs natively on Apple Silicon.

**No Intel software is needed or used.** The Intel compiler and Intel MKL that upstream requires are not
used here: MKL's random number generator is replaced by PCG, and the MKL vector helpers by standard C.
Any C compiler with OpenMP builds the code.

## Platforms

| Platform | Compiler | Status |
|---|---|---|
| **macOS on Apple Silicon (arm64)** | Apple clang + Homebrew `libomp` | Primary target; built and tested in CI on every pull request |
| macOS on Apple Silicon (arm64) | Homebrew GCC (`gcc-15`) | Builds and runs the smoke test; not run in CI |
| Linux x86-64 | GCC with OpenMP | Builds and runs (benchmarked on Ubuntu 24.04); not run in CI |
| Windows | — | Not supported by the portable build. Upstream's own `Makefile.bat` needs Intel's compiler |

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

Plain `make` builds upstream's Intel targets, so always name a target.

## Build on Linux (x86-64)

```
sudo apt install build-essential    # gcc with OpenMP
make MCsquare_portable
```

## Run

Run the binary directly with a configuration file. Upstream's `MCsquare` launcher script selects the
Intel builds and is not used here.

```
./MCsquare_arm64 Sample_input_data/config.txt
```

Materials are read from `./Materials` if present, otherwise from the directory in the environment
variable `MCsquare_Materials_Dir`. A quick check with 1000 primaries:

```
mkdir -p smoke_test_output
./MCsquare_arm64 Sample_input_data/smoke_test_config.txt
```

## Tests

The C regression tests build with clang + libomp and run with UndefinedBehaviorSanitizer, as CI does:

```
make test_transport test_remove_tmp test_material_labels test_angle
```

## Differences from upstream

Apart from portability, this fork fixes defects found while porting. Each is recorded in an issue here:

- **Secondary emission angles (#16).** Upstream interpolates the angular table of nuclear-inelastic
  secondaries at an energy in eV between brackets in MeV, so the sampled angles do not follow the
  ICRU data. It reads past the end of the angle table for about one sample in four. The measured
  effect is that dose just outside a large field is underestimated. The fix has been reported upstream.
- **Out-of-range emission angle (#24).** An angle index outside the table now aborts with a message
  instead of producing an angle.
- Other upstream defects fixed:
  - #5: a buffer overflow in beamlet mode, because path buffers were too short.
  - #8: a ray that misses the CT read uninitialised data in `Transport_to_CT`.
  - #12: temporary folders were removed through `system()` with a path from the configuration; they are
    now removed without a shell.
  - #20: an inner loop in `CT_Transport_SPR` reused the outer lane index. That code is not compiled in
    default builds.

The physics of the corrected sampling has not yet been validated against measurement or a TOPAS/Geant4
reference.
