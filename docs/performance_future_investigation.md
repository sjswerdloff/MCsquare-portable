# Performance: a possible future investigation (not started)

**Status: parked. Correctness comes first** (sjswerdloff, 2026-10-01). Nothing here has been profiled, built or
measured beyond the facts marked as measured. It records what we know and what we'd try, so the work can start from
evidence when it's wanted. The question asked was whether portable MCsquare could run 2× or faster, through compiler
flags or specialised code up to assembly, without losing accuracy. The honest answer is **possibly, and unknown until
profiled**.

## What we know (measured or read from source, 2026-10-01)

- **Build:** `MCsquare_arm64` is `clang -O3` plus `EXTRA_CFLAGS` (`-mcpu=native` on the Studio), with libomp. There's no
  link-time optimisation, no profile-guided optimisation and no vectorised maths library.
- **Precision:** already single precision (`VAR_COMPUTE`, `VAR_DATA` and `VAR_SCORING` are float, `src/include/define.h`),
  so there's no easy precision lever.
- **Written for Intel's compiler:** `VLENGTH 16` with `#pragma omp simd` loops, sized for AVX-512. NEON holds 4 floats.
  Upstream relied on icc and Intel's vector maths (SVML) for exp/log/sqrt/sin/cos inside those loops. **Unconfirmed
  suspicion:** on arm64, clang doesn't vectorise some of these loops, and the maths calls in them are scalar libm calls.
- **Throughput (planning runs, 200 MeV, 1e7 primaries):**
  - portable on the Mac Studio at 24 threads: 49–60 s in one set (29e970c2), 124.6 s in another (ba204d51);
  - upstream 85bf2911 + fix, icl on the Lenovo at 4 threads: about 223 s.
  - **The 2.5× spread between Studio sets is unexplained** (possibly other load), so there's no stable baseline yet.
- **Threading:** primaries go to threads through one shared atomic counter (#39), which may contend at 24 threads.

## Accuracy guardrails (must come before any speed claim)

1. **Changes that keep the arithmetic the same** (link-time and profile-guided optimisation, threading and
   work-distribution changes that don't alter per-history arithmetic): on **1 thread**, where runs are seed-reproducible
   (#39 affects >1 thread only), require **bit-identical** `Dose.raw` against the current build for several seeds and
   energies. Pass/fail, no judgement.
2. **Changes that alter floating-point results** (vectorised maths, NEON intrinsics with different rounding): bitwise
   identity isn't possible. They need the statistical comparison this repo already has: R80, σ and ring energy
   fractions against the current build over multiple seeds, reported against the reference bands in
   `validation/topas_design.md`. **These are physics-adjacent and need sjswerdloff's approval** under the
   no-physics-changes rule for portable.
3. **Excluded:** `-ffast-math` and similar. It licenses reordering and dropping operations in ways that can't be
   audited for medical-grade use.

## Candidate levers, cheapest first

| # | lever | scope | expected | accuracy class |
|---|---|---|---|---|
| 0 | clean baseline on an idle Studio, then profile (`sample`/Instruments, clang `-Rpass=loop-vectorize` / `-Rpass-missed`) | all | tells us where time goes | none, measurement only |
| 1 | link-time optimisation (`-flto`) and profile-guided optimisation (`-fprofile-generate` / `-fprofile-use`) | any compiler/host | typically 10–30% | 1: bit-identical at 1 thread |
| 2 | hand out primaries in blocks instead of one atomic increment each | all hosts | unknown until profiled | 1, if per-history arithmetic is unchanged |
| 3 | vectorised maths: SLEEF via clang `-fveclib`, or Apple Accelerate/vForce | SLEEF: any ARM or x86; Accelerate: macOS | the most plausible route to 2× **if** maths calls dominate | 2: statistical, needs approval |
| 4 | NEON intrinsics for one or two hot kernels the compiler fails on | any AArch64 | only if the profile shows it | 2 |
| 5 | Metal GPU port | Apple only | where 10×-class gains live for this kind of code | a port, not an optimisation; full validation |

Hand-written assembly isn't listed separately. Intrinsics get the same vector instructions and stay reviewable, and
assembly is only worth it if a profiled kernel resists both the compiler and intrinsics.

## Portability: ARM generally versus Apple only

- **Any AArch64** (Apple M-series, AWS Graviton, Ampere, Linux ARM): levers 1–4 with SLEEF, and host tuning via
  `-mcpu=native` as today.
- **Apple only:** Accelerate (also the only supported route to the M-series matrix units) and Metal.
- **ARM servers only:** SVE (Graviton/Neoverse). The M1–M3 chips don't have it, so it can't be developed or tested on
  the Studio.
- **Recommendation:** target ARM generally (NEON + SLEEF + LTO/PGO) and keep Apple-only paths optional, added only if a
  profile justifies them. The x86 builds benefit from levers 1 and 2 as well.

## Suggested first step, when this is picked up

A clean baseline and a profile on the Studio, then lever 1 behind the 1-thread bit-identical check. That shows within
about a day whether 2× is plausible without changing the arithmetic. It must not change the build used by an
in-progress validation study; the TOPAS comparison (#32/#33) pins its own build.
