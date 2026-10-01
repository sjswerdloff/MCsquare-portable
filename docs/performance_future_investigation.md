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

1. **Same-arithmetic check, for the deterministic path.** On **1 thread**, where runs are seed-reproducible (#39
   affects >1 thread only), require **bit-identical** `Dose.raw` against the current build for several seeds and
   energies. A change is in this class only if it actually MEETS that criterion. No optimisation is assumed to
   preserve IEEE results because of its name; LTO and PGO usually do, but the test decides.
2. **Production-thread checks, for anything touching threading or work distribution.** A 1-thread pass is not
   sufficient here, because the changed behaviour only exists with >1 thread. At the intended thread count require:
   - actual primary accounting: simulated count against requested, with no lost or duplicated work;
   - an account of how RNG streams are constructed and in what order draws happen, and what that changes;
   - a multi-run statistical comparison against the current build.
3. **Floating-point-changing work** (vectorised maths, NEON intrinsics with different rounding): bitwise identity isn't
   possible. Such work needs **change-specific accuracy requirements, justified and agreed BEFORE any speed claim**,
   with an uncertainty and acceptance decision for each quantity it could affect. **The reference bands in
   `validation/topas_design.md` are reporting references for the TOPAS comparison, not acceptance criteria for an
   optimisation**: a difference inside those bands does not establish that no accuracy was lost. Study plots against
   the bands stay descriptive. These changes are physics-adjacent and need sjswerdloff's approval under the
   no-physics-changes rule for portable. The review and evidence have to be for the production configuration (build,
   flags, threads), even when a 1-thread control passes.
4. **Excluded by policy:** `-ffast-math` and similar flags, because they relax floating-point semantics (ordering,
   NaN/inf handling, reassociation) that the transport code relies on.

**Timing baselines are recorded with their conditions:** compiler and version, flags, binary sha256, host, thread
count, inputs (case and config hashes), and the actual simulated history count. Otherwise an apparent speed-up can be
a workload change: #39's shared counter can overshoot N, and an idle host is not a loaded one.

## Candidate levers, cheapest first

| # | lever | scope | expected | accuracy class |
|---|---|---|---|---|
| 0 | clean baseline on an idle Studio, then profile (`sample`/Instruments, clang `-Rpass=loop-vectorize` / `-Rpass-missed`) | all | tells us where time goes | none, measurement only |
| 1 | link-time optimisation (`-flto`) and profile-guided optimisation (`-fprofile-generate` / `-fprofile-use`) | any compiler/host | heuristic only: often 10–30% in other codes, not measured here | 1, if the 1-thread bitwise criterion is actually met |
| 2 | hand out primaries in blocks instead of one atomic increment each | all hosts | unknown until profiled | 2: production-thread checks (1-thread bitwise is not enough) |
| 3 | vectorised maths: SLEEF via clang `-fveclib=SLEEFGNUABI` (needs libsleef at link time), or Apple Accelerate/vForce | SLEEF: any ARM or x86; Accelerate: macOS | the most plausible route to 2× **if** maths calls dominate (heuristic) | 3: change-specific criteria, needs approval |
| 4 | NEON intrinsics for one or two hot kernels the compiler fails on | any AArch64 | only if the profile shows it | 3 |
| 5 | Metal GPU port | Apple only | heuristic: other GPU Monte Carlo codes report 10×-class gains; nothing measured for this engine | a port, not an optimisation; full validation |

Hand-written assembly isn't listed separately. Intrinsics get the same vector instructions and stay reviewable, and
assembly is only worth it if a profiled kernel resists both the compiler and intrinsics.

## Portability: ARM generally versus Apple only

- **Any AArch64** (Apple M-series, AWS Graviton, Ampere, Linux ARM): levers 1–4 with SLEEF, and host tuning via
  `-mcpu=native` as today.
- **Apple only:** Accelerate and Metal. Accelerate is also Apple's only supported route to the CPU-side AMX matrix units. Those are distinct from the GPU, which is reached through Metal, and from the Neural Engine, which isn't relevant here.
- **ARM servers only:** SVE (Graviton/Neoverse). The M1–M3 chips don't have it, so it can't be developed or tested on
  the Studio.
- **Recommendation:** target ARM generally (NEON + SLEEF + LTO/PGO) and keep Apple-only paths optional, added only if a
  profile justifies them. The x86 builds benefit from levers 1 and 2 as well.

## Suggested first step, when this is picked up

A clean, recorded baseline and a profile on the Studio, then lever 1 checked against guardrail 1. That shows where
the time goes, and whether lever 1 meets the bitwise criterion and what it buys. It does not by itself say whether 2×
is reachable. None of this may change the build used by an in-progress validation study: the TOPAS comparison
(#32/#33) pins its own build.
