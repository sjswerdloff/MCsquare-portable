# Same-host comparison programme (PLAN, for sjswerdloff's approval)

**Status: proposed 2026-10-01. Nothing here has run, and nothing here authorises a run.** Written to execute
unattended while sjswerdloff travels (week of 5 October 2026), if he approves it. Once a part is approved, its run
list, seeds, histories, analysis script and margins are committed **before** its first run, and nothing is changed
after results are read. Results feed [`validation_report_draft.md`](validation_report_draft.md).

## Host rulings (sjswerdloff, 2026-10-01, in conversation with connor-227743e6)

- **Windows:** for the same source code, results from the HP and the Lenovo are interchangeable, so where a Windows
  run takes place doesn't matter.
- **Load:** the work is divided roughly evenly between the two Intel machines. Linux runs only on the HP (WSL2), so
  when Linux work is needed the Windows load goes to the Lenovo. This supersedes his earlier "no portable MCsquare arm
  on the HP" (relayed by clement-7074f29f).
- **So: part A (Windows) runs on the Lenovo, part B (Linux) runs on the HP.** Each part then has one fixed host, and
  host is not confounded within any part-A or part-B contrast.
- **Compute jobs under two hours go ahead without a further ask.** Four jobs below are estimated at over 2 h and are
  marked; they wait for his go-ahead.

## What "same-host" means here, and what it doesn't

Each contrast compares **configuration bundles**, not code alone. Holding the host fixed removes hardware and OS from
a contrast; it does not make the code the only difference, because the compiler, RNG implementation, runtime and math
libraries and build flags change with it. The question each part answers is whether the two configurations give
**equivalent expected values of the chosen endpoints** on this case. It is not a claim of code identity, and it doesn't
attribute a difference to any single factor.

| part | contrast | matched | differs (confounded within the contrast) | host |
|---|---|---|---|---|
| **A. Intel, Windows 10** | upstream 85bf2911 + fix vs portable | OS, thread count (4), cases, histories, endpoint script, seed policy; host (Lenovo) | source tree; compiler (icl 2021.1 vs MSYS2 gcc 16.2); RNG (MKL vs PCG); OpenMP runtime; math library; build flags | Lenovo |
| **B. Intel, Linux** | upstream + fix (icc) vs portable (gcc 13) vs portable (icc) | host (HP, WSL2), threads (3), cases, histories, script | upstream-icc vs portable-icc: source tree and RNG, same compiler. portable-icc vs portable-gcc: compiler, runtime and flags, same source. Not a full factorial: source and RNG stay together | HP |
| **C. Apple Silicon vs TOPAS** (confirmatory) | portable vs TOPAS opt0 (opt4 secondary) | host (Mac Studio), case, histories, endpoint definitions | everything about the codes, per the frozen design (TD) | Mac Studio |
| **D. Bridge** | portable on the Lenovo (from A) vs portable on the Studio | source commit, case, histories, endpoint script | host, CPU architecture, OS, compiler (MSYS2 gcc vs Apple clang), runtime, threads (4 vs 24). A platform-bundle contrast, not platform alone | reuses A's and C's portable runs |

## Cases

- **P (pencil beam):** `validation/mcsquare_pencil_case.py` at 100, 150 and 200 MeV, Gaussian σ = 3 mm, water
  400 × 400 × 350 mm, 1 mm grid, nozzle 1 mm (TD). 70 MeV stays out (grid-limited, #43).
- **F (field edge):** the #31 platform-study case, 200 MeV 15 × 15 cm in the default scanner's HU 0 →
  **Schneider_AT_AG_SI4**, endpoints exactly as #31. **Prerequisite: #31 merged.**

## Endpoints, per energy (fixed now, before any run)

The current endpoint slabs (3, 100, 200 mm) only suit 200 MeV. At 100 MeV the primary range is about 77.5 mm and at
150 MeV about 158 mm (alden-ec2221c7, #45 review), so a 200 mm slab lies beyond the primary core at both, and 100 mm
does at 100 MeV. Slab depths are therefore set per energy at about half and about 0.8 of the primary range:

| energy | R80 | σ and ring slabs (mm) | rings (radius, mm) |
|---|---|---|---|
| 100 MeV | yes | 40, 60 | 5–10, 10–20, 20–40, 40–80, 80–200 at each slab |
| 150 MeV | yes | 80, 125 | same |
| 200 MeV | yes | 100, 200 | same |

That is 13 endpoints per energy (R80, σ at 2 slabs, 5 rings at 2 slabs), 39 for case P. Case F keeps #31's 13.
R20 and the central-axis dose are **not** endpoints for case P: `pencil_endpoints.py` doesn't compute them.

**Prerequisite code change:** `pencil_endpoints.py` takes its slab depths from a `--depths` argument (default 3, 100,
200, so existing results are unchanged). It is reviewed and merged before any part-A run.

## Pre-specified analysis (A, B if authorised, D)

- **Estimand:** for R80 and σ, the run-level mean difference (arm − reference); for ring fractions, the geometric-mean
  ratio. Point estimate, Welch SE and Welch–Satterthwaite df, pointwise 90% and 95% intervals, computed by the same
  method as `validation/report_tables.py`. Reference: upstream + fix (A, B); portable on the Studio (D).
- **Margins, proposed for sjswerdloff's decision, every endpoint listed:**

  | endpoint | margin |
  |---|---|
  | R80 | ±0.05 mm |
  | σ, each slab | ±0.02 mm |
  | rings 5–10, 10–20, 20–40, 40–80 mm, each slab | ratio [0.98, 1.02] |
  | ring 80–200 mm, each slab | ratio [0.95, 1.05] (wider: see sizing; **decided** by sjswerdloff 2026-10-01) |
  | case F | the #31 margins |

- **Test:** per endpoint, TOST at one-sided α = 0.05, i.e. the 90% interval inside the margin.
- **Claims and multiplicity, kept separate:**
  1. **Joint claim** (primary): "equivalent on all endpoints of the part" holds only if every endpoint's TOST passes at
     unadjusted α = 0.05. This is an intersection–union test, so no multiplicity adjustment applies to it.
  2. **Individual claims** (secondary): each endpoint's equivalence, with Holm across the part's whole family, one
     family per part. Part A: 39 (P) + 13 (F) = 52 endpoints. D: 39 (P only).
  Holm is never applied to the joint claim.
- **Invalid values:** a failed σ fit, an R80 with multiple crossings, a non-finite value or a zero ring fraction in any
  run makes that endpoint **not established** (it fails the joint claim). It is reported as such, never imputed, and no
  run is excluded. A run whose simulated count is below N fails the run.
- **Three outcomes per endpoint:** equivalent (90% interval inside the margin); not equivalent (90% interval wholly
  outside); **inconclusive** (otherwise). Inconclusive is a legitimate result, not a failure to be rescued with more
  runs after the fact.

## Sizing (from planning variance, separate from timing)

Planning data exist only at 200 MeV (§5 of the report: 4 runs per arm, 1e7 histories). With B runs per arm the SE
scales as about √(4/B) of the planning SE. TOST power is about 80% when (margin − |true difference|) / SE ≥ about 2.6.
At B = 8, taking the true difference as the planning estimate (which still carries host confounding):

| endpoint (200 MeV) | planning SE (4 + 4) | projected SE (8 + 8) | margin | (margin − \|δ\|) / SE |
|---|---|---|---|---|
| R80 | 0.0038 mm | 0.0027 | 0.05 | 18 |
| σ@200 | 0.0014 mm | 0.0010 | 0.02 | 18 |
| ring 200 mm 40–80 (log) | 0.0021 | 0.0015 | log 1.02 = 0.0198 | 8.5 |
| ring 200 mm 80–200 (log) | 0.0052 | 0.0037 | log 1.05 = 0.0488 | 11 |
| ring 100 mm 80–200 (log) | 0.0188 | 0.0133 | log 1.05 = 0.0488 | 3.5 (1.4 at ±2%) |

So B = 8 at 1e7 is adequate for every 200 MeV endpoint except the 80–200 mm ring at the shallower slab, which is
underpowered at ±2% and about adequate at ±5%; hence the wider proposed margin. For the individual Holm-adjusted claims,
the smallest threshold in a 52-endpoint family needs about 3.9 rather than 2.6, so that ring may be inconclusive
individually. **100 and 150 MeV have no planning variance.** Their power is assumed comparable and not shown; rows whose
intervals are too wide are reported inconclusive. B = 8 per arm unless Alden's precision targets for part C say
otherwise.

**Timing (resource only; upper estimates, since 100 and 150 MeV run faster than the 200 MeV rates used).** One job per
arm × case, each 8 runs; pencil = 8 × 3 energies × 1e7, field edge = 8 × 3e7.

| host | job | rate used | estimate |
|---|---|---|---|
| Lenovo, 4 thr | A upstream icl, pencil / field edge | 223 s per 1e7 (task 10480) | about 1.5 h each |
| Lenovo, 4 thr | A portable MSYS2 gcc, pencil / field edge | about 400 s per 1e7 (estimated from the HP at 3 threads) | **about 2.7 h each: over 2 h** |
| HP, 3 thr | B upstream icc and portable icc, pencil / field edge | about 150 s per 1e7 (#3 c26744, sample plan) | about 1 h each |
| HP, 3 thr | B portable gcc 13, pencil / field edge | about 360 s per 1e7 (#3 c26744) | **about 2.4 h each: over 2 h** |

About 8.3 h on the Lenovo and 8.8 h on the HP in all, so the load is roughly even. The jobs under 2 h start once their
prerequisites merge. The four marked jobs wait for sjswerdloff's go-ahead; each job's first run checks its own time
against these estimates. Part C per the frozen design; TOPAS is the long pole, scheduled by clement-7074f29f.

## Seeds

Seeds are **non-duplicated assignments**, not disjoint random streams. Portable seeds thread `t` with PCG stream `t` and
state `RNG_Seed + 1e4·t` (report §7), so two portable runs reuse a state exactly when their seeds differ by a multiple
of 1e4 within the thread count. Rule: **no two portable seeds that can enter the same contrast differ by a multiple of 1e4** (the earlier what-if 9500xx and nuclear-off 9400xx blocks already do, harmlessly, because they are never compared; report §7). Blocks:
A upstream 960001+, A portable 961001+ (checked: no difference from any portable seed used so far is a multiple of 1e4); B, if
authorised, 962001+ / 963001+ / 964001+; C per TD. The run-list commit includes a check that enforces the rule.

## How it runs unattended

1. **Build once per host** (write-once directory, sha256 recorded); run jobs verify the hash and never compile
   (sjswerdloff's build/run separation). The workflow file is the run request; RUNS empty = build only.
2. **One workflow per part,** each with its own seed block.
3. **Outputs:** MCsquare on each host's system disk (sjswerdloff, #32 c27726): Lenovo under `C:\mcsq-win\ts\<sha12>` (short paths, #40); HP Linux under `~/fe-study/`;
   TOPAS on the T7.
4. **No human step mid-run.** Each job snapshots its scripts and inputs, writes `run.json`, refuses on any count, config
   or hash failure, and leaves the failure recorded. A failed part does not block the others.
5. **Watching:** clement-7074f29f and connor-227743e6 check task status and post per-part status to #32. Nobody reads
   endpoint results before every run of that part is complete.
6. **Order:** A and B start in parallel on their two machines once the `--depths` change (#47) merges (#31 is merged).
   C once the freeze PR is merged. D needs no new runs.

## Decisions needed from sjswerdloff before he leaves

1. Approve parts A, C and D, or name the parts to drop.
2. The margins above, every row (80–200 ring settled at ±5%; the others proposed).
3. ~~Merge #31~~: merged.
4. ~~Hosts~~: settled. A on the Lenovo, B on the HP.
5. Go-ahead for the four jobs estimated at over 2 h (timing table).

## Out of scope

Performance changes (`performance_future_investigation.md`, parked); any physics change to portable; validation of the
#16 fix against measurement.
