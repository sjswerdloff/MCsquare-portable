# Apples-to-apples comparison programme (PLAN, for sjswerdloff's approval)

**Status: proposed 2026-10-01. Nothing here has run.** Written to execute unattended while sjswerdloff travels (week
of 5 October 2026). Once approved, each part's run list, seeds, histories and analysis are committed **before** its
first run, and nothing is changed after results are read. Results feed [`validation_report_draft.md`](validation_report_draft.md).

## Why

The current evidence compares codes across hosts: upstream on the Lenovo against portable on the Mac Studio. Build,
compiler, RNG, host and OS are all confounded (validation report §5, §8). "Apples to apples" means: **same host, same
case files, same thread count, same histories, distinct seed streams, same endpoint script, varying one thing at a
time.**

## The comparisons

| part | question | varies | held fixed | host(s) |
|---|---|---|---|---|
| **A. Intel, Windows 10** | does portable reproduce upstream + fix? | code (upstream 85bf2911 + fix, icl 2021.1, vs portable, MSYS2 gcc) | host, OS, threads (4), case, H, seeds per arm | Lenovo (HP as its equivalent if the Lenovo is busy) |
| **B. Intel, Linux** | same question under Linux | code (upstream + fix, icc 2021.1.2, vs portable gcc 13; portable built with icc as a third arm to separate compiler from source) | host (HP, WSL2), threads (3), case, H | HP |
| **C. Apple Silicon vs TOPAS (confirmatory)** | portable vs TOPAS, per the frozen design | code (portable vs TOPAS opt0; opt4 secondary) | host (Mac Studio), case, H, endpoint definitions | Mac Studio |
| **D. Bridge** | does portable on Intel agree with portable on Apple for the pencil case? | platform | build source commit, case, H | Lenovo + Studio (reuses A's portable runs and C's portable runs) |

Cases, used identically in every part they appear in:
- **P (pencil beam):** `validation/mcsquare_pencil_case.py` at 100, 150 and 200 MeV, Gaussian σ = 3 mm, water
  400 × 400 × 350 mm, 1 mm grid, nozzle 1 mm (TD). 150 MeV is added because its falloff is wide enough for a 1 mm grid
  and it fills the energy gap; 70 MeV stays out (grid-limited, #43).
- **F (field edge):** the #31 platform-study case, 200 MeV 15 × 15 cm in the default scanner's HU 0 →
  **Schneider_AT_AG_SI4**, endpoints exactly as #31. **Prerequisite: #31 merged**, because its case files and analyser
  live on that branch.

## Pre-registered analysis (same for A, B and D)

- **Estimand:** the run-level mean difference (arm − reference) for R80, σ@100, σ@200 and R20; the geometric-mean
  ratio for ring fractions and central-axis dose. Field-edge lateral rows on their own scale. Report point estimates
  with pointwise 90% and 95% Welch intervals, as in the merged design (#41).
- **Reference** for A and B: upstream + fix. For D: portable on Apple.
- **Equivalence statement:** A and B are code-identity questions, so they get an equivalence claim with margins
  **set by sjswerdloff before any run**. Proposal for his decision:
  - field edge: the #31 margins;
  - pencil: R80 ±0.05 mm, σ ±0.02 mm, rings [0.98, 1.02].
  These are tighter than the TOPAS reference bands because they compare the same physics. TOST at one-sided
  α = 0.05 per endpoint, Holm across endpoints.
- **Invalid data stays visible:** failed fits and absent crossings are reported, never imputed; there's no selective
  run exclusion.
- **Counts:** each run records requested and simulated primaries (#39). An arm whose simulated count is below N
  fails the run.

## Sizing (to finish inside the week, measured throughputs)

| part | arm | per 1e7 at 200 MeV | runs (B per case × cases) | estimate |
|---|---|---|---|---|
| A | upstream icl, Lenovo 4 thr | about 223 s (task 10480) | 8 × 3 pencil + 8 field edge (3e7) | about 3 h |
| A | portable MSYS2 gcc, Lenovo 4 thr | about 400 s (estimated from HP 3-thread 1576 s / 3e7) | same | about 6 h |
| B | upstream icc / portable gcc / portable icc, HP 3 thr | 15 s / 36 s / 15 s per 1e6 (sample plan, #3 c26744) | same | about 1 day for the three |
| C | portable, Studio 24 thr | 49–125 s (spread unexplained) | per frozen B, H | hours |
| C | TOPAS opt0 (+opt4), Studio 16 thr | about 17 min | per frozen B, H | the long pole; scheduled by clement-7074f29f |

H = 1e7 for the pencil cases, matching planning; B = 8 per arm unless Alden's precision targets for C say otherwise.
The estimates are planning figures; the first run of each arm checks its own time against them.

## How it runs unattended

1. **Build once per host** (write-once directory, sha256 recorded); run jobs verify the hash and never compile
   (sjswerdloff's build/run separation). The workflow file is the run request; RUNS empty = build only.
2. **One workflow per part,** each with its own seed block, so no two arms ever share a stream:
   A 960001+ (upstream) / 961001+ (portable); B 962001+ / 963001+ / 964001+; C uses TD's confirmatory ranges.
3. **Outputs:** MCsquare on the system disk of each host (sjswerdloff, #32 c27726). Lenovo/HP Windows under
   `C:\mcsq-win\ts\<sha12>` (short paths, #40); HP Linux under `~/fe-study/`; TOPAS on the T7.
4. **No human step mid-run.** Each job snapshots its scripts and inputs, writes `run.json`, refuses on any count,
   config or hash failure, and leaves the failure recorded. A failed part does not block the others.
5. **Watching:** clement-7074f29f and connor-227743e6 check task status (the `/actions/tasks` listing is readable)
   and post per-part status to #32. Nobody reads endpoint results before every run of that part is complete.
6. **Order:** A and B start first (Intel hosts, independent of the Studio). C starts once the freeze PR is merged.
   D needs no new runs.

## Decisions needed from sjswerdloff before he leaves

1. Approve the programme, or name the parts to drop.
2. Equivalence margins for A and B (proposal above).
3. Merge #31 (prerequisite for the field-edge case files).
4. Whether the HP may stand in for the Lenovo in part A if the Lenovo is busy (proposed: yes, recorded per run;
   both treated as equivalent Windows 10 systems).

## Out of scope

Performance changes (`performance_future_investigation.md`, parked); any physics change to portable; validation of
the #16 fix against measurement.
