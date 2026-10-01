# Portable MCsquare: correctness, cross-platform agreement and comparison with upstream and TOPAS

**DRAFT, 2026-10-01. Not a publication. Every number cites its source. Results labelled PLANNING or DIAGNOSTIC come
from runs made BEFORE the TOPAS comparison design was frozen; they are not confirmatory and must not be quoted as
final.** The confirmatory acquisition and a same-host comparison programme are in
[`apples_to_apples_plan.md`](apples_to_apples_plan.md).

Source key: `#N cID` = comment ID on MCsquare-portable issue/PR #N; `rvID` = PR review ID; `UR` = the upstream report,
submitted as [OpenMCsquare work item 42](https://gitlab.com/openmcsquare/MCsquare/-/work_items/42); `TD` =
`validation/topas_design.md`.

## 1. What Portable MCsquare is

A fork of OpenMCsquare (upstream commit `85bf2911`) that builds with free compilers on Linux (gcc), macOS on Apple
Silicon (Apple clang + libomp) and Windows (MSYS2 UCRT64 gcc). It replaces Intel-specific dependencies: Cilk Plus →
OpenMP SIMD, the MKL random-number generator → PCG, and MKL vector helpers → standard C (README). The aim is that the
**source** is portable and each host builds it natively, tuned for that host (sjswerdloff).

## 2. Defects found and fixed relative to upstream 85bf2911

| issue / PR | defect | effect of the fix (as measured) |
|---|---|---|
| **#16 / #22** (merged `d661faf`), reported upstream as work item 42 | In the nuclear inelastic samplers (proton, deuteron, alpha), the 13-bin secondary emission-angle table is interpolated at the secondary energy **after** it has been rescaled to eV, between brackets in MeV. The cumulative goes negative or non-monotone, and `angle_index` reaches 13, past the end of `ICRU_angles` (#16 body; UR) | Scope: instrumented scratch copies of **portable** main, proton secondaries only, sample plan with range shifter, 1e5 primaries, arm64. As-is: cumulative negative in 79% of samplings, non-monotone in 89%, index 13 in 27%, 52% in the 0–5° bin. Fixed: 0, 0, 0, and the sampled bins match, to about 0.01 per bin, bars computed with **the code's own sampling convention** (per-bin weights, no solid-angle weighting), not an independent physical reference (#16 body, c26781). Two upstream-recipe executables were inspected: the Linux icc binary (UR §4) and an icl build made on the HP (#16 c27259, not run, and not compared bytewise with the Lenovo executable that produced the Windows dose rows). Both read `ICRU_angles[13..14]` as 0.0 padding, so in those two binaries index 13 samples isotropically over 0–90°. This is not established for upstream executables in general |
| **#24 / #25** (merged `ffe7441`) | The fork's earlier clamp of `angle_index` to [0,12] turned a broken table into a plausible 165–180° angle | Replaced by an explicit range check and `abort()` with a FATAL message, in all three species; tested per species with mutants (#25 c27100) |
| #8 | `Transport_to_CT`: a ray that misses the CT reads uninitialised data (README) | Only `dE` was garbage; counts not materially affected (#3 c26631) |
| #5, #12, #20 | Beamlet-mode path buffers too short; `system()` temp-folder removal; an inner loop reusing the outer lane index (not reached in default builds) | README |
| #10, #14, #18 (merged `3a84c5c`) | Per-lane recomputation of whole-vector cross sections (up to 32 evaluations where upstream does 2) | `Dose.raw` and `LET.raw` byte-identical at 1 thread (seed control differs); about 31% less MC time on arm64 (#18 body) |
| #28 / #29 (merged `b818b50`) | Windows: `_mkdir` arity; `rd /s /q` replaced by a Win32 walk that removes reparse points without following them | Native removal test passes on both Windows hosts (#29 c27328) |

**Effect of the #16 fix on dose (field-edge case).** One 200 MeV 15 × 15 cm field, mid-range depth 127 mm, lateral dose
as % of the central-axis dose at the same depth. ⚠️ **Phantom material:** the sources call it a "300 mm water cube"
(UR, #16), but the default scanner maps 0 HU to **Schneider_AT_AG_SI4**, not water (TD; #31 c27572). It is the same
phantom for every row below, so the comparisons stand; the label must be corrected in any publication.

| build, host | 5 mm out | 10 mm | 20 mm | 30 mm | source |
|---|---|---|---|---|---|
| upstream 85bf2911 icc, as-is, HP Linux | 6.02 | 3.82 | 1.90 | 1.07 | UR §4; #16 c27032 |
| upstream 85bf2911 icc **+ fix**, HP Linux | 6.97 | 4.64 | 2.33 | 1.23 | same |
| change | +0.95 | +0.82 | +0.44 | +0.16 | c27032 states noise ≤ 0.05 without naming the statistic |

In this case, upstream's **computed** dose next to the field edge is lower than the corrected implementation's by up
to about 0.9 percentage points of central-axis dose. This is a change in computed dose relative to the corrected code,
not evidence about physical dose or clinical accuracy: the fix has not been compared with TOPAS or measurement (§8).
The Windows icl build gives changes of +0.98 / +0.77 / +0.37 / +0.18; its two-seed differences after the fix are
0.13 / 0.07 / 0.04 / 0.01 at 5 / 10 / 20 / 30 mm (#16 c27255).

## 3. Upstream vs Portable on Intel (same field-edge case, both with the fix)

HP and Lenovo are treated as equivalent Windows 10 systems (Intel i5-6500T and i5-6400T, both Skylake; sjswerdloff,
2026-10-01). Mean of 2 seeds × 3e7, 3 threads. The sources state seed variation differently and the figures are not
interchangeable: c27032 gives "noise ≤ 0.05" (statistic not named); c27255 gives the two-seed difference (after the fix
0.13 / 0.07 / 0.04 / 0.01); c27347 gives the half-spread of the two seeds (at 127 mm, Windows ≤ 0.009, Linux up to
0.052).

| build | host, OS | 5 mm | 10 mm | 20 mm | 30 mm | source |
|---|---|---|---|---|---|---|
| upstream + fix, icc 2021.1.2 | HP, Linux (WSL2) | 6.97 | 4.64 | 2.33 | 1.23 | #16 c27032 |
| upstream + fix, icl 2021.1 | Lenovo, Windows 10 | 7.00 | 4.64 | 2.25 | 1.21 | #16 c27255 |
| portable, MSYS2 gcc 16.2.0 | HP, Windows 10 | 6.98 | 4.69 | 2.33 | 1.23 | #28 c27347 |
| portable, gcc 13 | HP, Linux (WSL2) | 6.99 | 4.69 | 2.30 | 1.21 | #28 c27347 |

Three of the four rows agree within about 0.05 points in every column. The Lenovo upstream 20 mm cell (2.25) is about
0.08 below the others, larger than that cell's own two-seed difference (0.04). **This is descriptive (n = 2 per build),
not an equivalence test, and it can't say whether that cell is noise.** Same-geometry identity
between the platform-study case and the upstream field-edge runs rests on #28 c27347's comparison and is to be
re-verified against the case files.

**Speed, same host and compiler.** HP, 3 threads, 1e6 primaries, sample plan with range shifter: upstream icc 15.1 s;
portable built with icc and the same flags 15.4 s, after #18 (#3 c26744). Portable with gcc 13 takes 35.7 s, and
29.1 s with `-march=native`. On this one sample, built with the same compiler, the two implementations run within 2%;
that doesn't establish that every remaining performance difference is the compiler's.

**Dose before the #16 fix** (sample plan, HP, 1e6): portable gcc vs upstream icc RMS difference 1.61–2.12% of max
(#3 c26577, c26582). MCsquare reported a global statistical uncertainty of about 5.9% for these runs; that is a different
quantity from an RMS percent-of-maximum difference and isn't used here as a calibrated bound on it. The range-shifter "outside geometry" count differed
(about 26k vs 22k per 1e6); the bisect traced it to the fork's emission-angle clamp, and with the energy fixed it is
about 30k (#3 c26699). Both pre-fix counts were wrong.

## 4. Portable across platforms (study pe1, #31)

200 MeV 15 × 15 cm field-edge case, 13 endpoints, Linux (HP) as reference, 16 seeds per platform, look 1 of a planned
two-look design. The margins were set by sjswerdloff and recorded before acquisition (#31 c27370). Analysis and design
amendments were made during wave 1, before the first inference, so the design itself calls this a pre-analysis
specification, not a wholly pre-run registration (#31). Verdicts (#28 c27572 = #31 c27571):
- **Windows vs Linux: EQUIVALENT, 13 of 13 endpoints.**
- **macOS (Apple Silicon) vs Linux: EQUIVALENT, 13 of 13 endpoints.**

Largest point estimates: lateral within ±0.03 points (margins ±0.2 / ±0.5), central axis within −0.08% (margin
±0.5%), R80/R20 within ±0.01 mm (margin ±0.5 mm). The claim is bounded: it's reference-based (each platform against
Linux, not Windows against macOS), for one field, with Windows pooled over two hosts (8 + 8). Wave 2 was not needed.

## 5. Upstream vs Portable, pencil beam, 200 MeV (PLANNING; different hosts)

Upstream 85bf2911 + fix, icl, Lenovo, 4 threads (task 10480, run root `C:/mcsq-win/ts/8d900ae9cbb9`, seeds
900201–4) against portable, Mac Studio M3 Ultra, 24 threads (commit 29e970c2, seeds 900101–4), 4 × 1e7 each.
Difference = upstream − portable; ratio = upstream / portable. Point estimate (SE), Welch df, pointwise 95% interval:

| endpoint | estimate (SE) | df | 95% |
|---|---|---|---|
| R80 | +0.0008 mm (0.0038) | 6.0 | [−0.0084, +0.0101] |
| σ@100 mm | +0.0000 mm (0.0007) | 5.4 | [−0.0018, +0.0018] |
| σ@200 mm | −0.0017 mm (0.0014) | 3.4 | [−0.0058, +0.0024] |
| rings at 100 mm, 20–40 / 40–80 / 80–200 | 0.9980 / 0.9978 / 1.0018 (log SE 0.0015 / 0.0022 / 0.0188) | | [0.9944, 1.0015] / [0.9917, 1.0040] / [0.9556, 1.0503] |
| rings at 200 mm, 20–40 / 40–80 / 80–200 | 0.9994 / **1.0070** / 1.0064 (log SE 0.0012 / 0.0021 / 0.0052) | | [0.9965, 1.0023] / **[1.0019, 1.0121]** / [0.9935, 1.0195] |

The 200 mm 40–80 mm ring differs by about 0.7%, and its 95% interval excludes 1; with build, compiler, RNG, host and
thread count all differing, nothing here says which of them causes it. Every other row's interval includes zero (or 1).
**This is a configuration comparison, not a build attribution (TD).**

**Reproduction.** The per-run endpoint records for this section and §6, each with its source path and sha256, are in
`validation/report_data/`; `validation/report_tables.py` recomputes both sections' tables, with 90% and 95% intervals,
and `validation/tests/test_report_tables.py` checks it against the values printed here. All three arms' endpoints are
from `pencil_endpoints.py` sha256 `711b8a4a…` (main at 172b84f). The portable runs were first scored with the earlier
29e970c2 version and re-scored with 711b8a4a for this draft; all 136 shared fields are identical.

## 6. Portable (Apple Silicon) vs TOPAS (PLANNING)

OpenTOPAS 4.3.0 / Geant4 11.3.2, EM option 0, against portable on the Mac Studio. Pencil beam, Gaussian σ = 3 mm,
400 × 400 × 350 mm water (`Scanners/Water_Phantom`), 1 mm grid (TD). 200 MeV, 4 × 1e7 per code, one endpoint script
for both codes. TOPAS `dose.bin` is the scorer that excludes neutrons, gammas and their descendants (TD); that removes
particles MCsquare doesn't transport, but it doesn't equalise charged-secondary yields or electron transport between
the codes. TOPAS runs: seeds 900001–4, `planning/E200_opt0_seed90000N_n10000000_th16_a1` on the T7.

Point estimate (SE), pointwise 95% interval, and the reporting band. The bands are reporting references, not
acceptance criteria (TD); "outside" below means the whole 95% interval lies outside the band.

| endpoint | portable − TOPAS (SE) | 95% | reference band (TD) |
|---|---|---|---|
| R80 | **−0.512 mm (0.003)** | [−0.520, −0.504] | ±0.3 mm, **outside** |
| σ at 3 mm depth | −0.0021 mm (0.0006) | [−0.0038, −0.0004] | — |
| σ at 100 mm | +0.0026 mm (0.0006) | [+0.0008, +0.0044] | ±0.1 mm, inside |
| σ at 200 mm | +0.0037 mm (0.0015) | [−0.0003, +0.0077] | ±0.15 mm, inside |

Reproduced by `validation/report_tables.py` (§5); R80 matches #32 c27679 (−0.51). Ring energy-fraction ratios, portable/TOPAS
(#32 c27757; SE of log in brackets):

| depth | 5–10 mm | 10–20 | 20–40 | 40–80 | 80–200 |
|---|---|---|---|---|---|
| 100 mm | 0.9956 (0.0004) | 1.0346 (0.0012) | 1.0671 (0.0018) | 0.9720 (0.0015) | 0.7947 (0.0179) |
| 200 mm | 1.0304 (0.0003) | 0.9404 (0.0011) | **0.9048** (0.0015) | **0.8887** (0.0032) | 0.8898 (0.0117) |

Bands for the rings are [0.90, 1.10] for 20–40 and 40–80 and [0.75, 1.25] for 80–200. The 200 mm 40–80 ring is
outside (95% interval [0.881, 0.896]); the 200 mm 20–40 ring's interval [0.9015, 0.9081] sits just inside. 95% intervals
for every ring are in the `report_tables.py` output.

**The R80 gap is open.** Depth shape at 200 MeV (#32 c27679): peak −0.46, R90 −0.46, R80 −0.51, R50 −0.56,
R20 −0.57 mm; 80–20 width −0.06 mm, so close to a shift but not rigid. At 100 MeV R80 is −0.43 mm (2 × 1e6,
PLANNING; #32 c27679).
Inputs checked so far that do not explain the gap: TOPAS geometry, water material and density, BDL energy and spread,
and the stopping table. That doesn't exhaust static or scoring explanations. MCsquare uses `SP_GEANT4`; its water table equals Geant4 option 0 in stopping power to a ratio of
1.00000–1.00006 over 50–400 MeV (TD stage 0). The integrated CSDA ranges agree within 0.008%, about 0.01 mm
(clement-7074f29f's stage 0 record, as stated in his review of this draft). Two independently written raw
readers agree.

**Nuclear-off diagnostic** (DIAGNOSTIC, #32 c27720; 4 × 1e6 per cell, nuclear processes off in both codes): the
gap persists, −0.52 mm with nuclear on and −0.56 mm off. The contrast is C = −0.041 mm (Welch SE 0.005 mm), a small non-zero
sensitivity, with no attribution. A 70 MeV point is grid-limited at 1 mm for both codes, and MCsquare's fine-grid
control failed its pre-stated criterion, so there is **no 70 MeV comparison** (#43).

**Secondary-energy what-if** (a throwaway sensitivity build authorised by sjswerdloff, never an arm; #32 c27733,
c27734, c27757). It treats the nuclear secondary-energy table rows as point values (piecewise-linear density)
instead of upper-node masses. With these tables and energy grids, that raises the mean secondary energies by half an
energy bin: at 200 MeV, about +5.7% for protons and +27% for alphas. 200 MeV, 4 × 1e7, against portable main:
- R80: −0.0007 mm (SE 0.0028), no change detectable at this precision.
- Total deposited energy ×1.0062.
- At 200 mm depth the outer rings rise: 20–40 / 40–80 / 80–200 mm ×1.0126 / ×1.0215 / ×1.0535. That's roughly a fifth
  of the 40–80 mm deficit against TOPAS (ratio 0.9078 against 0.8887 for main).
- At 100 mm depth the 20–40 and 40–80 mm rings **fall** (×0.9934, ×0.9837), contrary to the frozen direction-only
  expectation.

The pattern is mixed and depth-dependent. It doesn't explain the R80 gap, and it doesn't say which table convention
is correct (#32 c27757).

Stated hypotheses, not findings (TD): upper-node energy mass in the secondary-energy tables (convention unresolved),
local deposition of δ-electrons, and σ-window leakage. Solid-angle weighting is ruled out (applied at load).

## 7. Reproducibility

Portable is seed-reproducible at 1 thread and **not** at more than one: primaries are handed out through a shared
counter, and the simulated count can exceed N (#39). **How seeds map to random streams.** In `compute_simulation.c`
each thread `t` seeds PCG with stream selector `t` and initial state `RNG_Seed + 1e4·t + 1e5·call`. So runs with
different seeds use the same per-thread streams from different starting states; they are not disjoint streams, and two
runs collide exactly when their seeds differ by a multiple of 1e4 within the thread count. Within each contrast reported
here no two portable seeds differ by a multiple of 1e4, so no state is reused inside a contrast. Across contrasts
there is reuse: the what-if seeds 950001–4 and the nuclear-off seeds 940001–4 differ by exactly 1e4, so thread t + 1
of one run starts where thread t of the other does. Those runs are never compared with each other. Otherwise independence between runs is the usual assumption that
different PCG starting states give effectively independent sequences. Upstream's MKL seeding was not inspected. Upstream + fix on the Lenovo, with the same binary, seed and
4 threads, is not reproducible either (#39 c27678). TOPAS reproduced runs from their seeds in the configurations
tested (OpenTOPAS 4.3.0, seed 101 at 8 and 16 threads, seed 900013 across three grids; TD), which is not shown in
general. All comparisons here therefore use non-duplicated seed assignments per arm and statistics over runs, not
same-seed identity.

## 8. Limitations of what is reported here

1. Sections 5 and 6 are PLANNING data. The confirmatory design (`validation/topas_design.md`) is not frozen; run
   counts and histories (B, H) are still open.
2. No upstream-vs-portable comparison on the same Intel host for the pencil-beam case; the field-edge comparison
   (§3) is n = 2 per build.
3. No TOPAS (or measurement) validation of the #16 fix itself yet.
4. One field for the platform study; Windows pooled over two hosts.
5. Huang et al. 2018 (JACMP, doi:10.1002/acm2.12420) is cited in our design from an agent's page summary. Its figures
   need checking against the PDF before any comparison with published agreement.
6. The phantom label correction in §2 applies to the upstream report and #16 as well.

## 9. Contributors (from the record)

sjswerdloff (direction, margins, decisions); connor-227743e6 (portable fixes, MCsquare arms, upstream report);
clement-7074f29f (TOPAS install and arm, stage 0, diagnostics); alden-ec2221c7 (statistics, review);
cora-2f1e43dc and others (reviews). Authorship for any submission is sjswerdloff's decision, subject to each
contributor's consent, an accurate account of contributions, and the venue's criteria.
