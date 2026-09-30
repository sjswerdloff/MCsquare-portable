# Portable MCsquare across platforms (Windows, Linux, macOS arm64): pre-registered design

Margins decided by sjswerdloff, 2026-09-30. The macOS arm64 arm was added on his instruction the same day, before any study data.

Status: **wave 1 runs started; no study result has been read.** This is pre-analysis specification, which is not the
same claim as pre-run registration (see Order of events). The only data seen so far is the two-seed pilot on #28
(comment 27347), which is used for planning and is **excluded** from the study analysis. Changes after the first
study result would be visible in this file's git history.

**Order of events (sjswerdloff, 2026-09-30: "just run them all").** Wave 1 was started before this design and its
implementation were reviewed. The per-seed metric code (`platform_study_metrics.py`, with synthetic R80/R20 tests
in `test_platform_study_metrics.py`) and this design were committed at the commit the run uses. The analysis script
(`platform_study_analyse.py`) is committed before any `STUDY_RESULT` is read, and the git history shows that order.
Recorded timestamps (NZDT):
- 2026-09-30 17:04 — `7d07db0`, the commit wave 1 builds from (design, metric code, workflow).
- 2026-09-30 17:07 — `373961a`, the analysis script. Wave 1 was already running.
- 2026-09-30, during wave 1 — seeds windows 1008 and linux 2003 interrupted by a runner relabel; re-run at `600c66e`
  (rule in "Interrupted seeds").
- 2026-10-01, during wave 1, before any result was read — review-driven amendments (alden-ec2221c7, review 6821):
  the analyzer validates the dataset and fails closed before inference, and keeps a stop state between looks; the
  workflow refuses a failed build or stale output (prospective, for wave 2; no completed seed is re-run or dropped).
- First inspection of results: recorded here when it happens.

## Question

The same source revision is built three ways:
- `MCsquare_win_portable.exe`: Windows native, MSYS2 UCRT64 GCC;
- `MCsquare_portable`: Linux/WSL2, gcc;
- `MCsquare_arm64`: macOS arm64, Apple clang + Homebrew libomp.

Linux is the reference. There are two comparisons: **Windows vs Linux** and **macOS vs Linux**. For each, do the
dose distributions differ by less than a stated tolerance (**equivalence**, primary)? Is there **any detectable
difference** at all (secondary, reported but not used to stop)?

## Fixed conditions (identical on both platforms)

- **Machines (sjswerdloff, 2026-09-30, "just run them all"):**
  - **Linux** on the HP (DESKTOP-5H86O9N, i5-6500T, WSL2 Ubuntu), which has the only Linux runner.
  - **Windows native** on whichever Windows runner takes each seed: the HP or the Lenovo (DESKTOP-SR5GKKA, i5-6400T).
    There is no per-box label yet.
  - **macOS arm64** on the Mac Studio (M3 Ultra), pinned by its runner label `stuart-m3ultra-canary`.
  - All three arms run in parallel.
  - **Why mixing Windows hosts is acceptable:** both are Skylake. Builds target generic x86-64 (no `-march=native`),
    and both Windows boxes use the same MSYS2 GCC 16.2.0, so the same instruction stream should compute identically.
  - **Every seed records its host.** A Windows-by-host breakdown is reported as **exploratory** only. The confirmatory
    Windows arm pools all Windows seeds whichever host ran them, and no seed is excluded for its host.
- **Source revision:** one commit, recorded. Both builds come from it, at `-O3`, with their compiler versions logged.
- **Case:** 200 MeV, 15 × 15 cm field, 300 mm water cube with 2 mm voxels, `BDL_default_DN_RangeShifter.txt`, the
  default Scanners conversion files. `Num_Threads 3`, `Num_Primaries 3e7`, `Dose_MHD_Output True`.
  - **Correction, 2026-10-01, before any study result was read (description only; the input is unchanged):** the
    "water cube" is 0 HU through `Scanners/default`, which maps HU 0 to `Schneider_AT_AG_SI4` (label 44, density
    1.000), not to `Water` (label 17). Every arm uses the identical input, so the comparison is unaffected; the medium
    is a soft-tissue mixture at unit density, not water.
- **Seeds: distinct streams per platform, so the two samples are independent.**
  - Windows: wave 1 uses 1001–1016, wave 2 uses 1017–1032.
  - Linux: wave 1 uses 2001–2016, wave 2 uses 2017–2032.
  - macOS: wave 1 uses 3001–3016, wave 2 uses 3017–3032.
  - No seed runs on more than one platform. `Num_Threads 3` on all three.

## Metrics per seed (fixed voxel locations, computed identically on both platforms)

- **Lateral (8 endpoints).** Dose 5, 10, 20 and 30 mm outside the field edge at 127 mm and 201 mm depth, as % of
  the CAX dose at the same depth. These are exactly the rows of `validation/field_edge_profile.py`.
- **CAX (3 endpoints).** Absolute dose at 127 mm, at 201 mm, and at voxel index 23 (253 mm, the peak depth in the
  pilot).
  - **Analysed on the log scale:** each seed contributes ln(dose), and the CI is for mean ln W − mean ln L, which is
    the log of the geometric-mean ratio W/L.
  - **Margin, asymmetric:** a relative margin of ±δ% means the ratio must lie in [1 − δ/100, 1 + δ/100], applied
    exactly as ln(1 − δ/100) ≤ CI ≤ ln(1 + δ/100).
  - **No uncertain denominator:** the analysis never divides by an observed Linux mean.
  - **Invalid dose:** a CAX dose that is zero, negative or non-finite in any seed **fails** that endpoint for the
    study. It never enters the log, and the endpoint cannot pass.
- **Lateral rows need no such step.** Each is a ratio within one seed (lateral dose over that same seed's CAX dose at
  that depth), so the per-seed value is already on the percentage-point scale that the margin uses.
- **Distal edge (2 endpoints).** R80 and R20, compared in mm. The crossing rule, fixed now:
  - D(z) is the CAX depth-dose, sampled at voxel-centre depths z, averaged over the same central 10 × 10 voxels as
    the CAX rows.
  - Dmax is its maximum over all depths, and z_peak is the depth of that maximum.
  - For each level L = 0.8·Dmax and L = 0.2·Dmax, scan depth increasing from z_peak. **R_L is the first adjacent
    pair (z_i, z_i+1) with D(z_i) ≥ L > D(z_i+1)**, found by linear interpolation between those two voxels.
  - **Ambiguous:** more than one downward crossing of L distal to z_peak (the curve re-crosses). R_L is still the
    first crossing. That seed is flagged, and the number of flagged seeds per platform is reported.
  - **Absent:** no crossing before the end of the phantom. R_L is undefined for that seed. If any seed on either
    platform has an undefined R_L, that endpoint cannot pass, so joint equivalence cannot be claimed. The endpoint
    is reported as failed-to-measure, not dropped.
  - **Synthetic tests before wave 1:**
    - a smooth distal fall-off with analytically known R80 **and** R20;
    - a curve with a noise bump that re-crosses 80%, which must return the first crossing and flag it;
    - a curve that never falls below 20%, which must report R20 as undefined.

That makes 13 endpoints.

## Tolerances δ (equivalence margins), **decided by sjswerdloff 2026-09-30**

These are research thresholds for whether two builds of the same code agree, not clinical acceptance criteria.
The options were A (tighter) and B (looser); the decision per endpoint:

| endpoint | margin |
|---|---|
| lateral, 10/20/30 mm out (6 endpoints) | ±0.20 percentage points of CAX |
| lateral, 5 mm out, penumbra (2 endpoints) | ±0.50 pts (about a 0.2 mm edge shift at 201 mm, where the fall-off is ~2.6 pts/mm; about 1 mm at 127 mm, ~0.46 pts/mm) |
| CAX dose (3 endpoints) | ratio W/L in [0.995, 1.005] |
| R80, R20 (2 endpoints) | ±0.5 mm |

Pilot planning numbers, normal approximation, 90% power, one-sided α = 0.025 per look:
- **Most endpoints:** 16 seeds per platform suffices for option A if the true difference is near zero.
- **201 mm, 5 mm out:** per-seed SD ≈ 0.13 pts, observed difference +0.084. At option A (0.10) this is effectively
  unreachable (~1500 seeds) if that difference is real, and needs ~39 seeds if the true difference is 0. At option B
  (0.50) it needs about 3 seeds whether or not that difference is real. At 0.20 it needs about 29 seeds if it is real,
  and about 10 if not.
- **Caveat:** SDs from two seeds have one degree of freedom, and the true SD can be several times larger. These
  numbers are a guide, not a guarantee. The interim look (below) re-estimates them from the study's own data.

## Analysis (fixed now)

- **Per endpoint:** Welch two-sample comparison of the per-seed values, Windows versus Linux, with a confidence
  interval for the mean difference W − L (on the log scale for CAX). It is **unpaired**. The samples are
  independent by construction, because each platform uses its own seeds (see Fixed conditions). Whether same-seed
  pairing would have helped is **unestablished**; two pilot seeds cannot estimate that correlation.
- **Equivalence:** TOST for each endpoint. It passes if the (1 − 2α) CI for W − L lies inside ±δ.
  - Equivalence is claimed **per comparison** (Windows vs Linux, macOS vs Linux), only if all 13 of its endpoints
    pass. That is an intersection-union test, so no multiplicity adjustment is needed for either claim.
  - If both comparisons pass, the claim is: **both tested configurations are equivalent to the Linux reference, for
    these endpoints and margins.** It is not a claim that Windows and macOS are within the margins of each other
    (e.g. CAX ratios 1.004 and 0.996 both pass against Linux, but differ from each other by ~0.8%), nor that each
    Windows host individually agrees (the Windows arm pools the tested host mixture), nor that other fields, media or
    halo configurations (e.g. #33's pencil beams in water) behave the same. It can inform the choice of platform for
    #33, not certify it.
  - **Each comparison stops independently** under the stopping rule. A comparison that has stopped runs no more
    seeds; the Linux reference keeps running while either comparison needs wave 2.
- **Sequential looks:** at most two, with α split by Bonferroni: α = 0.025 at each look, so a 95% CI per look. That
  keeps the overall one-sided error at or below 0.05 whatever happens at the first look.
- **Difference detection (secondary, confirmatory at familywise 0.05):** two-sided Welch test per endpoint,
  Holm-adjusted across all **26** endpoint-comparisons at **α = 0.025 per look**, so two looks give familywise ≤ 0.05. A detected
  difference inside the margin is compatible with equivalence and is reported as such.
- **Non-equivalence (confirmatory, study-wide):** at each look, simultaneous two-sided CIs for all 26 endpoint-comparisons at
  confidence 1 − 0.025/26 (Bonferroni over endpoints and both comparisons; each two-sided CI covers both directions). Over two looks, the
  familywise chance of a false non-equivalence claim is ≤ 0.05. A per-look 95% CI lying outside ±δ is only an
  **exploratory flag**.

## Stopping rule (fixed now)

Applied to **each comparison separately**. After **wave 1** (16 seeds per platform):
1. **All 13 TOST CIs (95% at this look) inside ±δ:** equivalence claimed. Stop.
2. **Any endpoint's simultaneous CI (confidence 1 − 0.025/26) lies entirely outside ±δ:** non-equivalence is
   claimed for that endpoint under the stated error budget. Stop. This is a statistical decision, not an
   impossibility: more independent data could still shift that interval.
3. **Otherwise:** run **wave 2** (seeds 1017–1032) and analyse **all 32 seeds per platform pooled**. Wave 1 is
   never discarded. Then apply 1 and 2 with the look-2 α.
4. **After wave 2 and still neither:** report *inconclusive*, with the intervals. Any third wave needs a new,
   separately pre-registered α allocation; the error guarantee above covers only two looks.

## Interrupted seeds (added 2026-09-30 during wave 1, before any result was read)

A seed whose job is interrupted by infrastructure before it writes a record is **re-run with the same seed number
from the same source tree**. It produced no data, so re-running it cannot select on results. A seed that completed and wrote
a record is never re-run and never replaced.

Wave 1 had two interruptions when runner labels were reset: Windows 1008 (job 14114, Lenovo) and Linux 2003 (job
14077, HP). Both were killed mid-simulation. Neither printed its `RUN` line, produced a `STUDY_RESULT`, or created a
retained directory.

They were re-run at 600c66e (branch `connor/platform-study-rerun`) on sjswerdloff's instruction, started before the
main run finished. `git diff --stat 7d07db0 600c66e` shows one added file, the re-run workflow, and nothing else. So
the source, measurement code and data are identical to 7d07db0, and each re-run records its own commit in its record.

## Provenance (every seed, retained)

- **Kept per seed:**
  - `Dose.raw`, `Dose.mhd`, the config and the plan file;
  - sha256 of those and of the BDL, material and scanner files;
  - compiler version, source commit, host, start and end time;
  - the per-seed metrics as one JSON line in the job log.
- **Where it is kept:** under a per-study directory on the machine (`C:\mcsq-win\fe-study\<study-id>\` and
  `~/fe-study/<study-id>/` in WSL2 and on the Mac Studio). It is written once, and nothing in it is overwritten or deleted by any workflow.
  Artifact upload to Gitea currently fails (it advertises an upload URL on port 80), so the files stay on the HP
  unless that is fixed.

## Cost and timing

- **Per seed:** about 26 min on Windows and about 35–41 min in WSL2 in the pilot.
- **Wave 1:** about 7 h of Windows plus about 10 h of WSL2 on the same HP. If they run concurrently they contend for
  its four cores, so the runs are best scheduled overnight and with sjswerdloff. The wall-clock time does not affect
  the dose.
- **Wave 2:** about the same again.
