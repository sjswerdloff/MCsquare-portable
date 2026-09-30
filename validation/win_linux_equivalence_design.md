# Windows vs Linux dose: pre-registered design (DRAFT, margins awaiting sjswerdloff)

Status: **design only. No study data exists.** The analysis below is fixed before any study run. The only data seen
so far is the two-seed pilot on #28 (comment 27347), which is used here for planning and is **excluded** from the
study analysis. Changes after the first study result would be visible in this file's git history.

## Question

On the same machine, the same source revision gives dose distributions from `MCsquare_win_portable.exe` (Windows
native, MSYS2 UCRT64 GCC) and `MCsquare_portable` (Linux/WSL2, gcc). Do these differ from each other by less than a
stated tolerance (**equivalence**, primary)? Is there **any detectable difference** at all (secondary, reported but
not used to stop)?

## Fixed conditions (identical on both platforms)

- **Machine and runners:** the HP (DESKTOP-5H86O9N, i5-6500T), Windows native and WSL2 Ubuntu. This needs the
  Windows job pinned to the HP by a per-box runner label.
- **Source revision:** one commit, recorded. Both builds come from it, at `-O3`, with their compiler versions logged.
- **Case:** 200 MeV, 15 × 15 cm field, 300 mm water cube with 2 mm voxels, `BDL_default_DN_RangeShifter.txt`, the
  default Scanners conversion files. `Num_Threads 3`, `Num_Primaries 3e7`, `Dose_MHD_Output True`.
- **Seeds:** wave 1 uses 1001–1016; wave 2, if run, uses 1017–1032. The same seed list runs on both platforms.
  Seeds are not treated as pairing, because the pilot showed no correlation (see Analysis).

## Metrics per seed (fixed voxel locations, computed identically on both platforms)

- **Lateral (8 endpoints).** Dose 5, 10, 20 and 30 mm outside the field edge at 127 mm and 201 mm depth, as % of
  the CAX dose at the same depth. These are exactly the rows of `validation/field_edge_profile.py`.
- **CAX (3 endpoints).** Absolute dose at 127 mm, at 201 mm, and at voxel index 23 (253 mm, the peak depth in the
  pilot). Compared as % relative difference.
- **Distal edge (2 endpoints).** R80 and R20, the depths distal to the peak where the CAX depth-dose falls to 80% and
  20% of its maximum, found by linear interpolation between voxels. Compared in mm.

That makes 13 endpoints. The distal-edge metrics are new code and will be tested before wave 1 on a synthetic
depth-dose with a known R80.

## Tolerances δ (equivalence margins): **sjswerdloff to decide**

Proposed options; the margin decides the answer more than the number of seeds does:

| endpoint | option A (reproducibility-grade) | option B (clinical-grade) |
|---|---|---|
| lateral, 10/20/30 mm out | ±0.10 percentage points of CAX | ±0.20 pts |
| lateral, 5 mm out (penumbra, steep gradient) | ±0.10 pts | ±0.50 pts |
| CAX dose | ±0.5 % relative | ±1.0 % |
| R80, R20 | ±0.5 mm | ±1.0 mm |

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
  interval for the mean difference W − L. It is unpaired because the pilot's same-seed differences were as large as
  its between-seed spread; the three-thread OpenMP schedule and differing libm break the pairing.
- **Equivalence:** TOST for each endpoint. It passes if the (1 − 2α) CI for W − L lies inside ±δ.
  - Equivalence of the platforms is claimed only if **all 13 endpoints pass**. That is an intersection-union test,
    so no multiplicity adjustment is needed across endpoints for that joint claim.
- **Sequential looks:** at most two, with α split by Bonferroni: α = 0.025 at each look, so a 95% CI per look. That
  keeps the overall one-sided error at or below 0.05 whatever happens at the first look.
- **Difference detection (secondary):** two-sided Welch test per endpoint, Holm-adjusted across the 13 endpoints at
  each look, and reported. A detectable difference inside the margin is compatible with equivalence and is reported
  as such.

## Stopping rule (fixed now)

After **wave 1** (16 seeds per platform):
1. **All 13 CIs inside ±δ:** equivalence shown. Stop.
2. **Any endpoint's CI entirely outside ±δ:** non-equivalence shown for that endpoint. Stop, because more seeds
   cannot bring it inside.
3. **Otherwise:** run **wave 2** (seeds 1017–1032) and analyse **all 32 seeds per platform pooled**. Wave 1 is
   never discarded. Then apply 1 and 2 with the look-2 α.
4. **After wave 2 and still neither:** report *inconclusive*, with the intervals. Any third wave needs a new,
   separately pre-registered α allocation; the error guarantee above covers only two looks.

## Provenance (every seed, retained)

- **Kept per seed:**
  - `Dose.raw`, `Dose.mhd`, the config and the plan file;
  - sha256 of those and of the BDL, material and scanner files;
  - compiler version, source commit, host, start and end time;
  - the per-seed metrics as one JSON line in the job log.
- **Where it is kept:** under a per-study directory on the machine (`C:\mcsq-win\fe-study\<study-id>\` and
  `~/fe-study/<study-id>/` in WSL2). It is written once, and nothing in it is overwritten or deleted by any workflow.
  Artifact upload to Gitea currently fails (it advertises an upload URL on port 80), so the files stay on the HP
  unless that is fixed.

## Cost and timing

- **Per seed:** about 26 min on Windows and about 35–41 min in WSL2 in the pilot.
- **Wave 1:** about 7 h of Windows plus about 10 h of WSL2 on the same HP. If they run concurrently they contend for
  its four cores, so the runs are best scheduled overnight and with sjswerdloff. The wall-clock time does not affect
  the dose.
- **Wave 2:** about the same again.
