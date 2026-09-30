# Portable MCsquare vs TOPAS: pre-registered design (issue #32)

Status: **DRAFT, not frozen.** Revised for alden-ec2221c7's review 6822. Still to fix before freezing: B, H and the
margins, after the planning runs. No stage-1 data
exists. The TOPAS smoke tests (1e4 histories) checked the plumbing only and are not study data. After freezing,
changes are visible in this file's git history, and anything changed after data is read is marked as such.

Authors: connor-227743e6 (design, MCsquare arms, analysis), clement-7074f29f (install, TOPAS arm, stopping-power
and beam-model measurements), alden-ec2221c7 (statistics, review). An independent review by a Fable 5.1 agent was
requested by sjswerdloff; its source claims were re-read at `3807be0` before being adopted, and those not re-read
are marked.

## Question

How well does MCsquare **after the #16 fix** agree with an independent Geant4/TOPAS model, on depth dose, lateral
spread and the low-dose halo? Two MCsquare builds are compared with TOPAS: portable MCsquare, and OpenMCsquare
85bf2911 with the same fix applied (sjswerdloff, 2026-10-01). No before/after contrast is made here; the pre-fix
comparison is the upstream report's, and it is not rebuilt.

Not claimed: better agreement than Huang et al. 2018 (doi:10.1002/acm2.12420). They used a different Geant4
(10.x), different physics lists and a real nozzle model, so the two studies are not comparable on agreement.

## Stage 0: stopping-power tables, then MCsquare's R80 characterised (no TOPAS)

- **Which table MCsquare uses:** `Materials/Water/G4_Stop_Pow.dat` (mass stopping power, MeV cm²/g,
  scaled by the voxel density from the HU conversion; it matches Geant4's `ComputeTotalDEDX`, which equals the electronic value here), because `define.h` sets `DB_STOP_POW SP_GEANT4`, with linear
  interpolation between rows (`compute_EM_interaction.c`, `Total_Stop_Pow`). clement-7074f29f compared it with
  G4EmCalculator option 0 for G4_WATER (I = 78 eV) on Geant4 11.3.2 at all 800 rows: it **closely matches** option 0
  (ratio 1.00000–1.00006 over 50–400 MeV). That does not establish which Geant4 version generated the table; the
  residuals below 3 MeV (0.9919 at 0.5 MeV) and their provenance are unresolved.
- **0a, table integration only:** CSDA ranges integrated from G4_Stop_Pow.dat and from PSTAR (ICRU 49, I = 75 eV,
  https://physics.nist.gov/PhysRefData/Star/Text/appendix.html), 1/S trapezoid over **every row** from 0.5 MeV to
  E, with the residual range below 0.5 MeV stated separately rather than dropped. Reported: G4 longer by 0.48, 0.48
  and 0.45% at 100, 150 and 200 MeV (clement-7074f29f). Option 4's figures (+0.08/+0.26/+0.27 mm) are
  table-versus-option-4 **integral** differences, not measured R80 shifts. Any table comparison uses every row
  (sampling only at 100/150/200 MeV aliases option 4's oscillation into a constant).
- **0b, full-physics R80 characterisation:** MCsquare's R80 (the distal 80% dose crossing, rule as in stage 1) at
  100, 150 and 200 MeV in `Scanners/Water_Phantom`, reported next to the 0a CSDA ranges. R80 is a dose crossing,
  not a path-length integral: multiple scattering, range straggling, nuclear losses and the scoring definition all
  enter. No acceptance band is applied; a prior expectation (R80 slightly short of CSDA) is stated as a
  hypothesis, not a pass criterion.

## Stage 1: pencil beams

### Arms

| arm | build | where |
|---|---|---|
| **portable** | portable main (includes the fix) | platform decided after #31: if Windows/macOS/Linux are equivalent, the Studio |
| **OpenMCsquare + fix** | 85bf2911 with `fix_secondary_angle_energy.patch`, built with icc 2021.1 (the #17 recipe) | Linux runner (HP) |
| **TOPAS opt0** | OpenTOPAS 4.3.0 / Geant4 11.3.2, EM option 0: same stopping power as MCsquare to <0.01% above 50 MeV | Mac Studio |
| TOPAS opt4 | the same with EM option 4 (secondary; its stopping-power table integrates to CSDA ranges +0.08/+0.26/+0.27 mm different at 100/150/200 MeV, an integral difference, not a measured R80 shift) | Mac Studio |

TOPAS modules for both TOPAS arms, set explicitly:
- EM (opt0 or opt4);
- g4h-phy_QGSP_BIC_HP, g4h-elastic_HP, g4ion-binarycascade, g4decay, g4stopping;
- EM range 100 eV to 600 MeV, production cut 0.05 mm.

The physics that actually ran is taken from TOPAS's own G4EmParameters dump (the "Use ICRU90 data" line), not
from the input file.

### Source

- **Energies:** 100, 150 and 200 MeV. Monoenergetic, Gaussian σ = 3 mm in x and y, at the phantom surface.
- **TOPAS:** 10σ cutoff per axis, no angular distribution, energy spread 0, source 1 mm upstream of the surface in
  vacuum.
- **MCsquare BDL:**
  - **two rows** bracketing the energy, with identical parameters and Mean = Nominal. With one row,
    `Sequential_Search` returns −1 and `compute_beam_model.c:326` reads `Nominal_Energies[-1]`; two equal
    energies divide by zero.
  - single Gaussian (Weight2 = 0), energy spread 0.
  - **Divergence 1e-6 rad. This is a floor.** `diagonalize()` computes the small eigenvalue as a difference of
    values near spot². Below 1e-6 it quantises (~9e-16) and the divergence comes out silently wrong: 1–3% at
    1e-7, a factor 1000 at 1e-8. At exactly 0 it produces NaN. (clement-7074f29f, clang arm64; the gcc and icc builds
    are checked before acquisition, see "Source and geometry check".)
  - Correlation 1e-3, so the denominator is nonzero by construction.
  - **Nozzle-to-isocentre 1 mm, isocentre on the entrance surface** (amended 2026-10-01; was 0). With 0, a proton
    starts ON the face, `Transport_to_CT` treats it as inside (inclusive bounds), it sits at voxel index NY (one past
    the grid) and every dose lands in the first depth bin, with exit 0. 1 mm matches the TOPAS source 1 mm upstream;
    MCsquare's air polynomial then removes 0.47–0.77 keV (200–100 MeV), about 0.001 mm of water range.
    SMX/SMY nonzero. `validation/mcsquare_pencil_case.py` refuses a source plane on or inside the CT, checked on the
    values as written to the files.
  - The BDL's first line must be exactly `--UPenn beam model (double gaussian)--`: any other header selects a legacy
    parser silently (`data_beam_model.c:35`) and every primary misses the phantom, exit 0.
  - The deliberate differences from TOPAS, recorded as such: divergence 1e-6 rad, correlation 1e-3, and 1 mm of air
    (0.47–0.77 keV) in place of 1 mm of vacuum.

### Phantom and scoring

- **Phantom:** water, 400 × 400 mm laterally × 350 mm deep. MCsquare uses `Scanners/Water_Phantom`.
  `Scanners/default` maps 0 HU to Schneider_AT_AG_SI4, not water: this applied to all earlier field-edge work
  (#16, the upstream report, #31).
- **Scoring:** one Cartesian grid in both codes, 1 mm in every direction.
  - **Beam axis placement, the same in every arm:** on the corner shared by the four central voxels (x = y = 0 is
    the edge between lateral indices 199 and 200), as in the TOPAS files (smoke-test centroid −0.08 to −0.01 mm).
    Established by construction and checked in the planning runs (see "Source and geometry check"), never by
    repairing a confirmatory run. Reason: each arm is compared directly with TOPAS, so a half-voxel offset
    between codes would enter σ(z) and the core annuli against a 0.2 mm margin, and nothing downstream would flag it
    (clement-7074f29f, review 6814).
  - Layout: TOPAS binary is x-fastest, depth = 349.5 − k mm. MCsquare's axis order is to be confirmed.
  - The shared analysis script is tested on an asymmetric synthetic phantom in both layouts before any real data.
- **TOPAS scorers:**
  - PRIMARY `Dose` excludes neutron and gamma descendants.
  - SECONDARY `DoseAll` is unfiltered.
  - Why: MCsquare transports neither particle (`compute_nuclear_interaction.c:137-149`: inelastic events emit only
    p/d/α, recoils deposit locally). They carry 0.8% of the integral dose at 150 MeV (clement-7074f29f, measured).
- **MCsquare dose:** dose to medium, which is water here.

### Endpoints

All from the TOPAS `Dose` scorer (primary) and the MCsquare dose grid, by one shared script. Absolute dose per
primary (and deposited energy per primary in the phantom) is reported for every arm alongside the normalised
endpoints, exploratory only.

**Confirmatory, 200 MeV only, per MCsquare arm against TOPAS opt0:**

| endpoint | definition |
|---|---|
| R80 | IDD = dose summed over the **whole 400 × 400 mm plane** at each 1 mm depth bin. R80 is the first downward crossing of 0.8 × max(IDD) distal to the maximum, by linear interpolation between the two voxel-centre depths bracketing it (the #31 rule, `platform_study_metrics.distal_level`). No crossing: undefined, the endpoint cannot pass. More than one crossing: first one used, run flagged. |
| σ at 100 mm, σ at 200 mm | slab = the **six** 1 mm depth bins whose edges span [d − 3, d + 3) mm, i.e. bin centres d − 2.5 … d + 2.5 mm (TOPAS indices k = 349.5 − centre), midpoint exactly d. The same slabs are used for the rings. A synthetic boundary-value test pins the selected indices before any data. Profiles = the slab's dose projected onto x (summed over y) and onto y. Each is fitted over \|x\| ≤ 10 mm with a **voxel-integrated** Gaussian (erf differences over each 1 mm bin), free amplitude, centre and σ, no background term. σ = mean of the x and y fits (confirmatory); σ_x and σ_y are also reported separately as diagnostics, to expose a directional defect. A fit that does not converge, or gives σ outside [1, 20] mm, is failed: the endpoint cannot pass. |
| halo fraction, 2 rings × 2 depths | same slabs. Energy fraction in a ring = dose in voxels whose centre radius r is in the ring, divided by the dose over the whole scored plane in that slab. In water with equal voxel volumes this is the fraction of **deposited energy** in the slab, not a dose ratio. Confirmatory rings: 20 ≤ r < 40 mm and 40 ≤ r < 80 mm. Exploratory rings: [0,5), [5,10), [10,20), [80,200) mm. |

Everything else is **exploratory**: 100 and 150 MeV, R20, distal 80–20, IDD shape, the other rings, TOPAS opt4, and
`DoseAll`.

**Scope of the observable.** TOPAS `Dose` is TOPAS's dose with neutron and gamma descendants removed by ancestry.
It is not MCsquare's physics: removing those particles does not make TOPAS's remaining nuclear yields or electron
transport match MCsquare's. `DoseAll` and the `DoseAll − Dose` map are kept and reported. The 0.8% integral
fraction is not a measured halo effect.

### Estimands, invalid data, statistics

- **Run-level estimands.** Each independent run gives one value per endpoint. The estimand is the mean over runs of
  the per-run endpoint, **not** the endpoint of a pooled mean dose field. For R80 and σ it is the run-level mean
  difference (arm − TOPAS). For halo fractions it is the difference of mean ln(fraction), i.e. the log of the
  geometric-mean ratio arm/TOPAS, with asymmetric bounds ln(1 − δ) and ln(1 + δ) as in #31. Welch two-sample
  intervals on those per-run values. The planning runs check that per-run values are approximately normal and that
  the estimators have converged at the chosen histories per run. If they are not, the estimand is revisited
  **before** confirmatory acquisition, not after.
- **Invalid data never passes.** A fraction that is non-positive, non-finite or above 1, an undefined R80, or a failed fit makes that endpoint
  unable to pass. There are no pseudocounts, and no run or ring is dropped selectively. Histories per run are
  chosen from the planning runs so that the confirmatory rings are non-zero in every run. An all-zero result says
  nothing about rare events and is reported as such.
- **Planning runs (separately labelled, excluded from inference):** a few runs per arm at 200 MeV to estimate
  per-run variance, sparse-ring behaviour, and grid and estimator convergence. They also serve as the source and
  geometry check below.
- **Fixed sample, one look.** After the planning runs, B runs per arm and H histories per run are frozen in this file,
  so that the **joint** probability that every confirmatory endpoint of an arm passes is at least 0.9 under the
  assumed true discrepancies (0 and half-margin, both stated). It is estimated by simulation that resamples the
  planning runs' per-run endpoint vectors, preserving their dependence; per-endpoint power alone is not the target
  (seven endpoints at 0.9 each give only about 0.48 jointly if independent). One-sided α = 0.05 per TOST, i.e.
  90% CIs.
- **REFRAMED 2026-10-01 (sjswerdloff, relayed by clement-7074f29f):** *"I'm willing to look at the results and see
  where the boundaries are in terms of what we have demonstrated (i.e. let someone else draw the conclusion)."* So
  the tolerances below are **REFERENCE BANDS for reporting and interpretation, not acceptance margins**, and there is
  no global pass label. The report gives estimates with pointwise 90% and 95% Welch intervals against each band,
  per endpoint, per arm; readers draw the boundary. Acquisition does not wait on a margin decision.
  **The detailed statistics wording of this section, and whether any joint statement survives, is
  alden-ec2221c7's to write; the older TOST/intersection-union text in this section is superseded where it
  conflicts, pending his revision.**
- **Reference bands, 200 MeV, arm − TOPAS opt0** (agreed by clement-7074f29f, alden-ec2221c7 and connor-227743e6
  after a Fable oracle's suggestion, with Alden's amendments): R80 ±0.3 mm; σ@100 ±0.1 mm; σ@200 ±0.15 mm; ring
  energy-fraction ratios 20–40 and 40–80 mm within [0.90, 1.10], 80–200 mm within [0.75, 1.25], as asymmetric
  log bounds. A "small-discrepancy reference tier" is reported alongside (R80 ≤ 0.1 mm, σ ≤ 0.03/0.05 mm, inner
  rings ≤ 3%, far ring ≤ 10%). An R80 difference near the band is an **investigation flag**, never an attribution.
- **Scoring and estimator sensitivity, stated:** with the same TOPAS seed and build, R80 moves 0.12 mm between 1 mm
  and 0.5 mm depth bins. That is scoring-grid and estimator sensitivity combined, not crossing interpolation alone
  (clement-7074f29f, planning). It depends on falloff shape, so it
  need not cancel between codes, and it is material against a ±0.3 mm band.
- **Portable vs OpenMCsquare + fix:** reported as a build difference, exploratory. They differ in RNG (PCG vs MKL
  VSL), performance refactors and compiler.
- **Frozen before acquisition:** full source commits (portable main at the freeze commit; 85bf2911 plus the patch's
  sha256), the per-arm seed ranges (distinct, never shared between arms), thread counts, B and H, and the bands.
  Seed ranges: planning 900001+ (TOPAS) and 9001xx/9002xx (MCsquare); confirmatory TOPAS opt0 910001+, opt4
  911001+, portable 920001+, OpenMCsquare + fix 930001+.
- **Seeds and reproducibility, per code:** TOPAS reproduced runs from their seeds in the configurations tested
  (clement-7074f29f: OpenTOPAS 4.3.0, seed 101 at 8 and 16 threads, and seed 900013 across three grids). That is not
  shown in general. **MCsquare does not at more than one thread** (#39): primaries go to per-thread streams through
  a shared counter, so a seed identifies streams, not a realisation, and the simulated count can exceed N. Each
  MCsquare run records requested and simulated primaries; distinct seeds remain distinct streams.
- **Huang et al. 2018 is context, and its figures are TOPAS-against-MEASUREMENT** (corrected 2026-10-01; this line
  read as if they were MCsquare comparisons): R80 differences generally under 0.1 mm on average (§3.A.2), range
  within 0.6 mm (§3.D), spot size 0.1 ± 0.1 mm (§3.B). MCsquare appears only against TOPAS on patient plans (§3.F:
  99.2% of lung iCTV voxels within 3%). None of it is a band here. (Read through an agent's page summary; the two
  quoted sentences should be checked against the PDF before the report cites them.)

### Source and geometry check (before confirmatory acquisition)

- **Geometry by construction:** voxel edges at integer mm with x = y = 0 on the shared corner of the four central
  voxels in both codes; source position and direction as specified; MCsquare's axis order confirmed on an asymmetric
  synthetic phantom.
- **Diagnostic, from the planning runs:** the fitted centroid of each arm, with its uncertainty from the planning
  runs' spread. A centroid outside ±0.1 mm by more than that uncertainty stops the study for a geometry
  investigation **before** any confirmatory run. No confirmatory run is repaired or discarded after being looked at.
- **Beam model on the real builds:** the 1e-6 rad divergence / 1e-3 correlation BDL is checked on each actual build
  used (portable: gcc on Linux/Windows, Apple clang on macOS; OpenMCsquare: icc): sampled σ_x and σ_θ against the
  BDL, and no NaN. clement-7074f29f's check compiled `diagonalize()` alone with clang -O2 on arm64; it does not
  substitute for the real builds.

## Relation to #31

#31 tests selected endpoints of one 200 MeV field in Schneider_AT_AG_SI4, not pencil-beam halos in water. It can
inform which platform runs the portable arm, not certify equivalence for this benchmark. The portable configuration
used here (platform, compiler, commit, threads) is frozen and reported, and agreement is claimed for it only.

## Interpretive hypotheses, stated before confirmatory data (not findings)

- **Upper-node energy mass** (clement-7074f29f): the nuclear sampler puts each secondary-energy interval's mass at its
  UPPER node. If the ICRU 63 table rows are point values, mean secondary-proton energy is about 5 in 100 low, which
  predicts a narrower, shallower halo in MCsquare. If rows are interval averages, the sampler is exact. The table
  convention is UNRESOLVED; a trapezoid diagnostic build would be a separately authorised what-if (sjswerdloff),
  never an arm, and neither of its outcomes settles the convention.
- **δ-electrons deposited locally** in MCsquare (TOPAS transports them): a small, one-directional sharpening of the
  core, unbounded.
- **σ fit window leakage:** |x| ≤ 10 mm is about 2σ at 200 mm with no background term, so halo differences can leak
  into σ, in the direction under test. A diagnostic fit with a fixed halo term is reported alongside.
- Solid-angle weighting of the ICRU tables is applied at load (`data_nuclear.c`, `read_Nuclear_ICRU`) and was
  checked end to end on the real loader and sampler (clement-7074f29f): it is NOT a candidate.

## Unit-level check (cheap, before stage 1)

Histogram the fixed sampler's emission angles against the ICRU 63 input tables (extending
`tests/test_secondary_angle.c`).

## Provenance

- **TOPAS:** `validation/topas/` in this repository (`stage1_base.txt`, `make_run.sh`, `run_topas.sh`, from #35;
  cited by commit once #35 merges, not by the office tree). One directory per run, named for
  every choice, never reused. `provenance.txt` records host, UTC start and end, exit code, wall time, sha256 of the
  inputs and outputs, and the physics as run.
  Outputs go to `/Volumes/T7 Shield/MonteCarlo_runs/` on the Studio (sjswerdloff's call; the home volume is 94%
  full). `run_topas.sh` refuses to start below 100 GB free (office `9577c82`). Full 3-D is kept for every arm,
  ~0.86 GB per run (two double-precision scorers), about 51 GB for 10 seeds × 3 energies × 2 TOPAS arms.
- **MCsquare arms:** `.gitea/workflows/topas-mcsquare-planning.yml` (Studio) and `topas-openmcsquare-planning.yml`
  (Lenovo) from #40. Build and run are separate jobs: builds are write-once with a recorded sha256, and runs verify it
  and never compile. Each run root snapshots the generator, endpoint script, scanner tables, materials (copied or
  hashed) and build record; each run writes `run.json` (seed, requested and simulated primaries, overshoot, threads,
  energy, binary sha256, UTC start/end, `transport_status`) and `endpoints.json`. Studio outputs go to the T7, not
  the system disk. MCsquare's `Dose.mhd` is already per simulated primary; TOPAS's Sum is not divided by N, so raw
  absolute values are not compared across codes until a validated normalisation exists.
