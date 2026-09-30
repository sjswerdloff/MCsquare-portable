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
| TOPAS opt4 | the same with EM option 4 (secondary; range shifts by +0.08/+0.26/+0.27 mm at 100/150/200 MeV) | Mac Studio |

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
    1e-7, a factor 1000 at 1e-8. At exactly 0 it produces NaN. (clement-7074f29f, clang arm64; gcc/icc checked before acquisition,
    to match, not run.)
  - Correlation 1e-3, so the denominator is nonzero by construction.
  - Nozzle-to-isocentre 0, isocentre at the surface; SMX/SMY nonzero. (Review's citation `compute_beam_model.c:583-586`: otherwise an air energy-loss polynomial is applied that TOPAS does not model. Not re-read.)
  - The deliberate differences from TOPAS, recorded as such: divergence 1e-6 rad and correlation 1e-3.

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
| σ at 100 mm, σ at 200 mm | slab = the 5 depth bins centred on that depth (z ± 2.5 mm). Profiles = the slab's dose projected onto x (summed over y) and onto y. Each is fitted over \|x\| ≤ 10 mm with a **voxel-integrated** Gaussian (erf differences over each 1 mm bin), free amplitude, centre and σ, no background term. σ = mean of the x and y fits. A fit that does not converge, or gives σ outside [1, 20] mm, is failed: the endpoint cannot pass. |
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
- **Invalid data never passes.** A zero or non-finite fraction, an undefined R80, or a failed fit makes that endpoint
  unable to pass. There are no pseudocounts, and no run or ring is dropped selectively. Histories per run are
  chosen from the planning runs so that the confirmatory rings are non-zero in every run. An all-zero result says
  nothing about rare events and is reported as such.
- **Planning runs (separately labelled, excluded from inference):** a few runs per arm at 200 MeV to estimate
  per-run variance, sparse-ring behaviour, and grid and estimator convergence. They also serve as the source and
  geometry check below.
- **Fixed sample, one look.** After the planning runs, B runs per arm and H histories per run are frozen in this file,
  from sjswerdloff's margins, power 0.9, and an assumed true discrepancy (0 and half-margin both stated). One-sided
  α = 0.05 per TOST, i.e. 90% CIs.
- **Per-arm equivalence claim:** every confirmatory endpoint passes. That is an intersection-union test, so no
  multiplicity adjustment is needed. Each MCsquare arm's claim is made separately. Any individual endpoint
  advertised as a finding outside that joint claim is exploratory. Insufficient precision is reported as
  **inconclusive**, never as agreement.
- **Portable vs OpenMCsquare + fix:** reported as a build difference, exploratory. They differ in RNG (PCG vs MKL
  VSL), performance refactors and compiler.
- **Frozen before acquisition:** full source commits (portable main at the freeze commit; 85bf2911 plus the patch's
  sha256), the per-arm seed ranges (distinct, never shared between arms), thread counts, B and H, and the margins.
- **Margins:** sjswerdloff's decision, made after the planning runs show the variance. Huang's published agreement
  (R80 0.1–0.6 mm, penumbra 0.5 mm, spot σ 0.1 ± 0.1 mm) is context, not a margin.

### Source and geometry check (before confirmatory acquisition)

- **Geometry by construction:** voxel edges at integer mm with x = y = 0 on the shared corner of the four central
  voxels in both codes; source position and direction as specified; MCsquare's axis order confirmed on an asymmetric
  synthetic phantom.
- **Diagnostic, from the planning runs:** the fitted centroid of each arm, with its uncertainty from the planning
  runs' spread. A centroid outside ±0.1 mm by more than that uncertainty stops the study for a geometry
  investigation **before** any confirmatory run. No confirmatory run is repaired or discarded after being looked at.
- **Beam model on the real builds:** the 1e-6 rad divergence / 1e-3 correlation BDL is checked on the actual gcc
  (portable, Linux and macOS) and icc (OpenMCsquare) builds: sampled σ_x and σ_θ against the BDL, and no NaN.
  clement-7074f29f's check was clang on arm64.

## Relation to #31

#31 tests selected endpoints of one 200 MeV field in Schneider_AT_AG_SI4, not pencil-beam halos in water. It can
inform which platform runs the portable arm, not certify equivalence for this benchmark. The portable configuration
used here (platform, compiler, commit, threads) is frozen and reported, and agreement is claimed for it only.

## Unit-level check (cheap, before stage 1)

Histogram the fixed sampler's emission angles against the ICRU 63 input tables (extending
`tests/test_secondary_angle.c`).

## Provenance

- **TOPAS:** Clement's `make_run.sh` / `run_topas.sh` (office repo `6489450`). One directory per run, named for
  every choice, never reused. `provenance.txt` records host, UTC start and end, exit code, wall time, sha256 of the
  inputs and outputs, and the physics as run.
  Outputs go to `/Volumes/T7 Shield/MonteCarlo_runs/` on the Studio (sjswerdloff's call; the home volume is 94%
  full). `run_topas.sh` refuses to start below 100 GB free (office `9577c82`). Full 3-D is kept for every arm,
  ~0.86 GB per run (two double-precision scorers), about 51 GB for 10 seeds × 3 energies × 2 TOPAS arms.
- **MCsquare arms:** as in #31: per-run record.json with sha256s, compiler, commit and host.
