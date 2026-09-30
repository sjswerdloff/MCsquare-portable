# Portable MCsquare vs TOPAS: pre-registered design (issue #32)

Status: **DRAFT, not frozen.** Sections marked OPEN wait on alden-ec2221c7's statistics review. No stage-1 data
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

## Stage 0: range against the stopping-power tables (no TOPAS)

- **Beams:** monoenergetic pencil beams at 100, 150 and 200 MeV in `Scanners/Water_Phantom`.
- **Which table MCsquare uses:** `Materials/Water/G4_Stop_Pow.dat`, because `define.h` sets
  `DB_STOP_POW SP_GEANT4`. clement-7074f29f measured that this table is Geant4 11.3.2's **option 0** EM physics
  with I = 78 eV, using G4EmCalculator on G4_WATER at all 800 rows: the ratio is 1.00000–1.00006 over 50–400 MeV.
  The table drifts below 3 MeV (0.9919 at 0.5 MeV); the source Geant4 version for that range is unknown.
- **0a, transport:** R80 against the CSDA range integrated from G4_Stop_Pow.dat. The expected offset, written
  before the run, is R80 0.1–0.3% short of CSDA, from range straggling and nuclear loss.
- **0b, table:** the same comparison against PSTAR (ICRU 49, I = 75 eV). Measured from the tables alone, G4 CSDA
  ranges are longer by 0.48, 0.48 and 0.45% at 100, 150 and 200 MeV.
- **Rule:** any comparison of stopping-power tables uses **every row**. Sampling only at 100, 150 and 200 MeV
  aliases option 4's oscillation into an apparent constant offset (clement-7074f29f).

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
    1e-7, a factor 1000 at 1e-8. At exactly 0 it produces NaN. (clement-7074f29f, clang arm64; gcc/icc expected
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
    Verified per arm by the fitted centroid of its first run lying within ±0.1 mm of 0 in x and y; an arm that fails
    this is fixed before any further runs. Reason: each arm is compared directly with TOPAS, so a half-voxel offset
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

### Endpoints (OPEN: confirmatory set to be sized by alden-ec2221c7)

Proposed **confirmatory set, 200 MeV:**
- R80;
- σ(z) at two depths, from a Gaussian fit to the lateral core (not second moments, which neutron/gamma outliers
  dominate: 3.9–5.0 mm against 3.0 mm in the smoke test);
- ring-integrated dose fractions at two depths: the dose in an annulus of the 1 mm grid (by voxel-centre radius)
  divided by the dose integrated over the **whole scored 400 × 400 mm plane at the same depth**, both from the same
  scorer (TOPAS `Dose` for the primary comparison). Chosen over central-axis normalisation because it neither
  inherits the noise of the central voxels nor depends on the two codes' absolute dose normalisation. The annulus
  edges are fixed before data (OPEN with the depths).

Everything else is **exploratory**: 100 and 150 MeV, R20, distal 80–20, the IDD normalised at 20–30 mm, the
secondary TOPAS arm.

Huang's "dose beyond 40 mm > 1e-4 of CAX" is a threshold, not a quantity equivalence can be tested on, so it is
reported descriptively only.

### Statistics (OPEN: alden-ec2221c7)

- Independent seeds per arm, with n fixed before the run. Timing input (clement-7074f29f, 150 MeV, option 0, both
  scorers, 1e6 histories on the Studio): 120 s at 8 threads (8.7 GB), 68 s at 16 threads (15.9 GB). So about
  11 min per 1e7-history TOPAS seed at 16 threads.
- Per-seed uncertainty comes only from the spread between seeds. The review reports that MCsquare's own
  uncertainty covers voxels above 50% of max only (not re-read).
- **Agreement with TOPAS, per MCsquare arm:** TOST against sjswerdloff's margins. Proposed from Huang as a starting
  point: R80 0.5 mm, σ 0.2 mm. The halo ratio margin needs its own argument. Portable and OpenMCsquare + fix are
  each reported against TOPAS opt0, and against each other; they differ in RNG (PCG vs MKL VSL), performance
  refactors and compiler, so a difference between them is reported as a build difference, not attributed.
- A stopping rule, written before data, as in #31.

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
