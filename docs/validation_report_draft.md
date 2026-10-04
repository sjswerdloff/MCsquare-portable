# Portable MCsquare: depth-dose and off-axis comparison with upstream OpenMCsquare and TOPAS

**DRAFT, 2026-10-04, in the form of the published MCsquare validations (Souris et al. 2016; Huang et al. 2018). Not
a publication.** The confirmatory same-host acquisition (commit `2f9dab40`, §4) is complete and its endpoint results
are in §5.0–§5.2, with figures and gamma tables made from the dose files. The TOPAS comparison
(§5.4) is a pre-specified descriptive comparison at 100, 150 and 200 MeV; it is not confirmatory and carries no
margin or outcome. The broad field at 100 and 150 MeV (part E, §5.5) is a confirmatory follow-up acquisition and is
complete. One heterogeneous case with a range shifter (the Lungman phantom, §5.6) is reported descriptively.

Sources are cited compactly: `#N` = an issue or PR on MCsquare-portable; `UR` = the upstream report,
[OpenMCsquare work item 42](https://gitlab.com/openmcsquare/MCsquare/-/work_items/42); `TD` = `validation/topas_design.md`;
`AP` = `docs/apples_to_apples_plan.md`.

## 1. Introduction

Portable MCsquare is a fork of OpenMCsquare (upstream commit `85bf2911`) that builds with free compilers on Linux
(gcc), macOS on Apple Silicon (Apple clang + libomp) and Windows (MSYS2 UCRT64 gcc). It replaces Intel-specific
dependencies: Cilk Plus with OpenMP SIMD, the MKL random-number generator with PCG, and MKL vector helpers with standard C.
The physics models, data tables and beam model are upstream's. This report compares the two code lines' dose in
water: depth dose and off-axis dose at stated depths, by field-standard metrics and gamma analysis. It also compares
Portable MCsquare with TOPAS/Geant4.

## 2. Issue found in upstream OpenMCsquare

In the nuclear inelastic samplers (proton, deuteron and alpha secondaries), the 13-bin secondary emission-angle table
is interpolated at the secondary's energy **after** that energy has been rescaled to eV, between table brackets in
MeV. The interpolated cumulative distribution goes negative or non-monotone, and the sampled angle index can reach
13, one past the end of the 13-entry `ICRU_angles` array (#16; UR). In the two upstream-recipe executables inspected
(Linux icc and a Windows icl build), the read past the end lands on zero padding, so index 13 samples isotropically
over 0–90°. That is not established for upstream executables in general (UR §4).

Measured in an instrumented scratch build of Portable MCsquare, on the corrected code's sampling convention
(proton secondaries only, 1e5 primaries), not by instrumenting an upstream executable: as-is, the cumulative
is negative in 79% of samplings and non-monotone in 89%; the index reaches 13 in 27%; and 52% of samples fall in the
0–5° bin. With the fix, the first three are 0, and the sampled bins match the table to about 0.01 per bin (#16).

**Dosimetric effect.** One 200 MeV, 15 × 15 cm field (§3.2). Lateral dose at mid-range depth (127 mm), in % of the
central-axis dose at that depth:

| upstream `85bf2911`, icc, HP (Linux) | 5 mm outside the field edge | 10 mm | 20 mm | 30 mm |
|---|---|---|---|---|
| as distributed | 6.02 | 3.82 | 1.90 | 1.07 |
| with the fix | 6.97 | 4.64 | 2.33 | 1.23 |
| difference | +0.95 | +0.82 | +0.44 | +0.16 |

(UR §4; #16. The source states seed noise ≤ 0.05 without naming the statistic.) As distributed, upstream computes
less dose just outside the field edge than the corrected code, by up to about 0.9 percentage points of central-axis
dose, in the region next to organs at risk. That is a difference in computed dose; whether the corrected code is
closer to physical dose is the subject of the TOPAS comparison (§5.4) and has not yet been shown.

**Reported upstream** as OpenMCsquare work item 42, with the patch. **Every comparison in this report applies the
same fix to both code lines:** upstream `85bf2911` + `fix_secondary_angle_energy.patch`, and Portable MCsquare, which
carries the same correction.

## 3. Materials and methods

### 3.1 Pencil-beam model (case P)

- **Phantom:** water, 400 × 350 × 400 mm (lateral × depth × lateral), 1 mm voxels, 0 HU read with
  `Scanners/Water_Phantom` (the water conversion) (TD).
- **Source:** one spot on the central axis, monoenergetic (energy spread 0), single Gaussian with σ = 3 mm in x and y,
  divergence 1e-6 rad and correlation 1e-3 (numerical floors), starting 1 mm upstream of the phantom surface.
- **Energies:** 100, 150 and 200 MeV.
- **Beam-line devices:** none. This is an idealised source, chosen so that the two codes and TOPAS see an identical,
  fully specified beam (TD).
- **Histories:** 1e7 per run, 8 independent runs per energy per arm.
- **Inputs:** `validation/mcsquare_pencil_case.py`; TOPAS uses the same geometry and beam (TD).

### 3.2 Broad-field model (case F)

- **Phantom:** 300 mm cube, 2 mm voxels, 0 HU read with `Scanners/default`. ⚠️ That scanner maps 0 HU to
  **Schneider_AT_AG_SI4, not water**. It is the same phantom in every arm, so the comparisons stand; earlier
  descriptions of this case as a "water cube" (#16, UR) must be corrected (TD).
- **Beam model:** `BDL/BDL_default_DN_RangeShifter.txt`, upstream's IBA dedicated-nozzle double-Gaussian model, with a
  nozzle exit to isocentre distance of 410 mm. Virtual source to isocentre is 2234.8 mm in x and 1859.1 mm in y. At
  200 MeV the mean energy is 200.484 MeV, the energy spread 0.441, and the primary spot σ is 2.09 mm (x) and 2.89 mm (y).
- **Field:** one energy layer at 200 MeV, a 15 × 15 cm square of equally weighted spots on a 5 mm grid (900 spots),
  isocentre at the cube centre, gantry 0.
- **Beam-line devices:** the beam model defines a 74.1 mm-WET binary range shifter (material 64, density 1.20 g/cm³);
  it is set **OUT** in this plan.
- **Histories:** 3e7 per run, 8 independent runs per arm (16 in the platform study, §5.3).
- **Inputs:** `validation/field_edge_make_cases.py 200 150`.

### 3.3 Code lines, compilers and platforms

| arm | code | compiler | host, OS, threads |
|---|---|---|---|
| A-up | upstream `85bf2911` + fix | Intel icl 2021.1 + MKL 2021.1.1 | Lenovo M700 Tiny (i5-6400T), Windows 10, 4 |
| A-port | Portable MCsquare | MSYS2 UCRT64 gcc | same Lenovo, 4 |
| B-up | upstream `85bf2911` + fix | Intel icc 2021.1.2 + MKL 2021.1.1, upstream's AVX2 flags | HP EliteDesk 800 G2 (i5-6500T), Ubuntu on WSL2, 3 |
| B-pgcc | Portable MCsquare | gcc 13 | same HP, 3 |
| B-picc | Portable MCsquare | icc 2021.1.2, upstream's flags | same HP, 3 |
| — | Portable MCsquare | Apple clang + libomp | Mac Studio M3 Ultra, macOS: 3 in the platform study, 24 in the TOPAS-stage MCsquare runs |

Each contrast compares two arms on the **same host** (AP), so build, compiler and code line are not confounded with
hardware.

### 3.4 Scoring and comparison metrics

- **Depth dose.** *Pencil beam:* integrated depth dose (IDD) over the whole scored plane; R80 by linear
  interpolation at the first distal crossing (R90 and R20 the same way, descriptive); distal fall-off R20 − R80.
  *Broad field:* the depth dose is the mean over the central 20 × 20 mm patch, not a whole-plane integral; R80 and
  R20 are taken on that curve by the same crossing rule.
- **Off-axis.** *Pencil beam:* σ from a voxel-integrated Gaussian fit over |x| ≤ 10 mm, and the fraction of the
  slab's energy in annuli 5–10, 10–20, 20–40, 40–80 and 80–200 mm from the axis. Depths per energy: 40 and 60 mm at
  100 MeV, 80 and 125 mm at 150 MeV, 100 and 200 mm at 200 MeV (AP). *Broad field:* the dose at 5, 10, 20 and 30 mm
  outside the field edge, in % of the central dose of the same profile, at 127 mm and at 201 mm depth. The profile
  is the mean of the four sides of the field over a few depth rows around the stated depth, not a single line
  through the axis (`validation/platform_study_metrics.py`).
- **Gamma analysis:** 3D gamma between arms on the mean dose of each arm:
  - 2%/2 mm and 1%/1 mm, global (normalised to the maximum dose), 10% low-dose threshold, as in the literature;
  - in addition, a **low-dose gamma** (global, about 1% threshold, or local) on the broad-field planes, because the
    region next to the field edge where upstream's defect acts lies below a 10% threshold.

  Every gamma table states its criteria and threshold.
- **Equivalence statistics:** beneath each comparison, the difference of means with Welch intervals and a two
  one-sided test (TOST) against pre-stated margins: R80 ±0.05 mm, σ ±0.02 mm, annuli out to 80 mm ×[0.98, 1.02],
  the 80–200 mm annulus ×[0.95, 1.05]; for the broad field, dose outside the field edge ±0.5 points at 5 mm and ±0.2
  points at 10–30 mm, central-axis dose ×[0.995, 1.005], R80 and R20 ±0.5 mm (#31). An endpoint is equivalent when
  its 90% interval lies wholly inside the margin, not equivalent when wholly outside, and inconclusive otherwise. A
  joint claim per contrast is made by intersection–union over all 52 endpoints (39 pencil-beam, 13 broad-field), with
  Holm adjustment for the individual claims (AP; analysis reviewed in #48).
- **Random numbers:** Portable MCsquare uses PCG. Each thread `t` starts from state `RNG_Seed + 1e4·t + 1e5·c` on
  stream `t`, where `c` counts calls to the simulation loop in the process. With the default statistical-uncertainty
  batching, that loop is called once per batch, so `c` runs from 1 to exactly 10 (MIN_NUM_BATCH; no uncertainty target is set) in every run
  here, not just 1. Two runs start an identical generator only if, on the same thread, their seeds differ by a multiple of 1e5
  that some pair of their batch counts can absorb. Within the confirmatory acquisition (§4), the seed assignments have
  been screened for calls 1–10 at the threads used and share no start (#48). That holds for that population only:
  earlier comparisons reused seeds across platforms (12345 and 67890 in the preliminary field-edge runs of §5.2). No
  screen establishes statistical independence. Neither code is bitwise reproducible at more than one thread, so every
  comparison uses statistics over independent runs, not same-seed identity (#39).

## 4. Confirmatory acquisition

The same-host comparison of §3.3 ran at commit `2f9dab404cea02f5072352f7ae5773dca27eb2a1`: 8 runs per energy per arm
for case P and 8 runs per arm for case F, on the Lenovo (part A: A-up, A-port) and the HP (part B: B-up, B-pgcc,
B-picc). That is 160 runs (120 case P, 40 case F). All 160 completed and were collected; none failed and none is
missing.

- **Collection and analysis.** `validation/apples_collect.py` verified the run trees against their recorded hashes
  and copied the endpoint records; `validation/apples_analyse.py` produced every confirmatory number in §5.1 and
  §5.2. Both ran at commit `bad1419be597a8286379a5d0d8473acfd68edc0f`. The analysis document is committed as
  `validation/report_data/apples_analysis_2f9dab40.json`; its dataset fingerprint is sha256 `0d5be975…33ba` over 321
  input files.
- **What was fixed in the tools after the runs.** The collector could not read part A's logs, which Windows
  PowerShell wrote as UTF-16; it refused all of part A and wrote nothing. It was corrected at `bad1419` before any
  collection succeeded, so, as far as connor-227743e6 knows, before any full-run endpoint was read. After the results
  were read, a second review found two checks missing from the analysis (a contrast whose two arms carry one binary;
  a run that is not usable); they were added in #53 (merged as `d80871ff`) and change no result (AP, chronology).
- **Independent checks** (silas-397300f6; #48 comments 28187 and 28192). All 208 endpoint rows (4 contrasts × 52)
  were recomputed from the plan text with separate code. All 208 agree with the analysis document on outcome, and the
  200 that have estimates agree to a relative difference of at most 3.6e-10. The endpoint scripts of each arm's
  snapshot, the same scripts the runs used, were re-run on a different host on the Dose files of all 160 runs: the
  case-F values and the case-P R80 and annulus values are identical to the recorded ones, and the case-P σ values
  agree within 7.4e-6 mm (margin 0.02 mm; the fit is not bit-reproducible across hosts). That shows the recorded
  values are what those scripts give on those Dose files. It does not test the endpoint definitions. The collector
  itself does not tie endpoint values to the Dose files (#52); for this collection that recomputation does.

## 5. Results

**Contrasts (same host, both code lines with the upstream fix):** A = Portable (A-port) against upstream (A-up),
Lenovo, Windows; B1 = Portable gcc (B-pgcc) against upstream (B-up), HP, Linux; B2 = Portable icc (B-picc) against
upstream (B-up), HP, Linux. **Differences are Portable − upstream; ratios are Portable / upstream.** Each cell is the
point estimate with its 95% interval and the outcome of the equivalence test: **E** equivalent, **I** inconclusive,
**NE** not equivalent (§3.4). The test uses the 90% interval, so a printed 95% interval can extend past a margin in
a cell marked E. **†** marks an endpoint equivalent unadjusted but not after Holm adjustment.

The tables in §5.0–§5.2 and Appendix A are written by `validation/report_confirmatory_tables.py` from the analysis
document, and a test fails if the report and the document disagree. Tables 1b, 2 and 5 and the figures are
**descriptive**: they come from `validation/report_dose_descriptives.py`, which reads the Dose files, and they carry
no margin and no claim. That script first recomputed R80 for all 120 pencil-beam runs and all 13 broad-field
endpoints for all 40 broad-field runs from the Dose files and found them identical to the recorded values.

### 5.0 Summary of the confirmatory comparison

<!-- BEGIN GENERATED: confirmatory-summary -->
| contrast | joint claim (all 52 equivalent) | equivalent, unadjusted / Holm | inconclusive | not established | not equivalent |
|---|---|---|---|---|---|
| A | not established | 46 / 46 | 4 | 2 | 0 |
| B1 | not established | 46 / 46 | 4 | 2 | 0 |
| B2 | not established | 46 / 45 | 4 | 2 | 0 |
<!-- END GENERATED: confirmatory-summary -->

In each of the three contrasts, 46 of the 52 pre-specified endpoints are equivalent, 4 are inconclusive, 2 are not
established, and none is not equivalent. The joint claim, equivalence on all 52 endpoints, is therefore **not
established in any contrast**.

The endpoints that fall short are the same six each time. All six are energy fractions in the pencil beam's far halo
(40 mm or more from the axis) at 100 and 150 MeV:

<!-- BEGIN GENERATED: confirmatory-not-equivalent -->
| endpoint | A | B1 | B2 |
|---|---|---|---|
| P100/ring_40_40_80 | 1.0157 [0.9873, 1.0449] I | 1.0350 [1.0040, 1.0670] I | 0.9958 [0.9595, 1.0335] I |
| P100/ring_40_80_200 | not established | not established | not established |
| P100/ring_60_40_80 | 0.9833 [0.9377, 1.0311] I | 0.9749 [0.9232, 1.0296] I | 0.9784 [0.9216, 1.0388] I |
| P100/ring_60_80_200 | not established | not established | not established |
| P150/ring_80_40_80 | 1.0000 [0.9883, 1.0118] E | 1.0028 [0.9926, 1.0131] E | 1.0089 [0.9963, 1.0218] E† |
| P150/ring_80_80_200 | 0.9427 [0.8695, 1.0221] I | 1.0054 [0.9423, 1.0728] I | 1.0056 [0.9465, 1.0684] I |
| P150/ring_125_80_200 | 0.9696 [0.8877, 1.0590] I | 0.9814 [0.9073, 1.0616] I | 0.9541 [0.8754, 1.0398] I |
<!-- END GENERATED: confirmatory-not-equivalent -->

(The 150 MeV, 80 mm, 40–80 mm row is equivalent in all three contrasts; it is listed because in B2 it does not
survive Holm adjustment.)

How these are to be read:

- **The two endpoints that are not established** are the 80–200 mm annulus at 100 MeV, at both depths. The energy
  fraction there is exactly zero in all 8 runs of all 5 arms (80 values). The endpoint is a ratio analysed on the log
  scale, so the estimator is undefined in these observations. That is a statement about this estimator at 1e7
  histories per run. It does not show that no dose reaches that annulus, and it does not show a difference between
  the configurations.
- **All six endpoints stay in the pre-specified family of 52.** The joint claim is reported as not established. It is
  not recomputed over a smaller family chosen after the results were seen.
- **Inconclusive is a result, not a difference.** It means the 90% interval is neither wholly inside the margin nor
  wholly outside it. Whether an interval excludes 1 (or 0) is a separate question from the equivalence class. In B1,
  the 100 MeV, 40 mm, 40–80 mm annulus is 1.035 with a 95% interval that excludes 1, and it is inconclusive. In A,
  the dose 30 mm outside the field edge at 127 mm differs by +0.025 points with a 95% interval that excludes 0, and
  it is equivalent, because the whole interval lies far inside the ±0.2 margin (Table 4).
- **Why the far halo is inconclusive.** The number of histories was sized from 200 MeV planning data. At 200 MeV
  every annulus is equivalent in every contrast. At 100 and 150 MeV the outer annuli hold less energy, and their run-to-run
  spread is too large for the margins at 1e7 histories. From the observed standard errors, and assuming the true
  difference is zero, 80% power to show equivalence at the stated margins would need roughly 4 to 19 times the
  histories for the 40–80 mm annulus at 100 MeV and 3 to 7 times for the 80–200 mm annulus at 150 MeV, **in both
  arms** (margin / SE = t(0.95) + t(0.90): at a true difference of zero both one-sided tests must reject); adding histories to one arm alone reduces the standard error by at most a factor of 1.4. The 80–200 mm
  annulus at 100 MeV would need far more, since no run scored any energy there. That is design information for a
  separate, pre-specified acquisition. These results are not to be rescued by adding runs to this one (AP).
- **Intervals that exclude the null value.** Six of the 150 confirmatory 95% intervals exclude it. If every true
  difference were zero, 7.5 would be expected; that is a benchmark for calibrated intervals, not a test that there
  are no effects, and the intervals are not independent (endpoints of one run share its noise, and B1 and B2 share
  the B-up runs). Four of the six are in contrast A, all in the dose outside the field edge: 127 mm depth at 30 mm
  (+0.025 points), and 201 mm depth at 10, 20 and 30 mm (+0.070, +0.031, −0.028). One is B2's 200 MeV, 200 mm,
  10–20 mm annulus (0.9988). Those five are small differences that lie well inside their margins, so they are
  equivalent. The sixth is the B1 row above (1.035), which is inconclusive. Portable built with gcc is also above
  Portable built with icc on the same host at that endpoint (Appendix A), and A (gcc against icl) points the same
  way without excluding 1; at the other 100 MeV depth the same annulus points the other way in all three contrasts.
  That is pointwise evidence of a non-zero difference at one endpoint. It does not establish non-equivalence, and
  it does not isolate a compiler effect. A follow-up acquisition should name it in advance.

**Annular energy budget of the unresolved endpoints (descriptive, post hoc).** For each annulus endpoint that is
not equivalent, the table gives the share of the slab's energy that the annulus holds (the mean over the 16
upstream runs), the largest departure from 1 inside any of the three 95% intervals, and their product.

<!-- BEGIN GENERATED: halo-bound -->
| endpoint | share of the slab's energy (upstream mean) | largest \|ratio − 1\| within the three 95% intervals | product (illustrative) |
|---|---|---|---|
| P100/ring_40_40_80 | 0.0337% | 6.7% | 0.0023% |
| P100/ring_40_80_200 | no energy scored in any run | no interval | — |
| P100/ring_60_40_80 | 0.0110% | 7.8% | 0.0009% |
| P100/ring_60_80_200 | no energy scored in any run | no interval | — |
| P150/ring_80_80_200 | 0.0181% | 13.1% | 0.0024% |
| P150/ring_125_80_200 | 0.0040% | 12.5% | 0.0005% |
<!-- END GENERATED: halo-bound -->

The unresolved annuli hold between 0.004% and 0.034% of the energy at their depth. The product is an illustration
of scale, not a bound on dose, for these reasons:

- The endpoints are normalised fractions. Two configurations can have equal fractions and different absolute
  energy in the slab; nothing in this table constrains the slab total.
- The share is itself an estimate, and its uncertainty is not carried into the product.
- The two annuli that scored no energy have no interval and no upper limit here.
- Turning an annular energy share into local dose needs a model of the field. For a uniform broad field in a
  homogeneous medium, with a laterally invariant kernel, the far halo of the surrounding spots contributes about the
  same share of the local dose, which would put the products at the order of 0.002% of local dose. That
  approximation does not hold at a field edge, in heterogeneous anatomy or for an arbitrary plan.

**The clinical impact of the unresolved endpoints is therefore not established by this comparison.** The direct
test is the one made at 200 MeV: the broad-field case, whose out-of-field dose endpoints were all equivalent
(Table 4). Part E ran the same case at 100 and 150 MeV: every endpoint is equivalent in every contrast (§5.5).

### 5.1 Depth dose

![Fig 1. Integrated depth dose at 100, 150 and 200 MeV: mean over 8 runs per arm, all five arms overlaid, each
normalised to the maximum of its host's upstream arm; beneath, the difference Portable − upstream for each contrast
(% of the upstream maximum).](figures/fig1_idd.png)

**Table 1. Pencil beam: R80 difference (mm). Margin ±0.05 mm.**

<!-- BEGIN GENERATED: table-1-r80 -->
| energy | contrast | R80 difference (mm) |
|---|---|---|
| 100 MeV | A | 0.0000 [−0.0015, +0.0014] E |
|  | B1 | +0.0004 [−0.0010, +0.0017] E |
|  | B2 | −0.0002 [−0.0012, +0.0008] E |
| 150 MeV | A | −0.0004 [−0.0016, +0.0008] E |
|  | B1 | −0.0006 [−0.0021, +0.0009] E |
|  | B2 | −0.0005 [−0.0023, +0.0013] E |
| 200 MeV | A | +0.0020 [−0.0013, +0.0053] E |
|  | B1 | +0.0013 [−0.0019, +0.0045] E |
|  | B2 | +0.0032 [−0.0004, +0.0068] E |
<!-- END GENERATED: table-1-r80 -->

R80 is equivalent at every energy in every contrast. The largest estimate is 0.003 mm (B2, 200 MeV), and every 95%
interval lies within ±0.007 mm, against a margin of ±0.05 mm. There is no sign of an energy dependence. The
broad-field R80 and R20 (Table 4) agree within ±0.022 mm at the 95% level.

R90, R20 and the distal fall-off of the pencil beam were not pre-specified endpoints; Table 1b gives them as
descriptive values. Every estimate is within ±0.003 mm and every 95% interval within ±0.007 mm. One of the 27
intervals excludes 0 (the fall-off at
200 MeV in A, −0.003 mm), with no multiplicity adjustment. In Fig 1 the depth-dose curves of the arms differ by at
most 0.06% of the maximum.

**Table 1b. Pencil beam: R90, R20 and distal fall-off (R20 − R80) differences, Portable − upstream (mm), estimate [95% interval] from 8 runs per arm. Descriptive: no margin, no outcome.**

<!-- BEGIN GENERATED: table-1b-range-descriptive -->
| energy | contrast | R90 difference (mm) | R20 difference (mm) | distal fall-off, R20 − R80, difference (mm) |
|---|---|---|---|---|
| 100 MeV | A | 0.0000 [−0.0008, +0.0007] | −0.0001 [−0.0006, +0.0005] | 0.0000 [−0.0014, +0.0014] |
|  | B1 | +0.0002 [−0.0005, +0.0009] | −0.0001 [−0.0007, +0.0006] | −0.0004 [−0.0015, +0.0007] |
|  | B2 | −0.0001 [−0.0006, +0.0004] | −0.0002 [−0.0010, +0.0005] | 0.0000 [−0.0011, +0.0010] |
| 150 MeV | A | −0.0002 [−0.0014, +0.0010] | +0.0007 [−0.0006, +0.0020] | +0.0011 [−0.0009, +0.0030] |
|  | B1 | −0.0006 [−0.0023, +0.0011] | +0.0001 [−0.0012, +0.0015] | +0.0007 [−0.0011, +0.0025] |
|  | B2 | −0.0002 [−0.0021, +0.0016] | +0.0001 [−0.0014, +0.0017] | +0.0006 [−0.0014, +0.0027] |
| 200 MeV | A | +0.0022 [−0.0007, +0.0051] | −0.0007 [−0.0025, +0.0011] | −0.0027 [−0.0052, −0.0002] |
|  | B1 | +0.0010 [−0.0022, +0.0041] | +0.0012 [−0.0002, +0.0025] | −0.0002 [−0.0027, +0.0024] |
|  | B2 | +0.0020 [−0.0015, +0.0055] | +0.0007 [−0.0016, +0.0029] | −0.0025 [−0.0063, +0.0013] |
<!-- END GENERATED: table-1b-range-descriptive -->

**Table 2. Gamma pass rates on the arm mean doses, integrated depth dose and 3D (global, 10% low-dose cutoff): % of evaluated points passing (points evaluated), then the split-half noise control of the upstream arm. Descriptive.**

<!-- BEGIN GENERATED: table-2-gamma-pencil -->
| energy | contrast | IDD, 2%/2 mm | IDD, 1%/1 mm | 3D, 2%/2 mm | 3D, 1%/1 mm |
|---|---|---|---|---|---|
| 100 MeV | A | 100.000 (n = 79); control 100.000 (n = 79) | 100.000 (n = 79); control 100.000 (n = 79) | 100.000 (n = 6029); control 100.000 (n = 6035) | 100.000 (n = 6029); control 100.000 (n = 6035) |
|  | B1 | 100.000 (n = 79); control 100.000 (n = 79) | 100.000 (n = 79); control 100.000 (n = 79) | 100.000 (n = 6032); control 100.000 (n = 6030) | 100.000 (n = 6032); control 100.000 (n = 6030) |
|  | B2 | 100.000 (n = 79); control 100.000 (n = 79) | 100.000 (n = 79); control 100.000 (n = 79) | 100.000 (n = 6032); control 100.000 (n = 6030) | 100.000 (n = 6032); control 100.000 (n = 6030) |
| 150 MeV | A | 100.000 (n = 161); control 100.000 (n = 161) | 100.000 (n = 161); control 100.000 (n = 161) | 100.000 (n = 20289); control 100.000 (n = 20238) | 100.000 (n = 20289); control 100.000 (n = 20238) |
|  | B1 | 100.000 (n = 161); control 100.000 (n = 161) | 100.000 (n = 161); control 100.000 (n = 161) | 100.000 (n = 20263); control 100.000 (n = 20234) | 100.000 (n = 20263); control 100.000 (n = 20234) |
|  | B2 | 100.000 (n = 161); control 100.000 (n = 161) | 100.000 (n = 161); control 100.000 (n = 161) | 100.000 (n = 20263); control 100.000 (n = 20234) | 100.000 (n = 20263); control 100.000 (n = 20234) |
| 200 MeV | A | 100.000 (n = 265); control 100.000 (n = 265) | 100.000 (n = 265); control 100.000 (n = 265) | 100.000 (n = 55838); control 100.000 (n = 55768) | 100.000 (n = 55838); control 100.000 (n = 55768) |
|  | B1 | 100.000 (n = 265); control 100.000 (n = 265) | 100.000 (n = 265); control 100.000 (n = 265) | 100.000 (n = 55805); control 100.000 (n = 55757) | 100.000 (n = 55805); control 100.000 (n = 55757) |
|  | B2 | 100.000 (n = 265); control 100.000 (n = 265) | 100.000 (n = 265); control 100.000 (n = 265) | 100.000 (n = 55805); control 100.000 (n = 55757) | 100.000 (n = 55805); control 100.000 (n = 55757) |
<!-- END GENERATED: table-2-gamma-pencil -->

Every pencil-beam gamma passes at 100%, and so does the noise control (the upstream arm's four lowest seeds against
its four highest). **This says little.** As a check of the method, the same analysis was run on the A-up mean dose
against altered copies of itself. With a 3 mm shift in depth, the pass rate at 2%/2 mm falls to 78–91% for the
integrated depth dose but only to 95–97% in 3D. With a uniform 3% increase in dose, the 3D pass rate stays at 100%
at both criteria. For a single pencil beam, most voxels above the 10% cutoff lie on steep lateral gradients, where
the distance criterion absorbs a dose difference. Gamma at these criteria is therefore a weak test here; it is
reported for comparison with the literature, and the endpoint tables above are the evidence. The 10% cutoff also
restricts the analysis to the core of the beam: the far-halo annuli, whose equivalence is unresolved (§5.0), are
outside it. The gamma compares the raw arm means, normalised to the maximum of the upstream mean; the arms are not
rescaled to each other.

### 5.2 Off-axis dose

**Pencil beam.** ![Fig 2. Lateral profiles in x and y through the axis at each energy's two depths (100 MeV: 40 and
60 mm; 150 MeV: 80 and 125 mm; 200 MeV: 100 and 200 mm), log dose axis so that the halo is visible, arms overlaid;
first lateral axis.](figures/fig2_pencil_lateral.png) ![Fig 2b. The same for the second lateral axis.](figures/fig2b_pencil_lateral_y.png)

**Table 3. Pencil beam: spot σ difference (mm; margin ±0.02 mm) and annulus energy-fraction ratios (margin
×[0.98, 1.02] out to 80 mm, ×[0.95, 1.05] for 80–200 mm).**

<!-- BEGIN GENERATED: table-3-pencil -->
| energy, depth | contrast | σ difference (mm) | 5–10 mm | 10–20 mm | 20–40 mm | 40–80 mm | 80–200 mm |
|---|---|---|---|---|---|---|---|
| 100 MeV, 40 mm | A | −0.0001 [−0.0008, +0.0005] E | 0.9999 [0.9993, 1.0006] E | 1.0010 [0.9982, 1.0037] E | 1.0059 [0.9981, 1.0138] E | 1.0157 [0.9873, 1.0449] I | not established |
|  | B1 | +0.0003 [−0.0007, +0.0013] E | 0.9999 [0.9991, 1.0008] E | 0.9990 [0.9964, 1.0017] E | 0.9984 [0.9894, 1.0076] E | 1.0350 [1.0040, 1.0670] I | not established |
|  | B2 | +0.0002 [−0.0004, +0.0009] E | 1.0001 [0.9995, 1.0007] E | 0.9995 [0.9958, 1.0032] E | 1.0013 [0.9932, 1.0095] E | 0.9958 [0.9595, 1.0335] I | not established |
| 100 MeV, 60 mm | A | −0.0003 [−0.0011, +0.0005] E | 1.0000 [0.9998, 1.0003] E | 0.9993 [0.9975, 1.0010] E | 1.0022 [0.9962, 1.0081] E | 0.9833 [0.9377, 1.0311] I | not established |
|  | B1 | 0.0000 [−0.0013, +0.0013] E | 0.9997 [0.9989, 1.0005] E | 0.9992 [0.9973, 1.0011] E | 1.0031 [0.9985, 1.0078] E | 0.9749 [0.9232, 1.0296] I | not established |
|  | B2 | 0.0000 [−0.0011, +0.0010] E | 0.9999 [0.9993, 1.0006] E | 0.9988 [0.9968, 1.0008] E | 1.0041 [0.9992, 1.0090] E | 0.9784 [0.9216, 1.0388] I | not established |
| 150 MeV, 80 mm | A | −0.0007 [−0.0020, +0.0007] E | 0.9995 [0.9989, 1.0002] E | 0.9990 [0.9966, 1.0014] E | 0.9999 [0.9964, 1.0034] E | 1.0000 [0.9883, 1.0118] E | 0.9427 [0.8695, 1.0221] I |
|  | B1 | +0.0001 [−0.0011, +0.0012] E | 1.0002 [0.9996, 1.0008] E | 0.9995 [0.9969, 1.0020] E | 1.0012 [0.9986, 1.0039] E | 1.0028 [0.9926, 1.0131] E | 1.0054 [0.9423, 1.0728] I |
|  | B2 | +0.0006 [−0.0008, +0.0021] E | 1.0001 [0.9994, 1.0008] E | 1.0013 [0.9992, 1.0034] E | 1.0004 [0.9962, 1.0046] E | 1.0089 [0.9963, 1.0218] E† | 1.0056 [0.9465, 1.0684] I |
| 150 MeV, 125 mm | A | −0.0005 [−0.0017, +0.0007] E | 0.9999 [0.9993, 1.0005] E | 0.9999 [0.9987, 1.0011] E | 0.9995 [0.9968, 1.0023] E | 1.0042 [0.9940, 1.0145] E | 0.9696 [0.8877, 1.0590] I |
|  | B1 | −0.0004 [−0.0018, +0.0010] E | 1.0001 [0.9995, 1.0007] E | 0.9996 [0.9977, 1.0015] E | 1.0007 [0.9979, 1.0035] E | 0.9928 [0.9852, 1.0004] E | 0.9814 [0.9073, 1.0616] I |
|  | B2 | +0.0001 [−0.0015, +0.0016] E | 1.0003 [0.9996, 1.0009] E | 0.9991 [0.9976, 1.0005] E | 0.9994 [0.9961, 1.0027] E | 0.9931 [0.9862, 1.0000] E | 0.9541 [0.8754, 1.0398] I |
| 200 MeV, 100 mm | A | +0.0001 [−0.0011, +0.0012] E | 1.0001 [0.9993, 1.0009] E | 1.0001 [0.9976, 1.0025] E | 1.0009 [0.9979, 1.0040] E | 0.9983 [0.9942, 1.0024] E | 1.0160 [0.9953, 1.0372] E |
|  | B1 | −0.0003 [−0.0015, +0.0010] E | 0.9998 [0.9991, 1.0005] E | 0.9994 [0.9963, 1.0025] E | 1.0012 [0.9978, 1.0047] E | 0.9970 [0.9935, 1.0005] E | 0.9822 [0.9630, 1.0017] E |
|  | B2 | +0.0003 [−0.0010, +0.0015] E | 0.9998 [0.9992, 1.0004] E | 1.0015 [0.9994, 1.0036] E | 1.0013 [0.9980, 1.0045] E | 0.9963 [0.9902, 1.0025] E | 0.9874 [0.9663, 1.0091] E |
| 200 MeV, 200 mm | A | +0.0006 [−0.0015, +0.0028] E | 1.0002 [0.9995, 1.0008] E | 1.0005 [0.9988, 1.0023] E | 1.0001 [0.9972, 1.0030] E | 0.9995 [0.9946, 1.0044] E | 1.0034 [0.9860, 1.0211] E |
|  | B1 | +0.0003 [−0.0016, +0.0021] E | 0.9998 [0.9993, 1.0002] E | 0.9999 [0.9989, 1.0009] E | 1.0011 [0.9987, 1.0035] E | 1.0009 [0.9956, 1.0062] E | 1.0066 [0.9855, 1.0280] E |
|  | B2 | −0.0011 [−0.0030, +0.0008] E | 0.9998 [0.9994, 1.0003] E | 0.9988 [0.9977, 0.9998] E | 1.0006 [0.9993, 1.0019] E | 1.0029 [0.9972, 1.0087] E | 1.0060 [0.9853, 1.0271] E |
<!-- END GENERATED: table-3-pencil -->

- **Spot size.** σ is equivalent in all 18 cells. The largest estimate is 0.001 mm and every 95% interval lies within
  ±0.003 mm.
- **Core and near halo (out to 40 mm).** Every annulus is equivalent at every energy and depth in every contrast.
  The 95% intervals lie within 0.11% of 1 for 5–10 mm, 0.42% for 10–20 mm and 1.4% for 20–40 mm.
- **Far halo (40 mm and beyond).** At 200 MeV both outer annuli are equivalent at both depths in every contrast. At
  150 MeV the 40–80 mm annulus is equivalent and the 80–200 mm annulus is inconclusive. At 100 MeV the 40–80 mm
  annulus is inconclusive and the 80–200 mm annulus is not established (§5.0).

**Broad field.** ![Fig 3. Lateral profiles across the field edge at mid-range depth (127 mm), arms overlaid, log dose
axis; inset: dose 0–30 mm outside the edge.](figures/fig3_field_edge.png) ![Fig 4. Gamma maps on the plane at 127 mm depth for each contrast: 1%/1 mm global with a 10% cutoff, and 2%/2 mm
local with a 1% cutoff.](figures/fig4_field_gamma.png)

**Table 4. Broad field, 200 MeV: dose outside the field edge (difference in percentage points of the central-axis
dose at that depth; margin ±0.5 at 5 mm, ±0.2 at 10–30 mm), central-axis dose (ratio; margin ×[0.995, 1.005]) and
range (difference in mm; margin ±0.5 mm).**

<!-- BEGIN GENERATED: table-4-field -->
| endpoint | scale | A | B1 | B2 |
|---|---|---|---|---|
| 127 mm depth, 5 mm outside the edge | points | −0.031 [−0.089, +0.027] E | +0.009 [−0.027, +0.046] E | −0.024 [−0.068, +0.021] E |
| 127 mm depth, 10 mm outside | points | +0.017 [−0.030, +0.065] E | +0.006 [−0.032, +0.044] E | −0.008 [−0.049, +0.033] E |
| 127 mm depth, 20 mm outside | points | −0.004 [−0.044, +0.036] E | +0.019 [−0.009, +0.047] E | +0.017 [−0.016, +0.050] E |
| 127 mm depth, 30 mm outside | points | +0.025 [+0.008, +0.043] E | −0.004 [−0.027, +0.019] E | −0.008 [−0.033, +0.017] E |
| 201 mm depth, 5 mm outside the edge | points | −0.054 [−0.123, +0.015] E | −0.031 [−0.095, +0.033] E | −0.043 [−0.116, +0.030] E |
| 201 mm depth, 10 mm outside | points | +0.070 [+0.023, +0.117] E | +0.007 [−0.039, +0.053] E | +0.016 [−0.036, +0.069] E |
| 201 mm depth, 20 mm outside | points | +0.031 [+0.005, +0.057] E | −0.014 [−0.048, +0.020] E | −0.016 [−0.061, +0.028] E |
| 201 mm depth, 30 mm outside | points | −0.028 [−0.049, −0.007] E | −0.009 [−0.042, +0.023] E | −0.010 [−0.043, +0.023] E |
| central-axis dose at 127 mm | ratio | 1.0006 [0.9990, 1.0023] E | 0.9992 [0.9972, 1.0013] E | 0.9997 [0.9970, 1.0024] E |
| central-axis dose at 201 mm | ratio | 1.0010 [0.9986, 1.0034] E | 0.9999 [0.9975, 1.0023] E | 0.9999 [0.9981, 1.0016] E |
| central-axis dose at 253 mm | ratio | 1.0002 [0.9976, 1.0029] E | 1.0001 [0.9984, 1.0018] E | 1.0009 [0.9994, 1.0024] E |
| R80 | mm | +0.005 [−0.011, +0.021] E | −0.001 [−0.017, +0.015] E | −0.001 [−0.017, +0.014] E |
| R20 | mm | 0.000 [−0.012, +0.013] E | −0.008 [−0.020, +0.005] E | −0.007 [−0.021, +0.007] E |
<!-- END GENERATED: table-4-field -->

All 13 broad-field endpoints are equivalent in all three contrasts. Outside the field edge, where upstream as
distributed computes less dose (§2), the largest estimate is 0.07 points and every 95% interval lies within ±0.12
points of central-axis dose. Central-axis dose agrees within 0.10% (95% intervals within 0.34%).

**Table 5. Gamma pass rates on the broad-field arm mean doses, 3D (global unless stated): % of evaluated points passing (points evaluated), then the split-half noise control of the upstream arm. Descriptive.**

<!-- BEGIN GENERATED: table-5-gamma-field -->
| contrast | 2%/2 mm, 10% cutoff | 1%/1 mm, 10% cutoff | low-dose: 2%/2 mm, 1% cutoff | local 2%/2 mm, 1% cutoff |
|---|---|---|---|---|
| A | 99.9991 (n = 772072); control 99.9918 (n = 771394) | 99.9005 (n = 772072); control 99.1216 (n = 771394) | 99.9993 (n = 1013256); control 99.9938 (n = 1011009) | 99.7408 (n = 1013256); control 98.7499 (n = 1011009) |
| B1 | 99.9991 (n = 771760); control 99.9926 (n = 771311) | 99.9099 (n = 771760); control 99.0993 (n = 771311) | 99.9993 (n = 1012384); control 99.9944 (n = 1011189) | 99.7564 (n = 1012384); control 98.7198 (n = 1011189) |
| B2 | 99.9987 (n = 771760); control 99.9926 (n = 771311) | 99.9037 (n = 771760); control 99.0993 (n = 771311) | 99.9990 (n = 1012384); control 99.9944 (n = 1011189) | 99.7514 (n = 1012384); control 98.7198 (n = 1011189) |
<!-- END GENERATED: table-5-gamma-field -->

For the broad field, at least 99.90% of the evaluated points pass at 1%/1 mm with the 10% cutoff, and at least
99.74% pass the strictest analysis (2%/2 mm local with a 1% cutoff). In every cell the pass rate between code lines
is higher than the pass rate of the control within the upstream arm (99.10–99.12% and 98.72–98.75% for those two
analyses). The control compares means of four runs, which are noisier than the means of eight used for the
contrasts, so this is a descriptive comparison only: it does not show that the remaining failures are statistical
noise, which would need a noise-matched analysis.

**What the gamma analyses cover.** The cutoff is a percentage of the maximum of the upstream mean dose (at the
Bragg peak), for the local analysis as well as the global ones. With the 10% cutoff the evaluated region on the
127 mm plane is the field up to about its edge. With the 1% cutoff it reaches roughly 15 mm outside the edge
(Fig 4). The endpoints at 20 and 30 mm outside the edge are not covered by any gamma analysis here; Table 4 is the
evidence for them.

**Compiler (descriptive).** Portable built with gcc against Portable built with icc on the HP has no margin and
makes no claim (AP). Its 52 rows are in Appendix A. Fifty have estimates. Two 95% intervals exclude the null value:
the 100 MeV, 40 mm, 40–80 mm annulus discussed in §5.0, and the 200 MeV, 200 mm, 10–20 mm annulus, at 1.0011.

### 5.3 Portable MCsquare across platforms (study pe1)

Broad field (§3.2), 13 endpoints (central-axis dose, off-axis dose outside the field edge, R80 and R20), Linux on the HP
as reference, 16 runs per platform, margins set before acquisition (#31):

- Windows (HP and Lenovo pooled) against Linux: **equivalent on 13 of 13 endpoints.**
- macOS (Apple Silicon) against Linux: **equivalent on 13 of 13 endpoints.**

Largest point estimates: off-axis within ±0.03 points (margins ±0.2 / ±0.5), central axis within −0.08% (margin ±0.5%),
R80/R20 within ±0.01 mm (margin ±0.5 mm). The claim is against the Linux reference, for one field.

### 5.4 Portable MCsquare against TOPAS (descriptive)

OpenTOPAS 4.3.0 / Geant4 11.3.2, EM option 0; the TOPAS dose scorer excludes neutrons, gammas and their descendants.
Pencil beam (§3.1) at 100, 150 and 200 MeV. TOPAS: 8 runs × 1e7 per energy on the Mac Studio, acquisition frozen at
commit `5aedf3ac` (`validation/topas_halo_addendum.md`, #54). MCsquare: the three Portable arms of the same-host
acquisition (§4), 8 runs × 1e7 each, unchanged.

**This comparison is descriptive.** The MCsquare endpoints had been read before its design was fixed, and the
Portable arms ran on other hosts and with other compilers than TOPAS. There is no margin and no outcome. Intervals
are pointwise 95% (Welch). The three Portable columns share one TOPAS arm, so they are not three independent
comparisons. The analysis script was amended twice after the freeze and before any TOPAS dose file was opened; both
amendments are recorded in the addendum. TD's bands are shown at 200 MeV for reference only: they belong to TD's
confirmatory design, which has not been run.

**Table 6. Pencil beam, MCsquare against TOPAS: R80 and spot σ differences, MCsquare − TOPAS (mm), and annulus
energy-fraction ratios, MCsquare / TOPAS (fractions of each code's own slab energy). Estimate [95% interval], 8 runs
per arm. Descriptive: no margin, no outcome.**

<!-- BEGIN GENERATED: table-6-topas-halo -->
| energy | endpoint | A-port | B-pgcc | B-picc | TD band |
|---|---|---|---|---|---|
| 100 MeV | R80 (mm) | −0.4248 [−0.4262, −0.4234] | −0.4248 [−0.4260, −0.4236] | −0.4253 [−0.4260, −0.4247] |  |
|  | σ at 40 mm (mm) | −0.0107 [−0.0115, −0.0100] | −0.0105 [−0.0115, −0.0094] | −0.0105 [−0.0112, −0.0097] |  |
|  | 5–10 mm annulus at 40 mm | 0.994 [0.993, 0.994] | 0.994 [0.993, 0.995] | 0.994 [0.993, 0.994] |  |
|  | 10–20 mm annulus at 40 mm | 0.950 [0.947, 0.953] | 0.949 [0.947, 0.951] | 0.949 [0.946, 0.953] |  |
|  | 20–40 mm annulus at 40 mm | 0.926 [0.920, 0.933] | 0.927 [0.918, 0.935] | 0.929 [0.922, 0.937] |  |
|  | 40–80 mm annulus at 40 mm | 1.001 [0.966, 1.037] | 1.016 [0.979, 1.054] | 0.977 [0.937, 1.019] |  |
|  | 80–200 mm annulus at 40 mm | no ratio: non-zero in 0 of 8 MCsquare and 0 of 8 TOPAS runs | no ratio: non-zero in 0 of 8 MCsquare and 0 of 8 TOPAS runs | no ratio: non-zero in 0 of 8 MCsquare and 0 of 8 TOPAS runs |  |
|  | σ at 60 mm (mm) | −0.0352 [−0.0361, −0.0343] | −0.0347 [−0.0359, −0.0335] | −0.0347 [−0.0357, −0.0338] |  |
|  | 5–10 mm annulus at 60 mm | 0.991 [0.991, 0.991] | 0.991 [0.990, 0.992] | 0.991 [0.991, 0.992] |  |
|  | 10–20 mm annulus at 60 mm | 0.826 [0.825, 0.827] | 0.827 [0.825, 0.828] | 0.826 [0.825, 0.828] |  |
|  | 20–40 mm annulus at 60 mm | 0.793 [0.788, 0.797] | 0.794 [0.791, 0.797] | 0.795 [0.791, 0.798] |  |
|  | 40–80 mm annulus at 60 mm | 0.448 [0.428, 0.469] | 0.454 [0.435, 0.473] | 0.455 [0.433, 0.479] |  |
|  | 80–200 mm annulus at 60 mm | no ratio: non-zero in 0 of 8 MCsquare and 0 of 8 TOPAS runs | no ratio: non-zero in 0 of 8 MCsquare and 0 of 8 TOPAS runs | no ratio: non-zero in 0 of 8 MCsquare and 0 of 8 TOPAS runs |  |
| 150 MeV | R80 (mm) | −0.5397 [−0.5407, −0.5388] | −0.5401 [−0.5415, −0.5386] | −0.5399 [−0.5417, −0.5382] |  |
|  | σ at 80 mm (mm) | −0.0047 [−0.0058, −0.0035] | −0.0048 [−0.0055, −0.0041] | −0.0042 [−0.0053, −0.0030] |  |
|  | 5–10 mm annulus at 80 mm | 0.995 [0.995, 0.996] | 0.995 [0.995, 0.996] | 0.995 [0.995, 0.996] |  |
|  | 10–20 mm annulus at 80 mm | 1.031 [1.029, 1.034] | 1.032 [1.029, 1.035] | 1.034 [1.031, 1.036] |  |
|  | 20–40 mm annulus at 80 mm | 1.005 [1.001, 1.008] | 1.006 [1.002, 1.009] | 1.005 [1.000, 1.010] |  |
|  | 40–80 mm annulus at 80 mm | 0.784 [0.776, 0.792] | 0.784 [0.777, 0.790] | 0.788 [0.780, 0.797] |  |
|  | 80–200 mm annulus at 80 mm | 0.487 [0.462, 0.514] | 0.505 [0.483, 0.527] | 0.505 [0.488, 0.522] |  |
|  | σ at 125 mm (mm) | −0.0341 [−0.0354, −0.0328] | −0.0338 [−0.0350, −0.0326] | −0.0333 [−0.0348, −0.0319] |  |
|  | 5–10 mm annulus at 125 mm | 1.016 [1.015, 1.016] | 1.016 [1.015, 1.016] | 1.016 [1.015, 1.017] |  |
|  | 10–20 mm annulus at 125 mm | 0.850 [0.848, 0.851] | 0.851 [0.849, 0.852] | 0.850 [0.849, 0.851] |  |
|  | 20–40 mm annulus at 125 mm | 0.882 [0.879, 0.884] | 0.882 [0.880, 0.885] | 0.881 [0.878, 0.884] |  |
|  | 40–80 mm annulus at 125 mm | 0.855 [0.849, 0.862] | 0.853 [0.848, 0.858] | 0.853 [0.849, 0.857] |  |
|  | 80–200 mm annulus at 125 mm | 0.332 [0.313, 0.352] | 0.326 [0.304, 0.349] | 0.317 [0.293, 0.343] |  |
| 200 MeV | R80 (mm) | −0.5117 [−0.5147, −0.5087] | −0.5111 [−0.5139, −0.5083] | −0.5093 [−0.5126, −0.5060] | [−0.30, +0.30] |
|  | σ at 100 mm (mm) | +0.0020 [+0.0009, +0.0031] | +0.0015 [+0.0003, +0.0027] | +0.0020 [+0.0008, +0.0032] | [−0.10, +0.10] |
|  | 5–10 mm annulus at 100 mm | 0.995 [0.995, 0.996] | 0.995 [0.995, 0.996] | 0.995 [0.995, 0.996] |  |
|  | 10–20 mm annulus at 100 mm | 1.032 [1.030, 1.035] | 1.031 [1.027, 1.034] | 1.033 [1.030, 1.035] |  |
|  | 20–40 mm annulus at 100 mm | 1.063 [1.061, 1.065] | 1.064 [1.061, 1.066] | 1.064 [1.062, 1.065] | [0.90, 1.10] |
|  | 40–80 mm annulus at 100 mm | 0.969 [0.966, 0.973] | 0.969 [0.966, 0.973] | 0.969 [0.962, 0.975] | [0.90, 1.10] |
|  | 80–200 mm annulus at 100 mm | 0.812 [0.797, 0.828] | 0.797 [0.782, 0.812] | 0.801 [0.784, 0.818] | [0.75, 1.25] |
|  | σ at 200 mm (mm) | +0.0040 [+0.0020, +0.0060] | +0.0043 [+0.0026, +0.0060] | +0.0030 [+0.0012, +0.0047] | [−0.15, +0.15] |
|  | 5–10 mm annulus at 200 mm | 1.031 [1.030, 1.032] | 1.031 [1.030, 1.031] | 1.031 [1.030, 1.031] |  |
|  | 10–20 mm annulus at 200 mm | 0.939 [0.938, 0.941] | 0.939 [0.938, 0.940] | 0.938 [0.937, 0.939] |  |
|  | 20–40 mm annulus at 200 mm | 0.906 [0.903, 0.908] | 0.906 [0.903, 0.908] | 0.905 [0.903, 0.907] | [0.90, 1.10] |
|  | 40–80 mm annulus at 200 mm | 0.887 [0.883, 0.891] | 0.886 [0.883, 0.890] | 0.888 [0.884, 0.892] | [0.90, 1.10] |
|  | 80–200 mm annulus at 200 mm | 0.896 [0.884, 0.908] | 0.898 [0.884, 0.913] | 0.898 [0.884, 0.912] | [0.75, 1.25] |
<!-- END GENERATED: table-6-topas-halo -->

- **Range.** MCsquare's R80 is shorter than TOPAS's by 0.42 mm at 100 MeV, 0.54 mm at 150 MeV and 0.51 mm at
  200 MeV; the three arms agree within 0.003 mm at each energy. At 200 MeV the difference lies outside TD's ±0.30 mm
  band. In the planning runs at 200 MeV (4 × 1e7 per code, #32) the shift was close to rigid: peak −0.46, R90 −0.46,
  R80 −0.51, R50 −0.56, R20 −0.57 mm, and it persisted with nuclear interactions off in both codes (−0.56 mm). The
  water stopping power used by MCsquare matches Geant4's option 0 table to within 0.006% over 50–400 MeV, and the
  integrated ranges agree within about 0.01 mm (TD stage 0). **Most of the gap depends on TOPAS's step size**
  (below).
- **Spot size.** σ differs by at most 0.036 mm. At 200 MeV the differences are at most 0.005 mm and inside TD's
  bands; the largest are at the deeper slab at 100 and 150 MeV (about −0.035 and −0.034 mm).
- **Annuli at 200 MeV.** The ratios that have a TD band are inside it except the 40–80 mm annulus at 200 mm depth (0.886 to
  0.888, intervals entirely below 0.90). The 20–40 mm annulus at that depth is 0.905 to 0.906, and the 80–200 mm
  annulus is 0.80 to 0.81 at 100 mm and 0.90 at 200 mm. Each is within 0.02 of the planning value it replaces.
- **Annuli at 100 and 150 MeV.** At the deeper slab every ratio beyond 10 mm is below 0.9, and the lowest is in the
  outermost annulus that holds energy: 0.45 (100 MeV, 60 mm depth, 40–80 mm) and 0.32 to 0.33 (150 MeV, 125 mm, 80–200 mm).
  At the shallower slab the 80–200 mm annulus at 150 MeV is 0.49 to 0.51, and the 40–80 mm annulus at 100 MeV is
  consistent with 1 (0.98 to 1.02, intervals about ±4%). These annuli hold a small share of the slab energy (§5.0),
  and the ratios are of fractions of each code's own slab energy, not of absolute deposited energy.
- **No energy scored.** The 80–200 mm annulus at 100 MeV is zero in every run of both codes at both depths, so it
  has no ratio.

At 200 MeV, upstream with the fix and Portable are equivalent in every annulus, including the two outer ones
(Table 3). So, by inference from two separate comparisons, the lower far-halo fractions against TOPAS are common
to both MCsquare code lines and were not introduced by the port. That inference is not available for the 40–80 mm
annulus at 100 MeV and the 80–200 mm annulus at 150 MeV, where the same-host comparison is inconclusive (§5.0). None of this says anything about the
cause: physics models, scoring, estimators, and geometry or source conditions all remain possible.

**Diagnostic runs (clement-7074f29f; records in `validation/topas/stepsize/`).** These tests were not part of the
frozen addendum, and Table 6 does not use them. Each prediction was committed before its runs. All runs use the
TOPAS setup frozen for the addendum. The arms without a step limit reproduce the addendum's 1e7 means.

*Step size.* TOPAS was rerun with a maximum step size in the phantom (`d:Ge/Phantom/MaxStepSize`): 1e6 histories,
two seeds per arm, with seed pairs agreeing within 0.013 mm.

**Table 6b. TOPAS R80 against its step limit, and TOPAS − MCsquare (mm). Mean of two runs of 1e6 histories.
Diagnostic.**

| energy | step limit | TOPAS R80 | change from no limit | TOPAS − MCsquare |
|---|---|---|---|---|
| 100 MeV | none | 77.824 | — | +0.43 |
| 100 MeV | 0.1 mm | 77.556 | −0.27 | +0.16 |
| 100 MeV | 0.05 mm | 77.523 | −0.30 | +0.13 |
| 150 MeV | none | 158.745 | — | +0.54 |
| 150 MeV | 0.1 mm | 158.327 | −0.42 | +0.12 |
| 150 MeV | 0.05 mm | 158.272 | −0.47 | +0.07 |
| 200 MeV | none | 260.866 | — | +0.51 |
| 200 MeV | 0.5 mm | 260.694 | −0.17 | +0.34 |
| 200 MeV | 0.1 mm | 260.436 | −0.43 | +0.08 |
| 200 MeV | 0.05 mm | 260.392 | −0.47 | +0.04 |

- **Most of the R80 gap depends on the reference code's step size.** With a 0.1 mm step limit the gap falls from
  0.43–0.54 mm to 0.08–0.16 mm at all three energies. Going from 0.1 mm to 0.05 mm moves R80 by a further −0.033,
  −0.055 and −0.044 mm (100, 150 and 200 MeV), within the small change predicted for an effect that has nearly
  converged. Two step sizes cannot bound what remains below 0.05 mm. At 0.05 mm the gap is 0.04 to 0.13 mm.
- **No mechanism is claimed.** The Geant4 11.3.2 source applies a linear energy-loss approximation when a step loses
  less than 1% of the kinetic energy (`linLossLimit` 0.01, not changed by option 0 or TOPAS). In this phantom, steps
  end at 1 mm voxel boundaries, so that branch acts above about 83 MeV. An approximate model of it predicted a shift
  of −0.10 mm at 100 MeV, where −0.27 mm was measured, so that model does not account for the result quantitatively.
  The branch was not varied on its own, and other step-dependent processes are not excluded.
- **Where it was checked, the halo did not respond to the step limit.** At 200 MeV no appreciable change was
  observed in the spot σ or in the 40–80 mm annulus fraction at 200 mm depth, against a large change in R80. Other
  annuli and energies were not checked.
- **The deficit is in the scored charged-particle channel.** TOPAS's `Dose` scorer excludes neutrons, gammas and
  their descendants; its unfiltered `DoseAll`, run alongside, gives larger far-ring fractions of the slab energy in
  these runs (geometric-mean ratio over four runs: 1.09 for the 40–80 mm annulus at 200 MeV, 1.20–1.21 at 150 MeV,
  2.5–3.0 at 100 MeV, and 2.0–15 for the 80–200 mm annulus where `Dose` scored energy; at 100 MeV `Dose` scored none
  there). These are ratios of normalized fractions, not of absolute dose. They show that the deficit remains in the
  channel compared here; they do not validate the scorer's particle classification.

*Nuclear elastic.* At 200 MeV, TOPAS was run without hadronic elastic scattering (`g4h-elastic_HP` removed), 1e7
histories × 2, against the addendum's 8 runs. Without it the 20–40 and 40–80 mm annuli at 200 mm depth fall to 0.39
and 0.52 of their values with it, and the 40–80 mm annulus at 100 mm to 0.73. R80 moves +0.31 mm and σ at 200 mm
−0.18 mm. These ratios are of annulus fractions of the slab energy, not of absolute annular dose, and removing a
process changes transport as a whole, so they measure how sensitive the halo fractions are to elastic scattering, not
the share of the halo that elastic scattering deposits. That sensitivity is large next to MCsquare's 0.887 of TOPAS's
fraction in the 40–80 mm annulus at 200 mm (Table 6). A difference in how the two codes model elastic scattering is
therefore a candidate for the halo gap. This test removes elastic scattering rather than matching it to MCsquare's
model, so it does not show that elastic modelling is the cause.
It also does not explain the range gap: removing elastic makes TOPAS's range longer, and it is already the longer one.

### 5.5 Broad field at 100 and 150 MeV (part E, confirmatory)

§5.0 left open whether the field-edge agreement at 200 MeV holds at lower energies. Part E answers that with a
follow-up acquisition named before it ran (`docs/field_100_150_plan.md`).

**Design.** The broad-field model of §3.2 (15 × 15 cm, no range shifter) at 100 and 150 MeV. Each arm ran 8 runs per
energy at 6e7 histories, with the hosts, thread counts and code lines of §4: A-up and A-port on the Lenovo; B-up,
B-pgcc and B-picc on the HP. There are 17 endpoints per energy and 34 per contrast. Lateral dose is taken at 5, 10, 20,
30, 50 and 70 mm outside the field edge, at two depths per energy, in percentage points of the central dose of the same
profile. Central-axis dose is taken at those depths and at the maximum, as ratios, and R80 and R20 as differences in
mm. The margins are those of the 200 MeV field (#31). The analysis was written and reviewed before any part E run
existed, and frozen with the acquisition.

**Acquisition and amendments.** The full run was at commit `df2fabbd`. Two cells, A-port and B-pgcc at 150 MeV, were
stopped at 6 of 8 runs by the runners' default 3 h job timeout. With the timeout raised, both cells were rerun in full
at commit `b4b3bce3`, which differs from `df2fabbd` only in the two run requests and the plan (acquisition amendment
1). That amendment was recorded in the run request that made it (on the run branch, merged with this report), before any
dose value of the full run was read. The analysis takes those two cells from the rerun only and verifies each rerun
run against its own commit (analysis amendment 2, #66). It was reviewed and merged into main (`26a93235`) before any dose
value was read. Then all 80 runs verified (`--status`), and only then were endpoints computed, at `26a93235` (document
`validation/report_data/field_e_analysis_26a93235.json`). The order is recorded in the plan's chronology.

**Table 7. Part E, confirmatory contrasts: claims and counts.**

<!-- BEGIN GENERATED: table-7-part-e-summary -->
| contrast | joint claim (all 34 equivalent) | secondary joint claim | equivalent, unadjusted | equivalent, Holm |
|---|---|---|---|---|
| A: A-port vs A-up | established | established | 34 of 34 | 34 of 34 |
| B1: B-pgcc vs B-up | established | established | 34 of 34 | 34 of 34 |
| B2: B-picc vs B-up | established | established | 34 of 34 | 34 of 34 |
<!-- END GENERATED: table-7-part-e-summary -->

**Table 8. Part E, every endpoint: estimate [95% interval] and TOST outcome (90% interval). `lateral_D_X` is the dose X
mm outside the field edge at D mm depth, as a difference in percentage points of the central dose. `cax_D` and
`cax_max` are central-axis dose ratios, and `r80_mm` and `r20_mm` are differences in mm.**

<!-- BEGIN GENERATED: table-8-part-e-endpoints -->
| endpoint | scale | A | B1 | B2 |
|---|---|---|---|---|
| E100/lateral_39_5 | difference | −0.0202 [−0.0401, −0.0003] E | +0.0011 [−0.0239, +0.0261] E | −0.0072 [−0.0291, +0.0146] E |
| E100/lateral_39_10 | difference | +0.0175 [+0.0054, +0.0297] E | +0.0045 [−0.0125, +0.0215] E | +0.0037 [−0.0127, +0.0201] E |
| E100/lateral_39_20 | difference | +0.0013 [−0.0114, +0.0139] E | +0.0006 [−0.0121, +0.0132] E | +0.0090 [−0.0046, +0.0226] E |
| E100/lateral_39_30 | difference | +0.0010 [−0.0143, +0.0162] E | −0.0033 [−0.0097, +0.0030] E | −0.0008 [−0.0073, +0.0057] E |
| E100/lateral_39_50 | difference | −0.0023 [−0.0063, +0.0017] E | −0.0008 [−0.0056, +0.0040] E | +0.0008 [−0.0035, +0.0051] E |
| E100/lateral_39_70 | difference | +0.0002 [−0.0031, +0.0035] E | +0.0025 [−0.0015, +0.0066] E | +0.0008 [−0.0024, +0.0040] E |
| E100/lateral_61_5 | difference | −0.0274 [−0.0483, −0.0064] E | +0.0034 [−0.0283, +0.0351] E | −0.0022 [−0.0284, +0.0240] E |
| E100/lateral_61_10 | difference | +0.0021 [−0.0152, +0.0194] E | +0.0059 [−0.0126, +0.0244] E | +0.0039 [−0.0134, +0.0212] E |
| E100/lateral_61_20 | difference | −0.0039 [−0.0181, +0.0103] E | −0.0075 [−0.0183, +0.0034] E | +0.0043 [−0.0070, +0.0156] E |
| E100/lateral_61_30 | difference | −0.0003 [−0.0099, +0.0094] E | +0.0022 [−0.0068, +0.0113] E | +0.0128 [+0.0036, +0.0220] E |
| E100/lateral_61_50 | difference | −0.0024 [−0.0082, +0.0034] E | +0.0045 [+0.0002, +0.0088] E | +0.0037 [−0.0005, +0.0080] E |
| E100/lateral_61_70 | difference | −0.0010 [−0.0050, +0.0030] E | +0.0014 [−0.0021, +0.0049] E | +0.0017 [−0.0021, +0.0056] E |
| E100/cax_39 | ratio | 0.9990 [0.9981, 1.0000] E | 1.0001 [0.9989, 1.0014] E | 1.0002 [0.9991, 1.0014] E |
| E100/cax_61 | ratio | 1.0000 [0.9987, 1.0014] E | 1.0001 [0.9988, 1.0013] E | 1.0003 [0.9989, 1.0017] E |
| E100/cax_max | ratio | 1.0000 [0.9988, 1.0011] E | 1.0000 [0.9987, 1.0014] E | 1.0008 [0.9994, 1.0022] E |
| E100/r80_mm | difference | +0.0004 [−0.0010, +0.0017] E | −0.0001 [−0.0013, +0.0011] E | 0.0000 [−0.0011, +0.0010] E |
| E100/r20_mm | difference | 0.0000 [−0.0014, +0.0014] E | +0.0001 [−0.0013, +0.0015] E | +0.0002 [−0.0010, +0.0015] E |
| E150/lateral_79_5 | difference | +0.0184 [−0.0098, +0.0466] E | −0.0154 [−0.0416, +0.0108] E | −0.0113 [−0.0392, +0.0166] E |
| E150/lateral_79_10 | difference | +0.0007 [−0.0229, +0.0242] E | −0.0040 [−0.0267, +0.0186] E | −0.0224 [−0.0445, −0.0002] E |
| E150/lateral_79_20 | difference | +0.0012 [−0.0122, +0.0146] E | +0.0004 [−0.0145, +0.0153] E | −0.0038 [−0.0192, +0.0117] E |
| E150/lateral_79_30 | difference | +0.0027 [−0.0039, +0.0092] E | −0.0013 [−0.0101, +0.0074] E | +0.0046 [−0.0037, +0.0128] E |
| E150/lateral_79_50 | difference | −0.0008 [−0.0032, +0.0016] E | −0.0012 [−0.0047, +0.0024] E | +0.0009 [−0.0022, +0.0039] E |
| E150/lateral_79_70 | difference | +0.0001 [−0.0015, +0.0016] E | +0.0001 [−0.0009, +0.0011] E | +0.0001 [−0.0012, +0.0014] E |
| E150/lateral_125_5 | difference | −0.0058 [−0.0425, +0.0310] E | +0.0113 [−0.0198, +0.0424] E | −0.0111 [−0.0411, +0.0189] E |
| E150/lateral_125_10 | difference | −0.0206 [−0.0471, +0.0059] E | −0.0116 [−0.0364, +0.0132] E | −0.0109 [−0.0310, +0.0092] E |
| E150/lateral_125_20 | difference | −0.0021 [−0.0163, +0.0121] E | +0.0050 [−0.0073, +0.0174] E | −0.0007 [−0.0123, +0.0109] E |
| E150/lateral_125_30 | difference | −0.0038 [−0.0117, +0.0042] E | −0.0025 [−0.0124, +0.0075] E | −0.0034 [−0.0131, +0.0063] E |
| E150/lateral_125_50 | difference | −0.0008 [−0.0025, +0.0010] E | −0.0002 [−0.0023, +0.0019] E | +0.0007 [−0.0021, +0.0035] E |
| E150/lateral_125_70 | difference | −0.0001 [−0.0008, +0.0006] E | +0.0001 [−0.0008, +0.0009] E | +0.0001 [−0.0005, +0.0006] E |
| E150/cax_79 | ratio | 1.0000 [0.9987, 1.0013] E | 0.9996 [0.9982, 1.0010] E | 0.9997 [0.9981, 1.0014] E |
| E150/cax_125 | ratio | 1.0005 [0.9993, 1.0017] E | 1.0001 [0.9983, 1.0019] E | 0.9996 [0.9980, 1.0013] E |
| E150/cax_max | ratio | 1.0012 [0.9997, 1.0027] E | 1.0000 [0.9987, 1.0013] E | 1.0001 [0.9990, 1.0011] E |
| E150/r80_mm | difference | −0.0029 [−0.0067, +0.0009] E | +0.0004 [−0.0033, +0.0042] E | −0.0027 [−0.0078, +0.0025] E |
| E150/r20_mm | difference | −0.0054 [−0.0118, +0.0011] E | +0.0013 [−0.0033, +0.0059] E | +0.0012 [−0.0042, +0.0065] E |
<!-- END GENERATED: table-8-part-e-endpoints -->

- **Every endpoint is equivalent in every confirmatory contrast, before and after Holm adjustment,** so the joint
  claim holds for Portable MCsquare on the Lenovo (Windows) and with both compilers on the HP (Linux). This is the
  claim that §5.0 could not establish for the pencil beam's far halo at these energies. The broad field's endpoints
  are better determined at 6e7 histories, and they are a different quantity.
- The largest lateral difference is 0.027 percentage points (95% intervals within ±0.05), against margins of ±0.5
  points at 5 mm and ±0.2 points further out. Central-axis ratios are within 0.12% of 1 (intervals within 0.27%),
  and R80 and R20 differ by at most 0.003 and 0.005 mm.
- **Runs are not byte-reproducible for a seed.** Two runs with the same seed, binary and host gave different dose
  files (`docs/field_100_150_plan.md`, chronology). That is a fact about reproducibility, not about independence. The
  analysis assumes that runs are independent draws, as it does throughout; the run lists' seed screen
  (`field_followup_runs.py`) is a screen, not a proof of independence.

### 5.6 One heterogeneous case with a range shifter (Lungman, descriptive)

Every case above is a homogeneous phantom with no beam-line device (§7). This subsection reports one case that has
both. **It is descriptive.** It was not pre-specified and has no margin and no outcome. It is one beam and one run
per code, and the two runs are on different hosts with different compilers.

**Case.** The Lungman anthropomorphic chest phantom CT (512 × 512 × 426 voxels, 0.625 × 0.625 × 0.7 mm) with a 9 mm
simulated tumour (`tumours_100HU_1`); CT, mask and geometry by cora-2f1e43dc. One pencil beam at 127 MeV with the
range shifter in (`BDL_default_UN_RangeShifter`), gantry 135°, one spot aimed at the tumour centroid; plan and inputs
by clement-7074f29f. The HU-to-density curve (`Scanners/default`) and the beam model are generic. 1e7 primaries per
run; MCsquare reports a statistical uncertainty of 1.8% for the Studio run.

**Arms.** Upstream OpenMCsquare `85bf2911` with the fix, icl, on the Lenovo (4 threads; the A-up binary of §3.3),
against Portable MCsquare, Apple clang, on the Mac Studio (12 threads). The Lenovo run read the Studio package and
checked every file against its sha256. The package was assembled from the Studio run's inputs after that run, and its
config re-runs it (tumour mean within 0.3% of the stored dose at 1e6). Differences are upstream − Portable; ratios are
upstream / Portable.

**Table 10. Lungman pencil beam, upstream with the fix (Lenovo) against Portable (Mac Studio). Dose in Gy per MU
from the same normalisation in both. Descriptive: no margin, no outcome.**

| quantity | upstream + fix | Portable | comparison |
|---|---|---|---|
| tumour mean dose (Gy/MU) | 0.03348 | 0.03354 | ratio 0.998 |
| tumour minimum (Gy/MU) | 0.02386 | 0.02356 | ratio 1.013 |
| tumour maximum (Gy/MU) | 0.05164 | 0.05157 | ratio 1.002 |
| high-dose centroid (voxels above 50% of the body maximum) | | | within 0.05 mm on each axis |
| body voxels above 10% of the body maximum (n = 277,113) | | | mean difference −0.02% of the maximum; 95th percentile of the absolute difference 2.1% |
| gamma 2%/2 mm, global, 10% cutoff, body voxels | | | 99.996% of 277,113 points pass |

"Body" is CT above −900 HU. Air voxels are excluded because dose per unit mass in near-massless voxels is noise (the
hottest voxel of either run is in air at the beam entrance).

- **The two code lines agree on this case to within the statistical noise of the runs.** The body-voxel differences
  are about what two independent runs at about that uncertainty give.
- **Both codes place the beam identically.** That includes the convention that the plan's isocentre x is mirrored
  about the CT width, which the port did not change.
- **No range comparison is made.** Distal to the tumour the beam is in lung, where the dose along the axis sits on a
  plateau near 80% of its peak. Where a single-voxel axis profile crosses 80% is therefore set by noise: the two
  profiles cross 3.4 mm apart, and that is not a range difference. A range comparison would need a laterally
  integrated depth dose on a water-equivalent axis.
- **What this does not test.** A single pencil beam does not probe the dose outside a field edge, where the upstream
  defect acts, and both arms carry the fix. Nothing here is compared with TOPAS or with measurement, and one case
  does not establish agreement in heterogeneous anatomy in general.

## 6. Discussion

1. **Agreement between the code lines.** On the same host, and with the same fix in both, upstream OpenMCsquare and
   Portable MCsquare agree within the pre-stated margins on range, spot size, the core and near halo of the pencil
   beam at all three energies, the far halo at 200 MeV, and every broad-field endpoint. In the far halo at 100 and
   150 MeV the comparison is inconclusive or not established: the precision is not enough to show equivalence, and
   one pointwise interval there excludes 1 (§5.0). The joint claim over all 52 endpoints is therefore not
   established. The pattern of outcomes is the same in all three contrasts (icl against MSYS2 gcc on Windows; icc
   against gcc, and icc against icc, on Linux). No isolated compiler effect is established; the descriptive
   comparison of the two Portable builds has two pointwise intervals that exclude the null value (§5.2).
2. **The upstream defect in clinical terms.** As distributed, upstream computes 13–18% less dose than the corrected
   code between 5 and 30 mm outside the edge of the 200 MeV field (§2), up to 0.95 percentage points of central-axis
   dose. As an illustration, for a field that delivers 60 Gy on the axis that is about 0.6 Gy at 5 mm outside the
   edge, in the region where organs at risk sit; it is a difference between computed doses under that assumption, not
   an observed patient dose. With the fix in both code lines, the two differ there by at most 0.12 points at the 95%
   level (Table 4), about 0.07 Gy on the same scale. The defect is roughly ten times larger than any difference
   between the corrected code lines. Whether the corrected dose is closer to physical dose has not been shown (§5.4).
3. **Portability.** Portable MCsquare gives equivalent results on Linux, Windows and macOS (§5.3) and, with free
   compilers, matches the Intel-built upstream within the margins above. Treatment-planning research with MCsquare
   need not depend on Intel compilers.
4. **The TOPAS differences.** The pre-specified descriptive comparison at 100, 150 and 200 MeV (#54, Table 6) shows
   a range offset at all three energies (0.42 to 0.54 mm) and, at 100 and 150 MeV, outer-annulus fractions down to
   about a third of TOPAS's. By inference both are common to the two code lines and not introduced by the port.
   **Range:** most of the offset depends on TOPAS's step size. With a 0.1 mm step limit it falls to 0.08–0.16 mm,
   and to 0.04–0.13 mm at 0.05 mm (Table 6b, diagnostic); the change between those two step sizes is small, but two
   step sizes do not bound the remainder. The offset in Table 6 therefore says more
   about the reference calculation's default stepping in a 1 mm voxel phantom than about MCsquare. No mechanism is
   claimed. **Far halo:** the step limit did not change the one annulus checked. In TOPAS the 20–80 mm halo fractions
   at 200 MeV are strongly sensitive to hadronic elastic scattering, so elastic modelling is a candidate for the gap;
   its cause is not established. TD's confirmatory arms on one host have not been run.
5. **Speed.** Table 9 gives the wall time per 1e7 histories of every run in the confirmatory acquisitions, from the
   runs' own records. It is resource information, not a performance claim: one run at a time per host, each host's
   thread count, no tuning. **On these hosts, cases, energies and thread counts, the Portable builds made with the free
   compilers (gcc on Linux, MinGW-w64 gcc on Windows) took 2.4 to 2.7 times as long as the Intel-built codes.** The
   Portable build made with icc took about as long as upstream (B-picc against B-up). Each arm is a whole build: source,
   compiler, flags and libraries together. So this compares these builds as they ran, and does not isolate the
   compiler's share or rule out a contribution from the port's changes.

   **Table 9. Wall time per 1e7 histories, seconds: median (range) over 8 runs. Pencil beam and the 200 MeV field:
   acquisition `2f9dab40` (§4). Field at 100 and 150 MeV: part E (§5.5), 6e7 histories per run. Lenovo: Core
   i5-6400T (4 threads); HP: Core i5-6500T (3 threads).**

<!-- BEGIN GENERATED: table-9-run-times -->
| arm | host, threads | pencil 100 MeV | pencil 150 MeV | pencil 200 MeV | field 200 MeV | field 100 MeV (E) | field 150 MeV (E) |
|---|---|---|---|---|---|---|---|
| A-up | Lenovo, 4 | 83 (82–84) | 148 (147–148) | 222 (221–223) | 163 (160–165) | 60 (60–62) | 109 (107–110) |
| A-port | Lenovo, 4 | 206 (205–206) | 370 (369–370) | 554 (553–555) | 411 (409–419) | 154 (154–154) | 276 (272–277) |
| B-up | HP, 3 | 90 (87–91) | 162 (161–165) | 248 (248–252) | 179 (179–180) | 64 (64–65) | 118 (118–119) |
| B-pgcc | HP, 3 | 218 (214–222) | 392 (391–398) | 595 (592–609) | 441 (439–450) | 163 (160–165) | 296 (294–301) |
| B-picc | HP, 3 | 84 (82–85) | 152 (152–154) | 236 (233–237) | 171 (170–173) | 60 (60–60) | 112 (112–112) |
<!-- END GENERATED: table-9-run-times -->

## 7. Limitations

1. **Gamma is a weak test for the pencil beam.** At 2%/2 mm and 1%/1 mm with a 10% cutoff it does not detect a
   uniform 3% dose difference in 3D (§5.1). The pencil-beam comparison rests on the endpoint tables. All gamma
   values compare means of Monte Carlo runs, so they include statistical noise in both distributions.
2. **The joint equivalence claim is not established** in any contrast, because the far halo at 100 and 150 MeV lacks
   precision at 1e7 histories per run (§5.0). For the pencil beam, the clinical impact of those unresolved endpoints
   is not established. The broad field, which tests out-of-field dose directly, is equivalent at 100, 150 and 200 MeV
   (Table 4, §5.5).
3. **Homogeneous phantoms only.** Neither case has a density interface. Differences in lateral scattering and in the
   nuclear halo matter most clinically behind low-density tissue and at bone–air interfaces, and none of that is
   tested here. Planned as future work: a lung-density slab; a sinus-like cavity of air in bone; and the same cavity
   filled in steps, as in congestion. One anthropomorphic lung case has been run descriptively (§5.6); it is one beam
   and carries no claim.
4. **Neither confirmatory case has a beam-line device in the beam.** The upstream defect was found with a range
   shifter in, where nuclear secondaries from the shifter reach the phantom, so a range-shifter-in case would test it
   most directly. No aperture is modelled. The descriptive lung case of §5.6 has a range shifter in the beam, but
   it is a single pencil beam with the fix in both arms, so it does not test the defect.
5. **Energies and fields.** Pencil beams and one 15 × 15 cm field, each at 100, 150 and 200 MeV. Nothing outside
   that range is tested.
6. The TOPAS comparison is descriptive and covers the pencil beam only (100, 150 and 200 MeV); TOPAS and the
   Portable arms ran on different hosts, and the MCsquare endpoints had been read before its design was fixed.
   There is no comparison with measurement. The step-size diagnostic (Table 6b) is not pre-specified, uses 1e6
   histories and two runs per arm, and checks the halo in one annulus at one energy. The elastic test is at 200 MeV
   only, and it removes elastic scattering rather than matching it to MCsquare's model, so it cannot show the cause of
   the far-halo difference.
7. The broad-field phantom is Schneider_AT_AG_SI4 rather than water (§3.2).
8. **Provenance.** The collector verifies the Dose files against their recorded hashes but does not tie the endpoint
   values to them; for this collection a reviewer's recomputation on all 160 runs does (§4; #52). Each run's binary
   hash is checked against the snapshot of its own job, and the five arms carry five distinct binaries; nothing in
   the collection ties a hash to the build record.
9. Published agreement figures (Huang 2018): the four that the TOPAS design quotes were checked against the paper's
   full text (PMC6123159) on 2026-10-02 and are as quoted. Three are TOPAS against measurement (R80 differences
   generally under 0.1 mm on average; range within 0.6 mm; spot size 0.1 ± 0.1 mm). The fourth is MCsquare against
   TOPAS on one lung plan: 99.2% of iCTV voxels within 3% of the prescription dose. No comparison with them is
   made here.

## 8. Reproducibility and data

Sources, case generators, the workflows that produced every result, and the collection and analysis tools are in the
Portable MCsquare repository.

- **Run trees** (inputs, logs, records and Dose files, about 30 GB): two archives, part A sha256
  `ba6b8a3cdb06b6a0bc76d4ff514d3c750bdef917d05781e21ad145d9f8a74c7f` and part B sha256
  `33d042bb664277a15fcc14973ca71853b2e7202c2c02337c8a0dda3974e0560b`.
- **Collected endpoint records** (865 files with the collection manifest): archive sha256
  `e71859561b820bff7f887e8987c5c190542f45151ebc2f86cb9c77fecb8a7570`.
- **Analysis document:** `validation/report_data/apples_analysis_2f9dab40.json`, written by
  `validation/apples_analyse.py <collected> --parts A,B --json <file>` at commit `bad1419`.
- **Descriptive document and figures:** `validation/report_data/apples_descriptive_2f9dab40.json` and
  `docs/figures/`, written by `validation/report_dose_descriptives.py` from the run trees (pymedphys 0.41.0 for
  gamma). The annulus shares in §5.0 are in `validation/report_data/apples_ring_fractions_2f9dab40.json`, written by
  `validation/report_ring_fractions.py` from the collected records.
- **TOPAS halo document:** `validation/report_data/halo_868d1910.json` (and its table, `halo_868d1910.md`), written
  by `validation/topas_halo_compare.py` at commit `868d1910` from the 24 TOPAS runs frozen at `5aedf3ac` and the
  collected records above (321 files, dataset fingerprint sha256
  `0d5be975b294e0b77860e0ddbce5807caa526a0ca56481df5056effff40033ba`).
- **Part E** (§5.5): plan, amendments and chronology in `docs/field_100_150_plan.md`. Run trees: the full run of
  `df2fabbd`, archives part A sha256 `f6e42b7b8b47e52dba76a116373e26a393660472542051bc7b050f4eb787d94e` and part B
  sha256 `20a77f25651a6102758782a47be98af7622b845f1ec35fc9945cbd964bb95002`; the rerun of the two 150 MeV cells at
  `b4b3bce3`, part A sha256 `3bc98527d1cdde66242c6777bdb39ea450c72766d1902bb9a7b859913fa7d205` and part B sha256
  `da834d923894d6574fb95229b45aa8bcba5852b3ead85b9ee84b2c803fa37c16`. Analysis document
  `validation/report_data/field_e_analysis_26a93235.json` (and its table, `.md`), written by
  `validation/field_followup_analyse.py` at `26a93235` from all four (dataset fingerprint sha256
  `891933e22ce1e126076f22346a2bdf752e2c4008d60d6877710121d82fab0685`, 24,020 files).
- **Lungman case (§5.6, Table 10):** inputs in the Studio package `studio_E127_1e7` (each file with its sha256); the
  Lenovo run by `.gitea/workflows/lungman-lenovo.yml` at `cb740674` (archive sha256
  `3fcf204b46aef3067063c5c8bf13c2c96e5caec2865c0cd6ad5b2328aa817b9f`); the comparison in
  `validation/lungman/lenovo_cb740674c38f_vs_studio.json`, written by `validation/lungman_compare.py` (pymedphys 0.41.0
  for gamma). In that file the key "lenovo" is the upstream run and "studio" the Portable run.
- **Run times (Table 9):** `validation/report_data/run_times.json`, written by `validation/report_run_times.py` from the
  160 run records of the collected endpoint records above and the run records of part E's 80 analysed runs (the two
  rerun cells at `b4b3bce3`).
- **Diagnostic runs (Table 6b, elastic test):** predictions (committed before the runs), results, run scripts and the
  `DoseAll` ring data are in `validation/topas/stepsize/`. The TOPAS outputs were run at commit `6c4aaff81994`.
- **Tables:** `python validation/report_confirmatory_tables.py --check` confirms that the tables in this report are
  what those documents render.

The archives are held on the project's storage. A read-only link and a step-by-step reproduction guide will be added
before any wider circulation.

## 9. Contributors

sjswerdloff (direction, margins, decisions); connor-227743e6 (Portable MCsquare, MCsquare arms, upstream report);
clement-7074f29f (TOPAS installation, arm and diagnostics); alden-ec2221c7 (statistics and review); cora-2f1e43dc and
others (reviews). Authorship of any submission is sjswerdloff's decision, with each contributor's consent.

## Appendix A. Portable gcc against Portable icc (HP, Linux), descriptive

**Table A1.** B-pgcc against B-picc: estimate and 95% interval for each of the 52 endpoints. No margin and no claim
(AP). Differences are gcc − icc; ratios are gcc / icc.

<!-- BEGIN GENERATED: table-a1-compiler -->
| endpoint | Portable gcc − Portable icc, or ratio gcc / icc |
|---|---|
| P100/R80 | +0.0005 [−0.0007, +0.0018] |
| P100/sigma_40 | 0.0000 [−0.0009, +0.0010] |
| P100/sigma_60 | 0.0000 [−0.0012, +0.0013] |
| P100/ring_40_5_10 | 0.9999 [0.9991, 1.0007] |
| P100/ring_40_10_20 | 0.9996 [0.9960, 1.0031] |
| P100/ring_40_20_40 | 0.9972 [0.9870, 1.0074] |
| P100/ring_40_40_80 | 1.0394 [1.0030, 1.0772] |
| P100/ring_40_80_200 | not established |
| P100/ring_60_5_10 | 0.9997 [0.9990, 1.0005] |
| P100/ring_60_10_20 | 1.0004 [0.9983, 1.0026] |
| P100/ring_60_20_40 | 0.9991 [0.9952, 1.0029] |
| P100/ring_60_40_80 | 0.9964 [0.9470, 1.0484] |
| P100/ring_60_80_200 | not established |
| P150/R80 | −0.0001 [−0.0021, +0.0019] |
| P150/sigma_80 | −0.0006 [−0.0017, +0.0006] |
| P150/sigma_125 | −0.0005 [−0.0018, +0.0008] |
| P150/ring_80_5_10 | 1.0001 [0.9994, 1.0007] |
| P150/ring_80_10_20 | 0.9982 [0.9955, 1.0008] |
| P150/ring_80_20_40 | 1.0008 [0.9966, 1.0051] |
| P150/ring_80_40_80 | 0.9939 [0.9817, 1.0062] |
| P150/ring_80_80_200 | 0.9998 [0.9530, 1.0490] |
| P150/ring_125_5_10 | 0.9998 [0.9992, 1.0004] |
| P150/ring_125_10_20 | 1.0006 [0.9989, 1.0022] |
| P150/ring_125_20_40 | 1.0013 [0.9984, 1.0042] |
| P150/ring_125_40_80 | 0.9997 [0.9936, 1.0059] |
| P150/ring_125_80_200 | 1.0287 [0.9376, 1.1287] |
| P200/R80 | −0.0019 [−0.0056, +0.0019] |
| P200/sigma_100 | −0.0005 [−0.0019, +0.0009] |
| P200/sigma_200 | +0.0013 [−0.0006, +0.0033] |
| P200/ring_100_5_10 | 1.0000 [0.9993, 1.0007] |
| P200/ring_100_10_20 | 0.9979 [0.9950, 1.0009] |
| P200/ring_100_20_40 | 1.0000 [0.9975, 1.0024] |
| P200/ring_100_40_80 | 1.0007 [0.9944, 1.0071] |
| P200/ring_100_80_200 | 0.9947 [0.9727, 1.0171] |
| P200/ring_200_5_10 | 0.9999 [0.9995, 1.0003] |
| P200/ring_200_10_20 | 1.0011 [1.0001, 1.0021] |
| P200/ring_200_20_40 | 1.0005 [0.9981, 1.0030] |
| P200/ring_200_40_80 | 0.9980 [0.9937, 1.0023] |
| P200/ring_200_80_200 | 1.0006 [0.9831, 1.0183] |
| F/lateral_127_5 | +0.033 [−0.004, +0.071] |
| F/lateral_127_10 | +0.014 [−0.022, +0.050] |
| F/lateral_127_20 | +0.001 [−0.026, +0.028] |
| F/lateral_127_30 | +0.003 [−0.025, +0.032] |
| F/lateral_201_5 | +0.012 [−0.059, +0.084] |
| F/lateral_201_10 | −0.009 [−0.059, +0.041] |
| F/lateral_201_20 | +0.002 [−0.046, +0.051] |
| F/lateral_201_30 | +0.001 [−0.026, +0.028] |
| F/cax_127 | 0.9995 [0.9973, 1.0018] |
| F/cax_201 | 1.0000 [0.9974, 1.0025] |
| F/cax_i23 | 0.9992 [0.9974, 1.0010] |
| F/r80_mm | 0.000 [−0.017, +0.018] |
| F/r20_mm | −0.001 [−0.015, +0.013] |
<!-- END GENERATED: table-a1-compiler -->
