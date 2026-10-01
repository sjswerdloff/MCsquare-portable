# Portable MCsquare: depth-dose and off-axis comparison with upstream OpenMCsquare and TOPAS

**DRAFT, 2026-10-01, restructured to the form of the published MCsquare validations (Souris et al. 2016; Huang et al.
2018). Not a publication.** Results marked PRELIMINARY come from planning runs made before the confirmatory designs
were frozen; they are not confirmatory and must not be quoted as final. The confirmatory same-host acquisition
(commit `2f9dab40`, §4) is running and its results are not yet in this draft.

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

Measured on the corrected code's own sampling convention (proton secondaries, 1e5 primaries): as-is, the cumulative
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
| — | Portable MCsquare | Apple clang + libomp | Mac Studio M3 Ultra, macOS (platform study, TOPAS stage) |

Each contrast compares two arms on the **same host** (AP), so build, compiler and code line are not confounded with
hardware.

### 3.4 Scoring and comparison metrics

- **Depth dose:** integrated depth dose (IDD) over the whole scored plane; R90, R80 and R20 by linear interpolation
  at the first distal crossing; distal fall-off (R80–R20).
- **Off-axis:** lateral profiles in x and y through the central axis at stated depths; for the pencil beam, σ from a
  voxel-integrated Gaussian fit over |x| ≤ 10 mm and the fraction of energy in annuli 5–10, 10–20, 20–40, 40–80 and
  80–200 mm from the axis. Depths per energy: 40 and 60 mm at 100 MeV, 80 and 125 mm at 150 MeV, 100 and 200 mm at
  200 MeV (AP). For the broad field: the dose at 5, 10, 20 and 30 mm outside the field edge, in % of central-axis dose,
  at mid-range depth.
- **Gamma analysis:** 3D gamma between arms on the mean dose of each arm:
  - 2%/2 mm and 1%/1 mm, global (normalised to the maximum dose), 10% low-dose threshold, as in the literature;
  - in addition, a **low-dose gamma** (global, about 1% threshold, or local) on the broad-field planes, because the
    region next to the field edge where upstream's defect acts lies below a 10% threshold.

  Every gamma table states its criteria and threshold.
- **Equivalence statistics:** beneath each comparison, the difference of means with Welch intervals and a two
  one-sided test (TOST) against pre-stated margins: R80 ±0.05 mm, σ ±0.02 mm, annuli out to 80 mm ×[0.98, 1.02],
  the 80–200 mm annulus ×[0.95, 1.05]. A joint claim per contrast is made by intersection–union, with Holm adjustment
  across each contrast's endpoints (AP; analysis reviewed in #48).
- **Random numbers:** Portable MCsquare uses PCG. Each thread `t` starts from state `RNG_Seed + 1e4·t + 1e5·Num_call`
  on stream `t`, with `Num_call` = 1 here. No two runs share a seed, so no two runs share a generator. Neither code is
  bitwise reproducible at more than one thread, so every comparison uses statistics over independent runs, not
  same-seed identity (#39).

## 4. Confirmatory acquisition (in progress)

The same-host comparison of §3.3 is running at commit `2f9dab40`: 8 runs per energy per arm for case P and 8 runs per
arm for case F, on the Lenovo (part A) and the HP (part B). The collection and analysis code was reviewed and approved
before any result was read (#48). Its results will fill §5.1 and §5.2.

## 5. Results

### 5.1 Depth dose (PENDING §4)

Per energy: IDD curves of each pair of arms overlaid, with their difference; a table of R90, R80, R20 and fall-off
differences (mm) with intervals; gamma pass rates.

**PRELIMINARY (different hosts):** at 200 MeV, upstream + fix (icl, Lenovo, 4 threads) against Portable (Mac Studio,
24 threads), 4 × 1e7 each, R80 differs by +0.0008 mm (SE 0.0038; 95% [−0.0084, +0.0101]).

### 5.2 Off-axis dose (PENDING §4)

Per energy and depth: lateral profiles overlaid, σ differences and annulus ratios with intervals; gamma pass rates
on the planes; for the broad field, the profiles across the field edge.

**PRELIMINARY, broad field, both code lines with the fix** (mean of 2 runs × 3e7, 3 threads; dose in % of central
axis at 127 mm depth):

| code, compiler | host, OS | 5 mm outside | 10 mm | 20 mm | 30 mm |
|---|---|---|---|---|---|
| upstream + fix, icc 2021.1.2 | HP, Linux (WSL2) | 6.97 | 4.64 | 2.33 | 1.23 |
| upstream + fix, icl 2021.1 | Lenovo, Windows 10 | 7.00 | 4.64 | 2.25 | 1.21 |
| Portable, MSYS2 gcc 16.2.0 | HP, Windows 10 | 6.98 | 4.69 | 2.33 | 1.23 |
| Portable, gcc 13 | HP, Linux (WSL2) | 6.99 | 4.69 | 2.30 | 1.21 |

Three rows agree within about 0.05 points in every column; the Lenovo 20 mm cell is about 0.08 lower, more than its
own two-run difference (0.04). With n = 2 this is descriptive only (#16, #28).

**PRELIMINARY, pencil beam, 200 MeV, different hosts** (as in §5.1; difference upstream − Portable, ratio
upstream/Portable; 95% intervals):

| depth | σ difference | annulus 20–40 mm | 40–80 mm | 80–200 mm |
|---|---|---|---|---|
| 100 mm | +0.0000 mm [−0.0018, +0.0018] | 0.9980 [0.9944, 1.0015] | 0.9978 [0.9917, 1.0040] | 1.0018 [0.9556, 1.0503] |
| 200 mm | −0.0017 mm [−0.0058, +0.0024] | 0.9994 [0.9965, 1.0023] | **1.0070 [1.0019, 1.0121]** | 1.0064 [0.9935, 1.0195] |

The 200 mm, 40–80 mm annulus differs by about 0.7%, with an interval that excludes 1. Here host, compiler,
random-number generator and thread count all differ, which is why §4 compares on the same host.

### 5.3 Portable MCsquare across platforms (study pe1)

Broad field (§3.2), 13 endpoints (central-axis dose, off-axis dose outside the field edge, R80 and R20), Linux on the HP
as reference, 16 runs per platform, margins set before acquisition (#31):

- Windows (HP and Lenovo pooled) against Linux: **equivalent on 13 of 13 endpoints.**
- macOS (Apple Silicon) against Linux: **equivalent on 13 of 13 endpoints.**

Largest point estimates: off-axis within ±0.03 points (margins ±0.2 / ±0.5), central axis within −0.08% (margin ±0.5%),
R80/R20 within ±0.01 mm (margin ±0.5 mm). The claim is against the Linux reference, for one field.

### 5.4 Portable MCsquare against TOPAS (PRELIMINARY)

OpenTOPAS 4.3.0 / Geant4 11.3.2, EM option 0; the TOPAS dose scorer excludes neutrons, gammas and their descendants.
Pencil beam (§3.1), 200 MeV, 4 × 1e7 per code (TD; #32). Difference Portable − TOPAS, 95% intervals:

| endpoint | difference | 95% |
|---|---|---|
| R80 | **−0.512 mm** | [−0.520, −0.504] |
| σ at 3 mm depth | −0.0021 mm | [−0.0038, −0.0004] |
| σ at 100 mm | +0.0026 mm | [+0.0008, +0.0044] |
| σ at 200 mm | +0.0037 mm | [−0.0003, +0.0077] |

Annulus energy-fraction ratios, Portable/TOPAS:

| depth | 5–10 mm | 10–20 | 20–40 | 40–80 | 80–200 |
|---|---|---|---|---|---|
| 100 mm | 0.9956 | 1.0346 | 1.0671 | 0.9720 | 0.7947 |
| 200 mm | 1.0304 | 0.9404 | 0.9048 | 0.8887 | 0.8898 |

Portable's range is about 0.5 mm shorter than TOPAS's at 200 MeV, and about 0.43 mm at 100 MeV (2 × 1e6). The shift is
close to rigid: peak −0.46, R90 −0.46, R80 −0.51, R50 −0.56, R20 −0.57 mm (#32). The water stopping power used by
MCsquare matches Geant4's option 0 table to within 0.006% over 50–400 MeV, and the integrated ranges agree within
about 0.01 mm (TD stage 0), so the gap is **not yet explained**. It persists with nuclear interactions off in both
codes (−0.56 mm; #32). Annuli far from the axis carry 10–20% less energy in Portable at 200 mm depth.

## 6. Discussion

To be written with the confirmatory results. Points it will address: whether the two code lines, with the upstream
fix applied to both, agree within the stated margins on the same host; the size and clinical relevance of upstream's
field-edge underdose; and the open range and halo differences against TOPAS.

## 7. Limitations

1. §5.1 and §5.2 currently hold preliminary data only; the confirmatory results (§4) are pending.
2. The TOPAS comparison is preliminary and covers the pencil beam only; there is no comparison with measurement.
3. **Neither confirmatory case has a beam-line device in the beam.** The upstream defect was found with a range
   shifter in, where nuclear secondaries from the shifter reach the phantom, so a range-shifter-in case would test it
   most directly.
4. The broad-field phantom is Schneider_AT_AG_SI4 rather than water (§3.2).
5. Published agreement figures (Huang 2018) are to be checked against the paper itself before any comparison with
   them.

## 8. Reproducibility and data

Sources, case generators, the workflows that produced every result, and the collection and analysis tools are in the
Portable MCsquare repository. A step-by-step reproduction guide and the dose files (with a committed sha256 manifest
and a read-only link) will be added with the confirmatory results.

## 9. Contributors

sjswerdloff (direction, margins, decisions); connor-227743e6 (Portable MCsquare, MCsquare arms, upstream report);
clement-7074f29f (TOPAS installation, arm and diagnostics); alden-ec2221c7 (statistics and review); cora-2f1e43dc and
others (reviews). Authorship of any submission is sjswerdloff's decision, with each contributor's consent.
