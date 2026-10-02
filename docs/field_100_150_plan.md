# Same-host comparison, part E: field edge at 100 and 150 MeV (PLAN)

**Status: DRAFT for review. Nothing here has run.** sjswerdloff authorised the work on 2026-10-02 ("We aren't
really using the intel pcs for anything else right now, so go ahead"; "yes, go ahead with all of the work on the
intel pcs. just make sure we don't run out of disk space on those pcs (that's what the SMBWritable space is
for)"). The run lists, seeds, histories, endpoints, margins and analysis are frozen at the acquisition commit,
before the first full run, and nothing is changed after results are read. It follows
[`apples_to_apples_plan.md`](apples_to_apples_plan.md) (parts A and B); where this file is silent, that plan's
rules apply.

## Why

Parts A and B compared Portable MCsquare with upstream OpenMCsquare (with the fix) on the same hosts. In each
contrast 46 of 52 endpoints were equivalent and six fell short: pencil-beam annuli far from the axis at 100 and
150 MeV (report §5.0). Two of them could not be evaluated at all, because the 80–200 mm annulus at 100 MeV scored
exactly zero in every run. The clinical impact of those six is not established (#45 review 7080).

A pencil-beam annulus is not what a patient receives. What matters clinically is the dose outside the edge of a
field, which sums the far annuli of every spot in it. Parts A and B ran one broad field, at 200 MeV, and every
out-of-field endpoint was equivalent. This part runs the same field at 100 and 150 MeV, and looks further from the
edge than part A and B did.

## Question

At 100 and 150 MeV, for a 15 × 15 cm field: is the dose outside the field edge, the central-axis dose and the
range the same in Portable and in upstream with the fix, within the margins already used for the 200 MeV field?

It is one case, in a uniform medium, without a range shifter. It says nothing about heterogeneity or a range
shifter, and it does not resolve the pencil-beam annuli themselves.

## Arms, hosts and binaries

The arms, hosts and thread counts are those of parts A and B. **The binaries are the same files**: no arm is
rebuilt. Portable's source is unchanged since `2f9dab40` (nothing outside `validation/`, `docs/` and `.gitea/` differs
between `2f9dab404cea02f5072352f7ae5773dca27eb2a1` and this branch's base; the run job checks this against its own
commit and refuses otherwise).

| arm | build | host, threads | binary sha256 (as recorded by every part A or B run) |
|---|---|---|---|
| A-up | OpenMCsquare 85bf2911 + fix, icl 2021.1 | Lenovo, Windows 10, 4 | `f3a28398a399224c1e20d6e54e945571f6ee17693f4380ee7ad9c0e3ef84b966` |
| A-port | Portable `2f9dab40`, MSYS2 gcc 16.2 | Lenovo, Windows 10, 4 | `de93714ba856fc30605b1533247d5fa2be4e17c4c3e7cdc296bb5317c2aa60f4` |
| B-up | OpenMCsquare 85bf2911 + fix, icc 2021.1 | HP, Linux (WSL2), 3 | `b9106df25839a0e295565c0845dfbaff3ac3c80efa7f839d3f764eb9c87d768b` |
| B-pgcc | Portable `2f9dab40`, gcc 13.3 | HP, Linux (WSL2), 3 | `a21d406efac7fd67c8405df1bb86a86b3f419cb57a0ee82e12750e5824e16b7f` |
| B-picc | Portable `2f9dab40`, icc 2021.1 | HP, Linux (WSL2), 3 | `16e5d89911eba190cb9560be90c3e4b395c55399c5da6236362a629892e45ac1` |

The workflows have no build job. Each run job verifies its binary against the hash above and refuses on a
mismatch.

**Contrasts.** Confirmatory: A-port vs A-up; B-pgcc vs B-up; B-picc vs B-up. Descriptive: B-pgcc vs B-picc.

## Case

The #31 field-edge case, at two more energies: `validation/field_edge_make_cases.py 100,150 150`, which writes the
plans `E100_S150` and `E150_S150` and the same CT as at 200 MeV.

- Field 15 × 15 cm, spots on a 5 mm grid at equal weight, gantry 0, range shifter out.
- CT 150 × 150 × 150 voxels of 2 mm (300 mm cube), HU 0 throughout, which the default scanner maps to
  Schneider_AT_AG_SI4 (not water).
- Config as parts A and B, except the plan name and the number of histories.

## Endpoints: 17 per energy, 34 per contrast (fixed now)

All are computed by `validation/platform_study_metrics.py`, with the same estimators as at 200 MeV; the depths and
offsets become arguments.

| endpoints | 100 MeV | 150 MeV | definition |
|---|---|---|---|
| lateral dose, 12 | depths 39 and 61 mm | depths 79 and 125 mm | dose 5, 10, 20, 30, 50 and 70 mm outside the field edge, as a percentage of the central-axis dose at the same depth |
| central-axis dose, 3 | at 39 mm, at 61 mm, maximum | at 79 mm, at 125 mm, maximum | absolute dose in the central 20 × 20 mm; "maximum" is the largest value of that depth-dose curve |
| range, 2 | R80, R20 | R80, R20 | distal depths where the central-axis curve falls to 80% and 20% of its maximum |

- **Depths.** The voxel rows are centred on odd millimetres. 39 and 61 mm, and 79 and 125 mm, are the rows
  nearest the slab depths of the pencil-beam comparison (40 and 60 mm; 80 and 125 mm), which are about half and
  about 0.8 of the range. The lateral profile averages five rows (±4 mm), the four sides of the field and a 40 mm
  strip along each edge, as at 200 MeV.
- **Offsets.** 5 to 30 mm are the 200 MeV offsets. **50 and 70 mm are new.** The dose at a point outside the edge
  comes from annuli at least that far from each spot, so the offsets up to 30 mm are fed mainly by the annuli that
  were already equivalent in the pencil-beam comparison; 50 and 70 mm are fed by the 40–80 and 80–200 mm annuli,
  which are the ones that fell short. 70 mm is the furthest this phantom allows (the voxel centred 145 mm from the
  axis; the phantom ends at 150 mm). That voxel is two from the phantom's face, so it lacks scatter from beyond
  150 mm and is where a difference in boundary handling would land. Both arms have the same boundary, so the
  comparison is valid; it is not the dose 70 mm outside a field in a patient.
- **The maximum** replaces the 200 MeV case's fixed peak row (voxel index 23), which was chosen from a pilot. No
  pilot is run here.
- **No value is read before the freeze.** The smoke runs (below) are checked for completion only.

## Margins

Those of the 200 MeV field (#31, `validation/win_linux_equivalence_design.md`), unchanged:

| endpoint | margin |
|---|---|
| lateral at 5 mm | ±0.50 points of the central-axis dose |
| lateral at 10, 20, 30, 50 and 70 mm | ±0.20 points of the central-axis dose |
| central-axis dose (three) | ratio within [0.995, 1.005] |
| R80, R20 | ±0.5 mm |

The two new offsets take the ±0.20 margin of the other out-of-field points. The lateral margins are absolute, in
percent of the in-field dose, so they have a clinical reading that a ratio of two very small doses would not.

## Pre-specified analysis

As parts A and B: per endpoint, the difference of run-level means (arm − reference) for the lateral and range
endpoints and the geometric-mean ratio for the central-axis doses; Welch standard error and degrees of freedom;
TOST at one-sided α = 0.05; three outcomes (equivalent, not equivalent, inconclusive).

- **Joint claim, per confirmatory contrast:** equivalent on all 34 endpoints, each at unadjusted α = 0.05
  (intersection–union; no multiplicity adjustment).
- **Individual claims:** Holm within each confirmatory contrast, a family of 34.
- **A complete, verified population only.** A missing, duplicated, unexpected or unverifiable run makes the
  contrast PARTIAL: descriptive estimates, no joint claim, no Holm decisions.
- **Invalid values.** A non-finite or non-positive central-axis dose, an R80 or R20 without exactly one crossing,
  or a non-finite lateral value (it is a percentage of the same run's central-axis dose at that depth, so it is
  undefined when that dose is not finite and positive) in any run makes that endpoint not established. No run is
  excluded and nothing is imputed.
- **A lateral value of zero is a valid value.** The lateral endpoints are differences on an absolute scale, so a
  point where a run scores no dose contributes 0. A zero is informative here. In two 200 MeV field runs (one per
  Lenovo arm, 3e7 histories) the smallest non-zero voxel among the 400 that make a 50 or 70 mm value added between
  6e-7 and 4e-5 points to it, and the median non-zero voxel under 0.001 points; at 6e7 histories each is half
  that. So runs that are all zero bound both arms far inside ±0.20 points. A row where every run of both arms is exactly
  zero is marked "all runs zero in both arms": it says that neither code deposits dose there at this number of
  histories, not that a difference was measured and found small.
- **Zero standard error.** When both arms have zero variance (all-zero rows included), Welch's degrees of freedom
  are undefined. The outcome is then decided by the point difference against the margin: equivalent, with p = 0
  entered into Holm, if it is inside; not equivalent, with p = 1, if it is not.
- **What the lateral rows report.** Each arm's mean, in points, beside the difference and its intervals, so that
  a reader can see when a row is equivalent because both doses are very small.
- **Descriptive ratio at 50 and 70 mm.** For those rows the table also gives the ratio of the arm means (arm /
  reference) with a 95% interval, when both means are positive: the interval is on the log of each mean (delta
  method, standard error = SE of the mean / mean), with Welch–Satterthwaite degrees of freedom. It is descriptive
  only and enters no claim.
- **Named in advance.** In part B at 100 MeV, 40 mm, the 40–80 mm annulus, B-pgcc was 1.035 times B-up (95%
  [1.004, 1.067]). The ±0.20-point margin cannot show a difference of that size: 3.5% of the dose is 0.20 points
  only where the dose is about 5.7 points, and at 50 mm out it is far lower (about 0.3 points in those two 200 MeV
  runs). The equivalence rows at 50 and 70 mm answer whether any difference there is below 0.20% of the in-field
  dose. The descriptive ratio of the B-pgcc rows at 100 MeV, 50 and 70 mm, is the number that speaks to the 1.035.

## Sizing

There are no planning runs. The run-to-run spread of the 200 MeV field (8 runs per arm at 3e7 histories, parts A
and B) gives, for the worst of the three contrasts:

| endpoint (200 MeV) | SE at 8 + 8 runs, 3e7 | margin | margin / SE |
|---|---|---|---|
| lateral, 5 mm | 0.027 to 0.034 points | 0.50 | 15 to 19 |
| lateral, 10 to 30 mm | 0.011 to 0.025 points | 0.20 | 8 to 18 |
| central-axis dose (log) | 0.0011 to 0.0013 | 0.0050 | 4.0 to 4.6 |
| R80, R20 | 0.006 to 0.007 mm | 0.5 | 67 to 77 |

For 80% power at a true difference of zero, margin / SE has to be about 3.1 (t(0.95) + t(0.90) at 14 degrees of
freedom); for the smallest Holm threshold in a family of 34 it is about 4.9. The central-axis doses are the tight
ones. **Each run therefore uses 6e7 histories**,
twice the 200 MeV case, which takes them to about 5.7 to 6.5 if the spread at 100 and 150 MeV is like that at
200 MeV. That is an assumption: the same number of protons crosses the central patch at every energy, but nothing
has been measured at these energies. A row whose interval comes out too wide is reported inconclusive. No runs are
added after the results are read.

**The assumption is weakest for the maximum, above all at 100 MeV.** The 200 MeV figure was measured on a fixed
row. The maximum is a different estimator, and at 100 MeV the peak is narrow against the 2 mm rows, so its
run-to-run spread is not covered by that figure. It is also sensitive to range: a shift well inside the ±0.5 mm
range margin changes how the peak divides between rows and can move the maximum by an amount comparable to the
0.5% margin. So "maximum not equivalent, R80 equivalent" is a possible outcome that would be a range effect, not
a dose effect. Fixed now, for a maximum row that is inconclusive or not equivalent:

- it counts against the joint claim as it stands; the endpoint and its margin are not changed;
- the report gives, beside it, the R80 difference of the same contrast and, as a descriptive quantity only, the
  ratio (with its 95% interval) of the mean of the three central-axis rows centred on each run's maximum row,
  which depends less on how the peak divides between rows;
- if that three-row ratio lies within [0.995, 1.005] and the R80 difference is not zero, the row is described as
  consistent with a range effect on the 2 mm grid; otherwise as a difference in peak dose. Either way it is
  reported as not established or not equivalent, and the description is labelled post-estimation reading, not a
  test.

**8 runs per arm and energy, 6e7 histories each: 80 runs.**

## Seeds

In the existing blocks, new numbers, checked by `validation/apples_seeds_check.py` against every earlier seed of
the programme before anything runs:

| arm | 100 MeV | 150 MeV | smoke |
|---|---|---|---|
| A-up | 960041–960048 | 960051–960058 | 960941, 960951 |
| A-port | 961041–961048 | 961051–961058 | 961941, 961951 |
| B-up | 962041–962048 | 962051–962058 | 962941, 962951 |
| B-pgcc | 963041–963048 | 963051–963058 | 963941, 963951 |
| B-picc | 964041–964048 | 964051–964058 | 964941, 964951 |

## How it runs

1. **One workflow per host**, as parts A and B. The workflow file is the run request. One job per arm and energy,
   8 runs each: 4 jobs on the Lenovo, 6 on the HP.
2. **Smoke first** (1e5 histories, the smoke seeds), to show that each job completes on its host. Only completion
   and the transport checks are looked at.
3. **Each run** gets a new directory, writes `run.json` (seed, requested and simulated histories, threads, energy,
   binary sha256, commit, host, times, transport status) and the sha256 of its two dose files, and refuses on any
   count, config or hash failure. **No endpoint is computed on the run hosts.**
4. **Disk.** A run directory is about 41 MB and the whole part under 3 GB per host. Each job refuses to start with
   less than 10 GB free. When a host's jobs have all finished, its run tree is packed as one gzip tar, hashed, and
   copied to `\\192.168.1.212\SMBWritable\field_100_150\`; the copy is hashed on the share and compared. Only
   then is the tree removed from the PC, with the local tar (sjswerdloff, 2026-10-02: "move the run data off the
   pcs, on to the SMBWritable").
5. **Nobody reads an endpoint until every run of both hosts is complete and on the share.**

## Analysis inputs and what is verified

`validation/field_followup_analyse.py` reads the two extracted trees. For every run of the frozen lists it checks,
and refuses or marks the contrast PARTIAL otherwise:

- the run lists themselves, read from the workflow files at the acquisition commit;
- `run.json` against the frozen seed, energy, threads, histories, commit and the binary hash above, and against
  the config and the log's simulated count;
- the plan and CT against what `field_edge_make_cases.py` writes at the acquisition commit, and the materials,
  scanner and beam-model files against the committed ones;
- each dose file against the sha256 written beside it on the run host.

**It then computes every endpoint itself, from those dose files.** No endpoint value written on a run host is
used, which closes for this part the limit that #52 records for parts A and B.

## Compute

From the part A and B run times, scaled by the pencil-beam time ratios (100 and 150 MeV took 0.37 and 0.66 of the
200 MeV time) and by two for the histories: about 8 h on the Lenovo and 11 h on the HP, running at the same time.
An estimate; each job's first run shows the real rate.

## Chronology and amendments

Recorded here as they happen, not tidied: the acquisition commit, the smoke and full run identifiers, any change
made after the freeze and why, and the commit of the analysis that produced the reported numbers.
