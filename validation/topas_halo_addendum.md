# TOPAS halo comparison: addendum to the TOPAS design

**Status: DRAFT, complete for freezing. Nothing here is authorised to run.** It becomes FROZEN, at the commit that
changes this line, when sjswerdloff confirms the two decisions under "Decisions" in his own words; they reached
this file relayed (#54, comment 28305). Written at sjswerdloff's request (2026-10-02). It extends
`validation/topas_design.md` (TD), which stays the governing design; where this file is silent, TD applies.

## Why this is separate from the report's claim

sjswerdloff, 2026-10-02: the report has to demonstrate clinical equivalence between Portable MCsquare and upstream
OpenMCsquare. It does not have to address the difference in physics between MCsquare and TOPAS. This comparison is
therefore optional. It characterises a difference; it supports no equivalence claim and no clinical claim.

## What is already known (not confirmatory)

- At 200 MeV (4 runs × 1e7 per code, planning data), Portable's normalised energy fractions in the annuli far from
  the axis are about 10–20% lower than TOPAS opt0's (report §5.4). These are fractions of each code's own slab
  energy, not absolute deposited energy, and they are point estimates without intervals.
- On the same host, Portable and upstream with the fix are equivalent in every annulus at 200 MeV (report §5.0,
  Table 3). So, by inference, the lower fractions are common to both MCsquare code lines. Their cause (physics
  models, scoring, estimators, geometry or source conditions) is not established.
- At 100 and 150 MeV there is no comparison with TOPAS.

## Question

How large is the difference between MCsquare and TOPAS opt0 in the energy fraction of each annulus, at 100, 150 and
200 MeV, and how well is it determined? This is estimation only: no margin, no pass or fail, and no attribution of
cause (TD, "Primary reporting objective").

## Arms

| arm | build | host | runs |
|---|---|---|---|
| TOPAS opt0 | as TD (OpenTOPAS 4.3.0 / Geant4 11.3.2, EM option 0, TD's modules), from a `git archive` of the frozen commit | Mac Studio | 8 per energy × 1e7, new, seeds below |
| Portable MCsquare | the **existing** same-host runs at acquisition commit `2f9dab404cea02f5072352f7ae5773dca27eb2a1` | Lenovo (A-port), HP (B-pgcc, B-picc) | 8 per energy × 1e7 per arm, already made |

**TOPAS seeds** (a block of their own; TD reserves 910001+ and 911001+ for its confirmatory arms):

| energy | seeds |
|---|---|
| 100 MeV | 912001–912008 |
| 150 MeV | 912011–912018 |
| 200 MeV | 912021–912028 |

**The MCsquare arm is reused, and that limits what this comparison is.**

- Its endpoints **were read before this design was fixed** (report §5.0–§5.2). Every output of this comparison
  carries that statement.
- It ran on other hosts than TOPAS, with other compilers. Nothing yet relates those builds to a Portable build on
  the Studio: that is parts C and D of the same-host plan (TD's confirmatory arms and the bridge), which have not
  run and which this addendum does not replace.
- **All three Portable arms are compared with TOPAS, each separately.** None is selected and none are pooled: all
  three have been read, so a choice among them would be made with knowledge of the values, and they differ in host
  and compiler. They share one TOPAS arm, so the three sets of rows are not independent.
- The dataset is fixed by its fingerprint: sha256 `0d5be975b294e0b77860e0ddbce5807caa526a0ca56481df5056effff40033ba`
  over the 321 files `apples_analyse.py` reads from the collected tree. The analysis refuses any other.
- Named in advance, as the same-host plan requires: at 100 MeV, 40 mm, the 40–80 mm annulus, B-pgcc's fraction is
  1.035 times B-up's (95% [1.004, 1.067]). B-pgcc's ratio to TOPAS at that endpoint will differ from B-picc's by
  about that much for that reason.

Upstream with the fix is not compared with TOPAS. TOPAS opt4 is not included.

## Endpoints (all descriptive)

Per energy, at the same slab depths as the same-host comparison (100 MeV: 40 and 60 mm; 150 MeV: 80 and 125 mm;
200 MeV: 100 and 200 mm), by `validation/pencil_endpoints.py` with those depths, TD's slab rule and TD's `Dose`
scorer:

- the energy fraction in each annulus (5–10, 10–20, 20–40, 40–80 and 80–200 mm): geometric-mean ratio
  MCsquare / TOPAS, from the difference of mean logs of the per-run values;
- σ at each depth and R80: difference of run-level means, MCsquare − TOPAS.

Each with pointwise 90% and 95% Welch intervals: 13 endpoints × 3 energies × 3 Portable arms = 117 rows. TD's
reference bands are shown beside the 200 MeV rows only; they were set for 200 MeV and are not extended to the other
energies.

**An annulus with no scored energy.** In the same-host comparison the 80–200 mm annulus at 100 MeV was exactly zero
in all 40 MCsquare runs at both depths (80 values) at 1e7 histories. A log ratio is undefined there. For any
annulus in which any run of either code is zero, the row gives, per code, the number of runs with a non-zero value
and the arithmetic mean of the per-run fractions, and no ratio. No pseudocount is added. (The first draft said a
pooled fraction, energy summed over runs divided by slab energy summed over runs. The reused MCsquare records hold
fractions, not energies, so the mean of the per-run fractions is used for both codes.)

**An invalid value** (a failed σ fit, an absent or multiple R80 crossing) in any run makes that row "not computed",
with the count of invalid runs per code. The row stays in the table.

## Precision to expect (not an acceptance criterion)

From the same-host runs, the per-run standard deviation of the log fraction in MCsquare is about 0.03–0.06 for the
40–80 mm annulus at 100 MeV and about 0.06–0.08 for the 80–200 mm annulus at 150 MeV; at 200 MeV and in the inner
annuli it is smaller. If TOPAS varies similarly, 8 runs per code give 95% half-widths of roughly ±3–6% and ±6–9% on
those two ratios. That resolves a difference of the size seen at 200 MeV (10–20%) and does not resolve one under
about 5%. TOPAS's run-to-run spread at 100 and 150 MeV is not known; the intervals are reported as they come out,
and no runs are added after the results are read.

## Execution (clement-7074f29f)

- Runs are made by `validation/topas/make_run.sh` and `run_topas.sh` from a pinned `git archive` of the frozen
  commit, one after another, under `/Volumes/T7 Shield/SMBWritable/topas/halo_addendum/`.
- The Studio is shared with vivian-1a61bc9a's inference: `nice -n 19`, 12 threads to start, 16 if her throughput
  holds. The thread count is in each run's directory name and `run.txt`. Wall times depend on that load and are
  recorded, not compared. clement-7074f29f reports that TOPAS's results do not depend on the thread count in the
  configurations he tested; this design does not rely on it, since each seed is run once.
- A failed run is repeated with the same seed as a new attempt (`a2`, …); the failed attempt's directory stays.
- His estimates (#54 review 7063), at 16 threads: about 5.7 min per run at 100 MeV and 17 min at 200 MeV, scaled
  from measured 1e6 runs; 150 MeV is interpolated, not measured. About 0.90 GB per run, 21.5 GB in all.

## Analysis

`validation/topas_halo_compare.py`, committed with this file, with its tests.

- `--topas-root <dir> --status` reports whether the 24 runs are complete. It opens no dose file.
- The full analysis first requires exactly one run directory with `verdict: COMPLETE` for each of the 24 seeds, with
  `run.txt` and `stage1_base.txt` as recorded in its provenance, the base identical to the frozen commit's, one TOPAS
  executable for all runs, and no run directory outside the list. Otherwise it refuses, and no dose file is opened.
- It then checks each `dose.bin` against the sha256 in its provenance, computes its endpoints, verifies the MCsquare
  tree as `apples_analyse.py` does (frozen population, collection manifest, the fingerprint above), and writes one
  JSON document and one markdown table.
- **Nobody reads an endpoint until all 24 runs are complete.** `pencil_endpoints.py` is not run on a single TOPAS
  output of this set; the analysis above is the only reader.

## What it will not do

- Explain the 0.5 mm range offset or the halo difference.
- Change MCsquare's physics.
- Enter the report's equivalence results (§5.0–§5.2). It would replace the PRELIMINARY table in §5.4.

## Decisions

Relayed by clement-7074f29f (#54, comment 28305) as sjswerdloff's, 2026-10-02: "both, addendum first. reuse existing
Portable runs."

1. Both TOPAS programmes run, this addendum first, then TD's confirmatory arms.
2. The MCsquare arm is the existing same-host Portable runs; no new MCsquare runs are made.
