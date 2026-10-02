# TOPAS halo comparison: addendum to the TOPAS design

**Status: DRAFT. Nothing here is authorised to run.** Written at sjswerdloff's request (2026-10-02). It extends
`validation/topas_design.md` (TD), which stays the governing design; where this file is silent, TD applies.

## Why this is separate from the report's claim

sjswerdloff, 2026-10-02: the report has to demonstrate clinical equivalence between Portable MCsquare and upstream
OpenMCsquare. It does not have to address the difference in physics between MCsquare and TOPAS. This comparison is
therefore optional. It characterises a difference; it supports no equivalence claim and no clinical claim. The
compute is available at no marginal cost while sjswerdloff is away, so the runs can be made then or later, or not
at all.

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
| TOPAS opt0 | as TD (OpenTOPAS 4.3.0 / Geant4 11.3.2, EM option 0, TD's modules) | Mac Studio | 8 per energy × 1e7, new seeds |
| Portable MCsquare | a frozen commit of main | Mac Studio | 8 per energy × 1e7, new seeds |

- Upstream with the fix is not run again. The same-host comparison already relates it to Portable.
- The same-host runs at `2f9dab40` are **not** reused as the MCsquare arm. Their endpoints have been read, and they
  ran on other hosts. If sjswerdloff prefers to reuse them and run TOPAS only, the report says so and the comparison
  is labelled as using data read before this design was fixed.
- TOPAS opt4 is not included.

## Endpoints (all descriptive)

Per energy, at the same slab depths as the same-host comparison (100 MeV: 40 and 60 mm; 150 MeV: 80 and 125 mm;
200 MeV: 100 and 200 mm), by `validation/pencil_endpoints.py --depths`, with TD's slab rule and TD's `Dose` scorer:

- the energy fraction in each annulus (5–10, 10–20, 20–40, 40–80 and 80–200 mm): geometric-mean ratio
  MCsquare / TOPAS, from the difference of mean logs of the per-run values;
- σ at each depth and R80: difference of run-level means, MCsquare − TOPAS.

Each with pointwise 90% and 95% Welch intervals. TD's reference bands are shown beside the 200 MeV rows only; they
were set for 200 MeV and are not extended to the other energies.

**An annulus with no scored energy.** In the same-host comparison the 80–200 mm annulus at 100 MeV was exactly zero
in all 40 MCsquare runs at both depths (80 values) at 1e7 histories. A log ratio is undefined there. For any annulus in which either code has a
zero run, the report gives, per code, the number of runs with a non-zero value and the pooled fraction (energy in the
annulus summed over runs, divided by the slab energy summed over runs), and no ratio. No pseudocount is added.

## Precision to expect (not an acceptance criterion)

From the same-host runs, the per-run standard deviation of the log fraction in MCsquare is about 0.03–0.06 for the
40–80 mm annulus at 100 MeV and about 0.06–0.08 for the 80–200 mm annulus at 150 MeV; at 200 MeV and in the inner
annuli it is smaller. If TOPAS varies similarly, 8 runs per code give 95% half-widths of roughly ±3–6% and ±6–9% on
those two ratios. That resolves a difference of the size seen at 200 MeV (10–20%) and does not resolve one under
about 5%. TOPAS's run-to-run spread at 100 and 150 MeV is not known; the intervals are reported as they come out,
and no runs are added after the results are read.

## Before any run

1. The Portable commit, the TOPAS inputs, the run lists and seeds, and the analysis script are committed, and this
   file is marked frozen with those commits.
2. clement-7074f29f states the TOPAS run time per 1e7 histories at each energy and the disk needed (TD: about
   0.86 GB per TOPAS run, so about 21 GB for 24 runs).
3. Nobody reads an endpoint until every run of both arms is complete.

## What it will not do

- Explain the 0.5 mm range offset or the halo difference.
- Change MCsquare's physics.
- Enter the report's equivalence results (§5.0–§5.2). It would replace the PRELIMINARY table in §5.4.

## Decisions for sjswerdloff

1. Run it while away, later, or not at all.
2. New Portable runs on the Studio (as drafted), or TOPAS only against the existing same-host runs.
