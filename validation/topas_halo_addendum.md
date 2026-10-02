# TOPAS halo comparison: addendum to the TOPAS design

**Status: FROZEN at `5aedf3ac49d7a79567a7f413e95a244359600dd3` (2026-10-02); amended twice before the analysis, see
"Amendments" at the end.** The 24 TOPAS runs listed here are authorised, from a `git archive` of the frozen commit;
sjswerdloff confirmed the decisions under "Decisions" on 2026-10-02. The run list, seeds, endpoints and analysis
script do not change from the frozen commit; a change needed before the analysis is run is made as a recorded
amendment, in a new commit, before any endpoint is read. Written at
sjswerdloff's request (2026-10-02). It extends `validation/topas_design.md` (TD), which stays the governing design;
where this file is silent, TD applies.

## Why this is separate from the report's claim

sjswerdloff, 2026-10-02: the report has to demonstrate clinical equivalence between Portable MCsquare and upstream
OpenMCsquare. It does not have to address the difference in physics between MCsquare and TOPAS. This comparison is
therefore optional. It characterises a difference; it supports no equivalence claim and no clinical claim.

## What is already known (not confirmatory)

- At 200 MeV (4 runs × 1e7 per code, planning data), Portable's normalised energy fractions in the annuli far from
  the axis are about 10–20% lower than TOPAS opt0's (report §5.4). These are fractions of each code's own slab
  energy, not absolute deposited energy. The report tabulates these ratios as point estimates; its 4 runs per code
  would support intervals, as it gives for R80 and σ, but none are given for the annuli there.
- On the same host, Portable and upstream with the fix are equivalent in every annulus at 200 MeV (report §5.0,
  Table 3). So, by inference, the lower fractions are common to both MCsquare code lines. Their cause (physics
  models, scoring, estimators, geometry or source conditions) is not established.
- At 100 and 150 MeV there is no comparison with TOPAS in the halo (the annulus fractions). For range there is a
  planning value at 100 MeV (R80 about 0.43 mm shorter than TOPAS, 2 × 1e6 per code, report §5.4) and none at
  150 MeV.

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

13 endpoints × 3 energies × 3 Portable arms = 117 rows. A row carries pointwise 90% and 95% Welch intervals, except
a row under either of the two rules below (an annulus with a zero run; an invalid value): such a row has no ratio or
difference and no interval. A row whose values have zero sample variance in both codes keeps its ratio or
difference and carries no standard error and no interval (amendment 2). TD's
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
annuli it is smaller. **If** TOPAS's standard deviation of the log fraction equals MCsquare's, 8 + 8 runs give a 95%
interval on the ratio whose upper relative half-width is about 3.3% at a standard deviation of 0.03, 6.6% at 0.06
and 9.0% at 0.08 (t(0.975, 14) × SD / 2 on the log scale; the lower side is slightly narrower). These are the
widths to expect under that assumption. They are not a threshold: nothing here says which differences the
comparison will or will not detect. TOPAS's run-to-run spread at 100 and 150 MeV is not known; the intervals are
reported as they come out, and no runs are added after the results are read.

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
- The full analysis first requires exactly one run directory with `verdict: COMPLETE` for each of the 24 seeds, and
  no run directory outside the list. Each run must conform to the design and not only agree with its own record
  (amendment 1): `run.txt` byte for byte what `make_run.sh` writes for the directory's name; `stage1_base.txt` and
  the runner identical to the frozen commit's; a well-formed sha256 of the TOPAS executable, the same in all runs;
  no provenance key written twice; `dose.binheader` of the size and hash the provenance records and stating exactly the scorer,
  filter, component, grid, voxel widths, quantity and report that the base requests; `dose.bin` of the size that
  grid implies. An entry that looks like a run but is not named as `make_run.sh` names one is refused, not ignored.
  Otherwise the analysis refuses, and no dose file is opened.
- It then checks each `dose.bin` and its header against the sha256 in the provenance, checks the header again,
  computes the endpoints, verifies the MCsquare tree as `apples_analyse.py` does (frozen population, collection
  manifest, the fingerprint above), and writes one JSON document and one markdown table.
- **Nobody reads an endpoint until all 24 runs are complete.** `pencil_endpoints.py` is not run on a single TOPAS
  output of this set; the analysis above is the only reader.

## What it will not do

- Explain the 0.5 mm range offset or the halo difference.
- Change MCsquare's physics.
- Enter the report's equivalence results (§5.0–§5.2). It would replace the PRELIMINARY table in §5.4.

## Decisions

sjswerdloff, 2026-10-02: "both, addendum first. reuse existing Portable runs." (relayed by clement-7074f29f, #54
comment 28305), and then to connor-227743e6 directly: "I'm confirming the TOPAS halo comparison rulings."

1. Both TOPAS programmes run, this addendum first, then TD's confirmatory arms.
2. The MCsquare arm is the existing same-host Portable runs; no new MCsquare runs are made.

## Amendments

### Amendment 1, 2026-10-02 16:50 NZDT, before any endpoint was read (alden-ec2221c7, #54 review 7101)

**State of the data when it was made.** 10 of the 24 TOPAS runs had finished (the eight at 100 MeV and the first two
at 150 MeV); the rest were running or not started. No `dose.bin` of the set had been opened by connor-227743e6, and
clement-7074f29f reports opening none. The acquisition is untouched: `make_run.sh`, `run_topas.sh`,
`stage1_base.txt`, the run list and the seeds are those of the frozen commit.

**What the frozen analysis accepted and should not have** (each reproduced by alden-ec2221c7 on synthetic runs made
by the real wrappers):

1. a set in which no run records a sha256 of the TOPAS executable (24 absent values counted as one executable);
2. a `run.txt` that requests the EM option alone and not the six frozen physics modules, if its provenance recorded
   that file's hash: the check was of identity with the run's own record, not of conformity to the design;
3. a `dose.binheader` changed after the run (the X and Y bin counts exchanged): the header was never hashed, and
   it is what tells the reader the grid.

**What changes, in `validation/topas_halo_compare.py` only.** The gate now requires what the "Analysis" section
above lists. In addition a directory entry that looks like a run (a name starting `E<digit>`, or a directory
holding run files) but is not named as `make_run.sh` names one is a refusal; other entries are still listed as
ignored. The header and dose hashes and the header's content are checked again immediately before a dose file is
interpreted.

**What does not change.** The estimands, the endpoints, the intervals, the zero and invalid-value rules, the MCsquare
dataset and its fingerprint. Every new check reads provenance and configuration only; none depends on a dose value,
so none can be tuned to an outcome.

**Checked on the real runs.** The amended gate was run in `--status` mode on the run directory at the time above.
That mode reads `provenance.txt`, `run.txt`, `stage1_base.txt` and `dose.binheader` and takes the size of
`dose.bin`; it reads no dose value. The ten finished runs pass every new check; the refusal lists only the
fourteen runs not yet complete.

**Wording corrected in this file** (same review): rows under the zero-run or invalid-value rule carry no interval;
the precision paragraph gives interval widths under a stated assumption and no longer says what the comparison
"resolves"; "no comparison with TOPAS at 100 and 150 MeV" is limited to the halo, since a 100 MeV planning range
value exists; the 200 MeV planning ratios are described as tabulated without intervals, not as incapable of them.

### Amendment 2, 2026-10-02 19:10 NZDT, before any endpoint was read (alden-ec2221c7, #54 review 7114)

**State of the data when it was made.** 18 of the 24 TOPAS runs had finished (all at 100 and 150 MeV and the first
two at 200 MeV). No `dose.bin` of the set had been opened; the analysis had been run in `--status` mode only. The
acquisition is untouched.

**What amendment 1's analysis did and should not have.** For a row whose values are identical within each code
(zero sample variance in both), it gave a standard error of zero and 90% and 95% intervals of zero width: for
example a difference of -0.5 with the interval [-0.5, -0.5], or a ratio of 2 with [2, 2]. Zero sample variance in
eight Monte Carlo runs is not zero variance, so that interval states a precision that was not estimated. The
same-host field comparison has the same rule (#55).

**What changes, in `validation/topas_halo_compare.py` only.** Such a row keeps its ratio or difference and its run
counts, is marked `not estimated: zero sample variance in both codes`, and carries no standard error, no degrees of
freedom and no interval, in the JSON and in the table. Every other row is unchanged and is marked `welch`. A row
with zero sample variance in one code only is an ordinary Welch row, as before. In addition the byte count that the
provenance records for `dose.binheader` is now compared with the file's size (the hash already bound its content).

**What does not change.** The estimands, the endpoints, the intervals of every row with a non-zero standard error,
the zero-run and invalid-value rules, the gate, the MCsquare dataset and its fingerprint. The rule depends on
whether the standard error is exactly zero and on nothing else, and was fixed before any value was seen.

**Checked on the real runs.** The amended gate was run in `--status` mode at the time above: the eighteen finished
runs pass, including the new byte-count check, and the refusal lists only the six runs not yet complete.
