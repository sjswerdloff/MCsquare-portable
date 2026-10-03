# TOPAS step-size test: results (2026-10-03 20:48 NZDT)

Pinned code: MCsquare-portable `6c4aaff81994`. TOPAS 200 MeV, opt0, 1e6 histories per run, 12 threads,
two seeds per arm. Outputs: `/Volumes/T7 Shield/SMBWritable/topas/stepsize/6c4aaff81994/full/`.
Endpoints from `validation/pencil_endpoints.py`.

| arm | R80 seed 1 | R80 seed 2 | mean | shift vs ctrl | vs MCsquare (260.353) | wall per run |
|---|---|---|---|---|---|---|
| ctrl (no limit) | 260.8639 | 260.8691 | 260.8665 | — | +0.513 | 136 s |
| MaxStepSize 0.5 mm | 260.6980 | 260.6909 | 260.6944 | **−0.172** | +0.341 | 160 s (+18%) |
| MaxStepSize 0.1 mm | 260.4423 | 260.4295 | 260.4359 | **−0.431** | +0.083 | 273 s (+101%) |

- The control agrees with the 1e7 halo-addendum mean (260.865).
- Seed pairs agree within 0.013 mm.
- Spot σ at 100 and 200 mm and the 40–80 mm ring at 200 mm do not move with step size. So the effect is on range only, which fits an energy-loss mechanism rather than a scattering one.

**Against PREDICTION.md (committed 1b1fbe5 before any run):** the direction is confirmed. The size is LARGER than predicted: −0.17 vs ~−0.08 at 0.5 mm, and −0.43 vs ~−0.28 at 0.1 mm. H1's mechanism is supported; its size model under-estimates the effect.

**At 0.1 mm, TOPAS's R80 (260.436) is within 0.08 mm of MCsquare's (260.353)** and on the physics expectation from the Fable report (260.42 ± 0.05). The 0.5 → 0.1 mm trend has not been shown to have converged: an arm at 0.05 mm would show whether the residual shrinks further.

Bounds: 200 MeV only. Whether the same holds at 100 and 150 MeV (where H1 predicted a smaller share) is untested.

## Clarification (2026-10-04, after alden-ec2221c7's review)

"H1's mechanism is supported" (line 17) is too strong: the direction matched H1's prediction, which does not single out H1. Mechanism not identified; see the clarification in RESULTS_followup.md. Line 15 ("fits an energy-loss mechanism") is likewise an interpretation, not a test. The convergence question on line 19 is answered for 0.1 → 0.05 mm only (RESULTS_followup2.md).
