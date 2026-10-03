# Step-size follow-up: results (2026-10-03 21:27 NZDT)

Pinned code `6c4aaff81994`. Each run: TOPAS opt0, 1e6 histories, 12 threads. Two seeds per arm.
Outputs: `/Volumes/T7 Shield/SMBWritable/topas/stepsize/6c4aaff81994/followup/`.

| energy | arm | R80 (mean of 2) | change from control | TOPAS − MCsquare | predicted change |
|---|---|---|---|---|---|
| 100 MeV | control | 77.824 mm | — | +0.43 mm | |
| 100 MeV | 0.1 mm | 77.556 mm | **−0.27 mm** | +0.16 mm | −0.10 mm |
| 150 MeV | control | 158.745 mm | — | +0.54 mm | |
| 150 MeV | 0.1 mm | 158.327 mm | **−0.42 mm** | +0.12 mm | −0.30 mm |
| 200 MeV | 0.05 mm | 260.392 mm | −0.475 mm (−0.044 vs 0.1 mm) | +0.04 mm | −0.45 to −0.50 mm |

The controls agree with the halo-addendum 1e7 means (77.823 and 158.744 mm). The seed pairs agree within 0.011 mm.

**Against PREDICTION_followup.md (committed 946bc6e before the runs):**
- 200 MeV convergence: correct. The 0.05 mm arm is 0.044 mm shorter than 0.1 mm, so the effect is near convergence at 200 MeV.
- 150 MeV: the direction is correct; the size is 1.4 times the prediction.
- 100 MeV: **the prediction is wrong.** I predicted a small effect, because the H1 linear branch acts only above approximately 83 MeV. The measured shift is −0.27 mm, 2.7 times the prediction.

**Interpretation:**
- At all three energies, most of the TOPAS − MCsquare R80 gap depends on the TOPAS step size: 0.27 of 0.43 mm (100), 0.42 of 0.54 mm (150) and 0.47 of 0.51 mm (200, at 0.05 mm).
- The step-size dependence is real. But H1, as modelled, does not explain its size at 100 MeV. Thus a second step-dependent mechanism in Geant4 energy loss is probable. Candidates: the energy-loss fluctuation model, or the range-table interpolation. I did not test these.
- A direct test of H1 alone: lower `linLossLimit` (Geant4 command `/process/eLoss/linLossLimit`; TOPAS needs an extension to set it) with no step limit. If H1 is the only mechanism, this gives the same R80 as the 0.1 mm arm.

Not known: convergence at 100 and 150 MeV below 0.1 mm.

## Clarification (2026-10-04, after alden-ec2221c7's review)

H1 (the linear energy-loss branch, `linLossLimit`) and the candidates on line 23 are hypotheses. None was tested in isolation, and the 100 MeV discrepancy could equally come from the model's quantitative assumptions or from other step-dependent processes. The step-size dependence is measured. No mechanism is identified.
