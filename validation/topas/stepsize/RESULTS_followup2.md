# Step-size follow-up 2: results (2026-10-04 01:46 NZDT)

Pinned code `6c4aaff81994`. TOPAS opt0, 1e6 histories, 12 threads, two seeds per arm. Endpoints from `validation/pencil_endpoints.py` (`followup2_endpoints.jsonl`).
Outputs: `/Volumes/T7 Shield/SMBWritable/topas/stepsize/6c4aaff81994/followup2/`.

MCsquare R80 = halo-addendum TOPAS 1e7 mean + the addendum difference (report `halo_868d1910`): 77.398 mm (100 MeV), 158.204 mm (150 MeV), 260.353 mm (200 MeV).

| energy | step limit | R80 seeds | mean | change vs 0.1 mm | change vs no limit | TOPAS − MCsquare |
|---|---|---|---|---|---|---|
| 100 MeV | 0.1 mm (follow-up 1) | | 77.556 | — | −0.268 | +0.158 |
| 100 MeV | 0.05 mm | 77.5233, 77.5230 | 77.523 | **−0.033** | −0.301 | **+0.125** |
| 150 MeV | 0.1 mm (follow-up 1) | | 158.327 | — | −0.418 | +0.123 |
| 150 MeV | 0.05 mm | 158.2774, 158.2670 | 158.272 | **−0.055** | −0.473 | **+0.068** |

**Against PREDICTION_followup2.md (committed 3a21014 before the runs):** correct. Both changes are in the predicted "converged" range (−0.01 to −0.06 mm), and neither is larger than −0.10 mm.

**Result:** the step-size effect is near convergence at 0.1 mm at all three energies. At 0.05 mm the TOPAS − MCsquare R80 gap is +0.125 mm (100 MeV), +0.068 mm (150 MeV) and +0.04 mm (200 MeV).

Bounds: 1e6 histories, two seeds per arm (seed pairs agree within 0.011 mm). Diagnostic, not pre-specified.
