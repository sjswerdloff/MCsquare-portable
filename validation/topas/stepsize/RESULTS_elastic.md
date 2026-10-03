# Nuclear elastic off: results (2026-10-04 02:22 NZDT)

TOPAS opt0, 200 MeV, 1e7 histories × 2 (seeds 914201, 914202), physics modules without `g4h-elastic_HP` (the edit and the absence of `hElasticWEL_CHIPS_HP` / `hadElastic` were checked in a smoke run). Control: the 8 halo-addendum 200 MeV runs. Endpoints from `validation/pencil_endpoints.py`; comparison by `compare_elastic.py` (geometric-mean ratio for rings, mean difference for R80 and σ, 95% Welch intervals). Positive control: the control R80 mean is 260.866 mm, as in the addendum. The two elastic-off seeds agree within 0.0003 mm in R80.

| endpoint | elastic off / control | 95% |
|---|---|---|
| R80 | **+0.31 mm** (261.174 vs 260.866) | [+0.306, +0.309] |
| σ at 100 mm | −0.069 mm | [−0.073, −0.065] |
| σ at 200 mm | −0.181 mm | [−0.182, −0.179] |
| ring_200_0_5 | 1.117 | [1.116, 1.117] |
| ring_200_10_20 | 0.802 | [0.800, 0.804] |
| ring_200_20_40 | 0.391 | [0.381, 0.402] |
| ring_200_40_80 | 0.524 | [0.450, 0.611] |
| ring_100_10_20 | 0.678 | [0.676, 0.680] |
| ring_100_40_80 | 0.728 | [0.721, 0.734] |
| ring_100_80_200 | 0.990 | [0.978, 1.003] |
| raw IDD maximum | about 1.10 (3.76 vs 3.41) | |

Full table: `elastic_compare.txt`.

**Against PREDICTION_elastic.md (committed 10c261b before the runs): wrong, by a large margin.**
- R80: predicted a change below 0.02 mm. Measured +0.31 mm.
- Far rings: predicted 3 to 15% lower. Measured 25 to 61% lower (ring_200_40_80 0.52, ring_200_20_40 0.39). The 80–200 mm ring at 100 mm is unchanged.
- I under-estimated hadronic elastic scattering in water by an order of magnitude. Probably the cause is proton–hydrogen elastic scattering, which moves primary protons out of the core. I did not test this separately.

**What the result shows:**
- In TOPAS, hadronic elastic scattering gives about half of the dose in the 20–80 mm rings at 200 MeV, and it moves R80 by 0.31 mm and σ at 200 mm by 0.18 mm.
- MCsquare / TOPAS for ring_200_40_80 is 0.887 (addendum). Removing elastic gives 0.52. So the MCsquare halo deficit is much smaller than the whole elastic contribution, and a difference in how the two codes model elastic scattering is a candidate for it. This test does not show that it is the cause.
- R80 moves the wrong way to explain the R80 gap: removing elastic makes TOPAS longer, and TOPAS is already longer than MCsquare. The R80 gap is the step-size effect (RESULTS_followup*.md).

Bounds: 200 MeV only; 2 runs; elastic removed completely (not adjusted to MCsquare's model). Diagnostic, not pre-specified.

## Clarification (2026-10-04, after alden-ec2221c7's review)

The text above is kept as written. Where it differs from this section, this section applies.

- **The ring values are normalized fractions, not absolute dose.** Each ring endpoint is E_ring / E_slab within one configuration. An off/control ratio of 0.524 means the normalized 40–80 mm fraction at 200 mm fell by about 48%. The absolute ratio is 0.524 × (E_slab_off / E_slab_control), and the slab-total ratio is not in these records. The raw IDD maximum is a different statistic and does not supply it.
- **So the measured result is:** removing `g4h-elastic_HP` reduces the normalized 20–40 and 40–80 mm fractions at 200 mm to 0.391 [0.381, 0.402] and 0.524 [0.450, 0.611] of control, a 48–61% reduction. This is a removal sensitivity. Removing a process also changes trajectories, fluence and other interactions, so it is not a decomposition of dose by process. "Elastic gives about half of the dose in the 20–80 mm rings" (line 27) is not shown.
- **The MCsquare comparison** (line 28): the MCsquare/TOPAS fraction deficit (0.887) beside this sensitivity motivates elastic modelling as a candidate for the halo gap. It does not measure a shared absolute elastic share, and does not show the cause.
- **"The R80 gap is the step-size effect"** (line 29) should read: most of the observed R80 gap decreases under step refinement; a residual of +0.04 to +0.13 mm remains at 0.05 mm, and finer settings were not tested (RESULTS_followup2.md).
- **Line 24** (proton–hydrogen elastic) is a hypothesis, not tested.

An absolute claim would need the slab and annular totals, properly normalized, with their uncertainty.
