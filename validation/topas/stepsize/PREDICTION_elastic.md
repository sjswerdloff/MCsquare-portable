# Nuclear elastic off: prediction, written before any run (2026-10-04)

Question: how much of TOPAS's far-halo dose comes from hadronic elastic scattering (`g4h-elastic_HP`)? At 200 MeV the halo addendum (report `halo_868d1910`) gives MCsquare / TOPAS ring fractions below 1 in the far rings:
- ring_200_40_80: 0.887
- ring_100_80_200: 0.812
- ring_100_40_80: 0.969

So TOPAS has more far-halo dose than MCsquare.

Arm: TOPAS opt0, 200 MeV, 1e7 histories, 12 threads, seeds 914201 and 914202. Physics modules are the frozen six with `g4h-elastic_HP` removed. Everything else is the pinned archive `6c4aaff81994`, whose `validation/topas` is identical to the halo-addendum freeze `5aedf3ac49d7`.
Control: the 8 halo-addendum 200 MeV runs (unchanged physics). Endpoints are from `pencil_endpoints.py`, the same as in the addendum.

Prediction, TOPAS elastic-off / TOPAS control:
- R80: change smaller than 0.02 mm. Elastic scattering removes little energy from the primaries.
- Inner rings (5-10, 10-20 mm): within 2%.
- ring_200_40_80: down 3 to 10%.
- ring_100_80_200: down 5 to 15%.

Interpretation fixed now:
- If the far rings fall by about the MCsquare gap (about 11% at 200_40_80, about 19% at 100_80_200) or more, then the difference in elastic modelling can explain the gap.
- If they fall by less than 3%, then elastic is not the source, and the gap is in another process (non-elastic secondaries, or multiple scattering tails).

Confidence is low. I have no measured value for the elastic share of the far halo in water at 200 MeV.
