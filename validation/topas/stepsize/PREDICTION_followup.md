# Step-size follow-up: prediction, written before any run (2026-10-03)

Context: at 200 MeV the shifts were −0.17 mm (0.5 mm steps) and −0.43 mm (0.1 mm steps). The H1 model predicted −0.08 and −0.28, so the model under-estimates by a factor of approximately 1.5.

Runs: 1e6 histories, opt0, two seeds per arm.
- 200 MeV, MaxStepSize 0.05 mm (convergence check)
- 100 MeV and 150 MeV, control and 0.1 mm

Predictions, as R80 change from the control at the same energy:
- 200 MeV, 0.05 mm vs control: approximately −0.45 to −0.50 mm. That is at most 0.07 mm more than at 0.1 mm, if the effect is near convergence.
- 150 MeV, 0.1 mm: approximately −0.30 mm (model −0.20 × 1.5).
- 100 MeV, 0.1 mm: approximately −0.10 mm (model −0.06 × 1.5). The H1 linear branch acts only above approximately 83 MeV, so the effect at 100 MeV must be small.

If these predictions are correct, a gap to MCsquare remains at 100 MeV (approximately +0.3 mm) and at 150 MeV (approximately +0.2 mm). That residual then needs a second explanation.
