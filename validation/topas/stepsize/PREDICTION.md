# TOPAS step-size test: prediction, written before any run (2026-10-03)

Hypothesis H1 (Fable report; verified in the Geant4 11.3.2 source by a Sonnet agent): in
`G4VEnergyLossProcess::AlongStepDoIt`, when eloss < 0.01·E the step loss is `length * dEdx(E_pre)`
with no correction. Steps in our phantom are voxel-limited (1 mm), so TOPAS's range comes out too long.
The model gives +0.31 mm at 200 MeV for 1 mm steps and +0.23 mm for 0.5 mm steps.

Test: TOPAS 200 MeV, opt0, 1e6 histories, `d:Ge/Phantom/MaxStepSize` = none (control) / 0.5 mm / 0.1 mm.

Predicted R80 relative to the control:
- 0.5 mm: shorter by ~0.08 mm (model; the sign is the firm part, the size is soft)
- 0.1 mm: shorter by ~0.28 mm
- Per-run R80 noise at 1e6 is ~0.006 mm (estimated from the 1e7 intervals), so both are well above noise.

H1 is refuted for this phantom if the 0.1 mm arm moves less than ~0.05 mm or moves longer.
Control reference: TOPAS 1e7 mean R80 at 200 MeV = 260.865 mm (halo addendum); 1e6 runs should agree within ~0.02.
