# Step-size follow-up 2: prediction, written before any run (2026-10-04)

Question: is the step-size effect converged at 0.1 mm at 100 and 150 MeV? (At 200 MeV it is: 0.05 mm was 0.044 mm shorter than 0.1 mm.)

Runs: 1e6 histories, opt0, two seeds per arm, 12 threads, nice 19. Same pinned archive `6c4aaff81994`.
- 100 MeV, MaxStepSize 0.05 mm
- 150 MeV, MaxStepSize 0.05 mm

The controls and 0.1 mm arms are reused from follow-up 1.

Prediction, as R80 change from 0.1 mm to 0.05 mm at the same energy:
- If 100 and 150 MeV converge as 200 MeV does (about 10% of the 0.1 mm shift): 100 MeV about −0.03 mm, 150 MeV about −0.04 mm. Range for "converged": −0.01 to −0.06 mm.
- If the change is larger than −0.10 mm at either energy, the effect is not converged at 0.1 mm there. Then the 0.1 mm gaps in §5.4 are upper bounds, and a 0.025 mm arm is necessary.

My expectation: converged at both energies. Confidence is low at 100 MeV, because the H1 model did not predict the size of the shift there.
