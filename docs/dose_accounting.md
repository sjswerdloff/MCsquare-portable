# What MCsquare scores as dose, and what it drops

Portable MCsquare's physics is **unchanged** from upstream here. This page records where each part of a proton's
energy ends up, so a comparison with a full-physics code (TOPAS/Geant4, issue #32) can say which differences are
modelling choices rather than bugs. Read from the source at `3807be0`, with default `config.txt`.

## Where the energy goes

| Energy carrier | What MCsquare does | Scored as dose? | Where in the source |
|---|---|---|---|
| Primary proton, continuous loss (soft collisions, δ-electrons below `Te_Min`) | Condensed-history step with straggling | **Yes**, along the step | `compute_hadron.c`, `Compute_dE2` → `Energy_Scoring_*` |
| δ-electrons above `Te_Min` (0.05 MeV) | Sampled as a discrete ionisation; the electron is **not transported** | **Yes, locally** at the interaction point | `compute_hadron.c`, `Compute_Ionization_Energy` → `v_dE_hard`; `config.txt`: "currently locally absorbed" |
| Proton below `E_Cut_Pro` (0.5 MeV) | Transport stops | **Yes, locally**: all remaining energy | `compute_hadron.c`, the `Ecut_Pro` branches |
| Nuclear elastic (ICRU): nuclear recoil | Recoil not transported | **Yes, locally** | `compute_nuclear_interaction.c`, `Compute_Elastic_ICRU` |
| Elastic p–p: recoil proton | Transported as a secondary proton (local if below `E_Cut_Pro`) | **Yes** | `Compute_Elastic_PP` |
| Nuclear inelastic: secondary protons, deuterons and alphas | Transported (local if below `E_Cut_Pro`, or if `Simulate_Secondary_*` is False) | **Yes** | `Compute_Nuclear_Inelastic_proton/deuteron/alpha` |
| Nuclear inelastic: heavy recoils / fragments | Energy fraction from the nuclear table | **Yes, locally** | `Compute_Nuclear_Inelastic_recoils` |
| Nuclear inelastic: **neutrons** | Not produced | **No: dropped** | inelastic branch of `Compute_Nuclear_interaction` returns only the rows above |
| Nuclear inelastic: **gammas** | Not transported. With `Score_PromptGammas` True, prompt gammas go to their own tally (`PG_Scoring`) and return no dose | **No: dropped** | same; `Compute_PromptGamma` |
| Whatever else of the primary's energy the four rows above do not account for | Not tracked (I have not traced how the nuclear tables partition it) | **No** | same |

In short, the inelastic branch returns only the energy carried by p, d, α and recoils, and the primary is then
killed (`v_type = Unknown`). **Energy carried away by neutrons and gammas never reaches the dose grid.** Energy
below the electron and proton cutoffs is **deposited where it was produced**, not transported.

## How big is what is dropped

**Measured once, and only this:** in TOPAS (OpenTOPAS 4.3.0 / Geant4 11.3.2, EM option 0 plus the default
hadronic modules), for a 150 MeV Gaussian pencil beam (σ 3 mm) in a 400 × 400 × 350 mm water box, dose from
neutrons, gammas and **all their descendants** is **0.8% of the integral dose** in the box (filtered/unfiltered
= 0.9918, 0.9921, 0.9919 over three runs of 1e4 histories). clement-7074f29f, 2026-10-01, run files in the #32
TOPAS provenance.

What that does **not** tell us:
- **Where** the dropped dose lies. An integral fraction says nothing about whether it is spread thinly out of
  field or concentrated somewhere.
- Its size at other energies, or in anything other than a large water box. Out-of-box neutron dose is not in the
  box at all.
- Whether it is **clinically** significant: for example, in-field near the distal edge, low-dose bath at organs
  at risk, or peripheral dose.

## When to revisit (sjswerdloff, 2026-10-01)

No physics change to portable MCsquare for now. Stage 1 of #32 already scores TOPAS dose both with and without
neutron/gamma descendants (`Dose` and `DoseAll`) at 100, 150 and 200 MeV. So `DoseAll − Dose` gives the dropped
dose **as a 3-D map** at no extra cost. If that map shows differences that could be clinically significant, that
is the point to consider modelling neutron/gamma dose, and it is sjswerdloff's decision.

Not established here, and not claimed: any judgement of clinical significance. This page records only what is
scored and what is dropped.
