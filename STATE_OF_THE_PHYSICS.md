# Where we are on the physics — fully distilled (2026-07-07)

One page, no jargon walls: what petch is, what got proven, what each module does, and every
reference that carried weight, with what it was used for. Deep technical detail: PHYSICS.md.
Cycle-by-cycle log with gate numbers: FRONTIER_LOOP.md.

## The one-paragraph story

petch simulates how a plasma etches chip features, from first principles: particles (ions,
electrons, neutrals) fly in, surfaces charge up, chemistry removes atoms, the surface moves. This
campaign took the weakest pillar — surface charging — from a knob-tuned closure 36% off its
benchmark to a **zero-knob engine whose particle source carries a mathematical proof**, and in the
process found that the 30-year-old benchmark itself contains two defects. Along the way we added
three new chemistry/process modules (ALE, cryo, Bosch), each gated against published measurements,
one of them reproducing real SEM-measured silicon exactly. Every claim in the repo has a gate
script and a PASS/FAIL; every refuted idea is documented, not deleted.

## Module states (grade: what's derived vs calibrated)

| module | state | grade |
|---|---|---|
| Charging source (ions) | exact instantaneous-sheath derivation; bathtub IED verified to the analytic phase fraction | **derived** |
| Charging source (electrons) | Lambert×Maxwellian×collapse-bursts, backed by an exact invariance theorem | **derived + proven** |
| Charging dynamics | converged integrator, mirror BCs, quasi-Newton log balance, physical bounds | **derived** |
| HG 1997 benchmark | all scheme-independent observables match (flux 0.21/0.22, edge 7.8/7, foot E 28.0/28, horns 61/60); the voltage labels shown convention-dependent (same state reads 11–460 V through their under-specified kernel) and their published triple proven internally inconsistent | **closed at physics level** |
| ALE (Si/Cl₂/Ar⁺) | window, self-limitation dynamics, synergy all emerge from MD-fit yields; differentiable with exact gradients; inverse design demo | **derived (ROM)** |
| Cryo etch | Langmuir physisorption, E_ads measured (0.4 eV); 1.6× anchor reproduced | physics form + 2 calibrated params |
| Bosch DRIE | sequential punch-then-iso mechanics forced by SEM data; all 4 Ayon gates + smooth regime + 9.3× kill-test | mechanics derived; 2 rates from published endpoints |
| SF₆/O₂ etch + transport | ViennaPS-parity MC/Knudsen; de Boer knee via process params | physics + published params |
| Notching | mechanism gate passes under petch's own charging table; AR4 overetch trend = named frontier | first-principles path wired |

## What was PROVEN (not fitted) this campaign

1. **The invariance theorem**: for a Maxwellian plasma injecting electrons with the physical
   cosine-flux distribution, an oscillating RF sheath's barrier-selection and refraction cancel
   exactly — arrival angles are cos θ for ANY waveform. (One-line proof; MC-verified.)
2. **HG's cos^0.6 electron distribution is an artifact** of their unphysical uniform-in-angle
   injection (closed form for their convention: p ≈ 0.72–0.80; their 0.6 = histogram fit slop).
3. **HG's published state (floor 33 V + flux 0.22 + walls 7/39) is internally inconsistent**: a
   32 V cross-trench field zeroes the floor flux in any electrostatic field. Their voltage labels
   come from an under-specified charge→potential readout (same charge state reads 11–460 V through
   defensible variants of it); ours are the ion's-eye barrier by construction.
4. **SEE cannot fix a positive floor** (few-eV secondaries are recaptured by the +33 V well);
   **field reversal is absent at 400 kHz** (collapse-phase thermal flux covers the ion flux 2.6×);
   **the sheath is quasi-static to 10⁻³** (75 ps electron transit vs 2.5 µs period).
5. **Two of our own errors caught by mathematics**: a 1000× transit-time slip, and a v⊥ projection
   bug in the derivation (both found because derived numbers were checked against analytics).

## The remaining frontier (sequenced)

1. **ARDE × charging at deep AR** — charging halves the ballistic floor flux by AR 25 (measured);
   couple into the de Boer deep-trench run, gate the full collapse curve. (In progress.)
2. **Cross-tool parity** — same high-AR stack through petch + ViennaPS, profile overlay, charging
   prediction as the differentiator.
3. **Full-Poisson charging in physical units** — CHARGING_POISSON_PLAN.md (implementation-ready);
   gives absolute voltage labels and ends the convention question with units-honest electrostatics.
4. **Multi-line array geometry** — last HG fidelity item (the neighbor label).
5. **The paper** — theorem + artifact + inconsistency proof + derived source + gated validations.
6. **Then speed** (device-resident GPU) and the differentiable/inverse-design API layer.

## References, distilled (what each was actually used for)

**Charging benchmark + methods**
- Hwang & Giapis, J. Appl. Phys. 82, 566 (1997) — the benchmark (floor/edge/neighbor/flux vs AR);
  also the source of the stack dims (0.5 µm lines, 0.54 µm PR) that exposed our geometry inversion.
  Open: https://authors.library.caltech.edu/records/je8bd-j6v68
- Hwang & Giapis, J. Vac. Sci. Technol. B 15, 70 (1997) — their methods (Laplace-in-gas, charge
  update, mirror images, "isotropic flux" electron injection, 7000-step potential steady state).
  Open: https://authors.library.caltech.edu/records/ac5xn-zqb88
- Huang, Kushner et al., J. Vac. Sci. Technol. A 37, 031304 (2019) — MCFPM: variable-ε Poisson,
  capacitance-matched substrate, one-electron-per-ion; the template for the full-Poisson plan.
  Open: https://cpseg.eecs.umich.edu/pub/articles/JVSTA_37_031304_2019.pdf
- Ootera et al., Jpn. J. Appl. Phys. 33, 4276 (1994) — HG's cited Vlasov analysis (paywalled; the
  elementary Vlasov solution reconstructed and proven independently).
- Kawamura, Vahedi, Lieberman, Birdsall, PSST 8, R45 (1999) — RF-sheath IED theory (bathtub,
  ωτ_i scaling) grounding the derived ion source.
- Köhler et al., J. Appl. Phys. 57, 59 (1985) — the I₀(eV_rf/kTe) cycle-averaged electron flux
  (used in the flux-balance that rejected field reversal).
- Krüger, Wilczek, Mussenbrock, Schulze, PSST 28, 075017 (2019) + Hartmann et al., J. Phys. D 54,
  255202 (2021) — field reversal / electron-burst PIC physics (regime-bounded: matters at 13.56 MHz
  thin sheaths, not here). Open mirror: https://plasma.szfki.kfki.hu/~harti/resources/2021_Hartmann_JPD_54_255202.pdf
- Memos & Kokkoris, Micromachines 9, 415 (2018) — SEE modeling reference (SEE ruled out for the
  floor by recapture argument).

**ALE**
- Vella & Graves (2025), OSTI 2586627 — the three-layer site-balance ROM implemented verbatim
  (Eqs. 1–17 transcribed from the PDF after a scout summary proved wrong); window + dose gates.
  Open: https://www.osti.gov/pages/servlets/purl/2586627
- Kanarik et al., JVST A 33, 020802 (2015) + JVST A 35, 05C302 (2017) — the ALE synergy metric.

**Cryo**
- Antoun et al., Sci. Rep. 11, 357 (2021) — E_ads = 0.406 eV, residence-time law, −110/−120 °C
  cliff. Open: https://pmc.ncbi.nlm.nih.gov/articles/PMC7801591/
- Hsiao et al., Small Methods (2024) — the pseudo-wet CF₄/H₂ anchor (2.3 → 3.76 nm/s = 1.6×).
  Open: https://pmc.ncbi.nlm.nih.gov/articles/PMC11672179/ (the folk "2×" is actually 1.6×.)
- Lill & Berry, JVST A 41, 023005 (2023) — physisorption-etch theory framing.

**Bosch**
- Ayon et al., J. Electrochem. Soc. 146, 339 (1999) via the McVittie NNIN deck — Config R gates
  (28.2 µm / 434 nm / 140 nm / 250 nm). Deck: https://people.eecs.berkeley.edu/~pister/147fa14/Resources/BoschProc-STS.pdf
- Tillocher et al., Micromachines 12, 1143 (2021) — Config S (ultrafast, 60.8 nm/cycle, smooth).
  Open: https://pmc.ncbi.nlm.nih.gov/articles/PMC8537062/
- Park et al., Micro Nano Syst. Lett. 8, 14 (2020) — scallop measurement protocol.
- Ertl & Selberherr, Microelectron. Eng. 87, 20 (2010) — the academic 3D reference (validation-free;
  our bar). Open: https://www.iue.tuwien.ac.at/pdf/ib_2009/hashed_links/ep4PPErIJjnr4Y_us.pdf
- Laermer & Schilp, US 5,501,893 — the process definition.
- VLSet-AE, Microsyst. Nanoeng. (2026) — inverse SEM-measurement model (cited + distinguished:
  we are forward-predictive); its 16-run dataset = future Config T.

**Positioning / why-now**
- arXiv 2606.11231 / 2606.11247 (2026 perspectives) — name open differentiable fab-simulator
  infrastructure as the recognized unbuilt bottleneck.
- ViennaPS v3.6.0 (SoftwareX 2025) — the open baseline: no charging, no ALE gate, not differentiable.
- Coburn & Winters / Gottscho — ARDE conductance framing (petch transport is already Knudsen; the
  AR>20 collapse hunt is charging/IADF, not missing conductance).
