# C3 historical spec — directional plasma ALE (Si/Cl₂/Ar⁺)

## 2026-08-06 conservation and absolute-depth audit

This document began as an implementation spec. Its old “MD-validated” wording and acceptance anchors
must not be read as a current validation claim.

- The Vella–Graves model is an MD-regressed reduced-order model. It is useful for replay and gradients,
  but its printed equations are not elementally conservative: the coverage equations consume the
  SiCl₂ chlorine channel as `2 Y_SiCl2 θ²`, while the product equation emits SiCl₂ linearly in `θ`.
  At partial coverage the coupled equations create chlorine.
- At the first Vella–Hao experimental point (20.682 eV), the printed law gives 0.08134 nm/cycle for
  the reported 2 s Cl₂ + 3 s Ar-ion schedule. Replacing only the SiCl₂ product term by its
  atom-conservative `θ²` form gives 0.02913 nm/cycle, versus the digitized 0.11001 nm/cycle.
  Replaying the printed curve is therefore not an atomic-level validation.
- The current no-depth-fit result instead uses exact atom-counted endpoints from the released DeepMD
  `ALE/ALE.csv` trajectory, its independent bare-Si sputter table, and the experiment's measured
  positive-ion fluence. At 60/80/100 eV it predicts 0.4842/0.9453/1.4006 nm/cycle versus the
  digitized experimental curve's 0.5558/1.0453/1.5427 nm/cycle (12.9% maximum nominal error).
- This is a retrospective cross-source transfer, not a blind Tier-A prediction. The experiment did
  not measure a species-resolved IEAD; its plotted mean energy is inferred. The archived sputter
  uncertainty is reported separately, and missing energy/IEAD/model-form terms are not invented.

The executable record is
`results/curated/vella_hao_ale_depth/audit.json`; the dimensional implementation is
`src/petch/si_cl_ale_depth.py`.

Research scout verdict (2026-07-06): the highest-leverage next chemistry is **cyclic directional
Atomic Layer Etching**, implemented as a **three-layer site-balance ODE**. It is a cyclic wrapper
around chemistry petch already has (coupled coverages + IEDF + level set), specified in a 2025
MD-regressed paper, fast, trivially differentiable, and **unoccupied in open source** (ViennaPS: no
ALE, fixed chemistry; Kushner: closed/CPU; Neural Master Equation: needs MD training). CPU-buildable
now — no box.

## Primary benchmark
Vella & Graves (2025), "Si–Cl₂–Ar⁺ Atomic Layer Etching Window: A Fundamental Study Using MD and a
Reduced Order Model." OSTI: https://www.osti.gov/biblio/2586627 · PDF https://www.osti.gov/pages/servlets/purl/2586627
Synergy metric: Kanarik et al. JVST A 33, 020802 (2015) https://doi.org/10.1116/1.4913379 ;
predicting synergy JVST A 35, 05C302 (2017) https://www.osti.gov/servlets/purl/1376399

## The reduced-order model (historical source replay)
Three layers: top (θ₁, chlorinated surface), mixed (θ₂, Si-Cl mixed), crystalline. Per cycle:

**1. Modification step (Cl₂ dose, ANALYTIC — exact exponential, no integrator):**
  θ₁(t) = 1 − (1 − θ₁ⁱⁿⁱᵗ)·e^(−α t),   α = 2 J_Cl₂ γ_Cl / σ₁

**2. Bombardment step (Ar⁺, two coupled ODEs, mildly stiff):**
  dθ₁/dt = −(J_Ar/σ₁)( Y_Cl θ₁ + Y_SiCl θ₁ + 2 Y_SiCl₂ θ₁² + K_Cl^mix θ₁ )
  dθ₂/dt = −(J_Ar/σ₂)( Y_Cl θ₂ + Y_SiCl θ₂ + 2 Y_SiCl₂ θ₂² − K_Cl^mix θ₁ )

**3. Printed product/removal fluxes:**
  J_SiCl = J_Ar · Y_SiCl · (θ₁ + θ₂)
  J_SiCl₂ = J_Ar · Y_SiCl₂ · (θ₁ + θ₂)
  J_Si = J_Ar · Y_Si · (1 − θ₂)

The printed SiCl₂ product law is linear while its chlorine sink in step 2 is quadratic. It is retained
verbatim only to replay the published ROM; an atom-conservative product flux would be proportional to
`θ₁² + θ₂²` and does not reproduce the same depth.

**Energy-dependent parameters (all smooth → autodiff-clean), E in eV:**
  σ₂(E) = 0.77 (√E − √2.81) · σ₁            (mixed-layer depth; threshold ~2.81 eV)
  Y_Cl   = max(0, 3.2e-3 E − 4.25e-2)
  Y_SiCl = max(0, 1e-4  E − 1.4e-3)
  Y_SiCl₂= max(0, 2e-3  E − 2.9e-2)
  Y_Si   = max(0, 4.5e-5 E − 9e-4)

**Fixed constants:**
  σ₁ = 1e15 sites/cm² · γ_Cl = 0.25 · K_Cl^mix = 0.45 · J_Cl₂ = 9.8e17 cm⁻²s⁻¹ · J_Ar = 3.7e16 cm⁻²s⁻¹

## Historical replay gates (not experimental-validation gates)
1. **ALE window shape** (Fig 6, Si etched vs cycle, cyclic steady state after ~2 cycles):
   | Ar⁺ energy | regime | EPC (Å/cycle) |
   |---|---|---|
   | 15 eV | self-limiting (floor) | ~0.67 |
   | 17.5 eV | self-limiting | ~0.75 |
   | 20 eV | self-limiting (ceiling) | ~1.0 |
   | 22.5 eV | onset sputtering | ~1.7 |
   | 30 eV | pure sputtering | ~6 |
   → recover the **15–20 eV ALE window** and the non-self-limiting rise ≥22.5 eV.
2. **Cl uptake plateau** 0.8×10¹⁵ Cl/cm² per modification step (Fig 7).
3. **Experimental window-floor EPC** ≈ 0.68 Å/cycle (Matsuura, Kim — independent).
4. **Synergy** S = (EPC − α − β)/EPC (α = mod-only etch, β = removal-only etch): S→~100% inside
   the window, collapses outside.

## Notes / open parameters
- EPC[Å] = N_Si_removed[cm⁻²] / n_Si, n_Si ≈ 5e22 cm⁻³ · 1e8 Å/cm → /5e14. (0.67 Å ≈ 0.34 ML.)
- Bombardment dose/time not given explicitly in the scout report — Si-removed-per-cycle = ∫J_Si dt
  over the bombardment step; integrate to self-limiting saturation. May need the paper's bombardment
  time to nail the ABSOLUTE anchor; the WINDOW SHAPE (floor/ceiling/sputter-rise) is the real gate and
  should emerge from the yields alone. Fetch the OSTI PDF for the bombardment schedule before final gate.
- Numerics: modification step exact-exponential (zero cost); bombardment via exponential integrator or
  backward-Euler + implicit-function-theorem adjoint (differentiate the solution, not solver iters).

## Follow-on (same skeleton)
Bolt a temperature-dependent physisorption term onto the site-balance skeleton → **cryo SiO₂/Si₃N₄**
(the hot 2024-2026 unclaimed chemistry). Gate: CHF₃ cryo "rate doubles +20°C→−60°C" (JAP 133,113306,
2023); pseudo-wet HF/H₂O (Small Methods 2024, doi 10.1002/smtd.202400090). Thermal Al₂O₃ ALE (George,
Chem Mater 28, 2994) is the cleanest Arrhenius EPC-vs-T regression unit-test (0.14 Å @250°C → 0.75 @325°C).
