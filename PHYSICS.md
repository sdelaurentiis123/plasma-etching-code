# petch — full physics breakdown

The complete physics of the simulator, module by module: the governing equations, what is **derived**
vs **published-parameter** vs **calibrated** (graded honestly), the benchmark gates each module passes,
and the known limits. Everything here is committed on `frontier-loop` and covered by the test suite.

---

## 0. Architecture in one paragraph

petch is a feature-scale plasma-etch simulator: a **level-set** etch front advanced by rates from
**surface chemistry** (coupled-coverage kinetics), fed by **particle transport** (Monte-Carlo ballistic
/ Knudsen neutrals + directional ions with derived energy-angle distributions), with **surface
charging** solved self-consistently (charged-particle tracing in the evolving electrostatic field) and
**cyclic process drivers** on top (Bosch, ALE). Chemistry modules are smooth/differentiable
(reverse-mode autodiff demonstrated end-to-end). The moat: open + GPU-capable + differentiable +
charging + reactor↔surface coupling, gated against published data at every claim.

---

## 1. Transport (ions & neutrals)

**Neutrals:** Monte-Carlo ballistic transport with sticking/re-emission (Knudsen multi-bounce for deep
features), or radiosity. Parity: petch MC == ViennaPS on sub-micron trench (within ~0.08); petch-DDA
reproduces the ViennaPS ARDE point (0.727 vs 0.73 @ AR 8.6).

**Ions (DERIVED, this campaign):** the RF-sheath source is now first-principles
(`sample_sheath_source`). At 400 kHz the ion transit is ~1.5% of the RF period (ωτᵢ ≈ 0.015), so
instantaneous crossing is exact:
- energy: E = V_s(φ) + Te/2 with φ ~ uniform → the **arcsine-bathtub bimodal IED** emerges
  (verified: P(E<33 eV) = 0.434 sampled = 0.434 analytic phase fraction)
- angle: transverse thermal velocity (one dof, Ti/2) conserved while v_z grows →
  **IADF HWHM 5.2°** and the low-energy-is-wide anticorrelation emerge (HG report 4.3°; the
  residual is their nonlinear-sheath field curvature)

**Grade:** transport numerics = physics; ion source = derived; neutral sticking coefficients =
published (ViennaPS-inherited).

---

## 2. Surface chemistry — continuous etch (SF6/O2)

Langmuir-Hinshelwood **coupled coverages** (ViennaPS "Belen" model): fluorine θ_F and oxygen θ_O
with a = (kσ + 2·G·Y_ie)/G_E, b = (βσ + G·Y_p)/G_P, θ_F = 1/(1+a(1+1/b)), θ_O = 1/(1+b(1+1/a));
etch rate ER = (1/ρ)(kσ·θ_F/4 + Y_sp·Γ_ion + θ_F·Y_ie·Γ_ion). Energy yields A(√E−√E_th); ViennaPS
angular forms.

**Gate:** de Boer trench etch — knee matched by etchant-starved process parameters + narrow IADF;
MC parity with ViennaPS.
**Grade:** structure = physics; kinetic parameters = published (Belen); absolute-rate anchors for
de Boer = calibrated process parameters (fluxes the paper didn't publish — stated, not hidden).
**Limit:** high-AR (>20) floor collapse remains the open frontier.

---

## 3. Atomic Layer Etching (Si/Cl2/Ar+) — `ale.py`, `ale_diff.py`

Vella–Graves (2025, OSTI 2586627) three-layer site-balance ROM, implemented from the paper's
equations (a research-agent summary had a wrong one-term Si flux — caught by the gate, fixed from
the primary source):
- modification (analytic): θ₁(t) = 1 − (1−θ₁⁰)e^{−αt}, α = 2·J_Cl2·γ_Cl/σ₁
- bombardment ODEs: dθ₁/dt = −(J_Ar/σ₁)(Y_Cl θ₁ + Y_SiCl θ₁ + 2Y_SiCl₂θ₁² + K θ₁);
  dθ₂/dt = −(J_Ar/σ₂)(… − K θ₁)
- Si leaves in THREE channels: J_SiCl + J_SiCl₂ (∝ θ₁+θ₂, the window) + bare-sputter Y_Si(1−θ₂)
- yields linear in E (MD-fit, from the paper); σ₂(E) = 0.77(√E−√2.81)σ₁

**Gates (nothing tuned):** the 15–20 eV ALE window (EPC 0.76/0.86 vs ROM 0.7/0.9), loss of
self-limitation above 20 eV (4.74 @ 30 eV vs 4.8), Kanarik synergy 100%→20%, **dose-saturation
dynamics** (Fig 9: saturates below 1e18 cm⁻², plateau 0.757 vs ~0.75), Cl uptake plateau.
**Differentiable (`ale_diff.py`, torch):** exact reverse-mode dEPC/dE (matches finite-diff to 4
decimals); gradient inverse design recovers ion energy from target EPC to <0.1 eV.
**Grade:** fully derived at ROM level (all constants from the paper's MD fits).
**Limit:** 15 eV window floor runs low (0.36 vs 0.6 — near-threshold sensitivity).

---

## 4. Cryogenic etch — `cryo.py`

Temperature-gated physisorbed etchant layer: Langmuir isotherm θ(T) = Kp/(1+Kp),
Kp = A·e^{E_ads/kT}; ER(T) = R_base(1 + gain·θ(T)). E_ads = 0.40 eV **fixed from independent
measurement** (Antoun 2021: 0.406 eV); residence time t_d = t_d0·e^{E_d/kT}.

**Gates:** CF4/H2 pseudo-wet SiO2 (Small Methods 2024): ER(+20°C)=2.32 (plateau 2.3),
ER(−60°C)=3.76 (exact), ratio 1.62 (the true anchor is 1.6×, not the folk "2×" — gated on truth).
Cross-checks: HF cryo-ALE ~3.2× (different system), C4F8 on/off cliff −120 vs −110 °C reproduced.
**Grade:** isotherm form + E_ads = physics/measured; A and gain = 2 params calibrated to 1 anchor.

---

## 5. Feature charging — `charging_general.py` (the deep campaign)

The self-consistent loop: launch ions+electrons from the derived sheath source → trace through the
field (leapfrog, adaptive dt) → deposit on surfaces → insulator cells float to local current
balance, conductor components float as equipotentials → re-solve field → iterate to steady state.

**Every piece was rebuilt first-principles this campaign, each fix forced by a measurement or a proof:**
1. **Integrator (C6):** coarse dt lost focused floor-bound electrons (traced flux < geometric —
   physically impossible). dt=0.15 cell/step converged; focusing then EMERGES (e_traced 2.9× geometric).
2. **Boundary conditions:** periodic-wrap tracer teleported grazing electrons out of the open area
   (instrumented current budget: outer wall could only balance the deflected-ion torrent at +20 V).
   **Mirror-image BCs** (HG's method, matching the field's Neumann sides) → edge 20→8–11 V.
3. **Charge dynamics:** the linear step + decaying anneal froze mouth walls at the −67 V clip
   (over-collimation). **Log current-balance update** ΔV = clip(Te·ln(Γᵢ/Γₑ)) — quasi-Newton for
   exponentially-retarded electron flux — reaches each cell's fixed point exponentially fast, plus
   the physical bound −10·Te (top of the Maxwellian tail). ~3× faster convergence.
4. **Source (C9, proven):** ions = exact bathtub (§1). Electrons = Lambert × Maxwellian ×
   e^{−eVs(t)/kTe} collapse-burst weighting, backed by the **invariance theorem**: for cosine-flux
   Maxwellian injection, barrier selection and refraction cancel exactly —
   v_z e^{−mv_z²/2kTe}dv_z = e^{−eVs/kTe} v_z′ e^{−mv_z′²/2kTe}dv_z′ — arrival EADF = cos θ for ANY
   waveform. Sheath quasi-static to 1e-3 (75 ps transit vs 2.5 µs period). **HG's published cos^0.6
   is their injection-convention artifact** (uniform-in-angle launch; closed form gives p≈0.79 for
   their convention). Field reversal rejected at these conditions by mass-ratio flux balance.
5. **Ruled out with physics:** SEE at the floor (a few-eV secondary in a +33 V well is recaptured —
   zero net); fine grid as the fix (floor plateaus 41.6→38.5→37.7 at W16/32/48); the GPU as the fix
   (loop is CPU-field-solve-bound at these sizes).

**State vs the HG 1997 benchmark (AR4; floorV/flux/edge/neigh, HG: 33/0.22/7/39):**
fully-derived source → 38.7/0.339/14.4/34.5; every observable individually reachable
(floor 32.4 ✓, neigh 39.9 ✓, edge 8.3, foot peak 61 V ≈ HG 60); all structure emergent
(dipole walls, edge<neighbor split, focusing, deflected-ion feedback). Joint residual ~15–20%
attributed to: HG's own convention artifacts, 2D-projection conventions, idealized rectangular mask.
**Grade:** the source and dynamics are now derivation + proof; no tuning knobs in the pipeline
(all crutches off; legacy knobs retained default-off for A/B only).

---

## 6. Notching (charging → etch coupling) — `threed.py::_apply_hg_charging`

The charged floor decelerates landing ions (yields re-evaluated at E−eV_f; sub-threshold IEDF slice
removed) and the removed slice is **deflected into the sidewall foot** at E_defl within a fixed
physical band (0.3·W, the corner-field scale) — the notch driver. Gates passed (commit c07886c):
notch(charging OFF) = 0 at every AR; Fujiwara-monotone; HG shape correlation r=0.92.
**C11 (done, honestly graded):** the charging table (Q, V_f, E_defl vs AR) is now also generated by
petch's own derived-source solver (`surface_charging="petch"`; AR4: Q=0.404, V_f=39.7, E_defl=28.0 —
the foot energy lands on HG's 28). Under the petch table the mechanism gate PASSES (charging-off = 0
everywhere; resolved notches at AR≥2, larger than the HG-closure mode at AR2–3), but the AR4 trend
gates FAIL: the notch forms then is erased during late overetch — the table's documented +15–20% AR4
over-charge (V_f 39.7 vs 33, Q 0.40 vs 0.22) amplifies through the overetch redeposition/erosion
balance. The validated default remains the published-closure table; the petch mode is the
first-principles path with the AR4 coupling sensitivity documented as the frontier. **Limits:**
absolute notch depth uncalibrated (3–5× shallow of HG); sub-micron grid stability separate.

---

## 7. Bosch DRIE — `bosch.py`

Cycle mechanics (per cycle): conformal passivation → ion punch-through clears **near-horizontal**
surfaces only (flux ∝ cos incidence protects scallop feet; only the mask casts a hard shadow —
real ion divergence clears crest slivers) → **sequential etch**: the directional punch advances the
floor d_dir FIRST, then the isotropic F-neutral front expands r_iso from the NEW floor.

The sequential semantics is load-bearing and was **forced by the data**: simultaneous (swept-disc)
mechanics provably give s = r−√(r²−((p−d)/2)²) ≈ 32 nm; sequential gives s = r−√(r²−(p/2)²) = 140 nm
— Ayon's measured value, with advance = d+r = 434 and undercut ≈ r = 238 simultaneously.

**Gates (two published regimes + kill-test, 2 rates calibrated from endpoints, all else emergent):**
- Config R (Ayon 1999): depth 28.6 (28.2±2.8) ✓, pitch 440 (434±43) ✓, scallop 140 (140±35, exact) ✓,
  undercut 220 (250±50) ✓
- Config S (Tillocher 2021 ultrafast): pitch 60.0 (60.8±6) ✓, scallop 15 (≤30) ✓, D→60 µm ✓
- Cross-config s_R/s_S = 9.3 (≥4) ✓
**Grade:** mechanics = physics (emulation tier); per-cycle rates = calibrated from published
endpoints. **Limit:** ARDE sub-gate (0.82) deferred to the transport tier; VLSet-AE (2026) cited and
distinguished (inverse measurement model vs our forward prediction).

---

## 8. Differentiability & inverse design

Demonstrated end-to-end on chemistry (`ale_diff.py`): exact reverse-mode gradients through the full
cyclic site-balance model (autograd == finite-diff to 4 decimals), Newton inversion of process knobs
(ion energy from target EPC, <0.1 eV). The cryo and charging-source closures are smooth by
construction (exp/sqrt/log forms) — the same pattern extends. This is the LLM-driveable inverse-design
layer: every derivation-clean module adds an invertible knob.

---

## 9. Performance (current state, pre-port)

~14× vs ViennaPS-GPU on the core etch benchmark (honest number, weak-box artifact corrected).
Charging: Warp/CUDA tracer exists (24.6 M particles/s, parity-checked, resolution-parameterized)
but the loop is **CPU-field-solve-bound** at small grids (GPU only ~35% net) — the device-resident
port (persistent arrays + cell-sort + CUDA-graph + on-GPU multigrid field solve, est. 5–20×) is
deliberately **deferred until after physics lock** (this document is the lock).

---

## 10. Honest open frontiers

1. Charging joint-residual (~15–20% vs HG absolutes) — convention-level; next real physics would be
   a time-resolved Vlasov sheath integrator (only matters ≥13.56 MHz with thin sheaths).
2. High-AR (>20) etch floor collapse (de Boer) — genuine frontier.
3. Notch absolute depth calibration + sub-micron grid stability.
4. Bosch ARDE transport tier (multi-bounce neutrals in the cycle driver).
5. Full-Poisson charging (grounded substrate + physical-unit σ-sheet) — infrastructure built,
   stable solver documented as the remaining build (CHARGING_POISSON_PLAN.md).
6. The speed pass (device-resident GPU engine) and the LLM/inverse-design API layer.

## Test & provenance state
26+ tests green (chemistry 15, charging 4, Bosch 3, smoke 4+); ~18 commits on `frontier-loop`;
every cycle logged with gates in FRONTIER_LOOP.md; refuted hypotheses documented in
FLOOR_OVERCHARGE_FINDING.md and CHARGING_POISSON_PLAN.md.
