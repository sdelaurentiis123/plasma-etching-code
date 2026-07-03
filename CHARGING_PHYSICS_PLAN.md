# Charging physics: getting the mechanism config over the line

Deep plan, literature-first (researched 2026-07-03). Goal: the five simultaneous Workstream-A gates
from NEXT_ARC_PLAN.md — floor-flux 8-pt RMSE ≤ 0.05 *with the conductor on*, V_c(AR4) → 33 V ± 40%,
foot-ion energy RISING 15 → 27.5 eV, poly-line/foot-flux no-regress, Matsui open — so the notch can
run on OUR foot energies instead of the published E_defl table. Working rules unchanged: every step
gated, nothing tuned, honest fails documented.

## 0. Correction of the 2026-07-03 "framework limit" conclusion

The RF-phase experiment (RECONCILIATION.md) concluded the deep-AR electron deficit was a
"collisionless-2-D framework limit." **The literature refutes that.** Three anchors:

1. **Hwang & Giapis themselves** (JVST B **15**, 70 (1997) — the data we gate against) ran a 2-D
   collisionless MC with "realistic energy and angular distributions for ions and electrons" and got
   V_c = 33 V at AR 4. A 2-D collisionless model IS sufficient. Their electron/ion arrival
   distributions come from a **1d/2v Monte Carlo RF-sheath simulation** feeding a 2d/3v in-feature
   flux calculation (also: for a trench, v_y is conserved in-feature, so 2-D x–z dynamics is exact —
   the "3-D electron model" caveat applies to holes, not lines/trenches).
2. **Memos & Kokkoris** (Micromachines **9**, 415 (2018), open access) trace fully **isotropic**
   Maxwellian T_e = 4 eV electrons — wider-angle than anything we tried — with a variable-step stiff
   integrator (Matlab ode15s) and signed-distance-function termination. No orbit-dropping problem.
3. Same paper: **secondary electron-electron emission (SEEE) halves the charging potential**
   (45 → ~22 V) because secondaries/backscattered electrons "terminate preferentially at the
   [positive] valley region" — the exact redistribution channel our deep floor is missing.

So the deficit decomposes into (a) an **integrator artifact** (our fixed-step trace, cap 14·nz,
silently drops precisely the well-captured electrons that matter), (b) a **missing physics channel**
(SEE from the walls), and (c) an **arrival-distribution shortcut** (analytic weighting instead of a
sheath-MC-derived distribution). All three are fixable inside our 2-D collisionless framework.
RECONCILIATION.md gets a correction note when W1 lands numbers.

## W1 — Numerics: every electron trajectory terminates (do first, it re-baselines everything)

**Problem.** `charging2d.trace()`: fixed `dt = 0.45/max(|v|)`, hard cap `14*nz` steps, survivors
silently dropped. A 4 eV electron captured by the +40 V floor well orbits many times before landing;
the cap drops it → the deep floor is under-supplied *by the integrator*, and every electron-model
comparison made on top of this (first-order vs memoryless RF-phase) is confounded.

**Fix.**
1. Leapfrog (kick-drift-kick) stepping — symplectic, no numerical energy drift, so captured orbits
   neither artificially decay nor artificially escape. (Memos used adaptive ode15s; leapfrog + local
   dt is the vectorized equivalent for our many-particle trace.)
2. Adaptive local dt: `dt = min(0.45/|v|, 0.3/sqrt(|E_local|))` — resolve strong-field regions near
   the charged floor instead of blowing through them.
3. Step budget raised until the **survivor fraction gate** passes: alive-at-cap < 0.1% of launched
   electrons at AR 4 in the converged field. Instrument and PRINT the survivor fraction + a landing
   budget (floor/poly/PR/top) every gate run — no more silent drops.
4. Escape guard: particle above the mouth plane moving away with E_kin > |qV_local| → terminate as
   "returned to plasma" (physical, not a drop).

**Gates (W1).** (i) survivor fraction < 0.1% at AR 1–4; (ii) energy conservation on force-free paths
to machine precision, and < 1% drift per orbit in the well; (iii) re-run the 8-pt charging gate +
notching gates → new honest baseline for V_c, flux RMSE, foot energy. Expectation from the physics:
V_c(AR4) drops from 43.5 toward 33 as the dropped electrons land (they land mostly on the floor —
that is where the well pulls them). If V_c moves < 2 V, say so and move on — W2 is independent.
**Also re-test the memoryless RF-phase arrival model** (reverted 2026-07-03) on top of the fixed
integrator before discarding it permanently: its failure (V_c = 61) is confounded by the artifact.

**Effort.** 1–2 days. Pure numerics, no new constants. Keep the trace vectorized (the loop is
per-step over all alive particles already).

## W2 — Physics: secondary electron emission from the feature walls (the missing channel)

**Anchor.** Memos & Kokkoris 2018 (above): SEEE cut the charging potential ~50% in a directly
comparable geometry (polymer surface, T_e = 4 eV, Laplace + ballistic tracing). We need ~24%
(43.5 → 33). Kushner-lineage MCFPM work includes the same channel. Every constant below is
published — no tuning knobs.

**Model (implement verbatim from the paper).**
- Total electron yield σ_e(E) for polymer/PR: Dapor's Monte-Carlo curve (the only published set
  covering 0–50 eV; digitize from the Memos Fig.; cite both).
- Backscattering: Burke `η(E) = 0.115 (E/1000)^-0.223`; secondary yield `δ = σ_e − η`, clamped
  δ = 0 below 16 eV (Memos' constraint).
- Emitted secondaries: 1 eV (Seiler most-probable), cosine-isotropic from the local wall normal.
  Backscattered: elastic (retain E), cosine re-emission.
- Material-dependence: PR/polymer walls use the Dapor-PMMA curve; poly-Si line and oxide floor get
  their own published curves if the PR-only version under-delivers (poly-Si σ_max ≈ 1.1 @ ~250 eV —
  smaller at our energies; oxide σ higher). v1: walls-only (PR), floor absorbs — the mechanism is
  wall-emission feeding the floor.
- Each electron impact on an insulating wall: with probability η backscatter, else with probability
  δ emit a secondary (Monte-Carlo branch, weight-conserving); every emitted electron is traced like
  a primary through the same field (it gets pulled into the positive well — that is the mechanism).
  Cap the cascade generation at 3 (yields < 1 below 50 eV ⇒ geometric decay; assert total emitted
  weight per primary < 2).
- Charge bookkeeping: emission site gains +δ|e| per emitted electron (the wall charges POSITIVE
  where it emits) — this is also the published reason upper-PR walls don't pin hard negative, which
  our insul_vmin clip approximates today. SEE may let us RETIRE that clip: check whether the
  Vprwall floating potential lands in the −T_e band on its own.

**Gates (W2, on top of W1).**
1. V_c(AR4) → 33 V ± 40% (primary).
2. Floor-flux 8-pt RMSE — target ≤ 0.05; must not regress above the W1 baseline.
3. Foot-ion energy trend re-measured (needs V_c − V_poly to open with AR; SEE raising the deep-floor
   electron supply is exactly what lets the floor discharge relative to the poly line).
4. Matsui 300 eV unchanged; poly-line curve no-regress.
5. Sanity vs the literature effect size: SEEE ΔV_c should be a REDUCTION of roughly 20–50%
   (Memos: ~50%); if it overshoots below ~15 V at AR 4, the yield curve is being applied outside
   its validity — document, don't rescale.

**Effort.** 2–4 days (yield-curve digitization + MC branch + charge bookkeeping + gates).

## W3 — Source: 1-D RF-sheath MC for the joint (E, θ, φ) arrival distributions (the HG-faithful endgame)

**Anchor.** HG's own method: "surface flux calculated using the particle distribution function,
obtained from a Monte Carlo sheath simulation" (1d/2v, RF-modulated, collisionless), then traced
in-feature. The SCIRP 3-D implementation does the same with PIC (XPDC1). Nobody in the validated
literature uses an analytic electron arrival model — this is the structural upgrade that removes
our last unpublished approximation (the first-order burst weighting AND the memoryless variant are
both homemade).

**Build.** `src/petch/sheath1d.py`, ~200 lines, standalone + gated:
- Collisionless 1d3v RF sheath, Child-law profile, V_s(φ) = V_dc + V_rf·sin(φ) (Lieberman V_dc as
  now); ions enter at Bohm speed with T_i ≈ 0.05 eV transverse, integrate through the oscillating
  field (standard: gives the bimodal IEDF and the energy–angle anticorrelation we already impose
  analytically — now derived, one fewer assumption).
- Electrons: half-Maxwellian flux (T_e) injected at the sheath edge each phase bin; advance through
  the instantaneous field; record (E, θ, φ_RF) of those crossing the wafer plane. This produces the
  burst structure (arrivals cluster at sheath minima) with the correct residual-energy and angle
  correlations — no memoryless assumption, no cos^p guess.
- Feed both sampled distributions into `solve_trench_charging` as lookup tables
  (`arrival="sheath-mc"` flag; the analytic models stay reachable for A/B).

**Gates (W3).**
1. Sheath-MC IEDF reproduces the bimodal split ΔE and the θ(E) anticorrelation we currently impose
   (internal consistency — the ion side is already gated at the wafer level).
2. Electron arrival: burst fraction and EADF vs the published qualitative shape (quasi-isotropic,
   sheath-collapse-dominated); document the measured EADF width.
3. The five Workstream-A gates simultaneously — THE finish line.

**Effort.** 3–5 days. Independent of W1/W2 code paths but only meaningful measured on top of them.

## Order, budget, kill criteria

**W1 → W2 → W3**, full `charging_gate.py` + `notching_gate.py` after each, committed with measured
numbers (each full 8-pt gate ≈ 25–30 min CPU-local; PETCH_DEVICE=cpu, no box needed). W1 first
because it un-confounds every subsequent measurement; W2 is the highest-leverage physics; W3 removes
the last homemade approximation and is what "first-principles charging" means in print.

**Kill criteria.**
- W1: if survivor-fraction < 0.1% requires > 100× the step budget (CPU blow-up), switch the deep-well
  capture to an analytic orbit-average landing model (documented approximation), don't burn weeks.
- W2: if PR-only SEE moves V_c(AR4) < 3 V, add the poly/oxide curves once; if still < 3 V, the
  channel is not the lever in our geometry — document with the landing budget and stop.
- W3: if the sheath-MC electron EADF is statistically indistinguishable from the analytic memoryless
  model (K-S test on (E,θ)), the source was never the problem — record that as a (publishable)
  negative result and skip the integration.
- Global: if all three land and the five gates still don't hold simultaneously, the residual is the
  conductor/geometry model itself — bisect poly_um/pitch/boundary_um per the original plan, then
  stop. The honest plateau, precisely bounded, ships in the paper.

**What "over the line" buys.** The notching gate switches from HG's published E_defl(AR) table to
petch-computed foot energies → the notch becomes a fully first-principles prediction (charging solver
→ deflected energies → notch depth vs measured Nozawa/Fujiwara/HG), which no open code does end to
end. That plus the wafer-ARDE prediction is the paper.

## Sources

- Hwang & Giapis, JVST B **15**, 70 (1997) — method + V_c=33 V reference; 2-D MC, sheath-MC-derived
  distributions. https://pubs.aip.org/avs/jvb/article-abstract/15/1/70/470767
- Hwang & Giapis, JAP **82**, 566 (1997) — the gate curves (floor flux, potentials).
- Memos & Kokkoris, Micromachines **9**, 415 (2018), open access — SEEE halves charging potential;
  Dapor σ_e(E) for PMMA; Burke η(E); 1 eV isotropic secondaries; ode15s + SDF termination.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6187714/
- Zhang et al., Plasma Process. Polym. (2020; 2025) — EAD width controls electron penetration to the
  trench bottom. https://onlinelibrary.wiley.com/doi/10.1002/ppap.202000014 ;
  https://onlinelibrary.wiley.com/doi/10.1002/ppap.70037
- 3-D SiO2 charging implementation (open access) — PIC(XPDC1)-derived EEDF/IEDF as boundary input,
  FEM Laplace, absorb-on-hit. https://www.scirp.org/html/1-8102035_41882.htm
- Dapor, secondary-electron yield of PMMA (MC calculations, 0–1500 eV) — the σ_e(E) source cited by
  Memos; digitize via the Memos figure.
