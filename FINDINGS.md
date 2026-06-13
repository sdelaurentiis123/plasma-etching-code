# Phase 0 — Findings log

Running record of measured results as we retire the bias contributors. Numbers are on this
M1 Mac (CPU), 2D, vs the cached ViennaPS reference (width-8 depth 10.05 µm; normalized ARDE
`[0.932, 0.969, 0.985, 1.000]` for widths `[4, 6, 8, 12]`).

---

## Step 1 — Baseline parity
The `petch` package reproduces the original PoC exactly: width-8 center depth **9.18 µm** with
`rate_scale=0.29` (cached PoC value 9.40 — the 0.2 µm difference is MC run-to-run variance).

## Step 2 — Convergence harness
- **Grid:** depth plateaus by `dx=0.25` (0.25→0.125 changes depth 0.018 µm, below the
  seed-noise 0.063 µm). So the benchmark grid is essentially converged; contributor #5's
  *grid-resolution* part is small at `dx=0.25`. The −8.7% gap is therefore **not grid**.
- **Rays:** MC noise (across-seed std of depth) falls 0.25 → 0.12 µm over N = 5k → 40k,
  roughly `1/√N`. This is the noise floor QMC (Step 4) must beat.

## Step 3 — Contributor #1 (Belen coupled-coverage chemistry)
Implemented the exact ViennaPS Belen model (`src/petch/belen.py`): two coupled coverages with
`k_sigma=300`/`beta_sigma=0.04`, the standalone chemical-etch term `k_sigma·θ_F/4`, and
ViennaPS sticking `β_E=0.7 / β_O=1.0` (transport re-emission made consistent in the driver).

| model | width-8 depth | normalized ARDE | ARDE rmse |
|---|---|---|---|
| ViennaPS (ref) | 10.05 | 0.93, 0.97, 0.99, 1.00 | — |
| Langmuir @ rate_scale=0.29 | 9.18 (−8.7%) | 0.88, 0.93, 0.95, 1.00 | 0.038 |
| **Belen @ rate_scale=0.10** | **10.24 (+1.8%)** | 0.75, 0.86, 0.92, 1.00 | **0.110** |

**Result — two-sided:**
1. **Belen fixes the absolute depth** (+1.8% vs −8.7%) at its calibrated scale.
2. **Belen makes the ARDE lag *steeper* (worse)** — rmse 0.110 vs 0.038.

**Diagnosis (the valuable part):** Belen adopts ViennaPS sticking `β_E=0.7` (vs the PoC's 0.2).
Higher sticking consumes F near the trench top, so narrow/deep trenches get F-starved → steeper
lag. But **ViennaPS uses 0.7 too and still gets a flat ARDE**, which means our **neutral
transport over-depletes narrow trenches** (re-emission is insufficient — only 12 truncated
bounces, pseudo-random). The PoC's artificially-low `s_F=0.2` was *masking* this transport
deficiency by letting more F survive to the bottom.

**Conclusions:**
- Contributor #1 (chemistry) is **entangled with neutral transport** via the sticking
  coefficient; it cannot be cleanly retired in isolation.
- The real lever for ARDE is **better neutral re-emission** — the radiosity neutral solve
  (Step 4). That step now does double duty: variance-reduction speedup **and** the ARDE fix
  (restore flux to narrow-trench bottoms → flatten ARDE back toward ViennaPS).
- `rate_scale` (0.10 for Belen, 0.29 for Langmuir) is confirmed a **units / flux-normalization**
  issue (our `mc_flux` emits dimensionless open-field multipliers, not absolute ViennaPS fluxes
  in 1e15 cm⁻²s⁻¹), **separate from chemistry form**. Killing it needs flux-unit reconciliation,
  not a better coverage model.

## Step 4 — Eng speedups + transport diagnostics

**(a) Bounce-truncation hypothesis — REFUTED.** Belen ARDE at `n_reemit` 12 vs 40 is *identical*
(depths 8.38/10.24/11.16, same runtime) — the re-emission chain is already converged at 12
bounces. So the narrow-trench over-depletion is **not** bounce truncation. The cause is deeper
(re-emission angular law, flux normalization, or a saturation difference vs ViennaPS) — to be
resolved by the radiosity reformulation, not by adding bounces.

**(b) QMC (Sobol source launch) — WIN.** At fixed N=10000 (compact geom, 6 seeds):

| sampling | depth mean | across-seed noise (std) |
|---|---|---|
| pseudo | 6.409 | 0.224 µm |
| **sobol** | 6.443 | **0.118 µm** |

Same mean (no fidelity cost), **1.9× lower noise** → ~3.6× fewer rays for the same noise floor.
This is the "equal fidelity at far fewer rays" win, exact and measured.

**(c) Ions are over-sampled, but ions are NOT the bottleneck.** Width-8 belen depth vs ion ray
count: 20k→10.236, 5k→10.166 (unchanged), 1k→6.377 (collapses). So ions tolerate a ~4× ray cut
— BUT raytrace time barely moved (3.3s→2.2s) because the **many-bounce neutrals dominate the ray
cost**. Therefore the high-value speedup is the neutral solve, not the ion few-ray.

**Priorities converge:** the **neutral radiosity solve** is simultaneously (i) the dominant speed
lever (neutrals are the cost) and (ii) the ARDE fix (it recomputes the neutral equilibrium flux,
which is where the over-depletion lives). That is the next major build.

**(d) Neutral radiosity solve — built, and it OVERTURNS the hypothesis.** Implemented the exact
2D form-factor radiosity (`src/petch/radiosity.py`): `Γ = (I-(1-s)A)^-1 D`, all bounces, noise
free. It is consistent with MC (small-grid depth 6.67 vs 6.7; fully-exposed m_F → s). But on the
width-8 ARDE sweep it made ARDE **steeper**, not flatter:

| neutral model | norm ARDE (w4,w6,w8,w12) | ARDE rmse |
|---|---|---|
| belen + MC (12 bounces) | 0.751, –, 0.917, 1.00 | 0.110 |
| **belen + radiosity (exact)** | 0.697, 0.808, 0.889, 1.00 | **0.150** |
| ViennaPS | 0.932, 0.969, 0.985, 1.00 | — |

Since radiosity is the *exact* neutral physics and still gives steep ARDE, the steep ARDE is
**not** a neutral-transport-resolution problem. (Radiosity is currently O(M^3)-unoptimized →
~40 s/width vs MC ~5 s; its value here is determinism + differentiability, not speed-as-built.)

**(e) ROOT CAUSE FOUND — one cause for both symptoms.** Our model is in the **neutral-limited**
regime (θ_F flux-sensitive → F-depletion in narrow trenches drops the rate → steep ARDE);
ViennaPS's flat ARDE implies **F-saturated / ion-limited**. Confirmatory test — raising the F
flux toward saturation monotonically flattens ARDE toward ViennaPS:

| F flux | norm ARDE (w4,w8,w12) | rmse vs ViennaPS |
|---|---|---|
| ×1 | 0.751, 0.917, 1.00 | 0.112 |
| ×3 | 0.803, 0.935, 1.00 | 0.080 |
| ×6 | 0.815, 0.954, 1.00 | 0.070 |

**Both the absolute-depth gap (`rate_scale`) AND the ARDE-shape gap trace to the SAME root
cause: flux normalization.** Our `mc_flux` emits dimensionless open-field multipliers, not
ViennaPS's absolute fluxes (1e15 cm⁻²s⁻¹), so the F-to-ion ratio in the coverage sits in the
wrong regime. Reconciling the absolute flux units should fix ARDE and the absolute rate
**together**, likely retiring `rate_scale` as a byproduct.

**Phase-0 payoff:** the vague "~7% accuracy issue + ARDE lag" is now ONE localized,
evidence-backed root cause (flux normalization) with a demonstrated lever (the F/ion flux ratio).

**Next:** reconcile `mc_flux` flux normalization to ViennaPS absolute units. This needs ViennaPS's
flux convention + `unitConversion` (from source) and ideally one ViennaPS calibration run
(GPU/Linux box — not available on this M1). That single change is predicted to close both gaps.
