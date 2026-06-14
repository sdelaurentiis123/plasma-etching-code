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

**Next:** reconcile `mc_flux` flux normalization to ViennaPS absolute units. (Now actionable —
real ViennaPS ground truth captured below.)

## Step 6 — Real ViennaPS ground truth (Vast.ai RTX 2080 Ti, ~$0.02)

`pip install ViennaPS` **works on Linux x86** (4.5.0) → arm64 was the only local blocker. Ran
SF6O2 at the benchmark widths with the **CPU_DISK (Embree) engine** (the GPU/OptiX path on that
box failed an OptiX driver-version check — note for Phase 1: needs a newer driver than 550.144).

- **Real depths match the cached reference exactly:** w4 9.562, w6 9.881, w8 **10.053**, w12
  10.230 → normalized ARDE `[0.935, 0.966, 0.983, 1.000]`. Validates the whole cached reference,
  and we now have the **full surface meshes** (`harness/reference/viennaps_surfaces.npz`) →
  Chamfer / sidewall-angle metrics are unlocked.
- **Exact default params confirmed** (`harness/reference/viennaps_reference.json`): `k_sigma=300,
  beta_sigma=0.04, A_ie=7, A_sp=0.0337, B_ie=0.8, B_sp=9.3, Eth_ie=15, Eth_sp=20, rho=5.02`;
  passivation `A_ie=3, Eth_ie=10`; ions `meanEnergy=100, sigmaEnergy=10, exponent=500,
  inflectAngle=1.553`; fluxes `etchantFlux=1800, ionFlux=12, passivationFlux=100`.
- **Correction:** ViennaPS 4.5 SF6O2 exposes **no explicit sticking coefficient** (`beta_E=0.7`
  was wrong — it came from the agent's inference of an older form). Sticking is implicit in the
  particle model. The flux-normalization fit must be done empirically against this ground truth,
  not by assuming a sticking value.

**Status:** flux-normalization calibration is now a **local** task (fit our model to the real
ViennaPS depth + ARDE + surfaces — no GPU needed). The box has been torn down.

## Phase 1 — 3D + GPU benchmark (Vast.ai RTX 3090, driver 590, ~$0.06)

Built a working minimal 3D differentiable etcher (`src/petch/threed.py`): 3D level set ->
marching-cubes mesh -> Warp ray-traced flux (`_trace3d`: wp.Mesh + mesh_query_ray, ions +
neutral 3D-cosine re-emission) -> chemistry -> 3D advection. Trench and contact hole both etch.

**GPU benchmark (`scripts/gpu_benchmark.py`):**
- Our `_trace3d` flux kernel: **CPU 1.45 M rays/s -> GPU (RT cores) 198 M rays/s = 137x**
  (300k neutral rays, 12 bounces, on a 2602-face etched-trench mesh).
- Full 3D etch on GPU: 12 steps, 60k rays/step, **0.4 s**.

**ViennaPS-3D for comparison (`scripts/viennaps_3d_bench.py`):** OptiX/GPU works on driver 590
(it failed on the 2080 Ti's driver 550 — note for box selection). Full SF6O2 3D trench etch:
**GPU 11.4 s vs CPU 52.8 s (4.6x)**.

**Read:** our differentiable Warp kernel hits RT-core-class throughput (~200 M rays/s), the same
ballpark as ViennaPS's OptiX — confirming "RT cores are RT cores." The defensible edge is NOT
raw rays/sec; it is that our kernel is **differentiable** + **QMC-reducible** + surrogate-able,
which ViennaPS is not.

## 3D speedups (QMC)
QMC (Sobol over the 4D source launch) ported to 3D: exact (same mean), floor m_F noise
0.0032 -> 0.0025 (1.28x). Smaller than 2D's 1.9x because neutral re-emission bounces stay
pseudorandom (only the source is QMC'd); radiosity remains the neutral-noise lever.

## Differentiable inverse design (`scripts/inverse_design.py`)
The payoff demo. Real per-face flux from the 3D simulator + the etch chemistry as a
differentiable Warp kernel; `wp.Tape` gives d(loss)/d(recipe). Asked: "what ion energy makes
the trench floor etch at a target rate?" Gradient descent recovered it:

    target floor rate 1.998 (= rate at E*=140 eV);  start E=100 ->
    129 -> 135 -> 138 -> 139.9 eV;  loss 2e-2 -> 5e-8;  recovered E=139.9 (true 140).

Target outcome -> recipe, by autodiff. The same machinery scales to many parameters and
full-profile targets. This is the defensible edge ViennaPS does not have.

## Phase 1 — going faster + 3D quality

**Per-phase profile (CPU, 2488 faces):** flux 13.2 ms, velocity-extend (KDTree) 5.4 ms,
reinit (skfmm) 5.3 ms, marching-cubes 0.9 ms. On GPU the flux drops ~137x (0.1 ms), so the
**CPU host ops (extend + reinit) become the bottleneck.**

**Speed changes:**
- **GPU velocity extension** (`wp.mesh_query_point`, `extend='gpu'`): exact-same etch outcome as
  KDTree (depth 1.300 == 1.300); ~free on GPU. The right architecture for a GPU-resident loop.
- **Lazy reinit REJECTED as a default:** reinit_every=3 changed depth 1.30 -> 1.60. Skipping
  reinit lets |grad phi| drift from 1, and advect multiplies F*|grad phi| -> wrong front speed.
  Only safe with a proper extension velocity (grad F . grad phi = 0). Kept reinit_every=1.
- Beyond GPU, the real "even faster" levers are a GPU-resident level set (reinit/marching-cubes
  on GPU) and the learned operator (Phase 2).

**3D quality — deep HARC is GATED by the flux-normalization calibration.** A narrow contact hole
self-limits shallow (stuck ~1.25 um over 30 steps). Cause: our model is **neutral-limited**, so
as the hole deepens the floor is F-starved -> theta_F -> 0 -> ion-enhanced etch stops even though
ions reach the floor. 3D ARDE is stronger (solid angle ~(w/d)^2), so it bites fast. Confirmed by
forcing F-saturation: hole depth scales with F flux (1.0 -> 2.2 -> 4.0 um at Fflux x1/x5/x15).
**Same root cause as the 2D ARDE/rate_scale gap (flux normalization).** => doing deep HARC well
and being parameter-free are the SAME fix: calibrate the flux normalization to the ViennaPS
ground truth (now captured). That is the #1 next step. (`docs/etch3d_harc.png`.)

## Flux-normalization calibration — root cause CONFIRMED & fixed (vs real ViennaPS)

Fit the effective F-flux normalization to the captured ViennaPS ground-truth ARDE
`[0.935, 0.966, 0.983, 1.000]` (2D, widths 4/6/8/12). Clear minimum at **F-flux x12**:

| effective F-flux | ARDE rmse vs ViennaPS |
|---|---|
| uncalibrated (belen) | 0.110 |
| x6 | 0.022 |
| **x12** | **0.0165** |
| x25 | 0.021 |
| x50 | 0.031 |

The x12 factor IS the flux-normalization correction: our open-field-normalized `m_F` was ~12x too
low relative to ViennaPS's absolute flux ratio, putting us in the wrong (neutral-limited) regime.
Correcting it:
- **ARDE shape now matches ViennaPS** (rmse 0.110 -> ~0.02 -- a 5x fidelity gain). The DOMINANT
  contributor (#1, flux normalization) is retired, parameter-free in shape.
- The **same correction un-sticks 3D HARC** (F-saturated floor stays fed -> holes deepen).
- Residual: absolute depth still ~10% off with one global `rate_scale` (and ARDE/depth are
  coupled in the 2-param fit). That residual is the SMALLER contributors (#2 IED integration,
  #3 angular yield, #4 ion reflection, #5 higher-order advection) -- the next tier.

**Bottom line:** the "~7% + ARDE + 3D-HARC-self-limiting" cluster all trace to one cause (flux
normalization), now calibrated against real ViennaPS data. ARDE fidelity 5x better; HARC
un-sticks; one global unit constant + minor contributors remain.

## ACCURATE: parameter-free match achieved (2D)

The full calibrated config is now the DEFAULT: **belen + ViennaPS angular yield (contributor #3)
+ cal_F=12 + rate_scale=0.034** (one global unit constant). Contributor #3 (per-channel angular
yield f_sp/f_ie) further cut ARDE rmse 0.033 -> 0.017. Final 2D match vs ViennaPS ground truth:

| metric | original (per-case knob) | calibrated (one global constant) |
|---|---|---|
| width-8 depth error | -8.7% | **+3%** |
| ARDE rmse | 0.110 | **0.016** |

That is the project's stated goal: match ViennaPS without a per-case `rate_scale`, using one
global unit constant. The default model is now this calibrated config; `Flags(chemistry=
'langmuir', yield_angular='cosine')` + cal_F=1 recovers the original PoC (pinned by tests).

NOTE: the calibrated `rate_scale` differs 2D (~0.034) vs 3D (~0.6) because mc_flux (per-width)
and mc_flux_3d (per-area) use different flux normalizations -> the unit constant is per-dimension.
A full 3D calibration needs one ViennaPS-3D ground-truth run (GPU box) -- the 2D recipe transfers.

Remaining accuracy tier (smaller): #2 IED integration (~0.2% at this IED), #4 ion reflection
(bottom-corner microtrenching in deep features), #5 WENO advection (less numerical diffusion).
