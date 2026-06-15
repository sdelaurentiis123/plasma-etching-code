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

## 3D: contributor #4 + ViennaPS-3D ground truth + the honest 3D fidelity gap

- **Contributor #4 (ion specular reflection)** implemented in the 3D kernel (grazing reflect,
  energy retention) -> feeds the bottom corners; deepens a HARC hole 1.40 -> 1.70 um. Default off.
- **ViennaPS-3D ground truth captured** (Vast.ai A4000, `harness/reference/viennaps_3d_groundtruth
  .json`): SF6O2 contact HOLES depth 7.68/8.47/9.24 um at d=3/4/6 -> norm ARDE [0.832, 0.916, 1.0]
  (genuine 3D, steeper than the trench ARDE); trenches 9.29..10.09 at w=3..8.
- **Honest 3D fidelity gap.** Calibrating our 3D holes to that: cal_F~40 best fits the ARDE
  *shape at shallow depth* (rmse 0.08), but at ViennaPS-comparable DEPTH our narrow holes over-lag
  badly (d=3 norm 0.49 vs ViennaPS 0.83, rmse 0.26). **No single cal_F matches across depths** ->
  our neutral transport over-depletes deep narrow 3D holes more than ViennaPS (3D ARDE ~(w/d)^2
  bites harder). This is a model-form limit, not a constant: ViennaPS stays F-saturated at the
  deep narrow floor; our minimal MC neutral transport does not. Closing it needs a better neutral
  model (F-saturation / radiosity-with-correct-absolute-flux), not just calibration.

**Net 3D state:** the 3D etcher works (holes+trenches, GPU 137x, differentiable, un-sticks with
cal_F) and qualitatively does HARC; but **deep-HARC ARDE does not yet quantitatively match
ViennaPS** -- a named, real gap in the neutral transport regime. 2D is parameter-free matched;
3D needs the neutral-transport upgrade for deep features.

## Closing the 3D gap — investigation (web research + ViennaPS source)

**Found the missing mechanism.** Reading the ViennaPS source (`psPlasmaEtching.hpp`): the etchant
sticking is **coverage-dependent**, `S_eff = (1 - coverage) * beta` (Langmuir: radicals stick only
on BARE sites), and flux is recorded as ARRIVING (every hit), not stuck. On the saturated upper
sidewalls (high coverage -> low bare -> low sticking) radicals reflect and **penetrate deeper** to
the under-fed floor. Our model used CONSTANT sticking (0.7) -> radicals stick on the upper walls ->
floor starved. That is a real model deficiency we were missing.

**Implemented it** (`mc_flux_3d_coupled` + `_trace3d_cov` + a flux<->coverage fixed point, flag
`coverage_sticking`). It works mechanically: the etch goes **deeper** (more F reaches the floor).

**But it does NOT close the ARDE ratio.** Matched against ViennaPS-3D holes (d3/d6 depth ratio
0.832): both constant-sticking and coverage-sticking give **~0.62** -- and we run at FINER
resolution than ViennaPS (dx 0.25 vs 0.3), so it is not a grid artifact. Our narrow holes are
*steeper than the geometric shadowing ratio* (~0.74) while ViennaPS is *flatter*. So the dominant
cause is **not** sticking/re-emission/resolution.

**Narrowed residual:** the floor coverage balance. At the hole floor, vertical ions keep removing F
(ion-enhanced etch consumes the fluorinated layer), so theta_F stays low and the etch is amplified-
below-geometric; ViennaPS keeps the floor more F-saturated. Closing it needs the floor coverage/ion-
flux balance reconciled to ViennaPS (likely the absolute ion-vs-etchant flux ratio at the floor, or
a coverage-dependent ion yield), not just the neutral sticking. `coverage_sticking` is shipped
(default off, it is slower + changes normalization) as the correct first piece. **Honest: the 3D
deep-HARC ARDE gap is narrowed but still open.**

## Clean GPU rate-matched validation of coverage-sticking (A4000, `scripts/gpu_3d_validate.py`)

The CPU "~0.62 both, doesn't close it" above was **confounded** (compared at non-matched depth +
CPU MC noise). A clean GPU run that **rate-matches both configs to the same d6 hole depth** and
sweeps rate_scale tells a clearer story. Result (`harness/reference/gpu_3d_validate_result.json`):

| config | normalized ARDE [d3,d4,d6] | rmse vs ViennaPS [0.832,0.916,1.0] |
|---|---|---|
| coverage **OFF** (rate 0.50, d6=3.0) | [0.50, 0.667, 1.0] | **0.2397** |
| coverage **ON**  (rate 0.09, d6=2.75) | [0.727, 0.909, 1.0] | **0.0603** |

**Coverage-dependent sticking cuts the 3D ARDE shape error ~4x** (0.24 -> 0.06) at matched depth;
the d4/d6 ratio is near-exact (0.909 vs ViennaPS 0.916). It also etches **deeper at equal rate**
(more F reaches the floor: at rate 0.09, cov-ON d6=2.75 vs cov-OFF 1.5; cov-ON reaches 6.75 um by
rate 0.24 while cov-OFF plateaus ~3 um). So coverage-sticking IS the dominant missing 3D mechanism,
confirmed cleanly — supersedes the earlier "doesn't close it" once depth is matched.

**Two honest residuals remain.** (1) Matched depth here is ~3 um, not ViennaPS's 9.24 um — our
floor still starves before full HARC depth (cov-OFF plateaus ~3 um regardless of rate; the max
velocity stays high but sits on the re-pinned mask/sidewall, not the floor). (2) As cov-ON holes
deepen (rate 0.24 -> d6=6.75 um) the smallest hole over-lags again (d3/d6 0.727 -> 0.481), so at
full aspect ratio the gap reduces but does not vanish. Depth quantization (dx=0.25 -> +-1 cell)
also makes individual ratios noisy. Net: **coverage-sticking is validated as the right fix and
closes most of the gap; the residual is small-feature floor flux at the highest aspect ratios.**

### Production GPU-reinit NaN bug found (follow-up, not yet fixed)
While validating, the GPU Russo-Smereka reinit (`reinit_method="gpu"`) was found to develop a **NaN
instability in deep holes** (~step 16-20): phi corrupts, the `nan_to_num` velocity guard then zeroes
V (Vmax -> 0.000) and the depth runs negative. Isolated cleanly: **`reinit_method="skfmm"` (CPU) is
rock-stable** on the identical run (Vmax steady ~8, depth grows monotonically), so the bug is in the
GPU reinit kernel, not the flux/advection. The flux kernel (the expensive part) still runs on GPU;
skfmm reinit cost is negligible for this grid. The validation therefore uses skfmm. The earlier
"GPU reinit Russo-Smereka depth-exact" result held at shallower depth / fewer steps; it does NOT
survive deep-hole HARC runs. **Action item: debug the GPU reinit subcell/Godunov kernel for the
divide-by-zero / sqrt-of-negative that produces NaN at high curvature before using it in production.**

A defensive guard was added to `run_etch_3d` (`V = np.nan_to_num(...)`, finite-Vmax floor) so a NaN
can't propagate into the CFL substep count — but it masks rather than fixes the GPU-reinit root cause.

### GPU-reinit NaN crash — FIXED; a separate |grad| bias found and characterized
Reproduced the NaN locally (Warp runs on the M1 CPU — no GPU box needed) and fixed the crash:
1. **Russo-Smereka `D = s0/grad0` blowup** — a near-flat cell mis-flagged as an interface cell (deep-
   hole corner / re-pinned mask edge) drove `grad0 -> 1e-9` so `D -> inf -> NaN`. Fix: floor `grad0`
   at 0.5 and **clamp `D` to +-1.8*dx** (an interface cell is within a cell-diagonal of the contour).
2. **3D CFL** — the Godunov sweep is forward-Euler; `dtau=0.5*dx` is borderline in 3D (stable shallow,
   grows under stress). Lowered to `dtau=0.3*dx` (~dx/sqrt(3) margin).
Result: **no more NaN/blowup** (verified `scripts/repro_reinit_nan.py`: gpu now finishes with
`any_nan=False`, depth no longer runs negative). The reinit is also high-fidelity on smooth fronts —
`scripts/check_reinit_fidelity.py` (A): on a perturbed sphere it recovers the EXACT distance better
than skfmm (max|err| 0.036 vs 0.064); (B): on a real developed front it agrees with skfmm to mean
0.024 um with an **identical zero-contour**.

**BUT a separate, real bias remains on MASKED fronts:** measuring `|grad phi|` (upwind, mask cells
excluded) after one reinit of a real hole front gives gpu **mean 1.32 / std 0.43** vs skfmm 1.10 /
0.16 — and this is **iteration-converged** (identical at n_iter 24/48/96/160), so it is a wrong fixed
point, not under-convergence. Cause: the Russo-Smereka interface *freeze* (`out -> phi0/|grad phi0|`)
locks adjacent interface cells at mutually-inconsistent sub-cell distances when the input is not a
clean SDF — and the re-pinned mask boundary (`phi[mask]=mask_phi`) injects exactly such a kink each
step. The band then inherits `|grad|>1`, and since `advect = F*|grad phi|`, the front **over-moves ->
the full loop etches ~1-2 cells deeper than skfmm** (e.g. 2.0 vs 1.5 um wide-trench; 5.0 vs 2.75 um
on a stiff cov-ON hole where it compounds). skfmm (exact global fast-marching) does not have this.

**Status:** the crash is fixed and committed; **skfmm stays the trusted default**. For a fully
GPU-resident *and depth-accurate* loop the GPU reinit needs a true `|grad|=1` Eikonal solve on masked
geometry. Two candidate next steps: **(a)** a parallel **GPU fast-sweeping (FSM)** Eikonal solver
(gives `|grad|=1` like skfmm), or **(b)** the architecturally cleaner route — a proper **extension
velocity** (`grad F . grad phi = 0`, constant along normals) so `|grad phi|` stays ~1 under advection
and reinit becomes infrequent/amortized (this also retires the earlier "lazy reinit drifts" problem,
which was caused by the *non*-proper extension velocity). (b) is preferred: it removes per-step reinit
from the critical path instead of just making a subtly-biased reinit faster.

## Literature grounding — authoritative SF6/O2 physics (multi-agent dig, cited)

A deep primary-literature sweep (Belen/Ertl, Kushner, Graves, Steinbruchel/Gray, Gottscho, Coburn &
Winters, Donnelly, Flamm lineage; new refs added to `docs/references.bib` section E2) to (i) confirm
our constants are physical and (ii) find where ViennaPS approximates so we can exceed it. Headlines:

**Our model is physically sound — constants validated, approach is field-standard.**
- Yields (`A_ie=7, Eth_ie=15, A_sp=0.0337, Eth_sp=20, B_sp=9.3`), `k_sigma=300` (= Belen k*sigma_Si
  = 3e17 cm^-2 s^-1, *identical*), `beta_sigma=0.04`, coupled-coverage forms, and the `k_sigma*theta_F/4`
  chemical term all match the Belen 2005 SF6/O2 two-coverage model (JVST A 23(5), 1430; DOI
  10.1116/1.2013317) that ViennaPS implements. Our `betaE=0.7` IS correct for the SF6/O2 path (ViennaPS
  uses gamma_F=0.7, gamma_O=1.0; the generic 1.0 default is overridden there too).
- **F sticking is a *fitted/flux-dependent* parameter everywhere** (Belen 0.7, Marcos 0.1, Donnelly
  2017 flux curve 0.001-0.03, MD 0.98->0.23 as Si fluorinates). So our `cal_F` flux-normalization knob
  is the field-standard treatment, not a hack.
- **Our 3D-hole-steeper-than-2D-trench ARDE is CORRECT PHYSICS, not a bug** (Gottscho 1992, JVST B 10,
  2133): pure-AR scaling is hole R_bottom ~ (w/d)^2 vs trench ~ (w/d) — a hole confines the neutral
  acceptance solid angle in both lateral directions. Our holes [0.832,0.916,1.0] vs trenches [0.935,...]
  is exactly this.

**The ARDE master equation (Coburn & Winters 1989, APL 55, 2730; DOI 10.1063/1.101937):**
`R_bottom/R_top = K / (K + S - K*S)`, K = Clausing transmission (0.11 at depth/diameter=10), S = floor
reaction probability. **Neutral *delivery* limits, not product removal.** Crucially **as S->0,
R_b/R_t -> 1/K** (lower floor sticking => FLATTER ARDE). This is the quantitative form behind our
"neutral-limited (steep ARDE) vs F-saturated (flat ARDE)" diagnosis: ViennaPS's flat deep-ARDE is the
F-saturated regime, and our `cal_F=12` fix pushed us toward it (rmse 0.110->0.016). The 3D deep-HARC
parity job is therefore a **flux-normalization / floor-saturation** match, not new physics.

**Strategic reframe of "more accurate than ViennaPS" (decisive):** the most likely cause of the
*residual* deep-HARC floor over-lag is **surface charging** — differential in-feature potentials cut
the floor ion current ~60% by AR 4 (Hwang & Giapis 1997) — and **ViennaPS OMITS charging entirely**
(also omits volatile-product redeposition; Hoekstra & Kushner 1998 / Huard & Kushner 2017). Two
consequences:
1. **Matching ViennaPS-3D needs NO new physics** — it is the F-saturation/flux-normalization match above.
   Adding charging would *diverge* us from ViennaPS (our only non-experimental reference).
2. **Exceeding ViennaPS** = adding charging + redeposition + full-IEDF integration — but that is only
   *demonstrable* against EXPERIMENTAL wafer data, which we do not have. So per the agreed plan:
   match ViennaPS to parity + be faster now; stage the charging/redeposition physics as the
   "exceed once we have wafer data" track.

**Where a more rigorous model beats ViennaPS (ranked, for the exceed track):**
(1) surface charging (floor ion-current drop; fully omitted) — highest value/effort;
(2) volatile etch-product redeposition (narrows top, raises effective AR; omitted);
(3) flux/coverage-dependent F sticking (Donnelly 2017 curve vs constant gamma_F) — attacks floor
    starvation, cheap, keeps differentiability;
(4) full bimodal IEDF + energy-angle correlation vs mean+sigma Gaussian (only ~0.2-7% on total yield,
    but sharply nonlinear near the sqrt(E) threshold and for sidewall-vs-floor selectivity);
(5) our existing structural edge ViennaPS lacks: differentiability + QMC + (later) learned surrogate.

**Immediate next accuracy step (grounded):** reconcile the 3D absolute flux / floor-saturation regime
to match ViennaPS-3D deep-HARC, using the Coburn-Winters `K/(K+S-KS)` form as the target relation —
the physical version of extending `cal_F` to 3D, validated on a GPU box (deep holes, fresh ViennaPS-3D
surfaces). Coverage-dependent sticking (already validated to cut ARDE error ~4x) becomes the 3D default.

## GPU campaign — depth-resolved ViennaPS-3D parity (RTX 3090, scripts/*_3d_*.py + harness/reference)

Ran the full campaign on a fresh box: depth-RESOLVED ViennaPS-3D ground truth + a calibration sweep of
the two candidate levers (`cal_F`, then `betaE`) + a head-to-head speed benchmark. Results below;
ViennaPS used the CPU_TRIANGLE engine (OptiX/`libnvoptix.so` absent in the container — its `gpuAvailable()`
probe FATALLY hangs/crashes, so never call it; force the CPU engine).

**ViennaPS-3D depth-resolved ground truth** (`harness/reference/viennaps_3d_depth_resolved.json`,
holes d=3/4/6, durations 1/2/3/4.5 min). The ARDE steepens monotonically with depth exactly as Gottscho
1992 predicts (a hole confines neutral acceptance in both lateral dims):

| dur | d6 (um) | normalized ARDE [d3/d6, d4/d6, 1] |
|---|---|---|
| 1.0 | 3.29 | [0.899, 0.951, 1.0] |
| 2.0 | 6.37 | [0.862, 0.933, 1.0] |
| 3.0 | 9.24 | [0.831, 0.917, 1.0]  (reproduces the cached ground truth to 3 digits) |
| 4.5 | 13.21 | [0.793, 0.893, 1.0] |

**SPEED — we beat ViennaPS by ~23x (decisive WIN).** Same box, matched deep hole (d6~13 um): our
differentiable Warp etcher (GPU) **15.0 s** vs ViennaPS (CPU_TRIANGLE) **351.8 s** -> **23.4x faster** —
and that is with our reinit still on CPU (skfmm). At shallower depths the lead is ~10-16x. (Caveat:
ViennaPS GPU/OptiX was unavailable on this box, so it is our-GPU vs ViennaPS-CPU, not same-engine.)

**ACCURACY — found the right lever; substantial but not full parity (honest).**
- `rate_scale` does NOT fix deep-HARC: our d6 PLATEAUS at ~7 um regardless of rate -> hard floor
  starvation. (`petch_3d_match_result.json`, uncalibrated.)
- **`cal_F` does NOT fix it either** (`cal_F_sweep_3d_result.json`): swept 12->130 (10x the 2D fix);
  d6 stays ~7-7.5 um, ARDE stays steep. **The 2D flux-normalization fix does NOT transfer to 3D.**
  Reason: the F flux *reaching* the conductance-shadowed deep floor is ~0, so `Gb_E = Fflux*m_F*cal_F`
  ~ 0 * cal_F is still ~0. This REFUTES the earlier optimistic "matching ViennaPS needs only
  flux-normalization" reframe — see corrected [[accuracy-yardstick-vs-viennaps]].
- **`betaE` (F sticking) IS the lever** (`beta_sweep_3d_result.json`), as Coburn-Winters predicts
  (lower floor sticking S -> radicals reflect off saturated walls, penetrate deeper -> `R_b/R_t -> 1/K`,
  floor keeps etching, ARDE flattens). Lowering betaE 0.7 -> 0.08: max reachable d6 **7.0 -> 13.0 um**
  (un-starves; now matches/exceeds ViennaPS's deepest), and ARDE rmse-vs-VPS **0.207 -> 0.144**. betaE
  ~0.08 is the optimum (0.03 over-shoots). 0.08 sits within the literature F-sticking range (Donnelly
  2017 high-flux 0.001-0.03; feature-scale fits 0.1-0.7) — defensible, and physically sensible since
  deep 3D HARC needs more wall reflection to feed the floor.
- **Calibrated trajectory** (`petch_3d_match_calibrated_result.json`, betaE=0.08): our holes now reach
  14 um (ViennaPS 13.21) and track ViennaPS ARDE to **rmse ~0.11-0.12 across all four depths** (down from
  ~0.21). BUT a residual remains — the smallest, highest-AR hole over-lags: at d6=9.24 ours d3/d6=0.65 vs
  ViennaPS 0.831; at d6=13.21 ours 0.615 vs 0.793. Even at maximal wall reflection (betaE->0.03) our d3
  starves more than ViennaPS.

**Verdict.** Calibrated to a physically-grounded `betaE~0.08`, our 3D model now (i) **reaches full
deep-HARC depth** (floor starvation solved) and (ii) **matches ViennaPS-3D ARDE to rmse ~0.12 across the
depth trajectory** — close, parameter-light, but NOT exact. The residual d3 over-lag at the highest
aspect ratios is a genuine **transport model-form limit** (our bounce re-emission vs true Knudsen
molecular flow; plus the charging ViennaPS also omits) — exactly the physics flagged in the literature
dig as the "exceed" track. So: **we are decisively faster (~23x) and now close on 3D accuracy; closing
the last d3 residual needs the Knudsen-transport / charging upgrade, not another calibration knob.**

NOTE: betaE=0.08 is a 3D-HARC-specific effective sticking; the global default stays betaE=0.7 (correct
for SF6/O2 and pinned by the 2D tests). The campaign scripts set betaE per-run; do not change PAR globally.
