# Autonomous progress log (single writer)

Started 2026-07-12, continued from the Codex reconciliation handoff. Branch:
`codex/unified-engine-root-fixes`. **Local commits only, no push** (per user). One active writer.

Goal restated by the user: reach the OSTI/Oehrlein roadmap targets with a self-consistent engine
that **matches correct aspect ratios in a correct first-principles physics way** — no shortcuts, no
AR-shaped fudge; adaptive mesh/phase-space refinement (AMR); fast + GPU-accelerated at the end.

## Verified state (independently audited, not trusted from docs)

- **Branch is sound — KEEP it.** Forensic audit (Fable, read-only): 128 commits, one linear 24h
  chain, zero deletions, no concurrent-writer scars, docs *under*-claim. Rollback only re-introduces a
  known geometry bug and deletes real tested work. `pytest -q` = **296 passed, 1 skipped**.
- **The green tests are real conservation/analytic gates but contain ZERO aspect-ratio validation.**
  The one "aspect ratio ladder" test (`test_boundary_transport.py:89`) fires a *vertical monodirectional*
  beam into a straight trench and asserts floor flux = 1.0 at AR 1/4/16 — true by geometry, cannot show
  ARDE. All de Boer/ARDE work is loose scripts against the *legacy* monolith and was de-earned in the
  reconciliation. **"Matches aspect ratios first-principles" is unbuilt in the common engine.**

## Direction

1. **de Boer SF6/O2 ARDE through `feature-3d`** (radical-transport-limited → isolates AR transport from
   chemistry). First: a static neutral-flux-vs-AR gate that shows the floor-flux collapse and proves it
   CONVERGES under refinement. Then the rate closure with radical sticking as a *declared* input;
   calibrate low-AR, predict held-out AR40.
2. Jeon SiO2 depth-transfer (harness already built). 3. Charging only if it moves the profile above the
   error budget. (OSTI_OEHRLEIN_PROBLEM_MAP.md §6 ordering, adjusted so the pure-transport gate leads.)

## Work in progress

- **`scripts/deboer_arde_static.py`** (new): static ARDE floor-flux gate through the common engine —
  half-Maxwellian *flux* neutral source (analytic cosine law, no fitted angular closure) → ballistic
  face-gather (`gather_boundary_state_ballistic_3d`) → diffuse molecular-flow radiosity
  (`solve_diffuse_neutral_radiosity_3d`, one physical sticking `s`). Conservation exact (balance ~1e-13).

- **CONFIRMED: a FIXED angular quadrature flatlines the floor flux at high AR; the AMR requirement is
  real.** At s=1 (pure line-of-sight shadowing) a fixed 5-node quadrature flatlines transmission ~0.53
  for AR≥2 (~17x too high at AR16). Root cause, literature-confirmed (Coburn-Winters; JVST A 35
  05C301): the floor-reaching acceptance cone is ~arctan(1/A), so a fixed N-node quadrature aliases it
  once A≳N; angular samples must scale ∝A.

- **GPU-ready high-AR path found: QMC-sampled source + forward first-hit tracer.** The forward tracer
  batches every angular atom x source position into ONE Warp kernel; sampling the flux density with
  N=2^L Sobol points (`thermal_neutral_qmc_boundary_state`) concentrates rays by the physical cosine
  measure so the acceptance cone is resolved by raising N. It is **N-converged and fast** (AR16 in ~1s,
  stable across N=2^14..2^18): T(AR8)=0.065, T(AR16)=0.0355.

- **VALIDATED (three independent ways).** The common engine's ballistic neutral transport reproduces
  first-principles geometric shadowing (s=1). The two open items from the earlier correction are both
  closed:
  1. **Method reconciliation:** the adjoint gather CLIMBS to the forward+QMC value as it resolves
     (AR1: nt 24->48->72 gives 0.302->0.342->0.348 vs forward+QMC 0.345, ~1%). The earlier "8-14%
     disagreement" was the adjoint being angular-under-resolved (nt=24); properly refined, the two
     independent estimators agree.
  2. **Exact reference:** an INDEPENDENT pure-numpy analytic ray-trace of the same box
     (`reference_floor_transmission`, no engine code) gives 0.303/0.193/0.110/0.059 at AR1/2/4/8 —
     matching the opposed-strip view factor. The engine CONVERGES to it under GRID refinement
     (AR1: dx 0.02->0.0125->0.008 gives 0.345->0.295->0.310 vs ref 0.303). The coarse-grid dx=0.02
     over-prediction (~14% at AR1, ~3% at AR8; worse at low AR because 5 cells/opening staircases the
     walls) is grid error, not an engine bias.
  So getting a correct ARDE number requires BOTH adequate grid (>=~12 cells/opening) AND angular AMR
  (QMC-refined or nt ∝ A; a fixed quadrature flatlines/over-predicts the deep cone). Both are
  error-driven, no shortcuts. The earlier committed "premature 2% at AR1" (405e555) was the coincidence
  of nt=24 + dx=0.02 landing near the approximate analytic; corrected in 9eafb1f and closed here.

- **Committed gate:** `tests/test_arde_transport.py` — the first ARDE physics test in the suite:
  monotone ARDE collapse, engine == independent reference within grid tolerance, and the AMR-necessity
  regression (a fixed coarse quadrature over-predicts the deep flux; QMC does not).
  Analytic targets + citations: `ARDE_PHYSICS_REFERENCE.md`.

## Reactive s<1 family — VALIDATED vs independent particle Monte Carlo

The reactive case (walls+floor re-emit diffusely, reflection 1-s) is now validated by an INDEPENDENT
method, not just demonstrated. `scripts/arde_mc_reference.py::mc_reactive_transmission` is a stochastic
particle MC of the same box (region-aware straight-line tracing; react w.p. s or diffuse-reflect;
tally floor incidence). It self-checks to the s=1 geometric reference exactly (0.303/0.193/0.110 at
AR1/2/4), then across s the engine's deterministic radiosity agrees to ~1-3%:

| AR | s | particle MC | engine radiosity | ratio |
|----|-----|-------------|------------------|-------|
| 1 | 0.10 | 0.7666 | 0.7750 | 1.011 |
| 2 | 0.10 | 0.6264 | 0.6374 | 1.018 |
| 4 | 0.10 | 0.4214 | 0.4268 | 1.013 |
| 1 | 0.50 | 0.4185 | 0.4328 | 1.034 |
| 2 | 0.50 | 0.2737 | 0.2800 | 1.023 |
| 4 | 0.50 | 0.1485 | 0.1494 | 1.006 |

Behavior matches Coburn-Winters (monotone in AR, decreasing in s). Gated in
`tests/test_arde_transport.py` (now 8 tests). So the full ARDE TRANSPORT (s=1 geometric AND reactive
s<1) is validated by independent methods.

## Determinism

All estimators are bit-identical across repeated calls (seeded scrambled-Sobol + seeded default_rng;
deterministic Warp ray-mesh hits + numpy accumulation + GMRES): forward+QMC 0.3746781163693154,
adjoint 0.1434127263028523, particle MC 0.36712646484375, numpy ray-trace 0.19273757934570312 — each
identical on re-run. The gate is a stable deterministic regression.

## de Boer SF6/O2 ARDE — two-channel result

Experiment (de Boer/Blauw cryo, from the legacy npz): normalized rate [1.0, 0.43, 0.29, 0.20] at
AR [0, 10, 20, 40].

- **Radical-only (validated transport) is too steep, and this is DEFINITIVELY physics, not numerics.**
  Best radical-only normalized rate collapses to ~0.035 at AR40 (calibrated s~0.06 gives 0.25/0.09/0.035
  at 10/20/40). The MC floor is bit-identical across max_bounce 400..8000 at AR40 -> converged, no
  truncation. This RULES OUT the old "petch-MC under-samples the deep floor" suspicion (memory
  [[reconcile-craig-into-petch]]): the transport is correct Knudsen physics; the de Boer high-AR floor
  is a real second-channel effect.
- **Two-channel (radical + directional ion) reproduces the experiment.** de Boer SF6/O2 is
  ion-assisted; directional ions sustain the high-AR floor. `scripts/deboer_two_channel.py` adds a
  reduced directional-ion channel (Gaussian cross-slot angle, absorbing walls) and fits an additive
  radical+ion rate. Calibrated on the AR10,20 knee (RMSE 0.008), it PREDICTS the held-out AR40 floor:
  model NR [1.0, 0.43, 0.279, 0.169] vs experiment [1.0, 0.43, 0.29, 0.20], held-out AR40 error 0.031.
  Physical params: radical sticking s=0.06, ion IAD sigma~1 deg, ion/radical strength beta~0.4.
- **Caveats / remaining frontier:** the ion channel is a REDUCED analytic model (not yet the full
  common-engine ion transport), and the rate law is an additive assumption. The AR40 residual
  (0.169 vs 0.20) and the near-sub-degree ion IAD (sigma~1 deg) are the frontier the prior work
  flagged (sub-degree IADF and/or charging). Next: run the ion channel THROUGH the validated engine
  transport (narrow IonEnergyTransverseMaxwellianDensity) and, if the residual persists above the
  error budget, add charging.

## Phase 0 (approved plan) — de Boer ARDE from real coupled chemistry through the engine: WORKING

`scripts/deboer_feature3d.py` runs de Boer Si SF6/O2 through `feature_step_3d` with the VALIDATED
coupled kernel `ReducedSiO2FluorocarbonMechanism` reparameterized for Si-F (complex=SiF_x built by F
and removed by ion; polymer=O-passivation gating F access). The ARDE EMERGES from the coupling (not an
additive sum): normalized floor rate 1.0/0.87/0.74/0.53 at AR 0/2/4/8. The etch is genuinely
ion-assisted: ion_flux=0 -> rate exactly 0.0. This replaces the additive `deboer_two_channel.py`.

**Phase 0 COMPLETE — de Boer match through the coupled engine (`--deboer` mode).** Calibrating ONE
physical knob (F sticking s_F=0.08) on the de Boer knee (AR10,20): NR10/NR20 = 0.476/0.289 vs exp
0.43/0.29 (RMSE 0.033), then **predicting the held-out AR40 floor = 0.166 vs exp 0.20 (error 0.034)** —
through real coupled ion+neutral surface chemistry, calibrate-N/predict-N+1 style. Same accuracy as
the hand-built additive model (0.166 vs 0.169) but from the self-consistent engine. **The engine is
fast at high AR: AR40 runs in ~10 s at dx=0.01 — the feared "mesh wall" does not exist.** The AR40
under-prediction (0.166 vs 0.20) is the sub-degree-IAD / charging frontier (Phase 2). TODO: relabeled
Si product class; fast pytest gate; then Phase 1 (the root: differentiable + scalable transport).

## Roadmap (remaining)

1. de Boer: run the directional-ion channel THROUGH the validated engine transport (narrow
   IonEnergyTransverseMaxwellianDensity), replacing the reduced analytic ion model; sticking + ion IAD
   are DECLARED calibrated inputs with provenance/uncertainty; calibrate low-AR, predict held-out AR40;
   report grid/ray/digitization/model error separately. Add charging only if the AR40 residual exceeds
   the combined error budget.
2. GPU: run the forward+QMC path with device="cuda" (already threaded); accuracy-matched speed report.
3. Then Jeon SiO2 depth-transfer, then charging (only if it moves the profile above the error budget).

## Guardrails honored

Single writer. Local only. No benchmark/AR/region branch in governing physics. Sticking is the one
declared physical surface input with provenance. Refinement is error-driven. Device-agnostic transport
(Warp CUDA path ready for the GPU deployment).
