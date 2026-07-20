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

## Phase 1a spike (the root, differentiability) — CALIBRATION gradients are EXACT

`scripts/diff_calibration_gradient.py`. The near-term moat is gradient-based CALIBRATION (calibrate
declared surface params on structure N, predict N+1). Those params (sticking, yields) enter the diffuse
radiosity operator M(s)=I-(1-s)B SMOOTHLY, so d(floor flux)/ds is exact via implicit differentiation of
the linear fixed point -- one adjoint solve, `-(M^{-T}c)^T (B H)` -- matching central FD to ~1e-7:

```
  AR   s    analytic       central_FD      rel_err
   1 0.05  -2.00960e+00   -2.00960e+00    6.1e-08
   4 0.05  -4.72304e+00   -4.72305e+00    6.0e-07
```

Go/no-go: **GO** for the calibration half of the moat -- it is a solved problem (exact adjoint through
the radiosity fixed point; the same IFT structure will wrap the nonlinear charging fixed point later).
This cleanly separates the moat: CHEMISTRY-parameter (calibration) gradients are exact and cheap; only
GEOMETRY/shape gradients hit the discontinuous ray-hit boundary (the docs' "differentiability open" is
really this half -- inverse SHAPE design). Gated in `tests/test_arde_transport.py`.

## Phase 1 moat DEMO — data-efficient gradient calibration (`scripts/diff_calibration_demo.py`)

Recovering a K-band sticking map from floor fluxes, adjoint vs finite-difference L-BFGS-B:
```
   K  adjoint_solves  FD_solves  ratio   adj_err   fd_err
   2        81           114      1.4    3e-15    1.3e-5
   6       210           462      2.2    1e-13    3.4e-5
  12       210           897      4.3    8e-2     8e-2
```
Adjoint solve-count is flat in K (the full gradient is ONE solve); finite-difference grows O(K), so
the wedge ratio climbs 1.4->4.3. The exact adjoint recovers to machine precision (3e-15) vs FD 1e-5.
K=12 hits a physical IDENTIFIABILITY limit (12 depth bands from floor-only data at one AR is
under-determined) -- an honest result that quantifies how many structures must be observed.

## Self-consistency status (honest)

- Neutral radiosity re-emission IS a self-consistent (linear) fixed point (converged, conserved).
- **Self-consistent CHARGING loop is built, wired, AND tested**: `solve_dielectric_charging_steady_3d`
  converges the nonlinear current balance (I+ = I- per node to rtol 1e-14) to a self-consistent surface
  potential and reuses the field-converged ion events for chemistry; gate
  `tests/test_feature_step_3d.py::test_feature_step_solves_charge...` passes in the green suite.
- **DEMONSTRATED (`scripts/charging_selfconsistent_demo.py`):** the loop converges in 8 iterations
  (3 trust-region rejections) to a current-balance residual of EXACTLY 0.0 on all 32 active nodes
  (I+ = I- = floating condition) and a self-consistent surface potential -10.77..0 V (dielectric charges
  negative under net-electron arrival). So the self-consistent engine works, concretely. Known gap: a
  periodic-cell TRENCH trips a float32 cell-boundary tolerance in `lump_triangle_sheet_charge_3d` (verts
  exactly on the cell edge) -- a Phase-2 robustness fix, not physics. Phase 2 next: fix that, connect
  charging to notching + the deep-AR residual.

## Phase 2 robustness follow-up (2026-07-12)

- **Float32 endpoint projection is fixed and regression-gated.**
  `lump_triangle_sheet_charge_3d` now admits only a source-precision-scaled roundoff band at the nodal
  grid boundary and clamps admitted normalized coordinates to the exact endpoint.  A marching-cubes
  `float32(0.3)` endpoint previously normalized to `30.00000119` on a `dx=0.01` grid and was rejected;
  `tests/test_charging_poisson_3d.py` now reproduces that case while retaining charge conservation.
- **A real trench now passes sheet-charge projection, but production periodic-feature charging is not
  yet earned.**  The next failure is physical boundary handling, not the projection tolerance:
  `trace_boundary_state_field_3d` treats every lateral grid crossing as escape and exposes no periodic
  charged-trajectory option, unlike the periodic ballistic/radiosity paths.  Wide-angle half-Maxwellian
  electrons therefore leak from a periodic unit cell and the local current-balance solve does not
  converge.  Do not promote a trench charging gate by loosening its residual.  The immediate bounded
  task is a periodic, fixed-step charged trajectory map with crossing-segment hit handling.  For the
  centered translationally invariant slot, the existing natural lateral Poisson condition is the exact
  symmetry condition; general asymmetric periodic features still require explicit periodic nodal
  identification and retain that limitation.

### Phase 2b engine work after the projection fix

- The common fixed-step nodal-field tracer now wraps curved charged trajectories across lateral cell
  boundaries and checks every split crossing segment for a surface hit.  A zero-field oblique-ion gate
  proves open-cell escape becomes unit-probability periodic landing without changing kinetic energy.
  The option is threaded through physical charging and the public feature-step/solve APIs.
- The charging update now inverts the exact dense support-node response of the sparse Q1 Poisson
  factorization.  A manufactured gate requests four coupled surface-voltage increments and reproduces
  all four to ~3e-13; the old independent diagonal-capacitance approximation could not do this in a
  trench.  The common solver also exposes the already-proven type-II Anderson accelerator from the
  nodal lineage; it changes only nonlinear convergence, not the current-balance root.
- **The real-trench gate is still not promoted.**  With the physical continuous ion/electron densities,
  the total-current estimator is converged enough to expose the next architectural gap: forward QMC
  gives sparse/noisy *local* triangle currents, and the nonlinear solve stalls well above the per-node
  balance tolerance even though global current is resolved.  Finite-difference Newton was explicitly
  not retained: the earlier fixed-map campaign already showed the hard hit/escape Jacobian is ill
  conditioned.  Next is the known common-engine requirement, not parameter tuning: port the frozen,
  reversible, bidirectional per-surface current estimator to arbitrary 3-D triangles, then re-run the
  same trench gate with estimator uncertainty inside the current-balance budget.

### Phase 2c reversible 3-D current gather

- The arbitrary-triangle charged transport backend now has a reversible Liouville adjoint gather. It
  launches the declared numerical velocity proposal from triangle quadrature points, reverses the same
  fixed-step nodal Hamiltonian map used by forward transport, and scores the physical plasma-boundary
  density with the exact normal-velocity Jacobian. Proposal density changes variance only. Periodic and
  nonperiodic lateral source domains are classified explicitly.
- Independent gates establish the estimator measure before it enters the nonlinear solver: a flat
  zero-field Maxwellian returns unit landing and the exact `2*T_e` flux-weighted energy; a linear
  electron barrier returns `exp(-DeltaV/T_e)` and agrees with forward QMC; and a periodic trench agrees
  with independently launched forward QMC globally and across floor/lower-wall/upper-wall/top regions.
- Charging now selects `forward` or `adjoint` per charged species and merges their exact event measures
  before the unchanged signed-current projection and Poisson coupling. An end-to-end charging gate uses
  forward transport for a directional ion and adjoint transport for a Maxwellian electron. This is a
  reusable engine control, not a geometry or benchmark branch.
- **The production trench current-balance gate remains red and is not weakened.** Forward directional
  ions still sparsely resolve local wall support, while a broad local-frame adjoint ion proposal has
  excessive importance variance on beveled arbitrary faces. The demonstrated hybrid reduces the worst
  residual but does not meet the declared tolerance. The next bounded engine task is a face-oriented
  ion proposal (then adaptive bidirectional selection/error accounting), followed by the same unchanged
  convergence and charging-dipole gates.

### Phase 2d arbitrary-face directional-ion follow-up

- The reversible 3-D gather now supports a per-species `source_aligned` proposal frame in addition to
  the original `surface_local` frame. Source-aligned samples retain the global plasma-boundary
  direction while each arbitrary triangle supplies only the incident-normal projection. This closes
  the orientation/support defect for narrow ions without a face-, region-, or benchmark-specific law.
- A periodic-trench ion reciprocity gate now compares independent forward QMC with the source-aligned
  adjoint globally and over the same four depth regions used for electrons. Standard seven-point
  triangle quadrature reduces the global visibility-integration error from 0.94% to 0.48%; the gate
  requires both refinement improvement and a 0.6% absolute budget.
- The unchanged steady trench root was rerun with source-aligned adjoint ions, a stratified
  field-shift-support proposal, hybrid forward-ion/adjoint-electron transport, Picard, and Anderson.
  None met the 0.08 maximum local current-balance criterion (best audited maximum remained far above
  tolerance), so the charging dipole/profile gate is still not promoted. A worst-node coordinate
  relaxation experiment produced unphysical oscillatory voltages and was removed rather than added
  as speculative solver machinery.
- The remaining engine gap is now narrower: the estimator needs adaptive bidirectional local-current
  error control that follows the field-shifted narrow ion support per surface element. A single fixed
  broad proposal can have global support yet only a handful of effective samples in each shifted
  two-eV window. The next implementation must expose per-face variance/effective-sample diagnostics,
  refine or switch estimators from those diagnostics, and include that uncertainty in the existing
  current-balance confidence envelope before the nonlinear root can be earned.

### Phase 2e variance-certified bidirectional 3-D transport

- Primary-source review confirmed the estimator architecture. Veach's multiple-importance framework
  requires unbiased technique combination and Grittmann et al. (TOG 2019,
  doi:10.1145/3355089.3356515) show why density-only MIS can lose the benefit of stratification.
  Owen's scrambled-net analysis (SIAM JNA 1997, doi:10.1137/S0036142994277468) supports independent
  randomized-QMC replicates for error estimation. Backward particle-transport validation likewise
  treats independent forward/backward agreement as the judge rather than assuming the adjoint is
  correct (Niess et al., CPC 2018, arXiv:1705.05636).
- The common arbitrary-triangle engine now runs independent forward and reversible-adjoint ensembles,
  measures standard error per triangle, applies a direct-forward upper bound to unresolved adjoint
  zeros, requires cross-estimator consistency when both directions claim precision, and selects the
  lower-uncertainty event measure per face. Selected replicate events are averaged without discarding
  their energy-angle distribution, so chemistry consumes the certified measure rather than a scalar
  current reconstruction.
- Triangle position sampling is independently scrambled along with velocity sampling. On the real
  zero-field trench this exposed a deterministic visibility bias that fixed Gaussian points had hidden;
  with eight independent replicates both ion and electron maps certify on all 40 faces, global landing
  is within about 0.1% of conservation, 39/40 ion faces select adjoint, one selects forward, and every
  electron face selects adjoint.
- Charging accepts `bidirectional` per species and refuses uncertified current maps. A failed sampled
  nonlinear trial is rejected without advancing the charge state. A flat coupled gate proves the
  selected face-resolved energy-angle events pass through signed-current projection.
- A global sample-refinement experiment was deliberately removed: a few inconsistent triangles caused
  all faces and species to be retraced and made a six-state trench audit impractically slow. The next
  bounded port from the mature 2-D lineage is selective adjoint retracing on only inconsistent triangle
  indices (while retaining the full collision mesh), with global forward refinement only when its
  source-launched audit is unresolved. The real root and dipole gate remain unpromoted until that
  controller certifies every accepted field state and meets the unchanged 0.08 local-balance target.

### Phase 2f deterministic root epochs and selective retracing

- Adjoint gathering can now launch from an arbitrary subset of target triangles while retaining the
  full mesh for collision queries. Refined events replace only those faces; unaffected sparse
  energy-angle measures remain bitwise unchanged. This is phase-space estimator adaptation, not AMR:
  the geometry, Poisson grid, and level-set `dx` do not change.
- Charging discovery freezes the certified per-triangle forward/adjoint method map and fixed scrambled
  Sobol rules before nonlinear iterations. The inner current map is therefore deterministic for a
  given charge state. Certification failures reject sampled discovery trials; production root epochs
  must be followed by an independent replicate audit rather than rediscovering methods every step.
- On the coarse periodic trench, the frozen deterministic Picard map reduces RMS imbalance from 0.788
  to 0.627 in 30 accepted steps, but the worst local imbalance remains 0.978. Maximum-residual trust
  merit and inverse-Broyden secant experiments did not close the root and were removed. The remaining
  blocker is now the strongly nonlocal/nonsmooth current-response solve, not global landing accuracy,
  Poisson response, periodic trajectories, or ordinary Monte Carlo noise.

### Charging solver campaign Task 0 decision gate (2026-07-13)

- **Task 0A is complete but its precision gate is red.** Config
  `b0d263fbce2bf3f29946957b03145b4f204ef31cba09a39bcb886aa2745d2db7` ran 360/360 exact
  hard-visibility common-random-number pairs at nested levels 9/11/13, radii 0.1/0.05/0.025 V,
  seeds 401/409/419/421/431/433/439/443, and five stuck-state worst/dominant directions. At level 13
  the minimum signal/between-scramble-error is 0.512 (required 3); maximum nested-level difference is
  1.275 sigma. The paired restricted ensemble conditions are 22.01, 17.10, and 12.49, but are not
  decision-grade and are not dimensionally comparable to the archived full 47-by-47 Jacobian.
- **The corrected switch evidence is mixed on conditioning but clear on geometry.** Signed
  `log(cond_full/cond_no_switch)` replaces the invalid clamped attribution and is synthetic-gated in
  both signs. Switch removal is not consistently helpful. Level-13 switch-fraction slopes are
  0.851--0.950 (median 0.937), supporting a regular switching set over the measured radii. The result
  is invariant to active thresholds `1e-5`/`1e-4`/`1e-3`.
- **Task 0B rejects the proposed local planar electron split.** Config
  `28f4887dc11ee0e7a34291637ee68e3de6e9902753babf1f9964a986cd9f0b8f` passes the refined flat
  Maxwellian barrier gate (RMSE 0.000500, maximum error 0.001204, maximum unresolved trajectory
  fraction 0.000488). A region-calibrated -12..8 V trench sweep passes the declared derivative gate in
  0/4 regions, so Task 2b is not promoted.
- **Decision:** the handoff's `response below estimator uncertainty` branch applies. Tasks 1--4 were
  not coded or run; no convergence contract was changed. The unavailable RMS-0.627 coarse-3D state is
  not silently substituted: the preserved audit checkpoint is the same-lineage 47-DOF stuck AR4 state
  (RMS 0.597912), with input names and SHA-256 recorded in
  `CHARGING_SOLVER_CAMPAIGN_2026-07-13.md`. Python ThreadPool execution is permanently prohibited for
  these audits after the Numba workqueue abort; process workers pin all nested math runtimes to one
  thread.

### Charging solver campaign scope amendment and Task 1-pre (2026-07-13)

- **Physical time is unblocked; response-based accelerators are not.** The original decision table
  routes every resolved Task 0A outcome through the conservative derivative-free transient. The
  signal/error gate instead controls Tasks 3--4 and any measured response preconditioner. The local
  Maxwellian electron split remains rejected, and no replacement global preconditioner is promoted:
  none of the measured response radii is yet decision-grade.
- **The archived Jacobian was paired.** Its plus/minus launches reused the checkpoint, deterministic
  proposal seeds, frozen method map, and frozen levels. The large full-space frozen-map condition
  therefore must not be reinterpreted as an unpaired-difference artifact. The present restricted
  ensemble comparison remains underpowered and cannot attribute the discrepancy.
- **Task 1-pre passes on the actual stuck current map.** Config
  `21cdb54c47754ec652b424260684b8fc8e15bdabd2c76716df2cc65b9795c3bd` used common samples for the
  production horizon (`dt=0.01`, 4,000 steps), 4x and 8x horizons, and a half-timestep reference over
  eight level-9 scrambles. Every estimator reported zero unresolved trajectories. Horizon extension
  changed no current; timestep halving changed the ensemble ion/electron currents by 0.287%/0.170%
  and RMS/worst-node imbalance by 0.130%/0.095%. The historical stuck residual is not a finite-horizon
  artifact. Exact hard visibility is unchanged. Task 1 may enter with its original refinement and
  final-audit gates.

### Charging physical-time and PTC campaign (2026-07-13)

- **This is engine work, not a test-only workaround.** The 3-D charging engine now has a reusable,
  conservative multi-step physical-time API with separate ion/electron node and face currents,
  replayable charge histories, exact-throughput conservation, final current evaluation, frozen
  method-map propagation, and resumable integration failures. The established 0.08 convergence
  equation is nodal; triangle imbalance remains a separate discretization diagnostic.
- **Deterministic physical time is stable but does not close the root.** A 50-microsecond paired ladder
  at 250/125 ns agrees within 0.298% charge and 0.139% potential and reaches RMS 0.195/0.197, worst
  node 0.439/0.440. A 500 ns step is rejected by 31.6% potential disagreement. The physical dipole
  forms, conservation remains at roundoff, and longer exact-map trajectories fluctuate above the
  unchanged node contract.
- **The independent ion audit certifies the map, not convergence.** Sixteen bidirectional replicates,
  level-13 ceilings, and 32 face-position points certify all ion faces. With a level-12 electron
  proposal the unchanged endpoint is RMS 0.193551, worst node 0.445733, worst face 0.841180. Electron
  between-scramble uncertainty remains to be attached before any final-state claim.
- **Fresh scrambles do not remove the floor.** A 200-step stochastic physical trajectory using a
  separately frozen method map has final-window RMS 0.199220 +/- 0.006366 and worst node
  0.574845 +/- 0.016865; its confidence envelope is nowhere near 0.08. Conservation residual is
  1.40e-17 of absolute throughput.
- **Safeguarded current-direction PTC collapses.** Six accepted pseudo-steps and 14 rejections reduce
  the pseudo-step below 1e-11 s with best RMS 0.179424 and worst node 0.465949. No derivative,
  quasi-Newton, smoothing, or altered residual was used. Further frozen-map solver variants are
  closed; the next bounded path is the handoff's nodal/face/coarsened-patch equilibrium audit under
  grid refinement, with the 0.08 contract unchanged.
- **Estimator debugging hardened the engine.** Frozen bidirectional maps now skip discarded estimator
  directions/faces without changing selected currents. Surface-local folded grazing densities are
  explicitly rejected in source-aligned frames; this prevents the slow-ray horizon failure exposed
  during Task 1 entry.

### Charging equilibrium/discretization refinement (2026-07-13)

- **The residual is not only a raw-mesh statistic.** Exact hard-visibility config
  `d4ce71baa71c43f4d66e3c6af5533ddf58c43e6c2c1adf9a5c9b32a94b44905b` transfers the same voltage
  state from `dx=0.25` to `0.125 um` with `2.19e-14` Poisson reproduction error. Global imbalance is
  0.03862/0.03781, but fixed half-micron wall patches have RMS 0.162/0.173 and maximum 0.296/0.280.
  Four wall patches remain independently resolved outside 0.08 with the same sign on both grids.
- **Worst-node magnitude is sampling-sensitive; RMS and patches are much less so.** Forward level
  9 -> 11 changes fine node RMS 0.3465 -> 0.3440 but worst node 0.8924 -> 0.8000. Half-micron patch
  RMS changes 0.1648 -> 0.1735. A certified method map alone is therefore insufficient provenance on
  small faces; its scoring sample/position levels matter.
- **The refined transient is stable but still red.** Config
  `5258cc6c59942d2e4b02558b0bd49cff90e0636319ddfcba45eab1a5687a7cc9` extends the accepted refined
  trajectory to 15 microseconds. The final 10-microsecond timestep pair agrees to 0.0166% charge and
  0.00559% potential. Level-11/13 endpoint audits give node RMS 0.29975/0.29888 and worst node
  0.7717/0.7143. At level 13, half-micron patch RMS/max are 0.1287/0.2587 while global imbalance is
  0.00058: upper-wall redistribution, not missing global current, is the surviving mode.
- **Certified replay no longer drops back to base resolution.** Internally discovered estimator maps
  now freeze at their declared adaptive sample and face-position ceilings. Externally supplied maps
  keep their explicit scoring levels. The 0.08 contract and physical operator are unchanged.
- **Bounded next action:** do not launch another root solver. Attach certified per-face sampling levels
  to the method-map artifact or run the accuracy-matched fine transient on GPU before extending
  physical time. In parallel, scope the currently absent surface-conduction, bulk-leakage, secondary-
  electron, and reflection closures as physics-model questions, not numerical patches.

### Charging surface-physics closure audit (2026-07-13)

- **The surviving imbalance is a localized material-response question.** Exact absorbing-map preflight
  config `6fdaa46785f8f9f861a2b6aece83cb5d5282814839f0cdb570af13af6fc08e94` replays the refined
  15-microsecond state with the frozen ion method map and level-11 scoring. Integrated signed balance
  is +0.032 top, -0.126 upper wall, +0.003 lower wall, and -0.030 floor. The upper wall alone fails
  the four-region 0.08 screen; a global current multiplier would damage already-balanced regions.
- **Impact energies now constrain the missing physics.** Upper-wall electrons have mean/median impact
  energies 14.4/9.7 eV, while lower-wall and floor electrons are field-accelerated to means of
  49.8/41.1 eV. The physical response must consume impact energy, angle, material id, and surface
  state. Source temperature or local voltage is insufficient.
- **Ion SEE is sourced and bounded, but unlikely sufficient.** The event-weighted Sobolewski 2021
  kinetic Ar+-on-SiO2 fit gives yields 0.031--0.065 and 0.052 on the upper wall. This is recorded only
  as a diagnostic and was not applied. Electron backscatter/true SEE is the stronger literature-backed
  redistribution mechanism, but its SiO2 parameter table and full incident direction data are not yet
  in the engine.
- **Unified engine path is specified.** `CHARGING_PHYSICS_CLOSURE_AUDIT_2026-07-13.md` defines a richer
  charged-impact measure, material-tagged surface response, conservative full-field re-impact loop,
  and common charge-continuity operator consumed by physical time, PTC, and steady diagnostics. The
  existing neutral-emission machinery supplies the conservation pattern but cannot replace charged
  field trajectories. Bulk leakage is held by a Maxwell-timescale screen; surface conduction is held
  for sourced material/surface-state conductivity. No new governing physics was applied.
- **P0 impact preservation passes exact regression.** The shared sparse energetic event measure now
  optionally carries immutable impact position and full unit incident direction. Forward/adjoint
  field transport, bidirectional selection/replacement, and geometry filtering preserve it. The
  absorbing-map preflight summary, CSV, and PNG retain identical git object hashes after replay;
  57 targeted surface/transport/charging tests pass with one unavailable-CUDA skip. No response has
  yet been applied.
- **The charging current now uses one conservative surface-transfer contract.** Perfect absorption
  exactly reproduces the historical accumulate-flux-then-multiply current order, but now exposes
  incident, outgoing, and deposited signed charge rates. Step and physical-time results report the
  response balance separately from compatible-Q1 projection balance. Artificial reflection/emission
  accounting gates close at roundoff; no yield or reflection law has been activated.
- **Charged surface-origin re-impact reuses the production field tracer.** Sparse outgoing particle
  rates are launched from their preserved impact points, advanced by the same velocity-Verlet nodal
  field kernel, and converted to target-face flux exactly once. Every rate is classified as landed,
  escaped, or truncated; truncation raises by default. A periodic provenance defect was fixed so a
  wrapped hit stores its in-cell intersection rather than a covering-space coordinate. Synthetic
  landing/escape/horizon gates conserve particle rate exactly. The default governing operator remains
  the exact perfect absorber; non-absorbing responses require an explicit caller-supplied model.
- **P1 charged response/re-impact cascade passes.** Material ids, face normals, and material state now
  enter one response context. Physical time and the direct steady diagnostic call the same optional
  response/cascade path; default perfect absorption remains exact. Artificial perfect-absorber,
  specular closed-cavity, Lambertian-per-ion, explicit escape, Q1 projection, quadrature, timestep,
  launch-offset, and bounce-cap ledgers conserve signed charge. A capped cascade is returned with its
  unresolved particle charge for diagnostics and is rejected by a charging advance.
- **P2 sourced ion-SEE is real but insufficient at the audited state.** Config
  `603a7df78819f06ff380b2263674c2e8491dfddc61240362bf10c9329f8e94f7` applies only the Sobolewski
  2021 Eq. (8) kinetic Ar+-on-plasma-exposed-SiO2 yield, with deterministic Lambertian full-field
  re-impact and a declared 1/3/5 eV emitted-energy sensitivity. The central 3 eV case changes nodal
  RMS/worst from `0.29975/0.77171` to `0.29257/0.75842`; upper-wall signed imbalance changes only
  `-0.12574 -> -0.12243`. Angular and isolated flight refinement are stable, charge closes below
  `4.9e-15`, and roughly 64.5% of emitted electrons escape. This does not approach 0.08. No full
  transient or convergence claim is promoted while the SiO2 emitted-energy spectrum remains a
  bounded input and electron-impact response remains data-gated.
- **Charged-event lineage makes the nonlocal response visible.** Every emitted event now preserves
  landed/escaped/truncated termination and its landed face; source-region routing rows independently
  close to one. In the central 3 eV replay, top emission escapes entirely, upper-wall emission lands
  mostly on the floor (61.9%), lower-wall emission is 91.0% self-recaptured, and floor emission lands
  mostly on the lower wall (63.9%). The aggregate 64.5% escape statistic is therefore not a local
  loss law: it mixes open-top escape with strong deep-feature recapture. The routing matrix is
  committed beside its heatmap in `results/charging_ion_see_sensitivity_3d/`. This reinforces the
  engine contract: surface response produces particles, shared field transport selects their
  destinations, and only the central transfer ledger updates charge.
- **The physical-time final audit now uses exactly the integrated operator.** Review found that
  response-enabled updates used the charged cascade but the final current evaluation omitted the
  caller's response/material/cascade arguments and reverted to perfect absorption. It also used the
  original bidirectional options instead of the map frozen during the trajectory. The final audit now
  receives the same response state and frozen estimator configuration; a zero-update regression
  requires the returned final transfer to be the explicit cascade, not an absorber transfer.
- **Certified estimator provenance is now an engine result, not campaign folklore.** Every
  bidirectional species result retains its actual forward sample level, per-face adjoint sample
  level, per-face position-quadrature count, and replicate seeds. Physical time and the direct steady
  diagnostic freeze a discovered method map at those measured stopping levels. The older
  ceiling-based replay remains only as a conservative compatibility path for method-only archives.
  New equilibrium pilot and physical-time checkpoints serialize the provenance beside the method
  map, so later current audits can distinguish estimator selection from estimator accuracy.
- **Steady and transient diagnostics now return the same evidence classes.** The direct steady result
  retains final face currents, the charged-transfer/cascade ledger, its frozen estimator map, and
  sampling provenance instead of exposing only projected node currents. This does not promote the
  retired solver family; it makes any diagnostic use auditable against the same conservation and
  exact-operator checks as physical time.
- **P3 remains correctly data-blocked, not numerically blocked.** The 2026 paper's electron-TSE law
  refers its six material coefficients to NASA TM-79299 but does not publish the SiO2 values; the
  NASA memorandum derives and sensitivity-tests the NASCAP form while explicitly avoiding any
  particular-material prediction. The paper makes supporting data available only upon author
  request. Its backscatter form also leaves the effective material atomic-number choice unresolved
  for SiO2/conditioned layers. No curve digitization, PR constant, or guessed effective-Z model was
  admitted. Author-supplied coefficients or an independent plasma-exposed-SiO2 dataset are required.

### Charging co-evolution C1/C2 pre-sign-off work (2026-07-13)

- **Convergence-contract revision `CCA-2026-07-13-R2` is signed and active.** Stan's two binding
  riders require at least two physical patch scales (one no larger than any claimed feature) and cap
  every experimental-claim observable tolerance at the benchmark's experimental plus digitization
  uncertainty. The old per-node RMS/worst-node diagnostic remains mandatory on every run including
  failures. C3 is authorized; its task-entry gates and every handoff no-go remain binding.
- **C1 moving-surface charge remap passes.** The declared closure carries charge with retained or
  advancing material, removes and itemizes charge on an etched-away layer, and initializes newly
  exposed material uncharged. No-op is bitwise identical; translation and signed removal ledgers
  close; a 4/8/16/32-cell refinement ladder has observed RMS orders 1.124/1.204/1.187 with maximum
  relative charge-ledger error `1.48e-15`. Positive and negative inventories are conserved
  separately. The operator is solver-independent and was deliberately not wired into C3 before
  sign-off; revision R2 now permits that integration.
- **C2 bounded grazing reflection passes its common-engine certification.** Config
  `b3df9406a4ac848dafa410eee924ca587217a39fef22065dfc70e6b998c1ed5d` uses exact hard triangle
  visibility and deterministic manufactured 100 eV Ar+ events. The explicit literature-bounded
  sensitivity (`p_grazing=0.95`, exponent `3`, retention `0.90`) gives `P=0.94905` at cosine `0.1`.
  Particle and kinetic-energy ledgers close exactly, charge closes to `2.37e-15`, and a reflected
  path reverses to its source face within the launch-offset bound.
- **Reflection is now transport/chemistry physics, not a diagnostic-only bounce.** Every charged
  re-impact is concatenated into the exact chemistry-facing energetic event measure with face,
  position, energy, angle, and direction lineage. The manufactured straight-wall gate sends all
  128 grazing reflected events to the floor within `0.0101--0.0653` mesh units of the corner; the
  reflection-off case sends none. An incomplete bounce tail remains explicit and is refused. This
  earns a qualitative microtrench flux mechanism, not a calibrated Si/SiO2 reflection prediction.
- **Evidence is replayable.** `CHARGING_COEVOLUTION_C1_C2_AUDIT_2026-07-13.md` and
  `results/charging_coevolution_c1_c2/` contain the gate table, plot, CSV, JSON manifest, base git
  revision, implementation SHA-256 values, parameter sources/bounds, hardware, and numerical
  settings. No frozen-map root method, smoothing, volume Boltzmann term, altered visibility, or
  hidden bounce cutoff was introduced.

### Charging co-evolution C3 unified-engine integration (2026-07-13)

- **The signed C3 engine transaction is now implemented.** Authoritative face sheet charge projects
  through the compatible Q1 operator, advances by the conservative kinetic current ODE, and is
  checked against the independently updated nodal charge at roundoff. The final hard-visibility
  charged/reflected transport object is passed directly to common surface chemistry; the profile
  advances through the existing level-set/material path; C1 then remaps signed charge to the rebuilt
  surface and warm-starts the next geometry. No root solver, smoothed operator, or volume Boltzmann
  term was added.
- **Fixed physical time and safeguarded SER are one operator.** SER uses the residual-ratio update with
  declared timestep bounds, growth cap, and residual-growth rejection, and only activates in its
  declared residual region. The exact current/deposition operator, compatible Q1 projection, and
  final hard-visibility audit remain unchanged. Every evaluation records node RMS/worst, face
  RMS/worst, two-scale patch metrics, potential rate, incident/deposited charge, and both projection
  and surface-transfer conservation errors.
- **Both signed riders are executable.** Runs require at least two physical patch scales. An
  experimental profile-feature claim is refused unless one scale is no larger than its declared
  physical extent. Experimental claims require explicit observable tolerances no larger than the
  combined benchmark and digitization uncertainty. B2 gates the contract's exact
  `abs(Ji-Je)/Ji` statistic; the historical symmetric node/patch diagnostics remain present but do
  not gate acceptance. This distinction was corrected before any real-trench C3 claim.
- **Pulsed bias uses resolved co-simulation, not a quasi-static label.** Quasi-static mode refuses a
  waveform. `waveform_resolved` advances the same charge ODE once per explicit physical segment,
  drives chemistry/profile motion with that segment's exact endpoint transport, then uses the same
  signed remap/field rebuild. A two-segment 1 ns ion-rich/electron-rich smoke gate completes one
  physical update per segment and closes deposition charge to `8.09e-17` and `4.05e-17`; saturation
  is neither required nor claimed.
- **Manufactured C3 integration passes; scientific C3 closure remains pending.** Config
  `d3b5485aff03a950c82f6fb4a0161e76532b5120dc7f5a075bb340a7a4c444fc` gives planar B1 rate 0,
  node RMS/worst 0/0, patch maxima 0/0 at 0.25/1.0 micrometers, exact transport reuse, and zero
  deposition/remap ledger error. The full suite is 370 passed, 1 skipped. The real-trench cold/warm
  branch, timestep/grid/sample refinement, observable invariance, and independent high-sample B5
  audit remain required; C4 has not started. Evidence is in
  `CHARGING_COEVOLUTION_C3_AUDIT_2026-07-13.md` and `results/charging_coevolution_c3/`.
- **Legacy nodal charging checkpoints are not silently promoted to C3 face state.** On the identical
  archived 176-face refined mesh, the compatible Q1 face-to-node map has numerical rank 121,
  condition `1.41e18`, and a best reconstruction error of 0.584 relative L2 / 0.681 relative Linf.
  The inverse is neither unique nor accurate, so migration is refused. The real C3 trajectory must
  start from zero face charge or a face-authoritative checkpoint produced by C3 itself; no
  minimum-norm or regularized sheet-charge guess is admitted.
- **The first bounded real-trench C3 trajectory now writes a genuine face checkpoint.** Config
  `f188a7eb1eb6a7476313ffa44af810c47432fe71c3f940a8a8313f3478c7b96e` starts the coarse 40-face
  trench from zero sigma and reuses only the separately selected estimator map. Over two 0.125
  microsecond updates, node RMS/worst falls 0.760/0.995 -> 0.513/0.930 and exact ion-normalized B2
  patch max falls from 1111/663 to 0.963/0.963 at 0.25/0.50 micrometers. It remains far from 0.08;
  this is a cost/integration pilot, not convergence.
- **A real production blocker in weighted reflection was fixed with a declared error bound.** Strict
  zero-tolerance behavior refused the trench because fractional specular weight never becomes
  literally empty, even though the unresolved charge falls to `1e-43` by 16 bounces. A positive,
  explicit response-tail tolerance now closes only a decayed post-flight tail by absorbing it on its
  current impact faces; global charge remains exact and the spatial-current L1 error is bounded and
  reported as twice the tail fraction. The default remains zero, and a nondecaying closed specular
  cavity still refuses. Tightening `1e-10 -> 1e-12` changes face sigma by `1.44e-14` relative L2 and
  potential by `4.53e-15`, while the tighter maximum L1 bound is `8.52e-14`. Pilot evidence is in
  `results/charging_coevolution_c3_trench_pilot*/`, with the paired differences and artifact hashes
  in `results/charging_coevolution_c3_trench_pilot/tail_refinement_comparison.json`.
- **The first real C3 timestep-refinement sequence tightens, but is not declared closed.** At 2.5
  microseconds, fixed 125/62.5/31.25 ns runs give node RMS 0.3963/0.3956/0.3952 and B2 maxima
  19.81/19.68/19.67 at the 0.25 micrometer scale. Successive halving changes face sigma 3.76% then
  1.32%, nodal charge 2.12% then 0.843%, and potential 0.261% then 0.164%. The sequence is tightening
  (face-charge observed order about 1.50), but the finest face state is not yet invariant and the
  state is not saturated.
- **SER/PTC now performs real trial rollback on the same conservative operator.** The real trench
  exposed that the initial implementation (a) used denominator-sensitive B2 as a step-rejection
  merit, (b) halved a timestep without rolling back the already advanced charge, and (c) could use a
  new SER step for the face candidate while the nodal candidate retained the old duration. SER now
  safeguards the dimensional current residual, restores face/nodal state and both clocks, retries
  rejected trials, and preserves the compatible-Q1 mismatch below `2.46e-16`.
- **Signed charge cancellation no longer creates false conservation failures.** The engine records
  positive, negative, absolute-throughput, and signed-net inventories independently and normalizes
  roundoff by positive-plus-negative throughput. The failure that exposed this was only `9.12e-34`
  C of residual divided by a near-zero `9.42e-22` C signed net. The corrected real SER run's worst
  relative deposition-ledger error is `9.53e-17`.
- **Bounded PTC agrees with fixed physical time and accelerates it, but equilibrium remains open.**
  At 4.3549 versus 4.3750 microseconds, 80-step SER and 140-step fixed time differ by 0.339% in face
  sigma and 0.176% in potential. SER uses 42.9% fewer accepted steps and is 1.65x faster; its two
  growing-residual trials are explicitly rolled back. The endpoint still has node RMS/worst
  0.350/0.790, B2 11.50/10.59, and potential rate `1.46e6` V/s, so neither B1/B2 nor C3 closes.
  Paired hashes and differences are in
  `results/charging_coevolution_c3_trench_refinement/comparison.json`.
- **A hidden particle-trajectory resolution failure was caught before scientific promotion.** At a
  later charged checkpoint the former `0.005` trajectory step sent an Ar+ through the top mask and
  reported a solid-side hit: stored cosine `+0.998886`, gas-normal geometric cosine `-0.998886`.
  The common field transport now reconstructs incidence from terminal velocity and the explicitly
  declared gas normal, and refuses back-face or inconsistent lineage before current or response
  physics consumes it. Production callers pass level-set gas normals through primary, bidirectional,
  neutral-reuse, and re-impact paths; focused charging/transport tests pass.
- **The real-trench particle trajectory now has its own seven-level audit.** On the durable fixed
  checkpoint, `dt=0.005` refuses while `0.0025` through `0.000078125` pass without truncation. The
  final halving changes node RMS/worst by `0.0216%/0.0159%`, exact B2 at 0.25/0.50 micrometers by
  `0.0697%/0.0708%`, and potential rate by `1.44%`. The bounded pilot default is `0.0003125`; it is
  within `0.214%` in B2 of the finest replay but differs `4.05%` in potential rate, so final B1 claims
  require the two finer levels. All earlier `0.005` transient/PTC runs are downgraded to controller
  mechanics evidence and will not seed the resolved physical continuation. Evidence and the input
  face checkpoint are in `results/charging_coevolution_c3_trajectory_refinement/`.
- **Flight horizon and hit precision are now separate certified controls.** At the resolved
  `0.0003125` particle step, the old 50,000-step horizon (flight time 15.625) refuses unresolved slow
  adjoint electrons; 128,000 steps (time 40) passes and doubling to time 80 leaves reported currents
  identical. The bounded runner now defaults to the time-40 horizon. This is a trajectory-completion
  requirement, not a charging convergence tolerance.
- **Rare float32 shared-edge misses are repaired inside the unified exact operator.** A fresh
  zero-charge continuation reached a solid-side top-mask hit at 2.375 microseconds. Raw Warp replay
  was non-monotone under particle-step halving (`refuse / pass / refuse`), diagnosing triangle-edge
  floating-point degeneracy rather than ODE error. The fast path now certifies every hit from terminal
  velocity and the declared gas normal; only an invalid lineage is replayed from its original state
  with the same fixed-step Verlet path and edge-inclusive float64 hard triangles. Nothing is softened,
  dropped, or accepted incomplete, and primary plus charged re-impact paths share the fallback.
- **The hard-visibility repair is counted and refinement-tested.** At the exact failure checkpoint,
  particle steps `0.0003125 / 0.00015625 / 0.000078125` require `1 / 0 / 1` replays while node RMS,
  worst node, and both B2 scales change by at most `0.0975% / 0.0220% / 0.203%`; potential rate changes
  `2.89%`. The repaired 20-step physical trajectory completes to 2.5 microseconds with exactly one
  replay, closes deposition/transfer ledgers to `1.51e-16 / 4.51e-15`, and remains unconverged at
  node RMS/worst `0.3970/0.8870`, B2 `19.81/17.45`, and potential rate `1.50e6` V/s. Evidence is in
  `results/charging_coevolution_c3_lineage_replay/`; C3 remains open and C4 remains unauthorized.
- **Pre- and post-repair residual numbers are not compared across operators.** The historical
  `0.788/0.627/~0.30` values used different state authority, response content, trajectory controls,
  and hit certification. Old float32 histories are potentially vulnerable to the now-demonstrated
  rare edge miss, but no claim is made that every old sample contained one. The separate stuck-map
  audit also showed that its historical residual was unchanged by 4x/8x horizon extension; the new
  adjoint horizon failure resulted from reducing timestep without preserving total flight time.
  Post-repair baselines therefore start at commit `8317f07`, without contradicting that prior audit.
- **The repaired physical transient advances, but is not approaching every gate exponentially.** A
  60-step, 7.5-microsecond zero-charge reference reduces node RMS/worst from `0.760/0.995` to
  `0.300/0.723`, B2 from `1111/663` to `8.436/7.524`, and maximum potential rate from `2.22e8` to
  `1.03e6 V/s`. B2 is a ratio: `8.436` is 843.6%, not 8.436%, versus the `0.08` gate. A late-window
  exponential fits RMS well (`R^2=0.997`, tau 5.71 microseconds) but predicts a `0.240` floor; B2
  fits also remain far above the gate, while potential rate has only `R^2=0.300`. No saturation-time
  projection or fit-and-jump is earned.
- **Checkpoint/restart is bitwise invariant.** The uninterrupted 60-step reference and a 20-step
  checkpoint plus 40-step restart produce byte-identical checkpoint files and zero difference in
  every stored state array at 7.5 microseconds. This clears restart mechanics, not the still-pending
  cold-versus-remapped-warm stationary branch gate.
- **PTC is schedule-sensitive and only modestly helpful.** A 0.5% residual-growth safeguard suffers
  four rejections and collapses from about 142 ns to 8.77 ns, making it slower than fixed time. A 2%
  schedule advances 4.206 microseconds in 30 steps with no rejection; at the nearest fixed checkpoint,
  RMS and B2 agree within `0.059%` and `0.128%`, while potential rate differs `3.27%`. It saves about
  10.8% accepted steps, not orders of magnitude. The first schedule is rejected; the second remains a
  bounded same-operator accelerator pending matched endpoint and schedule refinement.
- **Replay fraction is now a permanent engine diagnostic.** The transport, reflection cascade,
  charging history, and run summary carry both replay count and eligible field-lineage denominator.
  The exact failure state needs one replay among 8653 eligible paths (`0.0116%`); the 7.5-microsecond
  endpoint needs zero among 9079. No monotone replay increase with accumulated charge is observed.
- **Out-of-the-box accelerators were screened against the actual state contract.** Potential-space
  physical evolution is already represented by the Poisson response applied to conservative current;
  making voltage authoritative would require inverting the known rank-deficient face projection.
  A matrix mobility remains a possible face-space pseudo-time preconditioner, not physical time, and
  is held until it can preserve the ledger and pass schedule refinement. Fit-and-jump is rejected by
  the decay evidence; ML proposals remain downstream of a converged exact-operator training target.
  Evidence is in `results/charging_coevolution_c3_decay_audit/`.
- **The repaired 7.5-microsecond state fails B1/B2 robustly but is not yet sample-refined.** Eight
  independent `11/9` scrambles give node RMS `0.29145 +/- 0.00247`, worst node
  `0.69934 +/- 0.01126`, B2 `7.228 +/- 0.511 / 6.481 +/- 0.527`, and maximum potential rate
  `(1.126 +/- 0.133)e6 V/s` (95% intervals on the mean). These are orders of magnitude outside the
  gates, but nested `11/9 -> 12/10` doubling still shifts B2 by 5.3--8.4% and potential rate by
  3.8--21.3%, so no precise plateau or relaxation-time claim is made.
- **B2's controlling region is preserved and localized, not thresholded away.** New checksummed
  current-audit artifacts retain every face current and both physical patch maps. The worst patches
  are the lower mask sidewall, carry about 3% of total throughput at 0.5 micrometers, and collect
  roughly 8--9 times more electron than ion current. The smaller ion denominator drives much of the
  B2 sampling variation. One exact replay occurs among 145,269 eligible `11/9` lineages, and charge
  ledgers remain at roundoff.
- **Independent-scramble provenance is now complete.** The C3 runner formerly refreshed forward
  samples with `--seed` while silently retaining adjoint proposal seeds 79/83. It now records and
  applies Ar+ seed `s` and electron seed `s+4`; default seed 79 is unchanged. Evidence is in
  `results/charging_coevolution_c3_sample_audit/` and the C3 audit report.
- **The apparent 7.5-microsecond floor was a short-window fit artifact.** Continuing the repaired
  fixed transient to 15 microseconds lowers node RMS/worst to `0.222/0.570` and B2 to `4.135/3.710`.
  A 125/62.5 ns comparison over the next 1.25 microseconds agrees within `0.031%` in face sigma and
  `0.057%` in potential. Maximum potential rate remains nonmonotone and above `1e6 V/s`.
- **A frozen low-sample map can steer physical time in the wrong mean direction.** From the common
  15-microsecond state, the fixed level-10/8 path raises floor potential `9.65 -> 11.28 V`, while
  independent level-11/9 audits predict a negative ensemble-mean floor rate. Eight fresh-scramble
  paths instead end at `9.208 +/- 0.331 V`; exact conservation alone does not remove finite-sample
  drift.
- **Fresh-scramble physical time is now a first-class unified-engine mode.** Every accepted update
  regenerates all forward/adjoint samples from a recorded seed epoch, and the final diagnostic uses
  the next unused epoch. Fresh-scramble SER is refused because noise cannot drive deterministic
  rollback. Halving 125 -> 62.5 ns changes the eight-path ensemble mean by `0.218%` in face sigma and
  `1.86%` in potential, with overlapping floor-potential/B2 intervals. This earns continued bounded
  stochastic integration, not C3 closure.
- **The proposed acceleration memo was filtered through the actual estimator/operator contracts.**
  `1/sqrt(N)` Richardson extrapolation is invalid for scrambled Sobol and nonlinear B2 maxima; the
  controlling ion faces use the frozen forward estimator, not adjoint; and a local voltage offset
  must be expressed as a globally consistent face-charge perturbation. Diagonal face pseudo-time is
  held until an ensemble-mean response and stable stationary state exist. Evidence is in
  `results/charging_coevolution_c3_stochastic_transient_audit/`.

- **The C3 late-time obstruction is a discrete state-space incompatibility, not evidence of an
  impossible kinetic equilibrium.** The declared 0.25 micrometer trench has 40 P0 face-charge
  degrees of freedom but a rank-34 face-to-Q1 nodal coupling, leaving six exact field-null modes.
  At segment 14, 96.8% of area-weighted stored face-charge norm and 95.3% of terminal-window
  net-current norm lie in that null component. Fifteen independent windows remain at least 0.99496
  cosine-aligned with their mean, so the component is systematic rather than sampling noise. Raw
  patch functionals have up to 0.426 dual sensitivity to the null space; identical Q1 fields can
  therefore report different raw B2 values. Evidence is in
  `results/charging_c3_q1_compatibility_audit/` and
  `CHARGING_C3_Q1_COMPATIBILITY_AUDIT_2026-07-14.md`.
- **A compatible-Q1 charge state is implemented without changing the resolved operator.** The Q1
  nodal load is authoritative and its unique area-weighted minimum-density-norm face
  representative is retained for remap/ledger operations. Zero-update projection preserves nodal
  charge to `4.46e-15` relative L1, potential to `4.29e-14` relative L2, and global charge to
  `2.01e-29 C`. A paired 500-step fork agrees with the legacy resolved trajectory within 1.75% in
  every audited field metric while reducing the face-state null fraction to `2.68e-15`. The runner,
  checkpoints, heartbeats, and unattended supervisor declare the state choice; raw B2 and retained
  node diagnostics remain reported.
- **Batch projective PTC is rejected at the compatible segment-16 state.** A fixed-state direction
  required 70 independent level-13 scrambles to reach signal/error 3.140. Exact paired scoring on
  unused epochs shows 2.5/5/10 microsecond candidates increase current-residual L2 by
  `2.058/5.161/10.898e-12 A`; physical-scale candidates are either inconclusive or worsen node
  metrics, and cost more to construct than direct physical steps. The reference continuation
  remains fresh-scramble compatible physical time. CCA-R2 is unchanged; the proposed separation of
  field-compatible balance from unresolved-current grid error is recorded, unsigned, in
  `CONVERGENCE_CONTRACT_R3_REVIEW_REQUEST.md`.
- **The long C3 continuation exposed and stopped on a particle/field topology mismatch.** Particle
  trajectories wrapped periodically in lateral x/y while the Q1 Poisson endpoints were independent
  natural-Neumann nodes. The stopped checkpoint has seam jumps of `20.413 V` in x and `24.150 V` in
  y; recomputing from its archived nodal charge reproduces those values to `2.90e-16` relative L2,
  proving the mismatch is in the declared operator rather than checkpoint corruption. The GPU was
  stopped after two completed compatible segments at `2.94875 ms`; it is idle and that archived time
  is not credited to the corrected model.
- **Periodic topology is now a shared engine invariant rather than a runner convention.** The Q1
  system algebraically identifies endpoint nodes before factorization, reduces charge to independent
  periodic DOFs, prolongs a bitwise-continuous voltage grid for transport, and constructs the
  compatible face state against that same reduced space. Coupled transport refuses a mismatch before
  tracing. The corrected seam is exactly zero in both axes; full tests pass (`401 passed, 1 skipped`).
- **The stopped inventory was converted only into a corrected warm proposal.** The intended periodic
  `72 x 40` face coupling has rank 22/nullity 18. Projection preserves its effective nodal load to
  `2.65e-15`, global charge to `1.06e-29 C`, and Poisson balance to `8.48e-30 C`, while changing the
  proposed field by 5.45% L2 (13.62 V maximum). Evidence and SHA-256 provenance are in
  `results/charging_c3_periodic_topology_audit/`. The next compute is bounded cold/warm corrected-
  operator verification, not a replay of the old millisecond march.
- **The repaired periodic physical reference is valid and still slowly charging.** Four consecutive
  400-step, 50 microsecond fresh-scramble CUDA windows complete with zero rejects, zero bounce or
  trajectory-horizon extensions, exactly zero voltage seam, and charge/transfer ledgers at
  `3.87e-15`--`4.95e-15` / `2.79e-15`--`2.84e-15`. Integrated B1 falls
  `60,510 -> 20,580 -> 17,389 -> 13,962 V/s`; Q1-resolved B2 falls
  `0.306 -> 0.200 -> 0.165 -> 0.151`, while raw B2 remains around 1.2--1.4 because it includes
  exact Q1-null face content. Every window passes independent integrity replay and fails signed R2.
- **The physical voltage motion separates into one collapsing coherent mode and late estimator
  jitter.** SVD of the four archived fixed-window displacement vectors assigns 96.597% of their
  energy to one direction. Its window projections are `-17.341/-2.990/-1.705/-0.374 V`, while the
  orthogonal L2 motion is `0.536/1.873/1.937/1.863 V`. The final small projection alone does not
  certify a stationary cloud because fixed-state endpoint audits at levels 11/12/13 cannot resolve
  the declared `1,000 V/s` B1: global signal/standard-error is at most 1.211 and the selected
  component intervals contain zero.
- **A decreasing-gain fresh-scramble tail is implemented as a bounded warm-start accelerator, not a
  convergence loophole.** The schedule `125 ns * (16/(16+k))^0.75` advances 400 conservative
  pseudo-steps (`10.1224 us` pseudo-time, zero physical time), preserves its gain age through
  restart, refuses frozen samples and stochastic SER, never self-certifies, and closes charge to
  `4.42e-14`. Config
  `2025866ccba91b4b8fc65af30b74b71b786f1056d6a3af535d82c21ce95c1e42` completes 400/400 with zero
  rejects; independent levels 11/12/13 endpoint audits reduce aggregate Q1/node diagnostics but do
  not establish a statistically decisive residual improvement. Full tests pass: 404 passed, 1
  skipped.
- **The fixed confirmation rejects overshoot and also rejects formal saturation.** From the tail
  checkpoint, config `c919399966fe14254d58bdcdfceacab5715240abccbdc0d03777dc87652ff35d`
  runs another 400 fixed physical steps/50 microseconds and moves `-2.272 V` along the archived
  dominant direction (same sign as the tail's `-2.200 V`), rather than returning toward the 200
  microsecond state. It is therefore the canonical warm checkpoint, SHA-256
  `2ed2c604728481fe217ccaa73d405905f79ccf24ec0cbb7ac0daf173ef51ce64`, but not an equilibrium:
  B1 is `17,748.6 V/s`, raw B2 `1.348/1.297`, Q1-resolved B2 `0.1681/0.1677`, and retained node
  RMS/worst `0.0385/0.0743`. Independent audit passes every integrity check and reports
  `contract_converged: false`.
- **The bounded July 14 compute campaign is closed without pretending the contract passed.** The
  exact operator, compatible state, fresh-scramble physical time, conservative restart, hard-hit
  replay, GPU path, and stochastic warm-start machinery are operational. Remaining blockers are
  continued resolved mean drift, an unaffordable pointwise stochastic B1 threshold, the raw/Q1
  mixed-space B2 contract, and missing physical-grid/B3 refinement. These require evidence and
  formal R3 review, not another frozen-map solver or unbounded march. Final evidence is in
  `CHARGING_C3_CLOSURE_AUDIT_2026-07-14.md` and
  `results/charging_c3_periodic_topology_audit/closure_mode_audit.json`. Vast instance `44895783`
  was destroyed after artifact retrieval; the final instance listing was empty.

## 2026-07-15 — unified charging product closure and finite-arrival ensemble wiring

- **The historical signed-R2 path remains the default and no convergence claim was rewritten.** A
  separate `CCA-PROFILE-STATIONARY-2026-07-15-DRAFT` contract now asks the profile-relevant question:
  whether a second fresh physical-time block changes potential, independently scored kinetic current,
  delivered species flux, profile velocity, or the predicted profile increment beyond declared
  tolerances. Endpoint scores use disjoint independent epochs, uncertainty expands every stochastic
  metric, hard visibility is mandatory, and the draft explicitly refuses experimental claims.
- **Independent transport ensembles are now averaged in the engine rather than in campaign scripts.**
  Sparse energy/angle/position/direction events are concatenated with replicate-scaled flux; neutral
  fields and hit/escape probabilities are averaged; incompatible operators/species/lineage schemas
  refuse. C3 can run two blocks, score both endpoint ensembles, reuse the accepted mean transport in
  chemistry, and still reports signed-R2 RMS/worst/B2 diagnostics.
- **The common charged path is public and restartable.** `PhysicalChargingProcess` invokes C3 directly;
  in-memory continuation preserves geometry, remapped charge, conservative surface state, and mesh
  fingerprint while deriving a new seed. `petch-charging-checkpoint-3d-v1` safely stores the same step
  boundary in `allow_pickle=False` NPZ, rejects undeclared arrays/types, and binds the source manifest
  by SHA-256. A zero retained-charge inventory also no longer fails C1 on an irrelevant geometric
  correspondence.
- **Physical shot noise is now distinct from quadrature noise.** For each sparse kinetic event the
  engine draws `N ~ Poisson(event_flux * physical_face_area * physical_duration)` and returns the
  exact compact flux `N/(area*duration)`. Perfect-absorber deposition is an integer multiple of
  elementary charge and the existing global ledger remains roundoff-closed. Poisson C3 runs require
  fresh sampling and cannot self-certify a stochastic snapshot. Surface-response branching is still
  conditional-mean unless an explicitly stochastic response law is supplied and is reported as such.
- **A true ensemble execution path and twist observables now exist.** Constant-boundary
  `physical_time_resolved` C3 co-evolves charge and geometry without a quasi-static gate;
  `PhysicalChargingEnsembleProcess` refuses non-Poisson/fake-shot-noise configurations and labels one
  realization non-predictive. Geometry-native centerline, lateral displacement, equivalent diameter,
  twist-onset depth/AR, ensemble variance, and systematic-direction z-score pass manufactured tilted
  hole and symmetric +/- tilt gates. N>=30, N/sample doubling, isotropy, AR sweep, and dense/sparse
  campaigns remain required before any twist claim.
- **Verification:** the complete local suite passes `425 passed, 1 skipped`; `git diff --check` and
  public-export smoke checks pass. Work is local only and nothing was pushed. Design and remaining
  gates are in `UNIFIED_ENGINE_WIRING_PLAN_2026-07-15.md`.

### Unified engine material/campaign completion increment

- **Profile claims now have geometry-native measurements rather than hand-read contours.** The common
  API measures left/right notch depth, notch asymmetry, bow width/expansion, hole centerline,
  displacement, equivalent diameter, and onset AR from declared physical bands/ROIs. Manufactured
  asymmetric trench and tilted-hole gates converge under refinement; ensembles carry componentwise
  standard errors/confidence bounds.
- **C4 and C5 are executable evidence protocols, not claims.** C4 now requires checksum-bound source
  imagery, replayable pixel coordinates/transforms, digitization uncertainty, a committed one-family
  calibration with at most two bounded parameters, exact held-out coverage, independent run-manifest
  hashes, charging-off causality, and a decomposed uncertainty envelope. The legacy hard-coded Hwang
  values are explicitly ineligible because they lack that provenance. C5 now executes nested N/2N and
  paired sample-doubling comparisons, an isotropy control, and systematic-direction scoring, while
  enforcing base N>=30 in a real campaign contract. Neither campaign has run, so no notch/twist
  validation claim is made.
- **Mask and substrate now share one engine without sharing one chemistry law.**
  `MaterialMechanismRouter3D` dispatches exposed face material IDs to independent mechanisms,
  namespaces and conservatively remaps their state, merges signed material/product ledgers, moves both
  material level sets, recomputes ownership after motion, and safely round-trips its state through the
  charged checkpoint. A non-routed mechanism is refused when multiple evolving materials are exposed.
- **The first honest redeposition feedback loop is wired.** Explicit transport-ready product
  populations can be assigned bounded, sourced sticking and density laws. Diffuse flight closes
  emitted = deposited + escaped material, captured same-material inventory becomes signed interface
  growth, and the next step can resputter that grown material through the same mechanism. Positive
  cross-material capture or pinned-surface growth refuses because either needs a new evolving material
  layer. Reactive SiO2 branching remains unresolved rather than guessed.
- **Run provenance now binds the actual operator.** Versioned charged/uncharged manifests include
  geometry/material array hashes, boundary distributions, species roles, mechanism/redeposition
  provenance, source and electrostatic coordinates, numerical settings, device, and exact-operator
  statement. C3 also itemizes float64 lineage replay, trajectory/bounce horizon recovery, the priced
  cascade-tail bound, charge-remap conservation, and redeposition conservation. The existing durable
  supervisor/heartbeat path remains the process-level recovery mechanism.
- **Verification:** focused unified-engine gates pass `97 passed`; the complete local suite passes
  `439 passed, 1 skipped in 129.94 s`; `git diff --check`, module compilation, and public-export/schema
  smoke checks pass. No GPU/remote resource was launched and nothing was pushed.
- **The assembled workflow was then executed outside pytest.**
  `scripts/unified_engine_smoke.py` ran the public two-material mask/substrate path with both level sets
  moving and same-material redeposition balance `1.16e-16`; ran public quasi-static C3 plus a safe
  disk checkpoint/resume with zero charge-remap error and exact transport reuse; and ran two distinct
  physical-Poisson finite-arrival realizations (seeds 81/182). The standalone C3 audit passed
  quasi-static, waveform, conservation and refusal gates; the independent planar charging demo
  converged in two iterations; a bounded unified Jeon step closed neutral balance to `1.14e-12` while
  correctly reporting nonpredictive inputs; and the separate legacy facade completed its 40-step hole
  example in `9.15 s`. These are operational smokes, not experimental C4/C5 validation.
- **The distributable package was tested rather than assuming the checkout represented it.** The
  first isolated-wheel import exposed a real release mismatch: wheel metadata declared `0.3.0` while
  `petch.__version__` declared `0.2.0`. Runtime version is now `0.3.0`, a regression test binds it to
  `pyproject.toml`, and rebuilt wheel SHA-256
  `7bcc0b74fa1b25c85b4c0d01385fedb78fc41daf0ee2574c08fb3c2f65cc3ed0` installs and imports from a
  fresh temporary environment with matching runtime/metadata versions. The installed wheel—not the
  source checkout—then passed `scripts/unified_engine_smoke.py`. Final verification passes
  `440 passed, 1 skipped in 129.91 s`, module compilation, and `git diff --check`.

## Roadmap (remaining)

1. Acquire and checksum the exact Nozawa/Hwang figure asset, record pixel transforms and digitization
   uncertainty, commit the C4 calibration split, then run the charging-off/on plus grid/time/sample
   ladders once. No held-out retuning.
2. Run C5 on a 3-D hole AR sweep with base N>=30, nested N doubling, paired transport-sample doubling,
   and the zero-direction isotropy control. Report distributions, never one deterministic twist.
3. Add production mask chemistry, reactive product branching, or a cross-material coating layer only
   when C4/C5 morphology provides causal evidence and sourced parameter bounds. The current refusals
   remain in force meanwhile.
4. Replay one versioned manifest in a clean environment and publish accuracy-matched CPU/CUDA event,
   conservation, and profile comparisons; batch independent realizations before optimizing kernels.
5. Run de Boer and Jeon held-out transfers through the same explicit common APIs, with grid/ray/
   digitization/model errors separated and charging enabled only by an observable causality test.
6. Extend adjoints through chemistry and moving boundaries only after the forward validation gates;
   add spatial AMR only when scale/refinement evidence earns its complexity.

## Guardrails honored

Single writer. Local only. No benchmark/AR/region branch in governing physics. Sticking is the one
declared physical surface input with provenance. Refinement is error-driven. Device-agnostic transport
(Warp CUDA path exercised in the bounded C3 campaign; device-independent audit artifacts retained).

## 2026-07-15 — ARDE history reconciliation and common moving-profile increment

- **The old deep-AR experimental claim is withdrawn at its source.** Git history and all reachable
  branches/reflogs were inspected. The normalized `1/.43/.29/.20` sequence came from an evaluated
  Blauw/Clausing curve, not direct de Boer Figure-9 pixels; no missing photon implementation was found.
  Direct Figure-9 pixels and their checksums are development data because they have now been scored.
- **Charged reflection now has one reusable engine bridge.** Face-resolved charged arrivals pass
  through the certified cascade, exact float64 shared-edge replay, lineage, and charge/energy/tail
  ledgers before the identical resulting flux reaches chemistry. Incomplete cascades and transport
  records without exact impact position are refused. A large-event reduction-order mismatch was fixed
  at the authoritative integrated-particle-rate ledger, without loosening conservation tolerance.
- **A real moving-profile defect was found and fixed.** Transport used periodic lateral topology while
  level-set advection/redistancing evolved duplicate endpoint planes independently. The seam could
  therefore open after motion and send a reflected particle into a shifted solid image. Godunov
  advection now wraps the unique periodic core, endpoint planes are reconstructed bitwise-identically,
  and redistancing uses wrapped padding. The feature and charging co-evolution paths share this rule.
- **The direct common-engine development score improved, but did not validate.** With the one legal
  `s_F=0.10` calibration unchanged, adding literature-bounded certified ion reflection improves the
  twelve transfer-point de Boer RMSE from `3.5446` to `2.5568 um`; the calibration-marker error is
  `0.0241 um`. The remaining miss is retained as evidence, not tuned away.
- **The operator was exercised on an actually moving trench.** A three-level nested-sampling gate ran
  common Belen chemistry + neutral radiosity + certified reflection for 150 s over five profile steps.
  Final centerline depths were `0.218634`, `0.217613`, and `0.217393 um`; the latest sampling delta is
  `-0.000220 um`, 4.6x smaller than the previous delta. The moving result differs from the identical
  frozen-rate counterfactual by `-0.001184 um`, about 5.4x the latest sampling delta. Charge and energy
  close at roundoff, the cascade tail remains bounded, and state remap conservation passes. This is a
  bounded mechanism/history demonstration, not experimental validation.
- **Verification:** the final complete suite passes `461 passed, 1 skipped in 134.60 s`; the
  regenerated moving-profile audit passes its
  conservation gate. Public evidence pages now label the legacy plot as a calculated-curve replay and
  state that independent feature-profile validation remains open. Work remains local; nothing pushed.

## 2026-07-15 — ViennaPS La Magna fluorocarbon chemistry in the common engine

- **The ViennaPS-class generic fluorocarbon mechanism is now selectable, not hand-waved.**
  `LaMagnaGarozzoFluorocarbonMechanism` implements the La Magna--Garozzo etchant, polymer, and
  etchant-on-polymer coverages; chemical, ion-enhanced, and physical removal; polymer saturation,
  deposition/removal, and finite film inventory; and signed material ledgers. Every transferred
  parameter names the exact ViennaPS 4.6.1 source commit
  `2956ed587984c6dc38be24c6e2390e10c9b2f0a7` and remains nonpredictive until experimental transfer
  evidence exists.
- **It uses the existing engine.** Mechanism-owned polymer growth is combined with recession in the
  common face velocity, survives material-ID routing, uses the same level-set/profile path, and can
  compose with the existing charged/reflected event measure. Safe checkpoints now register the state
  and preserve material-router remap bounds/modes.
- **A real state-remap type error was found by the end-to-end smoke.** Algebraic coverage fractions
  had been treated as conserved areal inventory, so a saturated coverage could exceed an artificial
  capacity when the triangulated surface contracted. The remapper now accepts explicit `intensive`
  versus `conservative` field declarations; all legacy states retain conservative defaults. The new
  model interpolates its three coverages while exactly conserving polymer-film and removed-material
  inventories.
- **Jeong radicals can now reach the chemistry without losing species roles.** The default archived
  `FC_total` closure is unchanged. The explicit `heavy_light` mode maps CF/CF2/CF3 to `FC_etchant`
  and C2F4/C3F6/C4F7 to `FC_polymer`; their one-way thermal fluxes sum back to `FC_total` at machine
  precision. `scripts/jeong_2023_transfer.py --chemistry-model lamagna_garozzo` selects the split and
  common neutral/surface fixed point; the old reduced law remains the default.
- **The first transferred-parameter result is correctly a failed validation.** A deliberately coarse
  80-step moving Jeong smoke completed in `19.66 s`, but predicted net deposition (`-82.0 nm` depth)
  at the `1223 nm` calibration anchor. At that source condition the grouped polymer thermal flux is
  about 59 times the grouped etchant flux, so raw Vienna defaults saturate polymer. No held-out datum
  was tuned; this result narrows the next scientific question to the radical role/weight and parameter
  transfer rather than engine execution.
- **Verification:** analytic coverage/yield parity, ion-only/deposition/film-depletion limits,
  conservative ledgers, growth advection, material routing, contracting-surface remap, checkpoints,
  boundary partition, public imports, compilation, and focused legacy gates pass. The full snapshot
  reports `514 passed, 1 skipped, 3 failed`; the same three unrelated failures reproduce at the exact
  pre-increment snapshot commit `56ac284` (one obsolete zero-field replay monkeypatch expectation and
  two Jeon-v1 data/split-count mismatches). They are not attributed to this increment and were not
  folded into the chemistry patch. Full evidence and remaining parity/validation gaps are recorded in
  `VIENNAPS_FLUOROCARBON_CHEMISTRY_AUDIT_2026-07-15.md`. Work remains local; nothing pushed.

## 2026-07-15 — Jeon evidence correction and nominal-width transfer audit

- The exact 20% CW development deck was recovered from the prior local Codex session and replayed
  exactly: complex-formation probability `6.309573444801936e-4`, substrate polymer-deposition
  probability `7.498942093324559e-5`, QMC level 10/seed 2017, and 4.25-degree component IADF. Every
  archived predicted depth replayed identically.
- The archived `open48` result did **not** score Jeon's nominal 60 nm marker. The 48 nm opening was an
  unsupported development hypothesis, and the scorer silently omitted it because target matching was
  exact. The earlier all-points calibration interpretation is withdrawn.
- The Jeon runner now writes a complete canonical input deck and SHA-256, serializes all surface
  probabilities and their bounds/sources, refuses unregistered simulated widths by default, and marks
  registered-width coverage explicitly. Diagnostic nonregistered openings require an opt-in and can
  never populate the calibration result.
- With the actual 60/80/100/150/180/200 nm geometries and the frozen deck, non-anchor interval coverage
  is `4/5` at 20% CW (`log-RMSE 0.14696`), `4/5` at the untouched 40% CW condition (`0.06281`), and
  `2/5` at the untouched 80% CW condition (`0.07963`). The nominal 60/200 predictions are
  `0.48175`, `0.48021`, and `0.48257`, versus experimental `0.34427`, `0.49578`, and `0.43000`.
- These runs advance only about `0.07 nm` in one step. They are useful initial transport/rate
  diagnostics, not replays of Jeon's roughly micrometre-scale evolved profiles. The experimental
  claim remains open; Jeong 2023's reported 1200 s duration is the absolute-time transfer path once
  its modeled radical densities are converted to reactive surface fluxes under an explicit closure.
- Verification now passes end to end: `517 passed, 1 skipped in 132.18 s`; `git diff --check` passes,
  the exact deck replay is deterministic, and the new geometry mismatch refusal fires before
  transport. The regression also removed an import-time Numba thread-environment side effect and made
  the legacy Jeon-v1 scorecard's common-wall-time development hypothesis explicit; neither change
  alters the physical operator. Work remains local; nothing pushed.

## 2026-07-16 — Jeong reactor-to-feature closure audit

- **The Jeong ion boundary is no longer hard-coded to one implicit Ar population.**
  `Jeong2023IonBoundaryClosure` accepts either a declared positive-ion density mixture with
  species-specific Bohm fluxes or direct species-resolved fluxes from diagnostics/a validated
  reactor model. The common feature engine already supported multiple energetic populations; the
  experiment adapter now preserves that capability and reports species flux, energy, hit, and
  escape diagnostics. Multispecies runs refuse the existing Ar-like reflection law until each ion
  has a declared material response.
- **A published Huang--Kushner channel missing from the reduced chemistry is now conservative.**
  Fluorocarbon ions in the 5--70 eV window can deposit polymer directly on oxide-fluorocarbon
  complexes with the Table-I `p0=0.1` law. Ar does not use this channel. Polymer deposition/removal
  remains ledger-closed, and species/energy/complex controls have manufactured tests.
- **The bounded test killed the channel as a Jeong slope closure.** Moving an intentionally generous
  20% of total ion flux to 15 eV `CF+` reduced the zero-D 200 nm density-sweep endpoint gain from
  `999.02` to `774.96 nm`; experiment is `451.47 nm`. It supplies only `40.9%` of the required
  reduction. Collisionless virtual-sheath broadening previously supplied at most `6.68%` yield
  flattening and no 5--70 eV population. No moving profile matrix was launched.
- **Historical Jeong scores are now automatically labeled stale when the operator changes.**
  `scripts/summarize_jeong_2023_validation.py` compares authoritative implementation SHA-256 values.
  The archived five held-out 200 nm runs used the old implicit polymer projection and now report
  `historical_stale_operator / historical_results_not_current_validation`; their numerical trends
  remain historical evidence, not current validation.
- **The missing response is quantified without pretending it is a model.** A diagnostic inversion
  of the three already-scored 200 nm flux depths requires effective response multipliers
  `1.2889 / 1.0543 / 0.7367`, equivalent to
  `Gamma_effective proportional to n_e^0.463`. Because the scored data determine that curve, it is
  explicitly development evidence and cannot certify prediction. The current zero-D anchor scale is
  `1.526120494`, but no profile scale is frozen.
- **Two cheap current-operator profile preflights refused safely.** At 90 and 180 steps the repaired
  operator moved the interface farther per step than the conservative surface-remap radius allowed
  (refusals at steps 3 and 6). This is a timestep/preflight result, not a physics failure. A finer
  profile schedule is not worth purchasing until a reactor boundary is evidenced.
- **Decision:** strict Jeong flux validation remains open on missing species-resolved ion flux/IEAD,
  complete radical wall fluxes, and hot-neutral production. The engine-side interfaces needed to
  consume those inputs now exist. Evidence is in
  `JEONG_2023_REACTOR_CLOSURE_AUDIT_2026-07-16.md`,
  `results/jeong_2023_reactor_closure_audit/audit.json`, and
  `results/jeong_2023_reactor_closure_audit/closure_audit.png`.
- **Verification:** focused boundary, surface, feature, and charged-cascade regressions pass
  `110 passed`; compilation and `git diff --check` pass for the increment. Work remains local;
  nothing pushed.

## 2026-07-17 — Krüger 5 nm fixed-pair terminal and periodic remap diagnosis

- **The 10 nm-calibrated pair did not transfer unchanged to the 5 nm operator.** The development
  trajectory reached a resolved internal pinch-off at `56.920695 s`, before the declared `60 s`
  endpoint. The pre-event state had depth `852.4593 nm`, top/mask opening `40.7589 nm`, maximum and
  top feature width `86.9869 nm`, remaining mask thickness `850.1051 nm`, and asymmetry
  `0.8006` cells. The topology gate observed one new enclosed gas cavity while solid and both
  material component counts remained unchanged. The candidate closing step was not remapped or
  accepted.
- **A numerical remap defect was isolated before that terminal.** At `56.482184 s`, reducing the
  candidate timestep by 8x left an apparent approximately `10.2 nm` remap jump unchanged. The old
  surface-state search used ordinary Euclidean coordinates even though x/y geometry and transport
  were periodic. A face adjacent through the `20 nm` periodic lateral seam therefore looked
  nonlocal. The repair certifies exact point-to-triangle distance after local candidate culling and
  uses periodic images for state interpolation. It still refuses a truly separated surface.
- **The repair was regression-tested before resume.** A cloned exact checkpoint accepted the next
  step with zero material-ledger residual and only one bounded subcell gas-node closure
  (`1.25e-25 m3`). The feature/remap/Krüger contract suite passes `84 passed`. The original run was
  then resumed from its exact checkpoint with an explicit operator-transition sidecar.
- **Claim restriction:** because steps `0--1814` used the old nonperiodic lookup, the terminal is a
  development-only mixed-operator result. It may diagnose the fine-grid response and propose a
  calibration trust-region direction, but it cannot support validation. An authoritative 5 nm/AMR
  confirmation must start at `t=0` with one checksum-bound operator.
- **The existing base-only response matrix gives a useful direction, not an endpoint solution.** At
  the pre-clog state, increasing the crosslinked-film fraction reduces mask-film growth and reducing
  the oxide-yield scale reduces excessive depth. A direct local solve suggests approximately
  `f=0.923`, `s=0.546`, but the censored `56.92 s` state is not a true `60 s` response, so this may
  seed only a bounded multi-fidelity trust-region step.
- **Evidence is preserved locally** under
  `results/krueger_2024_base_calibration_r17/fine_fixed_pair_development/`: complete `audit.json`,
  exact `checkpoint.npz`, profile/trajectory plots, and the numerical-operator transition receipt.
  SHA-256 for the audit is
  `18b5c81acbe69c95ba5dd276048641719b6e1e12ce12c036611fce7429199f81`.
- **The forward plan is now explicit.** `NUMERICS_CALIBRATION_AMR_PLAN_2026-07-17.md` defines
  triangle-overlap remap gates, initial/late short multiresolution bursts, sparse-narrow-band and AMR
  gates, and multi-fidelity calibration. `VIENNAPS_CODEBASE_AUDIT_2026-07-17.md` records the upstream
  4.6.1 architecture and a truth-first conformance suite. No held-out Krüger profile was opened.

## 2026-07-17 — Krüger gas-cavity continuation and complete 5 nm development response

- **The missing topology branch is implemented without weakening the default refusal.** The public
  feature engine defaults to `topology_change_policy="refuse"`. The explicit
  `continue_gas_cavity` mode accepts only periodic gas-cavity enclosure/opening when solid count,
  every per-material component count, and open-domain breakthrough are invariant. All other
  topology changes remain structured refusals. Accepted events and resume policy transitions are
  recorded in run history.
- **Manufactured and public-engine certification pass.** A genuine, non-injected keyhole evolved
  through closure, a positive-duration sealed interval, and reopening using the public transport,
  chemistry-result, level-set, CR2, topology, and conservative-remap path. Sealed faces received
  exactly zero external neutral and ion first-hit flux; access returned after reopening; remap
  residuals were at most `1e-14`; material ledgers were exact. The complete local suite passes
  `686 passed, 1 skipped in 143.98 s`.
- **The exact real checkpoint completed to 60 s.** Resume began from the accepted
  `56.92069539342471 s` checkpoint, accepted cavity enclosure at `56.9398350487125 s`, continued
  exposed-surface evolution, reopened the cavity at `58.66367366309934 s`, and reached exactly
  `60.0 s` in `828.75 s` wall time on Vast instance `45134431`. There was no second refusal.
- **The complete development endpoint is depth `900.8570 nm` and minimum mask opening
  `39.4664 nm`.** Relative to the two base calibration targets this is `+75.8570 nm` too deep and
  `-5.5336 nm` too narrow. Maximum/top feature width is `87.0097 nm`; remaining mask thickness is
  `850.0965 nm`; asymmetry is `0.9055` cells. This complete response supersedes the censored
  `56.92 s` value for proposal construction.
- **Conservation and solver evidence remained closed:** maximum material-ledger residual `0`,
  maximum neutral-radiosity balance error `2.91e-12`, maximum bounded unresolved material volume
  `1.25e-25 m3`. Six copied remote files match local SHA-256 exactly; `audit.json` is
  `98506d4b94fbdd48ed00baee254088c6e49fbd02321dcc2d799c2b0e7e5caa30`.
- **Claim restriction remains binding.** This path is mixed-operator development evidence because
  the periodic remap repair entered at `56.482184 s` and topology policy entered at resume. It can
  diagnose and propose, but cannot validate. The receipt and full artifacts are under
  `results/krueger_2024_base_calibration_r18/mixed_operator_topology_continuation/`.
- **Multiresolution preflight refused two invalid comparisons before compute.** Nominal 20 nm
  changes the published `130 x 20 nm` periodic cell into `120 x 40 nm`; aligned restriction of
  either available late 5 nm checkpoint to 10 nm changes physical topology. The valid next
  refinement experiment is the initial 10/5 nm paired short burst, followed by an earlier
  common-topology late checkpoint or verified AMR transfer. No held-out outcome was read.

## 2026-07-17 — Bounded 10/5 nm audit and calibration decision artifacts

- **The valid initial 10/5 nm pair completed without opening held-out data.** Both grids evolved the
  same physical `130 x 20 nm` periodic cell for `0.5 s`; the 10 nm case replayed the exact 20-step
  physical-time schedule accepted by the 5 nm case. On the RTX 4090 the 5 nm and 10 nm workers took
  `73.72 s` and `16.79 s` respectively and both returned `complete`.
- **The fitted early-time observables agree under refinement.** Depth rates are
  `13.66725/13.65726 nm/s` (10 nm minus 5 nm relative difference `-0.0730%`), opening rates are
  `-4.36593/-4.42851 nm/s` (`-1.4333%` relative to 5 nm), and remaining-mask-thickness rates are
  `0.49719/0.50233 nm/s` (`+1.0329%`). Integrated ion arrival is identical to roundoff; neutral
  arrival and area-weighted material velocity are close. A face-local maximum-width/maximum-velocity
  diagnostic is grid-sensitive and is not promoted to a calibration observable.
- **This passes shallow proposal fidelity, not late-profile equivalence.** The large 10-to-5 nm
  difference at 60 s develops after necking/accessibility/topology become nonlinear. The archived
  late checkpoint cannot be restricted to 10 nm without changing topology, so no late pair was
  fabricated. Continuous minimum-opening certification and an earlier common-topology checkpoint or
  verified AMR transfer remain required.
- **Artifacts were preserved and external compute stopped.** Ten remote files copied locally match
  SHA-256 byte-for-byte. Aggregate audit SHA-256 is
  `d211d8d2c24b94e5d65400c04e0d21203b8ae150cc0760428dc08d4fbfc5eeb4`; artifacts are under
  `results/krueger_2024_multiresolution_audit/initial_gpu_r18_20260717/`. Vast instance `45134431`
  is `exited/stopped`; the other instance was already stopped.
- **The next parameter direction is explicit and R1.9 authorizes exactly one bounded check.**
  `scripts/krueger_2024_development_trust_proposal.py` binds the complete mixed-operator response,
  base targets, old 10 nm Jacobian, R1.9 protocol, campaign, and continuation receipt. The full local
  candidate is `f=0.9355590084, s=0.5183583805`; a conservative last-step-scaled trust box yields
  `f=0.9004722560, s=0.5586489666`. Both remain development proposals. The ignored artifact SHA-256
  is `e1a7a97fd79054a3f01511ac4fa41db7fb78db2668bc0385da45e76713689ad7` after binding the
  current campaign revision. R1.9 permits only the
  safeguarded 10 nm, 60 s base response-model check and closes the 10 nm candidate sequence after
  that evaluation, regardless of outcome.
- **The paper-defined mask opening is now certified before the next run.** The final reopened R1.8
  checkpoint gives `w_m=39.46998 nm` from the mask material level set. A sealed throat now reports
  exactly zero while its internal pocket width remains a separate diagnostic. The wrong-reference
  `38.89991 nm` value came from applying a 1.0 um substrate default to the translated 1.8 um mesh;
  the production interface now requires the physical reference explicitly. The primary receipt is
  `results/krueger_2024_base_calibration_r18/mixed_operator_topology_continuation/opening_certification.json`,
  SHA-256 `9bdd11a0fd71e0dbe093f6a2076b9144e6e6c99cbe6a4cbb3a6f3eed48305956`.
- **The response model has been cross-examined rather than trusted.** The old R16-to-R17 direction
  was directionally useful but had `77.5%` vector error. A same-operator affine model through
  R15/R16/R17 is well enough conditioned for one local check and places the target inside their
  parameter triangle. Both models predict the safeguarded point near the base target
  (`43.774--44.151 nm` opening, `822.930--825.189 nm` depth). The exact acceptance envelope and
  stop rules are fixed in `KRUEGER_2024_RESPONSE_MODEL_AUDIT_2026-07-17.md` before execution.
- **The large-scale program is now consolidated.**
  `UNIFIED_ENGINE_VALIDATION_EXECUTION_PROGRAM_2026-07-17.md` organizes geometry/AMR, reactor
  boundaries, Krüger blind transfer, charged Nozawa/Hwang validation, independent HAR validation,
  Vienna conformance, speed/surrogates, and partner release into gated work packages. The immediate
  queue remains narrow: exactly one earned 10 nm response check, then a fine/AMR-versus-discrepancy
  decision before any long clean endpoint. Focused campaign regressions pass `38 passed`; the most
  recent complete suite passed `703 passed, 1 skipped` before the R1.9 identifier-only edits;
  no held-out outcome was read and nothing was pushed.

## 2026-07-17 — R1.9 response check closes 10 nm candidate sequence

- **The precommitted repaired-operator run completed cleanly.** The single R1.9 point
  `f=0.9004722559883319, s=0.5586489665864749` reached exactly `60 s` in `436` accepted adaptive
  steps and `467.55 s` wall time on the RTX 4090. It recorded no topology event or terminal
  refusal, exact material-ledger closure, and maximum neutral-radiosity balance error `2.75e-12`.
  Both Vast workers are stopped.
- **The coarse-grid mask opening is excellent; the depth response model is rejected.** Final
  opening is `45.0854 nm` against `45 nm`; depth is `853.2188 nm` against `825 nm`. The local
  response predicted `44.1514/822.9295 nm`, so actual-minus-predicted is `+0.9340/+30.2893 nm`.
  Scaled merit worsened from `2.6167` to `5.6438`, giving `rho=-1.3955`. The precommitted decision
  is `reject_response_model`; R1.9 closes the 10 nm candidate sequence and no second candidate is
  authorized.
- **Artifacts are local and byte-verified.** Remote audit SHA-256 is
  `b06c07aae21d26b52dd111866c0a381f8f61ad044e330e5ed826d88f3140765f`; launch manifest SHA-256
  is `821301496f133a16498e5051fe491827e5bdf2fd430f562cbf8e0ecf2d9e7aa5`; evaluation file SHA-256
  is `15047f554ca9d96a6efab4973e1433f552913a3656850ba08c26244b2c87dd77`. Evidence, profile, and
  comparison plots are under `results/krueger_2024_r19_response_check/`.
- **The failure mechanism is now evidence-backed.** The candidate begins shallower, reaches a
  `-4.952 nm` depth difference at `11.77 s`, crosses R17 at `14.89 s`, and ends `+15.314 nm`
  deeper. This is nonlinear geometry/access/state feedback, not a constant scale response. The
  stopped R17 worker also preserved the old source: chemistry, material mechanism, radiosity, and
  boundary transport are byte-identical, while `feature_step_3d.py` changed from `a1e8d02f...` to
  `61dbef70...` to repair periodic state remap and surface-distance certification. The old endpoint
  Jacobian therefore belonged to a pre-repair operator epoch.
- **The calibration data firewall is now executable.** Base Krüger transport uses the narrow
  `load_krueger_2024_base_boundary_fluxes` path and succeeds when held-out transfer observations are
  physically absent. Calibration reveals and held-out predictions now checksum-bind the same
  operator epoch. After the operator-epoch scoring contract and base-boundary data-firewall test,
  the full local suite passes `709 passed, 1 skipped`.

### Bounded frozen-checkpoint decomposition after R1.9

- **The no-evolution access screen completed under a real wall budget.** A first exact-CPU attempt
  was stopped when its post-return timer proved nonbinding. The replacement uses hard per-cell and
  total deadlines, atomic partial output, and an explicitly diagnostic reduced quadrature. The full
  R17/R19 checkpoint-by-parameter matrix completed in `49.37 s`; all four cells close the material
  ledger exactly and keep radiosity error below `2.62e-12`.
- **Accumulated geometry/state dominates the mask response.** At fixed mechanism parameters the net
  mask-rate magnitude drops by about `55--59%` between the R17 and R19 checkpoints. The parameter
  change itself moves that rate by only `7.45%` at R17 and `15.70%` at R19. The R19 state receives
  about `5.98%` more ion rate on SiO2 and `2.87%` less on the mask. Receipt SHA-256:
  `9eba876f153c9acc1881edc6d4894f6626682f8b251bac9b8bca975ced4e4237`.
- **`dt=0` cannot diagnose the oxide-yield direction.** A selected q3 evaluation also returns zero
  instantaneous SiO2 removal despite nonzero SiO2 ion arrival. Every saved oxide face is film-covered;
  the conservative mechanism exposes oxide only while advancing a positive duration. The next cheap
  test is therefore a frozen-geometry chemistry micro-step, bounded by the already-authorized next
  timestep, displacement, reaction-probability drift, and N/2N chemistry refinement—not a new profile.
- **The AMR path begins with a behavior-neutral seam.** `FeatureGeometry3D` and face ownership now
  live in dependency-neutral geometry state; a read-only uniform backend adds SI/mesh sampling,
  exact surface extraction, periodic metadata, and deterministic fingerprints without routing the
  engine through new storage. Backend/feature/checkpoint focused verification passes `86` tests.
- **Adaptive branch:** do not tune chemistry against this 10 nm miss. The unclosed late grid error
  is larger than the `28.2 nm` depth miss. Next are frozen no-evolution checkpoint diagnostics and
  sparse-band/AMR authority tests; held-out oxygen/power outcomes remain sealed.

### Positive-horizon chemistry exposes the next common-engine coupling error

- **The bounded positive-time diagnostic passed transport and conservation, then stopped at its
  predeclared chemistry gate.** Two q3 checkpoint transports completed in `62.54 s` total with exact
  material ledgers and maximum neutral-radiosity balance error `1.088e-12`. At
  `dt_next/16 = 0.00779536 s`, oxide removal is finite (`2796.75695` formula units,
  `1.26837e-25 m3`), N/2N removal agrees to `0.1312%`, and gross predicted motion is only
  `0.0210 dx`. The immutable audit is
  `results/krueger_2024_r19_response_check/frozen_surface_chemistry/audit.json`, SHA-256
  `694319374b4bd0067f88ca73678561120fed574e8c2fefac3fdc873f91bfa278`.
- **The failed quantity is surface-state/radiosity splitting, not a material ledger or geometry
  failure.** Flux-weighted reaction-probability drift is `0.2465%`, but the maximum local change is
  `5.961%`, above the declared `1%` cap. The production path currently estimates diffuse radiosity
  once and advances chemistry across a roughly `0.1247 s` profile step. Since polymer state controls
  neutral sticking, this can freeze the wrong redistributed neutral flux while the film evolves.
- **The lasting repair has been scoped without rewriting the engine.** On fixed geometry, direct
  transport and the diffuse form-factor operator are geometry/boundary quantities and may be cached;
  radiosity must be re-solved as state changes. `SURFACE_STATE_RADIOSITY_COUPLING_PLAN_2026-07-17.md`
  defines stable SC0--SC5 tasks. The required controller uses embedded step-doubling—one full
  chemistry step versus two half steps with a midpoint radiosity re-solve—and accepts the fine path
  only after every state increment, exchange inventory, per-face recession/growth, balance,
  conservation, and displacement gate passes. Local probability change remains a safety cap, not an
  error estimator.
- **Literature and independent review agree on the structure.** Established feature models describe
  coverage-dependent sticking and re-emitted neutral flux as a recursive problem and solve it
  iteratively or under a demonstrated timescale separation. Cached geometry is therefore sound;
  frozen state-dependent radiosity is not. A bounded diagnostic is authorized. Engine wiring, a long
  profile, GPU use, and held-out reveal remain held until cache identity, manufactured DAE, coupling-
  step refinement, and nested form-factor refinement pass.

### Coupled frozen-state result and first uniform-backend routing

- **The production cache identity gate passes exactly.** A first attempt correctly refused because
  disabling radiosity had also disabled periodic first-hit wrapping, losing `25--76%` of integrated
  incident rate. The repaired diagnostic captures direct transport inside the exact production
  periodic call. Its 1,833-face production and cached neutral fluxes are byte/hash identical, and
  the regenerated q3 receipt now records transport seed `241` and radiosity seed `10241`.
- **Two positive horizons give a stable physical direction.** At `7.795/15.591 ms`, the nominal
  R17 oxide removals are `2794.390/5405.351` formula units and R19 gives
  `2753.804/5325.748`, ratios `0.985476/0.985273`. All embedded, probability, ledger, radiosity,
  displacement, and tolerance-halving gates pass; gross motion is `0.02104/0.03991 dx` and
  tight/nominal disagreement is below `0.01%`. The `31.181 ms` cell hit the predeclared `480 s`
  total budget, classified as implementation/controller evidence only. Audit SHA-256:
  `0f337c8516542f4603926cc8b1d28d156abba19859d68ac2d6b6b6231a97706e`.
- **The first backend seam is now on the live uniform path without changing numerical values.** Both
  old and post-advection surface extraction plus face ownership route through
  `UniformFeatureGeometryBackend3D.extract_surface()`, then copy back into the legacy writable
  arrays. The legacy active-surface fingerprint remains unchanged; backend fingerprints are not
  substituted. Exact value/dtype/order/writeability tests and the feature/checkpoint/replay cluster
  pass `90` tests. This earns only a read-side seam; evolution remains dense and no AMR claim is made.
- **Periodic first-hit transport is no longer coupled to radiosity enablement.** The common feature
  API now has an explicit, backward-compatible `ballistic_periodic_lateral` declaration. A focused
  parity test proves periodic field-free first-hit transport runs without diffuse radiosity; an
  inconsistent request for periodic radiosity with nonperiodic first hit refuses. This closes the
  API defect that the first cache attempt exposed without changing old callers.
- **The bounded evidence is human-auditable.** The static figure
  `results/krueger_2024_r19_response_check/krueger_bounded_endpoint_diagnostics.png` is generated
  only from the three checksummed base-only JSON audits. Its sidecar records every plotted value and
  source hash; it labels the `31.181 ms` timeout as implementation/controller evidence only.
- **The reusable SC1 operator passes its standalone manufactured contract.** New module
  `src/petch/surface_radiosity_coupling_3d.py` binds immutable direct flux, form factors, face areas,
  material/species routing, solver controls, and an exact operator identity. Its embedded controller
  compares one full step with two half steps and a midpoint radiosity solve, accepts only the fine
  path, checks every state increment, all four exchange inventories, and per-face recession/growth,
  and keeps rejected trials transactionally invisible. Open, constant-probability, reflective-cavity,
  nonmonotone, replay, identity, minimum-step, incomplete-routing, product-refusal, and zero-duration
  tests pass; `22` SC1 plus adjacent radiosity/router tests complete in `5.01 s`. A `4x` accepted-step
  refinement contracts the manufactured cavity error by about `5.15x`. The module is not wired into
  profile evolution, and no combined neutral-to-material stoichiometric ledger is claimed because the
  mechanism API does not expose that mapping.
- **WP-AMR2 now has a shared exact surface primitive, still isolated from the live engine.**
  `src/petch/surface_mesh_3d.py` adds immutable `TriangleSurface3D` geometry with order-sensitive
  fingerprints, exact plane/edge/vertex nearest points, whole-triangle periodic images, strict
  material restriction, deterministic ties, and a certified bounded candidate path. Eleven focused
  tests pass, including bounded-versus-brute-force equivalence. It is not yet a BVH, does not yet
  compute overlap weights, and is not wired into state or charge remap; no AMR/performance claim is
  made.

### Full-state form-factor refinement preserves the global direction but blocks local profile authority

- **The audit receipt now contains the quantities its internal controller was already checking.**
  Every nominal/tight R17/R19 integration persists all seven final surface-state fields, every named
  per-face removed/outgoing/unresolved/deposited inventory, per-face recession/growth, and exact
  deterministic hashes. The old default seedless q3 certificate was replaced by the already-earned
  seed-241/10241 receipt, preventing an otherwise valid unattended 8-ray run from refusing at startup.
  Eleven focused audit tests and seven independent cross-ray receipt tests pass.
- **Three bounded levels separate global convergence from local ray-level sensitivity.** The 8/16/32-ray
  `dt_next/16 = 7.795 ms` runs all pass their coupling-step, probability, ledger, radiosity,
  displacement, and tight/nominal gates in `158.8/117.9/183.5 s`. From 16 -> 32 rays, integrated oxide
  changes by at most `0.024%` and the paired R19-R17 effect by only `0.0066%`; R19 remains about
  `1.45%` lower. But the maximum local fluorocarbon deposition/growth normalized L1 changes by about
  `18%` and the worst face by about `58%`. The global causal direction is earned; local moving-profile
  authority is not.
- **The response is a categorical view-factor sensitivity problem, not solver roulette.** One source
  triangle currently distributes 8--32 centroid-launched hard rays among 1,833 possible
  targets/escape. Spatially integrated oxide averages the observed level difference; local film
  growth does not. This is not yet a stochastic error bar: independent scrambles are still required,
  and the finite-triangle operator must integrate source area rather than use a centroid point source.
  No 60 s endpoint, held-out reveal, or SC3 wiring follows.
  `FORM_FACTOR_ACCURACY_AUDIT_2026-07-17.md` defines the lasting FF0--FF5 path:
  correctness/area-sampling/visibility/reciprocity gates, independent replicated scrambled-Sobol uncertainty,
  prospective integrated and fixed-physical-patch tolerances, adaptive allocation, and hierarchical
  refinement only if needed.

### Finite-area form factors and physical-patch scoring replace two mesh-sensitive instruments

- **The emitting-surface measure is now correct in the production-target mode.** The old Krueger
  receipts launched all diffuse rays from each triangle centroid; increasing angular rays could not
  remove that point-source bias. `triangle_area` now uses one nested four-dimensional scrambled Sobol
  rule for uniform source-triangle position plus cosine direction, while `legacy_centroid` preserves
  exact historical replay. On a close half-receiver, the old coarse centroid result differs from the
  finite-area result by more than `10%`; the finite-area result agrees with an `8x` surface-refined
  centroid reference within `1%`.
- **Independent scrambles are a first-class numerical object, not an afterthought.** The new immutable
  replicated form-factor controller requires distinct equal-level seeds, builds an exactly
  row-closing mean operator, hashes the mesh/domain/sampling/caller identity, measures raw area
  reciprocity without altering the operator, and reports Student-t downstream radiosity intervals.
  Its API minimum is four replicates; the first checkpoint coverage audit still requires eight.
- **Visibility ambiguity is explicit.** A float64 cell-by-cell periodic reference reports hard hit,
  open-top escape, or wrap-budget exhaustion. Manufactured Warp/reference parity passes exact and
  near shared edges, the preserved production grazing case, one periodic wrap, and open escape. The
  fast API's former `-1` ambiguity is closed by replaying every miss and gas-normal-invalid hit in
  float64; exhaustion refuses rather than becoming escape. Apparently valid fast hits still require
  one bounded full-event checkpoint parity audit before production promotion.
- **Physical patches no longer inherit centroid/remeshing noise.** The new convergence instrument
  clips every triangle against fixed SI patch boxes and transfers exact overlap area, split by
  material and oriented surface class. It conserves face and patch area, reproduces constant-field
  integrals under `1x -> 4x` retriangulation at two scales, uses a mixed absolute/relative norm near
  zero, and keeps worst-face values separate. The triangle operator is unchanged; patches score it,
  never smooth it.
- **Verification remains bounded.** The surface/form-factor/feature/coupling cluster passes
  `142 passed, 1 skipped`; charging plus the new patch utility passes `41` tests. No long profile,
  GPU worker, endpoint, held-out reveal, or chemistry retune was run. SC2 remains partial pending one
  real-checkpoint event-parity receipt and one replicated patch-scored shortest-horizon audit.

### Exact-periodic visibility closes FF1; replicated physical-patch scoring is now generic

- **The legacy periodic form-factor tracer failed its real-geometry cross-examination.** On the
  sealed base-only R1.9 checkpoint, eight finite-area rays from each of `1,833` triangles produced
  `14,664` events. The historical `_apply_bc` path disagreed with the cell-by-cell float64 authority
  on `991` events (`19` classifications and `972` target faces) and shifted `5.6463%` of the
  area-weighted source-row measure. Its seam inset accumulated a nonphysical translation on long
  grazing paths; this is a geometry-operator defect, not chemistry or sampling noise.
- **The repaired common path matches float64 exactly.** The candidate now advances one periodic cell
  at a time, maps seams to exact `0/domain` coordinates, and uses a tolerance only for float32
  distance comparisons. Finite-area source points use an owning-triangle interior limit; `911` of
  `14,664` points moved inward and `224` ultra-thin source events used the triangle centroid rather
  than being discarded. Fast and float64 receipts then match on every target, termination class,
  and wrap count, including rays with `320` wraps. There are zero false hits/escapes, wrap
  exhaustions, or solid-facing refusals. Audit SHA-256:
  `3c233b3c8861ece1c82dcf1edc4712ed8109d6cc1ec50bc9367371e78deb8d56`; identical event SHA-256:
  `6e357c6793dbf0d3e0300d785afb127101198c5846512caed48af0bb808054ad`.
- **The operator is promoted once, in the shared estimator.** `cellwise_certified` is now the default
  visibility mode for the ordinary and replicated form-factor APIs. Legacy float32, selective replay,
  and full-float64 modes remain named audit/replay options. The full local regression passes
  `820 passed, 1 skipped in 155.46 s`; no profile, held-out observation, parameter fit, GPU worker,
  or endpoint was run.
- **Replicate intervals now score physical fields, not triangle IDs.** The shared patch module maps
  an authority field and at least four independent replicate fields through one exact triangle/box
  overlap operator, reports Student-t intervals, separates mean-operator/nonlinear bias from sampling
  uncertainty, and gates their conservative sum at two or more fixed SI scales. Focused patch,
  form-factor, and boundary tests pass `56 passed, 1 skipped`. FF1 is complete; FF2/FF3 remain open
  only for the bounded eight-scramble, patch-scored shortest-horizon checkpoint receipt.

### The first replicated closure attempt stopped on compute authority, not physics

- **The physical-patch sliver instrument is repaired on the real mesh.** Boundary-plane coordinates
  retain input precision; only whole coplanar coordinates snap, intersected cells on both sides are
  enumerated, and roundoff-only pieces are removed before the retained pieces are renormalized to the
  authoritative triangle area. Exact-boundary, +/-1 ULP, float32-quantized, no-op, and inventory
  tests pass. On the `1,833`-face checkpoint, the smallest numerical patch area moved from about
  `1e-32 m2` to `6.11e-24 m2`; 20/40 nm scoring produces `272/151` physical patches and closes
  inventory to roundoff.
- **The Stage-A harness is scientifically ready but the direct preflight exceeded its first budget.**
  An independent review cleared the paired nested confidence gates, exact operator/artifact identity,
  physical normalizations, two-scale receipts, raw reciprocity, and paired R19-R17 sign contract.
  Two ordinary CLI launches then stopped at exactly the declared `60 s` zero-duration
  direct-transport limit, including one with warm JIT caches. No form-factor replicate ran and no
  scientific conclusion is permitted by either artifact. The active bounded repair is a persistent
  hash-bound direct-transport artifact with a one-time `<=120 s` construction allowance inside the
  unchanged `300 s` Stage-A ceiling; no ray escalation, held-out access, Stage B, or long profile is
  authorized.

### Replicated Stage A completes and localizes the remaining surface-transport error

- **The common hard-visibility operator now recovers new shared-edge misses inline.** The first
  replicate exposed a float32 event that looked like more than `1,024` periodic wraps; float64 hit the
  exact triangle after one wrap. `cellwise_certified` now replays only exhausted/solid-facing events
  with exact triangle data and retains a separately declared `4,096`-wrap exact emergency horizon.
  Across the complete 8/16-ray, eight-scramble screen, the operator replayed one/six events by level,
  recovered six shared-edge hits in total, and certified one real open-top escape after `1,077`
  wraps. Focused visibility/boundary tests pass `77 passed, 1 skipped`.
- **Boundary/mechanism species routing is unified rather than benchmark patched.** Mechanisms may
  advertise channels absent from a boundary deck; those zero-measure extras are projected out.
  Every supplied boundary species must still be covered, and declared-inert `C3F4` is now routed as
  an explicit zero probability instead of being indistinguishable from an omitted channel. The
  chemistry/router/coupling cluster passes `68` focused tests.
- **The sealed Stage-A receipt is a bounded precision hold.** It completed in `295.408 s`; the exact
  direct capture took `104.390 s` and is now stored in a checkpoint/operator/config/seed-bound
  artifact. Row closure, radiosity balance, material ledgers, nested samples, physical patch-scale
  eligibility, integrated SiO2/mask observables, and the paired R19-R17 direction pass. R19 remains
  lower by `5.35367e-14 m` at level 16 with a `1.22607e-17 m` paired 95% half-width. Local 20/40 nm
  flux/film receipts fail, with the worst level-16 patch about `48.5x` tolerance, and the tested
  horizon predicts `0.729--0.733 dx` motion versus the `0.05 dx` frozen limit. Stage B and global
  32-ray escalation remain held; the next task is a shorter frozen horizon and selected-source row
  allocation from the persisted artifacts. Audit SHA-256:
  `4e28066dcd0257feffac8983a462045b9fb5686eeaa762b3dce096422b3732b0`.
- **WP1-AMR2a's spatial index passes but is not silently wired.** Exact periodic nearest-triangle
  face/image/tie queries match brute force. The shared mixed-state transfer runs at `0.877x` the
  legacy time on `2,048` faces (gate `<=1.5x`); `30` focused geometry/remap tests pass. Live
  `feature_step_3d` promotion waits on one explicit regular-grid Kth-neighbor tie-authority test,
  because the legacy cKDTree chose such ties by traversal order while the new authority chooses the
  lower physical face id.

### Physical patch authority and the next estimator rung are now explicit

- **The patch statistic no longer lets an almost-empty Cartesian corner define a local mean.** The
  shared instrument records dominant-axis projected support and the represented nominal footprint,
  including the 20 nm periodic-y fundamental cell under a 40 nm requested patch. Every patch remains
  in the fixed-footprint integrated gate; a mean is gated only above the predeclared, geometry-wide
  `0.10` support fraction, and the all-patch mean remains a diagnostic. An all-ineligible comparison
  refuses. Stage-A schema is now v2, passes this contract through replicated, nested, and paired
  paths, serializes `0.05/0.075/0.10/0.25/0.50/0.75` sensitivity, and structurally prevents an old or
  incomplete receipt from authorizing Stage B. The focused controller/patch/charging cluster passes
  `76` tests.
- **The threshold has a durable geometry-only receipt, not a benchmark-score justification.** On the
  sealed 1,833-face base mesh, primary `0.10` support excludes `76/272` 20 nm patches but only
  `0.50425%` of actual area (`0.41677%` projected); at 40 nm it excludes `55/151` patches,
  `1.40056%` actual (`1.08456%` projected). Excluded area is still integrated and gated. The artifact
  read no transport, chemistry, profile, held-out observation, or GPU result:
  `results/curated/krueger_2024_physical_patch_support_sensitivity_2026-07-17.json`, SHA-256
  `6b61b10c74b7baf7280a7dfa03ec063d7a02ae6b19b46e6f832f3ebaec3ce148`.
- **The selected-source hypothesis failed its own cap.** Exact signed adjoint row decomposition closes
  to `2.96e-15`. All-patch scoring misleadingly placed 90% of its ranking mass in `234/1833` source
  rows (`12.77%`), while physically supported patches require `758/1833` (`41.35%`). The allowed
  458-row/25% plan captures only `78.99%`, so the persisted allocator reports
  `diffuse_source_error_blocker` and is historical diagnostic evidence only. Artifact SHA-256:
  `317cd1770b4abba8e36cc8d8ac98e137491e6a6414a06f055873a071824fe6f4`.
- **A reusable selective nested-RQMC primitive passes but is intentionally unpromoted.** It appends
  the identical Sobol indices to explicit source rows, merges integer hit/escape counts, leaves
  untouched rows bitwise identical, propagates certified-visibility refusals, and binds mesh,
  operator, seed, levels, and rows in an immutable receipt. The boundary/control/visibility cluster
  passes `75 passed, 1 skipped`. On a 200-face manufactured case it reduces the trace portion to
  `0.328x`, but tracing is already seconds while response scoring dominates; no Krueger wiring or run
  follows.
- **The simplest next scientific rung is uniform nested 16->32 at a verified short horizon.** A
  source-method audit found that all-row refinement adds only `234,624` events across eight scrambles,
  preserves paired uncertainty and equal reciprocal precision, and avoids same-data row-selection
  bias. Reciprocity/closure projection is retained only as a later manufactured control-variate
  experiment: fixed affine projection can remain unbiased, whereas nonnegative smoothing or
  constrained likelihood introduces finite-sample bias and its iid multinomial likelihood does not
  describe scrambled Sobol dependence. No Stage A, Stage B, profile, held-out, charging, or GPU run
  was launched by these increments.
- **Exact planar overlap remap is now a fast standalone authority, not yet a live-engine claim.** A
  new immutable sparse operator uses the deterministic periodic candidate index and exact float64
  triangle clipping, combines periodic images by physical face pair, never crosses materials or
  orientation, closes retained/removed/newly exposed extensive ledgers, and keeps intensive transfer
  convex with an explicit fill requirement for new area. Identity, retriangulation, subdivision,
  shifted-periodic, mixed-material, linear-convergence, and refusal gates pass (`33` focused tests).
  On a 2,048-face planar case it runs in `0.055661 s` versus `0.228537 s` for indexed KNN (`0.244x`).
  Curved, nonparallel, folded, and large-normal-motion surfaces refuse; no `feature_step_3d`, charging,
  AMR, or profile path was rewired. The next remap task is orientation-local patch grouping and
  corner/topology ownership, not pretending one plane is a trench.

### Matched ViennaLS topology conformance separates a bad reverse test from a real first-order lag

- **The original reverse-time comparison was ill-posed.** Reversing each engine from its own first
  grid-snapped closure compared different subcell cap thicknesses and produced nonconvergent petch
  reopening times (`0.125--3.125 s`). T1 now contains two independent analytic initial-value
  problems: close the same open rounded keyhole, then reopen the same prescribed 100 nm solid cap.
  The production directional hard-visibility audit is unchanged.
- **The pinned black-box receipt passes without claiming identity.** The adapter verifies clean exact
  ViennaLS/ViennaHRLE/ViennaCore/ViennaPS revisions, compiles an original API-only probe, applies hard
  timeouts, and runs 50/25/12.5 nm grids. At 12.5 nm, ViennaLS closes at `1.9375 s`, petch at
  `2.5625 s`, and the analytic event is `2.0 s`; petch's `0.5625 s` lag is exactly its declared one
  cell-crossing plus one-checkpoint localization budget. Both engines reopen the prescribed cap at
  the analytic `1.0 s` on all grids. The paired run passes in `39.47 s`; petch's conservative state
  ledger remains roundoff-exact. The event-time result is numerical conformance, not chemistry
  validation, and the visible petch closure bias remains a target for sparse/AMR comparison.
- **Verification is bounded and reproducible.** The topology/remap cluster passes `94` tests. Receipt
  and plot: `results/topology_petch_viennals_t1/audit.json` and
  `results/topology_petch_viennals_t1/petch_viennals_topology.png`. Profile-distance, multi-material,
  and periodic-translation extensions remain open; no Krueger endpoint or held-out data was run.

### Unified CUDA parity passes and the real Krüger step bottlenecks are localized

- **The device declaration is now unified and self-refusing.** The first bounded CUDA profile traced
  particles on `cuda:0` while the level-set module followed an unset `PETCH_DEVICE` onto CPU. The
  subprocess worker now binds both declarations to one input, and the profiler refuses a mixed
  device receipt. A corrected RTX 4090 run used one positive warm-up and one profiled step; it
  conserves material exactly, closes radiosity to `1.90e-12`, and has no topology event.
- **The correction changes no measured answer.** At identical seed, operator, grid, and 0.05 s
  physical duration, mixed versus unified CUDA differs by `0 nm` in depth, `8.59e-9 nm` in opening,
  and at most `1.02e-7 nm` over every common-unit endpoint measurement. This passes the explicit
  `1e-6 nm` bounded parity check.
- **The steady cost is not level-set reinitialization.** One warmed 5 nm step takes `7.780 s`:
  chemistry/material routing is `29.5%`, ballistic transport `28.8%`, diffuse exchange `15.3%`,
  legacy remap `13.7%`, other work `12.0%`, and redistance only `0.8%`. The Krüger worker inherited
  `legacy_knn`; the newer transfer operators were not silently credited. The naive fixed-step 60 s
  arithmetic is `5.19 h`, explicitly not an endpoint benchmark.
- **The paired 10/5 evidence supports calibration modeling, not local shape authority.** Initial
  depth/opening/mask-thickness rates differ by `0.073%/1.433%/1.033%`; maximum and top feature-width
  rates differ by `71.235%`. The latter is plotted rather than averaged away. No endpoint or held-out
  data was run. Report and machine receipt:
  `KRUEGER_2024_CUDA_PROFILE_REPORT_2026-07-17.md` and
  `results/krueger_2024_cuda_profile_summary/audit.json`. The postprocessor/profiler/multiresolution
  cluster passes `13` tests.

### Multi-fidelity model management is implemented and refuses stale Krüger evidence

- **The reusable controller keeps high fidelity in the loop.** It fits rank/condition-checked affine
  low responses, anchors high-minus-low discrepancy exactly at the current center, distinguishes a
  center-only receipt from an identified first-order correction, proposes physical-box trust steps,
  and accepts/grows or rejects/shrinks only from actual/predicted high-fidelity merit reduction.
  Synthetic value/gradient, rank-refusal, bounds, and trust-ratio tests pass.
- **The old grid-correction shortcut is not silently reused.** R1.9's proposed improvement produced
  `rho=-1.3955`, and its `feature_step_3d.py` hash differs from the current operator. The readiness
  scorer therefore treats every old endpoint derivative as historical, even though its parameter
  location remains the fixed location for the first clean fine anchor.
- **Short refinement evidence is classified correctly.** The `0.5 s` depth/opening rate differences
  (`0.073%/1.433%`) size shallow numerical error; the `71.235%` local width difference blocks shape
  authority. None is multiplied into a 60 s endpoint discrepancy. The real readiness status is
  `blocked_before_parameter_proposal` on five named blockers: remap backend, operator epoch, fine
  anchor, paired endpoint discrepancy, and current first-order response.
- **The next expensive action is one anchor, not an optimizer sweep.** After a bounded remap-backend
  comparison, run one clean current-epoch uniform-5-nm or certified-AMR base at the fixed R1.9 pair.
  A pass freezes calibration; a miss must first earn a current response direction before one
  safeguarded correction. No simulation, GPU, new pair, or held-out read occurred. Report and
  receipt: `KRUEGER_2024_MULTIFIDELITY_READINESS_2026-07-17.md` and
  `results/krueger_2024_multifidelity_readiness/audit.json`; the complete related cluster passes
  `31` tests.

### Real Krueger state-transfer gate selects common refinement at 10 nm

- **The remap choice is now an explicit part of the unified operator.** The Krueger multiresolution
  worker and CUDA profiler accept and record `legacy_knn`, `indexed_knn`, `partitioned_overlap`, or
  `common_refinement`; no campaign can silently inherit the old default.
- **A bounded same-state comparison ran exactly two 0.025 s steps per backend.** Every completing
  case began from identical analytic geometry/state and produced an identical first-step geometry,
  so the second step isolates the first physical consequence of remapped chemistry. Legacy and
  indexed KNN completed; planar overlap issued a certified nonparallel-face refusal; common
  refinement completed with no topology event, zero material-exchange residual, and maximum relative
  remap residual `1.04e-15`.
- **The global response is stable while the local shape remains appropriately diagnostic.** Relative
  to legacy after step two, common refinement changes depth by `+1.4729e-5 nm` (`21.6 ppm`), opening
  by `+2.1564e-4 nm` (`2.41 ppm`), remaining mask thickness by `-2.2756e-4 nm` (`-0.268 ppm`), and
  top width by `-0.19698 nm` (`-0.227%`). This promotes common refinement as the production candidate;
  it does not yet freeze it or authorize a long endpoint.
- **The next confirmation remains bounded.** A fresh RTX 5090 rental was destroyed after spending
  nearly eight minutes pulling its container without starting the simulation. No GPU remains
  running. The pending action is one paired 5 nm common-versus-indexed confirmation on a cached/ready
  GPU; local work proceeds to standalone sparse-band classification and halo packing rather than
  blind rerental. The remap cluster passes `93` tests. Report, plot, and machine receipts:
  `KRUEGER_2024_REMAP_BACKEND_AUDIT_2026-07-17.md` and
  `results/krueger_2024_remap_backend_audit/`.

### Fixed-dx sparse volume fails its Krueger cost gate before live integration

- **The standalone sparse authority is real and bounded.** It stores deterministic active/far blocks,
  partial boundaries, periodic neighbors, exact interface/material-junction payloads, and either
  replicated or sorted unique-node halos. All 72 analytic 10/5 nm shallow/mid/deep cases preserve the
  declared band and material ownership; the 5 nm endpoint canonicalization is only `5.724e-10` cells.
- **The geometry-specific performance hypothesis is rejected.** At the required safe 8-cell band and
  a 900 nm-deep 5 nm trench, best core storage is `1.415x` smaller, indexed halo storage is `1.148x`
  smaller, simultaneous storage is `0.634x`, and even the zero-overhead unique-node work ceiling is
  `1.921x`. These miss the predeclared 3x memory and 2x work gates.
- **No sparse solver rabbit hole follows.** The dense state is only 75,735 nodes/2.42 MB, while the
  live 5 nm profile spends 0.8% in redistance and most time in surface chemistry, ballistic transport,
  diffuse exchange, and remap. Fixed-dx sparse evolution is not wired; Krueger volume AMR is held
  pending a cost model. The earned next scaling work is symmetry-aware surface/transport reduction,
  with full 3-D retained for holes/twisting/asymmetry. Report, plot, and receipt:
  `BLOCK_SPARSE_KRUEGER_NO_GO_2026-07-17.md` and
  `results/block_levelset_manufactured_audit_3d/`.

### Surface-work preflight preserves full 3-D and removes duplicate event work

- **The existing exact extrusion path passes physics but does not earn a speed promotion.** On the
  zero-motion 10 nm Krueger base surface it certifies 724 triangles as 181 physical cross-section
  segments over two identical strips and closes diffuse balance to `5.73e-15`. Its local wall time is
  `30.262 s`, versus `28.085 s` for the current eight-ray 3-D path, so production remains unchanged.
  This bounded comparison read no held-out data and evolved no profile.
- **The accepted optimization changes bookkeeping, not the operator.** Already-validated immutable
  energetic events now survive active/material face selection without repeating direction and
  finiteness validation, and built-in frozen yield laws memoize the exact per-face reduction. Cache
  state is excluded from dataclass equality and serialization. A one-million-event manufactured
  comparison is bitwise identical, with `1.50x` subset and `4.03x` repeated-yield speedups.
- **The whole engine remains green.** The targeted surface/material/feature suite passes 119 tests;
  the repository-wide suite passes 955 with one expected skip. The next expensive action remains the
  bounded paired 5 nm remap confirmation, followed by one current-epoch fine base anchor—no endpoint
  sweep and no held-out reveal yet. Report and receipt:
  `KRUEGER_2024_SURFACE_WORK_PREFLIGHT_2026-07-17.md` and
  `results/krueger_2024_surface_work_preflight/audit.json`.

### Paired 5 nm CUDA gate selects and seals common-refinement surface transfer

- **The CPU slowness was measured rather than tolerated.** Identical local 5 nm indexed and common
  workers each reached the predeclared `600 s` timeout without completing a receipt. Both preserved
  Vast instances were destroyed. A fresh RTX 4090/32-effective-core instance (`45195298`) ran the
  same checksum-bound source archive (`12aee41`, archive SHA-256
  `21889b138f333c618f6fceb189d3618a84c4425b03f2630e9eed41879f30a795`).
- **The GPU comparison is genuinely paired.** Indexed and common refinement start from identical
  analytic 5 nm geometry and surface state, use seed `241`, and produce byte-identical first-step
  geometry before remap feedback can enter. Both complete two `0.025 s` steps without topology
  events, close the material ledger exactly, and keep maximum relative remap residual below
  `6.61e-16`.
- **Common refinement wins on physics and cost.** At step two, common minus indexed changes depth by
  `-4.48e-6 nm` (`-6.55 ppm`), opening by `+1.96e-4 nm` (`+2.19 ppm`), remaining mask thickness by
  `-1.60e-4 nm`, and maximum width by `-0.13485 nm` (`-0.153%`, 2.70% of one cell). Worker wall time
  is `14.56 s` versus `21.18 s`, a `31.3%` reduction. Common refinement is selected for the clean
  anchor; held-out profiles remain unread.
- **The selection now reaches the real campaign.** Commit `12aee41` records the backend in the
  pilot fingerprint, refuses backend changes on resume, requires one nonlegacy backend across the
  10/5 endpoints, and replays the frozen value into every blind transfer case. The aggregate bug
  that could treat a lone completing backend as paired is fixed; readiness accepts only two complete
  conservative 5 nm receipts. The remap/readiness/report cluster passes `27` focused tests, and the
  repository-wide suite passes `962` with one expected skip in `181.39 s`.
- **Next:** regenerate the checksum-bound readiness receipt with the paired audit, then run one clean
  current-epoch uniform-5-nm base anchor at the fixed R1.9 pair. No optimizer sweep and no held-out
  reveal precede that endpoint. Report, static plot, raw receipt, and environment manifest:
  `KRUEGER_2024_REMAP_SELECTION_2026-07-17.md` and
  `results/krueger_2024_remap_backend_audit_5nm_gpu/`.

### The 5 nm authority run exposes and closes a marching-cubes topology defect

- **The sealed `a15634f` run is a failed receipt, not an endpoint.** Seed `241`, config hash
  `32ca232d849e6e34f5378ddb2d5987c3ae826f8691b6ff3b8cb329fa366321fa`, hard visibility,
  common-refinement remap, and CUDA face-gather transport reached step 168 / 11.376815 s with
  155.697 nm depth and 76.709 nm opening. The material ledger remained exact and no topology event
  occurred, but one of 36,856 rays refused after a 1,024-wrap solid-facing exact replay. No held-out
  profile was read and the checkpoint is not resumable across the operator repair.
- **Exact replay found a mesh hole, not an undersized trajectory horizon.** The ray entered solid
  after 15.02 nm of path and reached the far surface with incoming cosine `-0.2699`. The old
  float32-scaled area cutoff had deleted 29 real sliver triangles: the raw 4,636-face contour had
  zero unmatched interior edges, while the filtered 4,607-face mesh had 56. Keeping every
  positive-area face makes the same ray hit its near-side face with cosine `+0.7698` and zero wraps.
- **Revision `418fb28` makes topology part of the operator contract.** CPU and GPU extraction retain
  all representable positive-area faces, cancel exact duplicate facets by orientation, and refuse
  unmatched interior or non-manifold edges. Manufactured sliver, missing-face, duplicate-facet,
  multi-material, transport, scattering, and Nozawa regressions pass; the repository-wide suite is
  green with 968 passed and one expected skip. Forensic receipt:
  `results/krueger_2024_base_authority_5nm_a15634f_topology_refusal/`.
- **Next:** restart the base authority case from zero on a clean descendant of `418fb28`; passing the
  old step-168 state is a required canary, not permission to reveal or tune against held-out data.

### Topology-complete meshes require separate quadrature support, and the real checkpoint passes

- **Retaining watertight slivers exposed a distinct source-discretization limit.** A clean descendant
  of `418fb28` reached step 113 / `8.961859 s` before 16 neutral quadrature rays launched from the two
  smallest retained faces into solid. Those faces are required for hard-visibility watertightness,
  but their minimum altitude is smaller than the declared launch inset and they cannot independently
  support a three-point source rule. Deleting them would reopen the mesh; moving every ray would
  change resolved transport.
- **Revision `f7f1df4` separates visibility geometry from quadrature support.** Every positive-area
  face remains in the collision mesh, material ledger, and state inventory. Only launch position and
  normal are inherited from the nearest edge-connected face that satisfies the declared inset,
  chosen by deterministic multi-source graph distance. At the failed state this affects 29 faces,
  `8.3764e-6` of total surface area, with maximum support distance `0.003336 um` (`0.667 dx`). All
  34,624 rays then classify with no replay or refusal. A bounded real-engine continuation advances
  steps 114--116 with exact material ledgers and radiosity error at approximately `1e-12`.
- **The continuation exposed one remap roundoff ghost, not lost inventory.** Two target triangles
  were exactly edge-on in the source-authoritative chart. Their true projected area was zero, while
  translated float64 clipping returned `2.7756e-17` mesh-area slivers and divided by a zero chart
  factor. Revision `efaa070` rejects charts below the standard `sqrt(eps)` numerical-rank threshold;
  the unmatched area remains in the existing removed/newly-exposed ledgers. The exact checkpoint
  replay records two rejected chart pairs among 439,498 aligned candidates, zero material-ledger
  residual, no visibility recovery, and unchanged step-114 profile observables.
- **This is preflight verification, not held-out validation.** The forensic replay explicitly records
  `held_out_profile_data_read=false` and `scientific_result=false`. Receipt:
  `results/krueger_2024_base_authority_5nm_efaa070_remap_forensic/audit.json`. The next action is the
  single preregistered clean `t=0 -> 60 s` base anchor under commit `efaa070`; only a passing base
  endpoint can freeze the operator and authorize the untouched oxygen/power predictions.

### Local level-set gradients close both late Krüger visibility refusals

- **The apparent launch failure was a target-normal defect.** At the preserved 5 nm and 10 nm
  checkpoints, the failing origins were in gas and exact rays reached real gas-to-solid crossings.
  Five triangles in each mesh were oriented backward because the old quarter-cell level-set probe
  crossed a neighboring interface in a thin fold. This produced 12 of 36,064 and 20 of 12,928
  solid-facing classifications, respectively.
- **Revision `c986d2a` makes orientation local and first-principles.** The engine retains each
  geometric triangle normal and aligns it with `-grad(phi)`, using the exact derivative of the
  trilinear nodal level set at the triangle centroid. A manufactured thin-fold regression binds the
  failure; a true level-set critical point refuses rather than guessing.
- **Both real failure states now pass and advance.** Full corrected visibility passes in 0.436 s
  (5 nm) and 2.047 s (10 nm) on CPU. One exact next engine step then completes at each resolution
  with zero material-ledger residual and bounded motion: step 204 at 5 nm (0.0556 cells) and step
  316 at 10 nm (0.2933 cells). The repository suite passes 974 tests with one expected skip.
- **Claim boundary:** no held-out profile was read, and the CPU candidates are not replacement
  authority endpoints. Compute stops here. The next expensive action remains one checksum-bound
  base completion, followed by the preregistered freeze-or-one-correction gate. Report and receipt:
  `KRUEGER_2024_NORMAL_ORIENTATION_CLOSURE_2026-07-18.md` and
  `results/curated/krueger_2024_normal_orientation_closure_2026-07-18.json`.

### The Krüger reveal gate now binds the actual R1.9 proposal and clean launch epoch

- **A stale handoff seam was closed before it could contaminate validation.** The freeze utility
  understood the superseded axisymmetric-secant schema but not the current R1.9 trust proposal, and
  it did not require the 10/5 launch manifests. The archived proposal also carried a transient
  protocol hash. It was regenerated from the same sealed base-only inputs against the committed
  protocol; the physical pair is byte-for-byte unchanged.
- **Revision `28ffd45` makes source lineage executable.** Freeze now requires paired checksum-valid
  10/5 launch manifests, one clean revision/archive, the declared 5/10 nm refinement relationship,
  common-refinement remap, exact R1.9 parameters, and launched executable hashes matching the code
  being frozen. These receipts are embedded in the reveal.
- **The old partial launches refuse for the right reason.** Their `feature_step_3d.py` hash predates
  the local-gradient normal repair, so they remain diagnostics and cannot authorize held-out work.
  A new clean `t=0` authority epoch is mechanically required. No held-out data, simulation, or GPU
  was used in this preflight. The suite passes 976 tests with one expected skip. Report:
  `KRUEGER_2024_REVEAL_GATE_PREFLIGHT_2026-07-18.md`.

### The clean-epoch base authority pair is launched (2026-07-20)

- **Re-orientation after handoff.** The uncommitted visibility source-relaunch instrumentation
  found on top of `5ac9458` was judged speculative: no documented open refusal motivates it, it
  softens the hard-visibility refuse-rather-than-guess contract, and the pilot summary canary test
  fails with its new keys. It is preserved unmodified on the side branch
  `wip/visibility-source-relaunch` (`9a59b0c`) and removed from the campaign branch. The campaign
  tree is clean at `5ac9458`; the full suite passes 976 with one expected skip on that tree.
- **New-epoch launch manifests authored and gate-validated before spending compute.** Paired
  10/5 nm manifests for source epoch `5ac9458bc419cff48ff95b35a11fd3f90b7eee1f` (archive SHA-256
  `75f25d599fe219ba278270c981ad2972d40219191968c9030d7244788cc8e462`) carry the same R1.9 fixed
  pair, operator, and data firewall as the sealed efaa070 manifests; only `feature_step_3d.py`
  changed (`c986d2a` local-gradient normal orientation). The pair was validated against
  `krueger_2024_freeze.py::_bind_launch_manifests` (including executable-hash-vs-current-source)
  before commit (`e3322dc`). The freeze gate binds azimuth/charging evidence by content, not by
  source epoch, so the sealed receipts remain valid and were not regenerated.
- **Both runs are live on checksum-verified boxes.** The archive SHA was re-verified on each box
  before extraction, and the four executable source hashes were recomputed on-box and match the
  manifests. 5 nm authority: Vast instance 45402811 (RTX 4090, 48 effective cores). 10 nm
  companion: Vast instance 45402814 (RTX 4090, 16 cores). Each runs under an unattended
  supervisor that resumes across `wall_budget_checkpoint` exits and stops after three consecutive
  process failures so a deterministic refusal is preserved as evidence rather than retried
  indefinitely. Missing wheels on the boxes (scikit-fmm built from source with g++, scikit-image,
  numba) were installed before launch; no repository source was modified.
- **Claim boundary:** nothing scientific is earned by launching. No held-out observation has been
  read. The endpoint alone decides freeze-or-one-correction per R1.8/R1.9.

### The fifth latent defect: launches inside the piecewise-linear artifact zone (2026-07-20)

- **Both clean-epoch authority runs refused deterministically** (5 nm step 289, 16/37920 rays;
  10 nm step 328, 9/13344): exact replay met solid-facing hard intersections at the mask necking
  region. The half-finished uncommitted relaunch found at re-entry was completed, gated, and
  proved insufficient on the preserved state (2/16 recovered) -- the pathology is deeper than a
  mis-launched quadrature point.
- **Forensics, in order, each on the real failing checkpoint:** an exhaustive brute-force
  crossing audit proved the float64 tracer sound (zero missed front crossings, so no seam bug);
  per-face authority sampling showed refusing launches lie within 0.32 dx of the trilinear
  surface while sitting inside neighboring piecewise-linear solid; grazing rays traverse
  extended contact slabs (up to ~90 dx of path) where exactly coincident facets correctly cancel
  as interior surface. Path-length and crossed-plane-perpendicular admission rules were each
  implemented, measured, and refuted before the final criterion was selected from data.
- **The operator recovery (`6789dc2`):** a ray whose first intersection is a back face never
  crossed a front face, so its launch is inside the discretization artifact zone; if the launch
  lies within 0.5 dx of the authority (trilinear) surface -- sampled by a vectorized closure the
  engine now passes into the visibility estimate -- it may exit through consecutive back-face
  crossings. A launch farther from the authority surface still refuses hard. Receipts count
  overlap exits and record exit depth through the estimate receipt, engine diagnostics, and
  pilot summary; manufactured regressions bind recovery, grazing admission, authority admission,
  and the still-fatal deep refusal. All four seeds on the failing state classify completely.
  Suite: 981 passed, 1 expected skip.
- **Claim boundary:** flux carried by recovered rays is <=0.1% of the quadrature with row weight
  1/8 -- below every stated error budget -- and every recovery is receipted, never silent.
