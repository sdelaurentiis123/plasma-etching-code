# Charging solver campaign: Task 0 decision record

Date: 2026-07-13

## Outcome

The Task 0 diagnostic gate is **closed**. Its response-precision result continues to hold Tasks 3--4
and any response-derived accelerator, but it no longer blocks the derivative-free physical-time
transient in Task 1. This is a scope correction: every resolved Task 0A branch in the original handoff
sent the campaign through Task 1, which consumes conservative current evaluations rather than a
finite-difference response.

- Task 0A produced a complete, paired, exact-operator audit, but the smallest level-13
  signal-to-between-scramble-error ratio is 0.512. The required response is therefore still below
  estimator uncertainty. The paired restricted condition numbers are diagnostics, not a basis for
  promoting a solver.
- Face-switch fractions scale almost linearly with perturbation radius (level-13 fitted slopes
  0.851--0.950; median 0.937). This supports a regular switching set at the measured scales, but it
  does not establish that switch noise caused the archived full-Jacobian condition number.
- Task 0B passes the planar Maxwellian barrier gate after resolving slow zero-field trajectories, but
  the local barrier term captures the dominant trench-region derivative in 0/4 regions. The proposed
  Task 2b planar boundary-current preconditioner is rejected.

The handoff's `response below estimator uncertainty` branch therefore applies: increase paired
precision before any solver conclusion. No frozen-map solver variant, transient, PTC, smoothing, or
stochastic fixed-point implementation was added.

## Preserved inputs

| Input | SHA-256 | Role |
| --- | --- | --- |
| `results/charging_task0_inputs/stuck_ar4_checkpoint.npz` | `6e092bf7704d74e5da03a42493559a74f08de086f9c83e18d254b7a7f93fe65d` | canonical available stuck state |
| `results/charging_task0_inputs/relaxed_ar4_checkpoint.npz` | `5964778f9576674d2ab7b98d97331bfe0bbb5acf612138a848f8396f7f93211f` | superseded pilot provenance only |
| `results/charging_task0_inputs/archived_fd_jacobian_0p025v.npz` | `0bcc115ab62d3960e74d0d30ed5a10a2c7788b17d3e7506c417eca9d08be6479` | dominant directions and archived full condition |

The serialized checkpoint matching the handoff's coarse-3D RMS 0.627 / maximum 0.978 state was not
present in the repository or temporary artifacts. The canonical available failure checkpoint is the
47-DOF `ar4_rule2_baseline` state, with raw RMS 0.597912 and maximum 1.37957, and is in the same lineage
as the archived 47-by-47 finite-difference Jacobian. The audit uses this checkpoint rather than the
partially relaxed `eval_epoch5` state (RMS 0.209756). Results must not be represented as measurements
of the unavailable coarse-3D checkpoint.

## Task 0A: paired ensemble response

Configuration hash: `b0d263fbce2bf3f29946957b03145b4f204ef31cba09a39bcb886aa2745d2db7`

- Exact hard-visibility kinetic operator; no smoothing.
- Common random numbers for `Q-delta`, `Q`, and `Q+delta` within every scramble.
- Eight scrambles: 401, 409, 419, 421, 431, 433, 439, 443.
- Nested Sobol levels 9, 11, and 13.
- Radii 0.1, 0.05, and 0.025 V.
- Five directions: the two worst checkpoint coordinates and three dominant archived-Jacobian
  directions.
- 360/360 paired evaluations completed. Every process pins OMP, OpenBLAS, MKL, NumExpr, and Numba to
  one thread. Python thread pools are prohibited because Numba's workqueue backend aborted under
  concurrent thread use.

Highest-level results (five-direction subspace):

| Radius (V) | Median frozen response | Median ensemble response | Minimum signal/error | Frozen condition | Ensemble condition (95% bootstrap CI) | No-switch condition |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.100 | 4.140 | 2.328 | 1.197 | 9.98 | 22.01 (9.49--26.55) | 12.69 |
| 0.050 | 4.573 | 2.675 | 0.718 | 8.66 | 17.10 (7.62--26.61) | 13.27 |
| 0.025 | 5.315 | 3.102 | 0.512 | 7.09 | 12.49 (7.86--19.26) | 21.81 |

These are identical-direction, identical-checkpoint, first-scramble frozen versus eight-scramble
ensemble comparisons. They are not comparable in dimension to the archived full 47-by-47 condition
number (~72,000 at 0.025 V). The earlier apparent two-orders-of-magnitude claim is not established by
this paired experiment.

The signed attribution is no longer clamped. It reports switch-component energy, signed projection,
full/no-switch condition numbers, and `log(cond_full/cond_no_switch)` with replicate-bootstrap errors.
At level 13 the signed log changes are +0.551 +/- 0.316, +0.253 +/- 0.373, and -0.557 +/- 0.351 from
large to small radius: switch removal is not a consistently beneficial conditioning operation. A
synthetic injected-switch test verifies both positive (harmful) and negative (beneficial) signs. The
results are unchanged when the active-node threshold is swept over `1e-5`, `1e-4`, and `1e-3`.

The switch-fraction scaling is much better resolved than the response conditioning. At level 13 all
five log-log slopes lie between 0.851 and 0.950 (median 0.937), close to the regular-boundary prediction
of one. See:

- `results/charging_task0a_stuck/ensemble_response.png`
- `results/charging_task0a_stuck/switch_scaling.png`
- `results/charging_task0a_stuck/condition_decomposition.png`
- `results/charging_task0a_stuck/summary.json`

## Task 0B: analytic electron boundary current

Configuration hash: `28f4887dc11ee0e7a34291637ee68e3de6e9902753babf1f9964a986cd9f0b8f`

The first flat run exposed a finite-trajectory-horizon artifact at exactly 0 V: 4.675% of slow grazing
reverse trajectories remained unresolved after 1,600 steps, and the flux was 0.95325 rather than 1.
Increasing the audited horizon to 25,600 steps reduced the unresolved fraction to 0.000488 and moved
the flux to 0.999512. This is resolution of the unchanged kinetic operator, not a relaxed gate.

The refined flat gate passes: RMSE 0.000500, maximum absolute error 0.001204. The repelling branch
matches `exp(V/T_e)` to approximately 1e-7 and the attracting branch saturates at unity.

The first trench pilot used charge-coordinate shifts that moved local voltage by tens to hundreds of
volts and could not resolve a 4 eV electron scale. The reported campaign instead uses Poisson linearity
to calibrate a separate charge shift per region and samples each region at mean local voltages
-12, -8, -4, -2, 0, 2, 4, and 8 V.

| Region | Derivative capture | Derivative correlation | Gate |
| --- | ---: | ---: | --- |
| top | -109.893 | -0.716 | fail |
| upper wall | -6.037 | 0.157 | fail |
| lower wall | -2.081 | 0.621 | fail |
| floor | -3.824 | -0.504 | fail |

Pass required capture >= 0.5 and correlation >= 0.8 on most regions. Result: 0/4. The analytic curve
is diagnostic only; every score and final statement uses the hard-visibility kinetic current. See
`results/charging_task0b/electron_boundary_audit.png` and
`results/charging_task0b/summary.json`.

## Decision amendment and bounded next action

1. Run Task 1 after a paired production trajectory-horizon/timestep audit. Task 1 remains subject to
   all conservation, timestep-refinement, independent-final-audit, branch, and exact-operator gates.
2. Continue to hold Tasks 3 and 4. Task 0A has not cleared the response-precision gate needed to
   select those accelerators.
3. Drop Task 2b's local planar Maxwellian split for this model; it passed flat physics but failed the
   actual trench stiffness criterion.
4. Hold Task 2a and any replacement preconditioner until Task 1 establishes a timestep-refined steady
   state. A measured nonlocal response is the only replacement currently consistent with Task 0B,
   and it must not be built until its paired basis clears the signal/error gate. No present radius does.
5. In parallel, acquire a higher-precision paired ensemble response at the stuck state only if it is
   needed to choose an accelerator.
   The level-13 minimum signal/error of 0.512 implies that brute-force replicate scaling alone would be
   expensive; use an accuracy-matched GPU run or a targeted estimator-variance improvement rather than
   another broad CPU sweep.
6. Preserve the unchanged 0.08 contract and exact hard-visibility final audit.

The archived 47-by-47 finite-difference campaign used the same deterministic checkpoint, frozen
proposal rules, and seeds for each plus/minus pair. It must not be described as an unpaired Jacobian.
Its large full-space condition number remains real evidence about that frozen discrete operator, while
the underpowered five-direction ensemble audit cannot yet attribute the difference.

## Task 1-pre: production trajectory horizon

Configuration hash: `21cdb54c47754ec652b424260684b8fc8e15bdabd2c76716df2cc65b9795c3bd`

The canonical stuck current map was evaluated with common samples across the production horizon
(`dt=0.01`, 4,000 steps), 4x and 8x horizons, and a half-timestep reference (`dt=0.005`, 64,000
steps). Eight scrambles at level 9 found zero unresolved ion-forward, ion-adjoint, or
electron-adjoint trajectories at every horizon. The 1x, 4x, and 8x current maps are identical.

Against the half-timestep reference, the ensemble ion and electron current changes are 0.287% and
0.170%; RMS and worst-node imbalance changes are 0.130% and 0.095%. The paired per-scramble RMS
change is noisier (3.23% +/- 0.86%), but the ensemble-mean current and residual gates pass. Therefore
the historical stuck residual is not explained by the 4.7% unresolved-electron artifact seen only in
the earlier flat zero-field audit, and Task 1's horizon entry gate passes. See
`results/charging_task1pre_horizon/trajectory_horizon_audit.png` and `summary.json`.

## Task 1 and 2a execution: engine transient works; equilibrium contract remains red

The engine now exposes a reusable `integrate_dielectric_charging_transient_3d` API rather than a
campaign-only loop. It records separate positive/negative face and compatible Q1-node currents,
worst/RMS node and face imbalance, exact per-step conservation, replayable nodal charge histories,
and a resumable checkpoint on transport failure. The convergence equation remains the established
active-node balance. Triangle balance is reported separately because triangles are integration
elements, not independent stored-charge degrees of freedom.

Two estimator bugs/inefficiencies were fixed during entry:

- A folded surface-local grazing proposal was briefly combined with a source-aligned ion frame. That
  swaps the ion's large source-normal speed into a horizontal component and creates arbitrarily slow
  reverse rays. The invalid combination is removed and now rejected by the engine.
- Once a separately certified forward/adjoint face map is frozen, the engine now traces adjoint rays
  only on adjoint-selected faces and skips a direction entirely for all-forward/all-adjoint maps. A
  regression gate proves the selected current is bitwise unchanged.

The deterministic exact-operator branch used paired timestep ladders from 1 ns through 250 ns and
warm continuation to about 160 microseconds. Key accepted comparison: over a 50-microsecond interval,
250 ns versus 125 ns steps agree to 0.298% in charge and 0.139% in potential, ending at RMS 0.195/0.197
and worst node 0.439/0.440. A 500 ns step is rejected (31.6% potential disagreement). The trajectory
forms the nonlocal dipole expected from the physics: upper regions become strongly negative while
lower-wall/floor regions become positive. Longer accepted-step trajectories fluctuate instead of
closing the 0.08 node contract; this is not a conservation or timestep blow-up.

An independent exact hard-visibility endpoint audit used 16 ion bidirectional replicates through
level 13, up to 32 face-position points, and an electron proposal at level 12. Every ion face
certified; the unchanged state scored RMS 0.193551, worst node 0.445733, and worst face 0.841180.
Electron between-scramble uncertainty is not yet attached, so this is not represented as a complete
final audit—and the residual is far above contract regardless.

The fresh-scramble Task 1 variant froze the method map from that separate pilot and used independent
scrambles at every physical step. Over the final 40 of 200 steps its statistically stationary values
are RMS 0.199220 +/- 0.006366 and worst node 0.574845 +/- 0.016865. Conservation residual is
1.40e-17 relative to absolute charge throughput. Unfreezing samples therefore changes the path but
does not remove the residual floor.

The bounded current-direction PTC used no derivative or quasi-Newton update. It accepted only
safeguarded physical-current directions and halved pseudo-time after residual worsening. Result: 6
accepted steps, 14 rejections, pseudo-step collapse below 1e-11 s, best RMS 0.179424 and worst node
0.465949. It does not converge and the PTC branch is closed.

Decision: no additional frozen-map solver variant. Promote the handoff's equilibrium/discretization
audit (compatible nodal balance versus area-weighted faces and coarsened patches under grid
refinement). Preserve the 0.08 contract and flag any demonstrated discretization-scale mismatch for
human review rather than changing it. The local planar electron preconditioner remains rejected;
Tasks 3--4 remain held by the Task 0A precision gate.

## Equilibrium/discretization audit: wall-scale imbalance survives refinement

The engine now exposes `current_balance_metrics_3d` for integrated currents on raw elements or an
explicit integer patch map. It reports the unchanged unweighted local RMS/maximum plus
throughput-weighted RMS and global balance. Face current density must be multiplied by physical area
before aggregation; this prevents unequal marching-cubes triangles from being averaged as peers.

Configuration `d4ce71baa71c43f4d66e3c6af5533ddf58c43e6c2c1adf9a5c9b32a94b44905b`
evaluated the same archived voltage field on `dx=0.25` and `0.125 um` grids. Trilinear voltage
transfer is reproduced by the refined Poisson solve to `2.19e-14` relative L2. A separate pilot chose
the ion forward/adjoint map (39/1 coarse, 164/12 fine); four independent scoring groups could not
reselect it. Every result uses hard visibility and full kinetic currents.

At forward level 11, coarse/fine global imbalance is 0.03862/0.03781. Raw-face RMS is 0.427/0.463;
compatible-node RMS is 0.206/0.344. Fixed 0.5 um patches are substantially less sensitive: RMS
0.162/0.173 and maximum 0.296/0.280. All four half-micron wall patches remain independently resolved
outside 0.08 with the same sign on both grids. Their signed coarse/fine RMS difference is 0.064.
Thus the current state contains a physical wall-scale redistribution mode; the floor is not solely a
single-triangle or nodal-projection artifact. This is an operator-refinement statement, not a claim
that the transferred voltage is a refined-grid equilibrium.

The forward sample-level comparison is also diagnostic. Raising level 9 to 11 changes fine node RMS
only 0.3465 to 0.3440 and half-micron patch RMS 0.1648 to 0.1735, but reduces the worst node from
0.8924 to 0.8000. Small fine-grid faces therefore make the worst-node statistic more estimator-level
sensitive than the RMS or physical patches. The contract is unchanged; the audit records both levels.

The refined physical transient then advanced the mapped state for 15 microseconds. Its paired
125/62.5 ns schedules agree to 0.0166% in charge and 0.00559% in potential over the final 10
microseconds. At the integration level, node RMS falls from 0.3379 to 0.3020 while the worst node
remains 0.8688. Independent frozen-map endpoint scoring at forward levels 11 and 13 gives stable RMS
0.29975/0.29888; the worst node continues to move with sampling, 0.7717/0.7143, but remains far above
0.08. At level 13, the half-micron patch RMS is 0.1287 and its maximum is 0.2587: the upper walls
remain electron dominated while the floor, top, and lower-wall half patches are inside 0.08. Global
imbalance is only 0.00058, demonstrating redistribution rather than missing total current.

One engine replay defect was fixed from this evidence. After adaptive certification, the transient
previously froze only the estimator method and silently returned to base sample/position levels.
Internally discovered maps now replay at the declared certification ceilings; externally supplied
maps retain their explicitly audited scoring levels. No solver or tolerance changed.

Decision: do not reinterpret the present per-node failure as pure discretization, and do not loosen
0.08. A longer refined transient is computationally justified only after attaching per-face certified
sampling levels (rather than a method-only map) or moving this accuracy-matched campaign to the GPU.
The missing physical closures—surface conduction, bulk leakage, secondary emission, and reflection—
should be assessed as model-scope candidates before another root algorithm.
