# Charging solver campaign: Task 0 decision record

Date: 2026-07-13

## Outcome

The Task 0 decision gate is **closed**. Tasks 1--4 were not implemented or run.

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

## Decision and bounded next action

1. Hold Tasks 1, 2a, 3, and 4. Task 0A has not cleared its precision gate.
2. Drop Task 2b's local planar Maxwellian split for this model; it passed flat physics but failed the
   actual trench stiffness criterion.
3. If the campaign resumes, acquire a higher-precision paired ensemble response at the stuck state.
   The level-13 minimum signal/error of 0.512 implies that brute-force replicate scaling alone would be
   expensive; use an accuracy-matched GPU run or a targeted estimator-variance improvement rather than
   another broad CPU sweep.
4. Preserve the unchanged 0.08 contract and exact hard-visibility final audit.
