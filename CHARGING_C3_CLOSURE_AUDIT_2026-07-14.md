# C3 charging closure audit

Date: 2026-07-14 (America/New_York)

Status: **the corrected unified charging engine is operational and integrity-verified; the declared
CCA-2026-07-13-R2 steady-state contract is not converged; C3 remains open and C4 remains blocked**

This is the terminal report for the bounded July 14 charging campaign. The paid GPU was destroyed
after every referenced artifact was copied and independently audited. No further continuation is
authorized by this report.

## Verdict

| Question | Answer | Evidence |
| --- | --- | --- |
| Does the unified engine advance the exact hard-visibility kinetic charging operator without losing charge? | **Yes** | Six terminal trajectories complete with zero rejected steps; charge and surface-transfer ledgers close at `3.87e-15`--`4.95e-15` and `2.79e-15`--`2.84e-15`. |
| Are particle wrapping, electrostatic topology, and the authoritative charge space now consistent? | **Yes** | The periodic Q1 seams are exactly `0 V`; particle/field topology mismatch is refused; the periodic face coupling is rank 22/nullity 18 and compatible updates retain no measurable null inventory. |
| Did the decreasing-gain tail overshoot or jump to a different branch? | **No** | Its independent 50 microsecond fixed-physical confirmation continued in the same dominant voltage direction. |
| Has the corrected state reached a statistically stationary physical equilibrium? | **Not yet** | The confirmation moved another `2.272 V` along the archived dominant mode and had integrated B1 `17,748.6 V/s`. |
| Does the state pass signed R2? | **No** | B1 is 17.75 times its `1,000 V/s` gate; raw B2 is `1.348/1.297` versus `0.08`; even the non-governing Q1-resolved diagnostic is `0.1681/0.1677`. |
| Is charging now usable as an engine capability? | **Yes, provisionally** | The exact operator, conservative transient, restart, GPU path, compatible state, fresh-scramble mode, and bounded stochastic warm-start path all work. The result may seed engineering/demo co-evolution, but not a C4 validation claim. |

The canonical warm checkpoint is
`results/charging_c3_periodic_topology_audit/c3_periodic_decreasing_gain_fixed_confirmation_50us_l11/face_checkpoint.npz`
(SHA-256 `2ed2c604728481fe217ccaa73d405905f79ccf24ec0cbb7ac0daf173ef51ce64`).
It is a valid advanced charging state, not a certified equilibrium.

## What the campaign actually fixed

The original failure was not one impossible nonlinear equation. It was a stack of independently
testable engine defects and estimator mistakes:

1. A frozen finite Sobol set injected persistent current bias into physical time. Fresh independent
   scrambles make the step error fluctuate instead of accumulating coherently.
2. Slow adjoint electrons exceeded a shortened total flight horizon. The production horizon is now
   refinement-certified and incomplete trajectories are refused.
3. Float32 ray tracing occasionally missed shared triangle edges. Every hit is certified and only
   invalid lineages are replayed against edge-inclusive float64 hard triangles.
4. P0 triangle charge contained modes that the Q1 Poisson field could not see. The authoritative
   state now lives in the compatible Q1 charge space while the raw diagnostic remains reported.
5. Particle trajectories wrapped laterally while Poisson did not. Periodic endpoint identification
   is now built into the common Poisson system, and mismatched topology refuses before transport.
6. A fixed low-sample late transient could not distinguish mean drift from estimator jitter. The
   engine now supports a declared decreasing-gain, fresh-scramble pseudo-time tail that preserves
   the same zeros, never credits physical time, never self-certifies, and is followed by an
   independent fixed-physical confirmation.

These are engine repairs, not test edits. Each one is now an invariant or an explicit run mode in
the unified transport/Poisson/charging path.

## Corrected physical reference

All rows use the repaired periodic Q1 operator, compatible face state, fresh scrambles, exact hard
visibility and reflection, CUDA transport, 400 accepted fixed steps of 125 ns, and an independent
terminal-window audit.

| Corrected physical window | Integrated B1 (V/s) | Raw B2, 0.25/0.50 micrometers | Q1-resolved B2, 0.25/0.50 micrometers | Node RMS/worst | Integrity |
| ---: | ---: | ---: | ---: | ---: | --- |
| 0--50 microseconds | 60,510.1 | 1.2104 / 1.2095 | 0.3060 / 0.3036 | 0.0285 / 0.0768 | pass |
| 50--100 microseconds | 20,579.6 | 1.1819 / 1.1434 | 0.1995 / 0.1989 | 0.0382 / 0.0959 | pass |
| 100--150 microseconds | 17,389.1 | 1.3854 / 1.3422 | 0.1646 / 0.1642 | 0.0379 / 0.0814 | pass |
| 150--200 microseconds | 13,962.3 | 1.3743 / 1.3394 | 0.1512 / 0.1508 | 0.0493 / 0.1202 | pass |
| 50 microseconds after stochastic tail | 17,748.6 | 1.3484 / 1.2974 | 0.1681 / 0.1677 | 0.0385 / 0.0743 | pass |

All five rows have zero rejected steps, zero bounce-budget extensions, zero trajectory-horizon
extensions, exactly zero periodic voltage seam, and independently reproduced endpoint Poisson
solutions. The maximum exact-lineage replay fraction is `8.10e-4`; replay is a counted repair, not
a discarded particle.

## Physical drift versus estimator jitter

An SVD was fit only to the four consecutive fixed-physical reference displacements. Its dominant
direction contains 96.597% of their displacement energy. The confirmation was excluded from this
fit and then scored against the archived direction.

| Motion | Dominant-mode projection (V) | Orthogonal L2 (V) | Interpretation |
| --- | ---: | ---: | --- |
| 0--50 microseconds | -17.3410 | 0.5358 | initial coherent charging |
| 50--100 microseconds | -2.9899 | 1.8728 | mean drift smaller; sampling motion visible |
| 100--150 microseconds | -1.7053 | 1.9368 | continued same-direction drift |
| 150--200 microseconds | -0.3742 | 1.8628 | one window is mostly orthogonal jitter |
| decreasing-gain tail from 200 microseconds | -2.2001 | 0.9633 | reduced-gain proposal remains on the branch |
| fixed confirmation from the tail | -2.2715 | 2.4415 | no reversal; physical stationarity is not established |

The last fixed confirmation is the deciding observation. Had the tail overshot, physical time would
have moved back toward the 200 microsecond state. It instead moved farther in the same dominant
direction (`cosine = -0.681`). The tail is therefore a safe warm-start accelerator for this branch,
but neither it nor the following fixed window is an equilibrium certificate.

The reproducible calculation and checkpoint hashes are in
`results/charging_c3_periodic_topology_audit/closure_mode_audit.json`, generated by
`scripts/charging_c3_closure_mode_audit.py`.

## Why pointwise B1 and raw B2 no longer answer the scientific question alone

Eight independent endpoint audits at three nested sample levels cannot resolve a pointwise B1 of
`1,000 V/s`: global signal-to-standard-error ratios are `0.546/0.253/1.211` at the pre-tail state and
`0.278/0.390/1.468` at the tail, and every selected largest-component 95% interval contains zero.
This is not evidence that B1 passes. It is evidence that the current estimator cannot certify such
a small instantaneous derivative at affordable sample counts.

Raw B2 has a different defect. On the intended periodic grid the conservative `72 x 40` face
coupling has rank 22 and nullity 18. Raw face-patch sums can see those 18 exact field-null modes,
while Poisson and particle trajectories cannot respond to them. The final confirmation reports raw
B2 `1.348/1.297`, of which `1.181/1.130` is explicitly attributed to the unresolved face component;
the Q1-resolved part is `0.1681/0.1677`. CCA-R2 nevertheless remains in force, so this report records
failure rather than silently replacing its gate.

## Decreasing-gain result

The tail uses

    dt_k = 125 ns * (16 / (16 + k))^0.75

for 400 fresh-scramble pseudo-steps. It advances `10.1224 microseconds` of declared pseudo-time and
zero physical time; the final gain is about `10.88 ns`. It completes 400/400 with zero rejects,
conserves charge to `4.42e-14`, and reduces the late orthogonal displacement from about `1.86 V` to
`0.963 V`. It correctly writes `converged: false` and cannot self-certify. Full local tests after
this implementation pass: **404 passed, 1 skipped**.

The fixed confirmation uses config hash
`c919399966fe14254d58bdcdfceacab5715240abccbdc0d03777dc87652ff35d`, summary SHA-256
`eccbc052d3b7398705ecc05cb1f7757741dd98ccbdc5a72e73be50bfc6f06ed4`, and source hashes recorded
inside its manifest. Its independent audit has `integrity_pass: true` and
`contract_converged: false`.

## Closure decision

The July 14 compute campaign stops here for a principled reason: another long run on the same grid
and sample budget can move the state, but it cannot by itself distinguish a `1,000 V/s` B1 signal
from estimator variation, cannot make the Q1 field control its 18 raw-face null modes, and cannot
supply the missing physical-grid/B3 refinement. More runtime is therefore not a convergence proof.

The lasting path is bounded and ordered:

1. keep the final confirmation checkpoint as the warm state for unified co-evolution and demos;
2. obtain one physical grid refinement and score field-compatible balance plus the raw-compatible
   discrepancy, without changing R2 unless R3 is explicitly signed;
3. replace single-scramble pointwise B1 certification with a preregistered, independent statistical
   stationarity test only through formal contract change control;
4. use multi-scramble/current-block averaging or decreasing-gain warm starts where wall-clock matters,
   while retaining one fixed-physical, exact-operator confirmation; and
5. start C4 validation only after the signed contract, B3 observable refinement, and B5 audit pass.

This is not an impossible charging problem. It is a working stochastic multiphysics engine whose
remaining blockers are a still-resolved slow physical drift and two evidence-contract/discretization
questions. Those are now isolated; no further frozen-map solver work is warranted.

## Evidence

- `results/charging_c3_periodic_topology_audit/periodic_terminal_audit.json`
- `results/charging_c3_periodic_topology_audit/window_50to100_audit/periodic_terminal_audit.json`
- `results/charging_c3_periodic_topology_audit/window_100to150_audit/periodic_terminal_audit.json`
- `results/charging_c3_periodic_topology_audit/window_150to200_audit/periodic_terminal_audit.json`
- `results/charging_c3_periodic_topology_audit/c3_periodic_endpoint_sampling_audit/audit.json`
- `results/charging_c3_periodic_topology_audit/c3_periodic_decreasing_gain_endpoint_sampling/audit.json`
- `results/charging_c3_periodic_topology_audit/decreasing_gain_fixed_confirmation_audit/periodic_terminal_audit.json`
- `results/charging_c3_periodic_topology_audit/closure_mode_audit.json`
- `results/charging_c3_periodic_topology_audit/canonical_warm_state.json`

GPU instance `44895783` was destroyed after retrieval; the final Vast instance listing was empty.
