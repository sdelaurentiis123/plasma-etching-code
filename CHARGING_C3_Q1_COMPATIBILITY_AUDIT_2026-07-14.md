# C3 charging: Q1 state-compatibility and acceleration audit

Date: 2026-07-14 (America/New_York)

Status: **two coupled discretization defects confirmed and repaired; corrected fresh-scramble
physical reference, decreasing-gain tail, and fixed confirmation complete; engine integrity passes;
CCA-R2 does not converge; C3 remains open and C4 remains blocked**

This report does not declare charging convergence and does not change CCA-2026-07-13-R2. It
documents why the legacy face-state trajectory could not satisfy its raw patch diagnostic, the
bounded engine change that removes the inconsistency, and the exact-operator tests that rejected
an attempted projective acceleration.

The terminal campaign verdict and canonical warm checkpoint are recorded in
`CHARGING_C3_CLOSURE_AUDIT_2026-07-14.md`. In short: four corrected 50 microsecond physical windows
reduced integrated B1 from `60,510` to `13,962 V/s` and Q1-resolved B2 from `0.306` to `0.151`, but
an independent fixed 50 microsecond confirmation after the decreasing-gain tail still moved
`-2.272 V` along the archived dominant charging direction. It passed every integrity check and
failed R2 at B1 `17,748.6 V/s`, raw B2 `1.348/1.297`, and Q1-resolved B2 `0.1681/0.1677`. The GPU was
destroyed after artifact retrieval. The engine is operational; formal steady-state convergence is
not claimed.

## Critical topology correction (added after the compatibility continuation)

The compatible-state continuation exposed a second, more fundamental wiring defect: particle
transport wrapped trajectories periodically in lateral x/y, but the Q1 Poisson system treated the
opposite endpoint planes as independent natural-Neumann boundaries. The saved field therefore had
discontinuous periodic seams which grew to **20.413 V in x and 24.150 V in y**. A particle crossing
the cell boundary saw a different potential after wrapping. This is not a slow physical mode and no
amount of integration can repair it.

The GPU run was stopped after two completed compatible segments at 2.94875 ms cumulative archived
time. No charging process or GPU compute client remains active. The common Poisson engine now
identifies periodic endpoint nodes before factorization, prolongs the independent solution back to
the closed interpolation grid, and refuses any particle/Poisson topology mismatch before transport.
The endpoint voltage differences are now bitwise **0 V** in x and y.

The intended periodic 0.25 micrometer operator has 72 independent volume nodes and a `72 x 40`
surface coupling of rank 22 and nullity 18, rather than the historical nonperiodic `135 x 40`, rank
34, nullity 6 coupling reported below. Projecting the stopped checkpoint into the corrected
periodic-compatible state:

- preserves the effective periodic nodal load to `2.65e-15` relative L1;
- preserves global face charge to `1.06e-29 C`;
- reduces periodic-null inventory from 60.15% to `1.65e-15`;
- closes the periodic Poisson residual to `7.62e-31 C`; and
- changes the proposed potential by 5.45% relative L2, with a 13.62 V maximum local change.

That last number makes the consequence explicit: the old trajectory is not a corrected physical
history. Its final inventory is retained only as a warm *proposal*. The earned next experiment is a
short, bounded, timestep-refined cold/warm comparison under the corrected periodic operator; the
archived 2.94875 ms will neither be replayed nor credited toward saturation.

## Executive finding

The long fresh-scramble transient was integrating two different state spaces:

1. kinetic deposition was stored as one piecewise-constant (P0) charge per surface triangle; and
2. electrostatic feedback consumed only the conservative Q1 nodal projection of that charge.

For the historical **nonperiodic** 0.25 micrometer trench field, the face-to-node map had shape
`135 x 40`, rank 34, and nullity 6. Six independent face-charge patterns therefore produced exactly
zero nodal load and
exactly zero change in potential, field, or particle trajectory. The physical-time update could
accumulate those patterns forever, while the Q1 field had no degree of freedom with which to oppose
them. Raw patch sums can also see those six modes. Consequently two face states with identical
electrostatics can report different B2 values.

This is a mixed-discretization compatibility defect, not evidence that the kinetic charging physics
lacks an equilibrium. It explains the otherwise contradictory late trajectory:

- raw B2 remained `1.56 / 0.813` at the two declared patch scales;
- Q1 terminal-window node RMS/worst were already `0.0101 / 0.0334`;
- 96.8% of the area-weighted face-charge norm and 95.3% of the terminal-window net-current norm lay
  in the Q1-null component; and
- B1 remained field-visible and nonzero at `3.66e4 V/s`, so the resolved state was not yet saturated.

The repair makes the Q1 nodal load authoritative and stores its unique area-weighted
minimum-density-norm face representative for remap and ledger operations. It removes only exact
field-null content. A zero-update projection preserves nodal charge to `4.46e-15` relative L1,
potential to `4.29e-14` relative L2, and total charge to `2.01e-29 C` while reducing the null
fraction from 96.75% to `2.72e-15`.

## Why the diagnosis is decisive

### The null component is systematic, not sampling noise

Across 15 independent 50 microsecond terminal windows:

- mean null-current L2 norm: `4.165e-12 A`;
- minimum cosine alignment of an individual null vector with the ensemble mean: `0.99496`; and
- averaging all 15 windows leaves `4.165e-12 A` and a 90.23% null fraction.

Zero-mean sampling noise would decorrelate and shrink. This component remains nearly collinear and
does not shrink, which identifies a deterministic basis mismatch at this grid.

### A patch sum is not automatically a function of the field state

The maximum dual sensitivity of the declared patch-sum functionals to Q1-null modes is 0.370 at
0.25 micrometers and 0.426 at 0.50 micrometers. Those are structural properties of the bases, not
current-estimator realizations. A raw patch numerator can therefore change without any change in
the state that Poisson and transport consume.

### Simple refinement does not make the algebraic mismatch disappear

| Grid spacing (micrometers) | Faces | Rank | Nullity | Maximum patch-functional null sensitivity |
| ---: | ---: | ---: | ---: | ---: |
| 0.500 | 16 | 16 | 0 | 0.000 |
| 0.250 | 40 | 34 | 6 | 0.426 |
| 0.125 | 176 | 121 | 55 | 0.362 |

The 0.5 micrometer grid is injective only because it is severely under-resolved. Refining the
surface creates more P0 triangle modes than the Q1 trace can represent. What must tighten under
physical refinement is the magnitude and observable effect of the unresolved current, not the
mere existence of the algebraic null space.

## Compatible-state validation

The engine now constructs the exact conservative face-to-node coupling once, exposes rank,
nullity, and condition number, and projects each accepted face update into the Q1-visible subspace.
The legacy raw-face path remains available only for historical replay.

A paired 500-step continuation from the same segment-13 state and the same fresh sampling epochs
gave:

| Quantity | Legacy | Compatible | Relative difference |
| --- | ---: | ---: | ---: |
| B1 (V/s) | 36,608.5 | 36,343.5 | -0.724% |
| Node RMS | 0.04815 | 0.04773 | -0.861% |
| Worst node | 0.13794 | 0.13552 | -1.75% |
| B2, 0.25 micrometers | 1.56124 | 1.56355 | +0.148% |
| B2, 0.50 micrometers | 0.81334 | 0.81027 | -0.378% |

The final compatible and legacy potentials differ by 0.0250% relative L2 and nodal charges by
0.0753% relative L1, while the stored compatible face state contains no measurable Q1-null
inventory. This is the expected result: the resolved trajectory is preserved while the
uncontrolled bookkeeping direction is removed.

The production continuation initially resumed from the projected segment-16 checkpoint at
cumulative physical time `2.82375 ms`, then was stopped after two 500-step audited segments when the
particle/Poisson periodic-topology mismatch was diagnosed. No earlier trajectory was restarted or
recomputed, and no further GPU time is being spent on the invalid operator.

## Projective/PTC acceleration result

Two projective directions were tested; neither was accepted.

1. The old 50 microsecond terminal-window current was evaluated along moving states. Positive
   jumps of 31.25--100 microseconds decisively increased the exact field-resolved residual. That
   lagged direction is rejected.
2. A current was then estimated at one fixed compatible segment-16 state. The estimator required
   70 independent level-13 scrambled-Sobol replicas to reach the preregistered signal/error ratio:
   signal `7.390e-13 A`, standard-error norm `2.353e-13 A`, ratio `3.140`.

The fixed-state direction was scored with exact hard visibility and unused common-random-number
epochs. Large candidates were rejected decisively:

| Pseudo-step | Maximum predicted voltage move | Paired change in current-residual L2 |
| ---: | ---: | ---: |
| 2.5 microseconds | about 0.80 V | `+2.058e-12 +/- 0.424e-12 A` |
| 5 microseconds | about 1.60 V | `+5.161e-12 +/- 0.483e-12 A` |
| 10 microseconds | about 3.20 V | `+10.898e-12 +/- 0.466e-12 A` |

At 0.125--0.5 microseconds the paired changes were smaller than their 95% confidence intervals;
at 1 microsecond node RMS and worst-node imbalance worsened significantly. Even a marginal small
candidate would cost 70 high-sample current evaluations to construct and more evaluations to score,
so it cannot accelerate a one-sample physical step economically. Batch projective PTC is therefore
closed for this late state. The reference remains fresh-scramble compatible physical time; a
predeclared decreasing-gain stochastic tail is the next accelerator only if the continued reference
demonstrates a stationary distribution rather than continued mean drift.

The apparently improving raw B2 of some rejected candidates is not an acceptance signal: those
same candidates worsened the authoritative Q1 residual by factors of two to seven. This is a direct
experimental demonstration of the mixed-space diagnostic conflict.

## Engine changes

- `CompatibleQ1SurfaceChargeProjector3D` exposes the conservative coupling, rank/nullity,
  area-weighted projection, reconstruction, and dual compatibility of face statistics.
- `NodalPoissonSystem3D` now supports periodic axes by algebraic endpoint identification before
  factorization. Charge reduction, equal-share canonical representatives, voltage response, and
  compatible face projection all use the same independent periodic space.
- The coupled transport path refuses periodic particle wrapping unless the Poisson system declares
  exactly the matching lateral periodic axes. This turns the campaign defect into an engine
  invariant rather than a configuration convention.
- `integrate_surface_charging_to_saturation_3d` has an explicit
  `compatible_q1_charge_state` mode, projects initial and accepted states, and records null-current
  and removed-inventory ledgers without changing the kinetic operator.
- The C3 runner and unattended supervisor persist the state-space choice in configs, checkpoints,
  heartbeats, summaries, and recovery records.
- A proposal-only compatible PTC layer enforces nodal/global conservation, a maximum voltage jump,
  zero physical-time credit, unused scoring epochs, and exact-audit-before-acceptance.
- Reproducible fixed-state ensemble-current and paired-candidate audit scripts record hashes and
  confidence intervals.

Focused tests cover manufactured null modes, exact field preservation, injective-grid identity,
periodic endpoint equality and representative invariance, particle/field topology refusal,
functional compatibility, conservative physical updates, and the PTC voltage safeguard. The full
suite after the decreasing-gain implementation passes: **404 passed, 1 skipped**.

## Contract consequence

CCA-R2 remains in force. Raw B2 and every per-node diagnostic continue to be reported on every run,
including failures. This audit does not authorize replacing raw B2 with a friendlier number.

It does establish that raw B2 is not a single-valued functional of the authoritative Q1 state on
the current grid. A draft R3 change-control request therefore separates:

1. field-compatible patch balance, which can define stationarity of the discretized coupled
   operator; and
2. raw-versus-compatible patch discrepancy, which becomes a mandatory spatial-discretization
   error and must tighten under grid refinement before an experimental claim.

Until that revision is explicitly approved and its refinement evidence exists, C3 cannot be
declared converged even if B1 closes. C4 remains blocked.

## Scientific basis

The repair follows the same compatibility principle used by charge-conserving finite-element PIC:
deposition must live in the discrete charge space consumed by the field equation, and null-space
content must not be allowed to corrupt the coupled dynamics.

- Wang et al., geometric electrostatic PIC and Whitney 0-form charge deposition:
  <https://arxiv.org/abs/2012.08587>
- O'Connor et al., quasi-Helmholtz FEM-PIC and spurious null-space charge:
  <https://www.sciencedirect.com/science/article/abs/pii/S0010465522000637>
- Crawford et al., exact charge-conserving scatter/gather:
  <https://arxiv.org/abs/1409.0854>
- Owen, randomized quasi-Monte Carlo with independent scramblings:
  <https://epubs.siam.org/doi/abs/10.1137/S0036142994277468>
- Polyak and Juditsky, averaging for stochastic approximation:
  <https://epubs.siam.org/doi/10.1137/0330046>

Projective integration literature permits acceleration only when the chosen coarse variables are
complete enough to close the coarse dynamics. The discovered six-dimensional hidden face space is
exactly the failure mode that must be repaired before such acceleration is meaningful:

- Equation-free projective integration for plasma PIC:
  <https://www.sciencedirect.com/science/article/abs/pii/S0021999107001726>

## Evidence artifacts

- `results/charging_c3_q1_compatibility_audit/audit.json`
- `results/charging_c3_q1_compatibility_audit/trajectory.csv`
- `results/charging_c3_q1_compatibility_audit/q1-compatibility-trajectory.png`
- `results/charging_c3_q1_compatibility_audit/fixed_state_ptc/`
- `results/charging_c3_periodic_topology_audit/audit.json`
- `results/charging_c3_periodic_topology_audit/periodicized_warm_proposal.npz`
- `results/charging_c3_periodic_topology_audit/periodic_topology_repair.png`
- `results/charging_c3_periodic_topology_audit/decreasing_gain_fixed_confirmation_audit/`
- `results/charging_c3_periodic_topology_audit/closure_mode_audit.json`
- `CHARGING_C3_CLOSURE_AUDIT_2026-07-14.md`

Every convergence number in this report is diagnostic evidence only. The final C3 claim still
requires B1, the approved B2 contract, B3 refinement, B4 ledgers, and an independent B5 exact-operator
audit.
