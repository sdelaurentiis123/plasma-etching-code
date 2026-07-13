# Charging convergence-contract amendment

Date drafted: 2026-07-13  
Status: **DRAFT — NOT IN FORCE**  
Required approver: Stan  
Approval: **not yet provided**

This document proposes a refinement-based replacement for the present frozen-trench charging
convergence contract. It does not modify `END_STATE_VERIFICATION.md`, does not loosen a runtime
threshold, and does not authorize the C3--C6 co-evolution campaign. The current per-node RMS and
worst-node criterion remains controlling until Stan explicitly signs this amendment.

## Governing target

Surface charge remains governed by

    d sigma / dt = J_i(sigma, Gamma) - J_e(sigma, Gamma),

where `Gamma` is the current material geometry and both currents are evaluated by the exact
hard-visibility kinetic transport operator in the self-consistent Q1 Poisson field. The proposed
contract identifies equilibrium with a timestep- and grid-refined stationary state of this
conservative ODE. Pseudo-transient continuation may accelerate approach to that state but may not
define a different residual or final operator.

## Grounds for amendment

All three grounds are required and have been recorded:

1. **Estimator/refinement evidence.** Task 0A measured a minimum small-radius
   signal-to-between-scramble-error ratio of 0.51 at level 13, below the required value 3. The raw
   per-node response therefore contains unresolved quadrature and mesh-scale hit structure at the
   scale used by the old diagnostic.
2. **Validation-practice evidence.** Published feature-charging/profile models reach charging
   saturation through transient accumulation and validate emergent profile or potential observables;
   they do not use machine-scale balance on every frozen-geometry node as the experimental claim.
3. **Strictly checkable replacement.** The proposed gates require timestep, sample, and physical-grid
   refinement, exact ledgers, and an independent exact-operator audit. They replace an unresolved
   discretization statistic with quantities that must tighten under refinement rather than merely
   removing a difficult test.

## Proposed amended contract

Every charging checkpoint used for a co-evolved profile claim must pass all five gates.

### B1 — potential saturation

- Report `max_n |dV_n/dt|` over a declared terminal time window.
- State the dimensional tolerance and time window in the run manifest before the run.
- Timestep halving must preserve the saturation decision and terminal potential within the declared
  observable tolerance.

### B2 — fixed-physical-patch current balance

- On patches of fixed physical size, report

      |integral_patch (J_i - J_e)| / |integral_patch J_i| <= 0.08.

- Patch boundaries and physical dimensions are manifest inputs, not grid-index selections.
- The patch statistic must tighten under at least one spatial-grid refinement. Report the measured
  convergence order; non-tightening results fail this gate and trigger a discretization audit.

### B3 — observable invariance

At minimum report floor potential and every profile observable used by the claim, such as notch
depth or bow width. Each must change by less than a preregistered dimensional tolerance under:

1. charged-transport sample-level doubling;
2. charging-timestep halving; and
3. one spatial-grid refinement.

The tolerance, extraction algorithm, and comparison checkpoint are declared in the run manifest.

### B4 — exact conservation

- The signed surface-transfer ledger closes at every charged-particle cascade.
- Compatible Q1 projection closes at every charging step.
- Moving-surface remap reports retained charge and charge removed with etched material separately.
- No bounce cap, topology event, or remap failure may silently discard charge.

### B5 — independent exact-operator final audit

- Use independent scrambles and a higher declared sample level than the evolution run.
- Use exact hard visibility and the full declared kinetic material-response operator.
- Report RMS and worst-node imbalance, fixed-patch imbalance, potential-saturation measures,
  estimator uncertainty, seeds, sample levels, and all conservation ledgers.

## Retained legacy diagnostic

Per-node RMS and worst-node relative current imbalance are computed and reported forever. After this
amendment is signed they become mesh-scale diagnostics rather than the sole convergence decision.
They may not be hidden, renamed, smoothed, or omitted from a failed run. Historical results remain
evaluated under the contract active when they were produced.

## Validity boundary

The quasi-static saturation contract applies only when the charge-relaxation time is well separated
from profile motion and the applied bias waveform. If a pulsed or chopped waveform is not slow relative
to charging saturation, the validity layer must refuse quasi-static output or switch to the same
physical-time operator in waveform-resolved co-simulation mode.

## Required run record

Every run records engine version and git revision, input checksums, geometry and material manifests,
scramble mode and seeds, estimator method and sampling provenance, timestep/grid/sample settings,
charge and mass ledgers, B1--B5 values, retained node diagnostics, wall time, and hardware.

## Sign-off

By signing below, Stan authorizes this amendment to replace the frozen-trench per-node criterion as
the charging convergence decision for C3 and later tasks. It does not waive any B1--B5 gate.

Approved by: ____________________  
Date: ____________________  
Revision approved: ____________________

