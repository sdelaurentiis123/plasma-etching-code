# Charging convergence-contract amendment

Date drafted: 2026-07-13

Revision identifier: **CCA-2026-07-13-R2**

Status: **APPROVED AND IN FORCE**

Required approver: Stan

Approval: **signed by Stan on 2026-07-13 after incorporation of both binding riders**

This document establishes a refinement-based replacement for the former frozen-trench charging
convergence decision. It does not silently modify `END_STATE_VERIFICATION.md` or loosen a runtime
threshold. Its signature authorizes C3 subject to the entry gates and no-gos below; it does not
authorize C4--C6 out of order. Per-node RMS and worst-node imbalance remain mandatory reported
diagnostics but are no longer the sole convergence decision.

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

## Amended contract

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
- The patch statistic must be reported at no fewer than two patch scales. For any claim concerning
  a specific profile feature, at least one reported patch scale must not exceed the physical extent
  of that feature.

### B3 — observable invariance

At minimum report floor potential and every profile observable used by the claim, such as notch
depth or bow width. Each must change by less than a preregistered dimensional tolerance under:

1. charged-transport sample-level doubling;
2. charging-timestep halving; and
3. one spatial-grid refinement.

The tolerance, extraction algorithm, and comparison checkpoint are declared in the run manifest.
For any checkpoint supporting an experimental claim, the declared observable tolerance must not
exceed the corresponding experimental uncertainty (including digitization uncertainty) of the
benchmark being claimed.

### B4 — exact conservation

- The signed surface-transfer ledger closes at every charged-particle cascade.
- Compatible Q1 projection closes at every charging step.
- Moving-surface remap conserves retained positive and negative charge inventories separately and
  reports charge removed with etched material separately; net-charge cancellation is not a valid
  conservation scale.
- No bounce cap, topology event, or remap failure may silently discard charge.

### B5 — independent exact-operator final audit

- Use independent scrambles and a higher declared sample level than the evolution run.
- Use exact hard visibility and the full declared kinetic material-response operator.
- Report RMS and worst-node imbalance, fixed-patch imbalance, potential-saturation measures,
  estimator uncertainty, seeds, sample levels, and all conservation ledgers.

## Retained legacy diagnostic

Per-node RMS and worst-node relative current imbalance are computed and reported forever. Under this
signed amendment they are mesh-scale diagnostics rather than the sole convergence decision. They may
not be hidden, renamed, smoothed, or omitted from a failed run. Historical results remain evaluated
under the contract active when they were produced.

## Validity boundary

The quasi-static saturation contract applies only when the charge-relaxation time is well separated
from profile motion and the applied bias waveform. If a pulsed or chopped waveform is not slow relative
to charging saturation, the validity layer must refuse quasi-static output or switch to the same
physical-time operator in waveform-resolved co-simulation mode.

## Required run record

Every run records engine version and git revision, input checksums, geometry and material manifests,
scramble mode and seeds, estimator method and sampling provenance, timestep/grid/sample settings,
charge and mass ledgers, B1--B5 values, retained node diagnostics, wall time, and hardware. Every
reflection, SEE, and conduction parameter is recorded with its value, source, and declared bounds.

## Authorization boundary under signature

Signature of this revision authorizes C3 according to the charging co-evolution handoff, including
wiring `remap_surface_charge_3d` and the certified reflection channel into the shared co-evolution
driver. It does not relax task-entry gates or any handoff no-go. In particular:

- per-node RMS and worst-node diagnostics are reported on every run, including failures;
- no parameter may be tuned on a Nozawa held-out target under any circumstances;
- quasi-static output is refused for pulsed-bias regimes unless timescale separation is demonstrated;
  waveform-resolved co-simulation is the permitted path otherwise;
- every reflection, SEE, and conduction parameter appears in the run manifest with its source and
  bounds;
- PTC may accelerate only the same physical-time charge ODE and must pass the same final exact-operator
  audit;
- positive and negative retained-charge inventories remain separately conserved and reported; and
- incomplete charged-particle bounce cascades are refused and may not be treated as closed.

## Sign-off

By signing below, Stan authorizes revision `CCA-2026-07-13-R2` to replace the frozen-trench per-node
criterion as the charging convergence decision for C3 and later tasks, subject to every task-entry
gate and authorization boundary above. It does not waive any B1--B5 gate.

Approved by: **Stan**

Date: **2026-07-13**

Revision approved: **CCA-2026-07-13-R2**
