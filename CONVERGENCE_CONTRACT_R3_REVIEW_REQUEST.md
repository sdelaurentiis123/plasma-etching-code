# Charging convergence-contract R3 review request

Date drafted: 2026-07-14

Proposed revision identifier: **CCA-2026-07-14-R3-DRAFT**

Status: **NOT IN FORCE — EXPLICIT STAN SIGN-OFF AND REFINEMENT EVIDENCE REQUIRED**

CCA-2026-07-13-R2 remains the governing contract. This document proposes a narrow B2 clarification
after the Q1 compatibility audit demonstrated that raw patch current is not a single-valued
functional of the discrete electrostatic state.

## New evidence requiring review

The initial audit used the historical nonperiodic Q1 operator, whose conservative P0-face to
Q1-node coupling had 40 columns, rank 34, and nullity 6. That audit found up to 0.426 dual
sensitivity of a raw patch sum to exact null modes and a coherent null-current component across 15
independent terminal windows. It also led to discovery of a more fundamental topology defect:
particle transport wrapped laterally while Poisson did not.

On the repaired intended periodic grid, the common operator has 72 independent volume nodes and a
`72 x 40` surface coupling of rank 22 and nullity 18. The declared 0.25/0.50 micrometer patch
functionals have maximum null sensitivities 0.790/0.498. Compatible projection preserves periodic
nodal charge, potential, and total charge to roundoff and leaves exactly zero voltage seam. Four
corrected 50 microsecond physical windows and a post-tail fixed confirmation retain null inventory
near `1.6e-15` while the raw/compatible B2 discrepancy remains large.

Therefore the current R2 raw B2 numerator can vary while the field, trajectories, and all future
electrostatic feedback remain identical. Requiring the field to drive that component to zero is not
a well-posed equilibrium condition for this discretization.

This evidence does **not** establish stationarity under the repaired operator. The post-tail fixed
confirmation still has integrated B1 `17,748.6 V/s`, raw B2 `1.348/1.297`, and field-compatible B2
`0.1681/0.1677`, and it continues in the same dominant charging direction. R3 therefore remains a
review request rather than a route for declaring the current checkpoint converged.

## Proposed B2 replacement

Replace only the B2 decision text with all of the following. Every other R2 clause remains binding.

### B2a — field-compatible fixed-patch balance

For each declared physical patch, conservatively project the net kinetic face current through the
same face-to-Q1 coupling used by Poisson and back to its unique area-weighted compatible face
representative. Report

    |integral_patch P_Q1(J_i - J_e)| / |integral_patch J_i| <= 0.08,

at no fewer than two patch scales. Patch construction, ion denominator, hard-visibility kinetic
operator, fixed physical sizes, and feature-scale rider remain exactly as in R2. The projection is
not smoothing and may not change the nodal current or potential response.

### B2b — unresolved-current discretization budget

Also report the original raw R2 statistic and the difference between raw and compatible patch
numerators at both patch scales. For any experimental claim:

- the raw-versus-compatible difference must tighten under one physical grid refinement;
- its effect on every B3 observable must lie inside the declared grid-discretization contribution
  to the uncertainty budget; and
- no patch with structural null sensitivity may be hidden, merged, or excluded after seeing the
  result.

Failure to tighten blocks the claim and triggers the discrete-equilibrium audit. It does not permit
loosening 0.08, adding conductivity/SEE without causality evidence, or dropping the raw diagnostic.

## Unchanged clauses

- B1, B3, B4, and B5 are unchanged.
- Raw B2, node RMS, and worst-node imbalance are reported forever, including failures.
- Exact hard visibility, independent final scoring, timestep/sample/grid refinement, and all charge
  ledgers remain mandatory.
- C4 remains blocked until C3 passes the signed contract.

## Evidence still needed before signature

1. **Complete:** field-compatible patch metrics are present in every C3 history, heartbeat, summary,
   and current audit while raw B2 remains reported.
2. Run the continued compatible physical-time trajectory to a statistically stationary B1 window;
   the July 14 post-tail confirmation explicitly did not satisfy this item.
3. Repeat the final current audit on one finer physical grid and quantify both B2a and B2b.
4. Confirm B3 floor potential and any claimed profile observable within the preregistered grid and
   experimental uncertainty.

## Sign-off

No signature is requested until the four evidence items above exist.

Approved by: **PENDING**

Date: **PENDING**

Revision approved: **PENDING**
