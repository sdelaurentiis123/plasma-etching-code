# Charging discretization design: geometry, field, and particles as one system

Date: 2026-07-11. Status: design gate for the next experimental implementation. The production tracer
and self-consistent solver are unchanged.

## Problem isolated by the reciprocity ladder

The frozen nonuniform AR4 audit separates source scoring from spatial discretization. With the derived
Liouville exit-state score and support-complete ion proposal, W32 ion backward/forward reciprocity closes
to +2.7% at adequate statistics (electron +0.6%). W16 electron error remains roughly +9--15%, depending
on the frozen field. Timestep reduction from 0.15 to 0.04 does not remove the error, while W16 to W32
spatial refinement does. The remaining leading error is therefore spatial/interface consistency.

Two attractive local fixes were implemented experimentally and rejected before commit:

1. A symmetric discrete-gradient mover preserved `v^2 + qV` to roundoff for a manufactured orbit, but
   bilinear interpolation through solid cell centers worsened W16 electron reciprocity from +8.6% to
   +18.4%. Exact conservation of the wrong reconstructed Hamiltonian is not physical accuracy.
2. A finite-volume `h/2` Dirichlet-face stencil, followed by its matching ghost-cell gradient, worsened
   the low-stat W16 electron audit to -17.9%. A face-correct field sampled as one nearest-cell vector is
   still not a compatible particle interpolation.

These failures rule out further correction factors on the current cell-centred mover.

## Required discrete physics

For the current rectilinear 2-D feature geometry, build a boundary-fitted nodal electrostatic system:

- gas-domain vertices carry `V`; material faces coincide exactly with grid edges;
- floating-surface voltage is imposed at the actual boundary vertices/edges, never at a covered-cell
  centre;
- the Laplace/Poisson weak form and particle field use one Q1 (bilinear quadrilateral) or P1 (triangular)
  basis;
- `E_p = -sum_a V_a grad(N_a)(x_p)` is evaluated only from the element containing the particle;
- a particle step is split at every element or material-face crossing, preventing tunnelling;
- absorption and deposited current are recorded at the exact face intersection;
- the adjoint uses the same orbit map in reverse and the already-derived exit-state source score;
- a later self-consistent charge deposition must use the same basis/face measure as interpolation.

This matches the central conservation result in modern PIC: field interpolation and deposition/field
discretization must be compatible, and discrete integration-by-parts identities matter. Energy
conservation alone is insufficient. Relevant primary sources:

- Markidis & Lapenta, *The Energy Conserving Particle-in-Cell Method*, arXiv:1108.1959.
- Chen, Chacon & Barnes, *An Energy- and Charge-conserving, Implicit, Electrostatic PIC Algorithm*,
  arXiv:1101.3701. Orbit substepping prevents tunnelling.
- Fichtl, Finn & Cartwright, *An Arbitrary Curvilinear Coordinate Method for PIC Modeling*,
  arXiv:1201.1476. Boundary-fitted coordinates and a compatible particle mover are treated together.
- Wang et al., *A parallel electrostatic PIC method on unstructured tetrahedral grids*, JCP 363,
  178--199 (2018), DOI 10.1016/j.jcp.2018.02.011. Vertex potential, finite-volume Gauss law, and
  arbitrary bounded geometry are verified as one algorithm.
- Wang et al., *Improved C1 shape functions for simplex meshes*, JCP 418, 109632 (2020). A potential
  basis with continuous electric field improves orbit fidelity on engineering geometries.
- AMReX embedded-boundary methods represent regular, cut, and covered cells separately and modify the
  elliptic stencil at the actual boundary; covered solid data are not treated as ordinary interpolation
  nodes.

## Implementation sequence and gates

1. Implement a standalone nodal Laplace solve for the existing axis-aligned AR geometry. Do not change
   `self_consistent_backward`.
2. Gate constant, linear, and separable harmonic potentials; require second-order potential convergence
   and the exact imposed face voltage.
3. Implement element-local field evaluation and exact DDA-style face crossing.
4. Gate Hamiltonian error, time reversal without absorption, and no tunnelling for grazing trajectories.
5. Run independent forward/backward electron and ion scoring in uniform fields.
6. Run the frozen AR4 W16/W32 audit. Promotion requires both species within 4%, monotone spatial
   convergence, and no regression when the timestep is halved.
7. Only then wire the nodal field through the deterministic fixed point and compare converged current
   residuals. HG remains a simulation reference, not a tuning target.

### First nodal implementation result (2026-07-11)

The standalone nodal Laplace solver and Q1 element-local tracer now pass parallel-plate potential/field,
uniform-field impact energy, one-cell no-tunnelling, diagonal two-face crossing, and manufactured linear-
field electron/ion reciprocity gates. The tracer resolves absorption at the first material-face event.

The frozen AR4 audit is not yet a joint pass. At W32 and the required high ion statistics, corrected ion
reciprocity is +1.6%, while electron reciprocity is -7.6%. A natural surface-Maxwellian proposal was an
unusable rare-event estimator on a +39 V floor. Replacing it with an exact support-complete mixture of
natural and barrier-shifted normal energies resolves the population and agrees with the legacy electron
gather (W16: -4.7% versus -4.2%; W32: -7.4% versus -7.6%). Thus proposal support is no longer the leading
electron suspect. Exact remainder-preserving lateral reflection fixes a real orbit bug but only improves
W16 by about 0.3 percentage point. The nodal field/tracer remains experimental and is not wired into the
charging fixed point. Adding four independent forward scrambles then exposed the apparent W16 electron
residual as reference-estimator uncertainty: the forward mean shifted from 0.0754 to 0.0733 with standard
error 0.00144, leaving legacy and shifted-Liouville adjoints at -1.4% and -1.9%. Corrected ion is +4.6%
with forward standard error 0.00109 and remains the near-gate statistics frontier. W32 must be repeated
with the same multi-scramble forward protocol before promotion.

The structured nodal implementation should remain embarrassingly parallel over particles and amenable
to Warp/CUDA. General 3-D curved geometry can later replace Q1 rectangles with AMReX-style embedded
boundaries or boundary-fitted simplex elements without changing the phase-space scoring derivation.

### Arbitrary-face support and estimator certification (2026-07-11 continuation)

Wiring the nodal tracer into a material-grid fixed point exposed a support error that floor-only gates
could not see. The adjoint surface proposal was expressed in global plasma-boundary coordinates on every
face. That is complete on a horizontal floor, but on a vertical wall global vertical velocity is
tangential and must span both signs. The old rule omitted valid upward-moving wall impacts and produced
forward/adjoint ion discrepancies as large as 8--11 replicate standard errors even when the timestep and
sample count were refined.

The proposal is now defined in each face's local tangent/inward-normal frame and rotated into global
coordinates with a unit Jacobian before time reversal. The physical density is still evaluated only at
the traced plasma exit. A zero-field vertical-wall gate checks the analytic `E[vx/vz]` Liouville factor
using an explicit local surface proposal. On the frozen nonuniform trench field, the former systematic
wall undercount disappears; remaining wall-adjoint uncertainty is proposal variance, and the independent
forward estimator is selected when it is better resolved.

Two statistical false-certification paths were also closed:

- method hysteresis may no longer retain an uncertified estimator when the complementary direction meets
  tolerance;
- when both directions claim precision, their cell currents must agree within a declared combined-error
  threshold, and bidirectional certification requires at least four independent replicates;
- direct equal-weight forward sampling carries a Bernoulli hit-count standard-error floor plus the
  existing zero-hit upper bound, preventing a few accidentally similar QMC scrambles from reporting
  implausible precision for rare deep-cell hits.

These are numerical support and evidence rules, not fitted physics. A corrected ten-step continuation
from the formerly failing state reduced certified current-balance RMS from 0.462 to 0.382 without an
estimator refusal. Full fixed-point convergence, initialization/damping invariance, and the AR/profile
ladder remain open and must precede any charging-validation claim.

### Compatible boundary-current unknowns (experimental, not promoted)

Longer continuation showed that independent covered-cell voltages averaged onto shared boundary nodes
still contain weak checkerboard modes. For the diagnostic trench, the map from cell voltages to actual
face-average boundary voltage has condition number about 183; directly inverting it would amplify current
estimator noise by more than two orders of magnitude. A separate nodal fixed-point candidate therefore
places dielectric voltage unknowns directly on the physical material-boundary vertices used by the Q1
field solve.

The first mass-lumped prototype assigned half of each face-total current to each endpoint. That is also
insufficient: on an open boundary chain, face-constant deposition has one alternating nodal null mode.
The tracer now returns the exact DDA hit face and intersection position; forward histories deposit through
the two linear endpoint shape functions, and adjoint face-position quadrature returns the same two endpoint
moments. The resulting endpoint-current basis has full row rank on the diagnostic geometry. Selected
replicate ensembles are deposited through this basis, so face covariance is measured rather than assumed.

This path is isolated in `charging_nodal_fixed_point.py`; the established solver has not been replaced.
On the same freestanding-wall numerical stress geometry, switching from face-constant to endpoint-resolved
deposition lowered the certified RMS residual from roughly 0.5 to 0.35 and the maximum from about 1.5 to
1.08 at the saved state. This is progress, not convergence. The filled-material trench is the physically
relevant geometry; checkpoints from the freestanding-wall and filled-material topologies must never be
interchanged, and the diagnostic harness now refuses such a mismatch.

A subsequent restart audit found that endpoint replicate arrays were cached before forward/adjoint cross
refinement. Face totals used the refined ensemble while nodal deposition used stale endpoint moments, so
an apparent pass was history dependent. Endpoint moments are now extracted only after all refinement; a
regression requires the selected endpoint mean to reflect the final refined face ensemble.

With that fix, the filled-material 20x18 trench (gas-only plasma row, ten-cell opening, fourteen-cell
depth), finite-transit 40+/-10 V RF sheath, 4 eV electrons, four QMC replicates, and level-14 adjoint ceiling
reached certified nodal current balance: max `|log(Gi/Ge)|` outside 2-sigma current intervals = 0.1428 and
RMS = 0.0382. Restarting the same evaluated state with beta changed from 0.25 to 0.1 re-certified at max
0.1053 and RMS 0.0266. The two results differ by at most 0.052 V on boundary nodes (0.0039 V RMS) against
surface potentials of order 27 V. This closes one nontrivial bulk-trench convergence/restart gate. It does
not close sample/grid/AR/initialization ladders, dielectric-volume permittivity/storage, SEE, or experiment.
The warm CPU path still costs roughly ten seconds per nonlinear iteration in this environment; CUDA
acceleration and proposal-variance reduction remain required product work.

### Exact Warp orbit backend (experimental)

The compatible Q1/midpoint/DDA orbit map now has a float64 Warp implementation. It preserves the four
midpoint iterations, adaptive step, first crossed-face ordering, lateral/bottom remainder reflection,
plasma exit state, impact energy, oriented hit normal, and exact hit position of the Numba reference.
Three parity gates cover wall/floor hits, exits/reflections, and nonuniform-field ion/electron trajectories:
all discrete outputs agree exactly and floating outputs agree to 2e-10 or tighter. The complete suite is
116 passing tests with Warp enabled. `PETCH_DEVICE=cpu` keeps the established Numba backend;
`PETCH_DEVICE=cuda` or `cuda:N` selects Warp. CUDA performance and full-solver CPU/GPU result parity are
still open and must be measured on an actual accelerator before any speed claim.

### Physical surface-charge / variable-permittivity Poisson mode (experimental)

The converged boundary-voltage candidate is a useful current-balance root solver, but dielectric voltage
is not an independent physical state. A dielectric stores free surface charge; its voltage follows from
the material permittivity, grounded/floating conductors, and Poisson's equation. The nodal candidate now
has an opt-in physical mode built on a reusable Q1 weak-form system:

- cellwise positive relative permittivity enters one compatible nodal stiffness matrix;
- dielectric state is nodal line charge in C/m (the 2-D per-unit-depth representation), with face sheet
  charge mass-lumped to the same endpoint basis as particle current;
- top plasma, grounded material, and floating-conductor nodes are explicit Dirichlet sets;
- one sparse factorization serves every field solve and supplies the exact diagonal Green response, in
  F/m, used to convert a resolved voltage correction into a physical charge correction; and
- every solve reports the free-node Poisson residual and electrostatic energy.

Manufactured gates reproduce a uniform parallel-plate voltage and a two-dielectric series capacitance to
machine precision, conserve lumped edge charge, and verify positive response capacitance. On a narrow AR1
filled trench with a three-cell SiO2 layer over a grounded bottom (a numerical gate, not an experimental
stack), twelve accepted high-resolution evaluations reduced certified current-balance RMS from 1.071 to
0.459 while the Poisson residual stayed below 8e-15 V. A longer pre-restart-contract continuation reached
0.326, but is not used as restart evidence because estimator state was not then serialized.

Checkpoint state now includes physical surface charge, adaptive forward/adjoint levels, estimator method
hints, and accepted gain age. A six-evaluation checkpoint and a fresh one-evaluation process replayed max
and RMS residual, currents, standard errors, method counts, potential, and electrostatic energy exactly.
This closes the restart invariant; it does not yet close Poisson-mode AR/grid/sample convergence.

On the filled-material trench, the endpoint-resolved candidate accepted 20 consecutive fixed-point
evaluations. Its certified current-balance RMS fell from 2.31 to a best value of 0.305 before fluctuating
at 0.35 as rare-hit sampling error became comparable to the remaining imbalance; the certified maximum
fell from 5.87 to a best value of 0.668. This establishes a stable descent path on the physical topology,
not final convergence. A higher-sample independent evaluation and AR/grid/initialization invariance remain
required before promotion.

The direct forward estimator is the physical support audit; the adjoint estimator is an importance-sampled
accelerator. When both are statistically admissible they must agree. When the adjoint misses a rare support
mode but the independent forward estimate satisfies its declared uncertainty tolerance, the solver now uses
the forward estimate and retains `adjoint_support_unresolved` as a diagnostic instead of rejecting valid
physics. It still refuses the cell when neither direction is certified or when two admissible estimators
disagree beyond the declared consistency threshold.

A zero-hit direct estimate carries a finite rule-of-three upper bound. The selector must preserve that bound
even when the adjoint also reports zero: an exact adjoint zero cannot erase finite forward uncertainty. This
case first became visible on the deepest electron-shadowed nodes of a narrow AR16 trench. After the fix, an
identical six-evaluation diagnostic (same finite-transit ion/electron boundary state and numerics; only trench
depth changed) reduced certified RMS current imbalance from 1.50 to 0.52 at AR1, 1.44 to 0.69 at AR4, and
0.82 to 0.53 at AR16. The AR16 sequence contained a recoverable sampling transient to 1.91, so this is a
cross-AR descent gate, not convergence or experimental validation. It motivates uncertainty-aware local step
control; it does not motivate an AR branch or a fabricated nonzero electron current.

## Signed error budget: what else can be wrong and how errors stack

The charging fixed point balances `Gi(V) = Ge(V)`. For small relative transport errors, the raw voltage
residual changes by

`delta r = delta log(Gi/Ge) ~= delta Gi/Gi - delta Ge/Ge`.

Thus equal-sign species errors cancel in current balance while opposite-sign errors add. Cancellation can
make voltage look correct while both fluxes are wrong. The fixed-point response further multiplies this by
the inverse slope/Jacobian of the coupled charging map; weakly observable or shadowed cells can turn a
small flux bias into a large voltage displacement. Every item below therefore needs a signed measurement,
not only a final floor-voltage comparison.

### Transport discretization errors

- **Boundary location:** cell-centred solid voltage displaces the electrical surface relative to the
  geometric hit face. This is currently the leading spatial suspect.
- **Field/particle incompatibility:** nearest-cell `E`, the Laplace stencil, and face absorption do not
  derive from one potential basis. A trajectory can conserve energy in a reconstructed potential yet be
  the wrong physical trajectory.
- **Post-step collision detection:** the tracer records the first solid cell after a step rather than the
  exact ray/face intersection. This perturbs impact energy, face identity near corners, and reversibility.
- **Element tunnelling and corner ordering:** one step can cross more than one grid face; the selected hit
  can depend on update order. Grazing wall trajectories are especially sensitive.
- **Time integration:** the production mover has bounded uniform-field energy error, but nonuniform-field
  time reversibility is not gated. The observed timestep invariance says this is not the dominant current
  error, not that it is zero.
- **Geometry staircasing:** W16 and W32 represent the same nominal dimensions with different discrete
  corners, conductor thicknesses, and sampled face sets. This changes both numerics and the represented
  physical geometry unless a common analytic boundary is used.

### Adjoint/source-estimator errors

- **Phase-space score:** the old ion score evaluated a 1-D implied source state instead of the actual exit
  state. The Liouville correction removes most of this bias.
- **Proposal support:** a 1-D shifted RF-energy proposal omitted valid 2-D surface states. The broad
  multiple-importance stratum closes support but converges slowly near RF-arcsine horns.
- **Quadrature variance mistaken for bias:** W32 ion error moved from +10% to +2.7% with one additional
  Sobol level. Scramble-to-scramble dispersion and a sample ladder must accompany every result.
- **Electron projection:** the source is a 3-D flux Maxwellian projected into 2-D dynamics. The discarded
  out-of-plane velocity is harmless only if its analytic marginalization and Jacobian are carried through
  the adjoint derivation. A naive Maxwellian exit-energy ratio was tested and rejected because it amplified
  the current orbit-energy/interface error.
- **Finite face quadrature:** face-centre launch produced 6--16% coarse-grid bias; Sobol integration across
  the physical face fixed that gate. Corners still contribute multiple oriented faces and require correct
  area pooling.
- **Finite source plane/domain:** reflecting lateral boundaries, source height, and simulated open-area
  width can redirect grazing particles. Existing source-height tests were inert, but boundary/domain
  convergence must be repeated after a new field discretization.

### Field and fixed-point errors

- **PDE residual versus transport consistency:** a small five-point Laplace residual only proves the
  current algebraic stencil is solved. It does not prove the stencil places material boundaries correctly
  or that particle interpolation is compatible.
- **Unconverged charging:** iteration 10 is a transient. Matching a reference there is not validation.
- **Damping and estimator noise:** nonlinear `log(Gi/Ge)` updates turn rare low-flux samples into large
  voltage steps. Clipping can stabilize the iteration while biasing the fixed point; it must not define the
  physics.
- **Corner pooling/conductor pooling:** a physical cell has one insulator voltage and a connected conductor
  has one equipotential. Treating faces independently creates contradictory currents; this bug was fixed,
  but the new discretization must preserve the invariant.
- **Inactive-cell observability:** when both species fluxes are below estimator resolution, local floating
  voltage is underdetermined. Arbitrary voltages there can still alter the global field and neighboring
  trajectories. This needs uncertainty-aware regularization derived from charge transport/capacitance, not
  a hidden clamp.

### Missing physical mechanisms that can masquerade as numerical error

- time-dependent RF sheath transit rather than instantaneous uniform phase;
- ion-neutral collisions and charge exchange between sheath edge and feature mouth;
- secondary-electron emission, electron-induced emission, and their energy/angular spectra;
- dielectric polarization, finite oxide thickness, trapped-charge depth, and material-dependent
  permittivity;
- surface/bulk charge conduction, leakage to substrate, and temperature-dependent relaxation;
- electron/ion reflection, neutralization, re-emission, and sputtered/redeposited charged products;
- evolving geometry and surface composition feeding back on charging during the etch cycle;
- reactor-scale nonuniform incident distributions and correlations between energy, angle, RF phase, and
  species flux.

These mechanisms must not be introduced to repair a failed reciprocity invariant. Reciprocity, geometry,
and conservation are numerical gates first. Only after they close should experimental discrepancies be
used to identify missing physics or infer declared material parameters.

### Required anti-cancellation reporting

For every charging validation point, report at least: forward and backward `Gi` and `Ge` separately;
signed reciprocity error for each species; Sobol scramble uncertainty; W and timestep ladders; PDE and
current-balance residuals; iteration count; floor/sidewall/conductor potentials; and whether the observable
is analytic, another simulation, calibrated experiment, or held-out experiment. A scalar floor voltage or
notch depth alone is insufficient evidence.
