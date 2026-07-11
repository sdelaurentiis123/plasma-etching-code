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
