# petch end-state verification contract

This file defines the slice that must be true before petch is described as a validated product
foundation. It is a completion contract, not roadmap copy. A green unit-test suite alone does not close
it, and an entry remains open unless the cited artifact measures the stated quantity on the current code.

## Product slice and boundaries

The first shippable slice is a **SiO2-first, feature-scale forward solver with physical reactor boundary
conditions**:

- Inputs: a 2-D/3-D material geometry; dimensional species fluxes and joint energy-angle-phase
  distributions; material permittivity/conductivity; stateful surface-reaction mechanisms with units,
  provenance, and uncertainty; and numerical accuracy targets.
- Governing core: conservative species transport, compatible electrostatics and dielectric/conductor
  charging, stateful surface balances, and conservative interface motion.
- Outputs: evolved profile; charge/potential/field; species-resolved surface flux, energy and angle;
  surface state; numerical/parameter/model uncertainty; convergence diagnostics; and an explicit validity
  decision that refuses unsupported extrapolation.
- Product wedge: calibrate uncertain physical inputs on structure N and predict held-out structure N+1.
  Verified gradients should make this data-efficient; they do not by themselves constitute inverse design.

The slice does **not** claim a predictive reactor model, arbitrary-material chemistry from recipe names,
closed-loop target-profile-to-recipe optimization, or learned replacement of known conservation laws.
Those are later layers connected through the same boundary and surface-mechanism contracts.

## Meaning of “general”

There is one transport/electrostatics/interface engine. Aspect ratio, feature name, benchmark identity,
surface region, or expected answer may not select a physical formula or estimator. Geometry and physical
material/source/mechanism inputs may change. Adaptive numerical work may change only from measured error.

Chemistry is a physical input module, not a universal hard-coded reaction law. Generality is earned when a
second chemistry/material mechanism runs through the unchanged transport, electrostatic, state, geometry,
UQ, and product contracts. Until that replay passes, “chemistry-extensible” is accurate and
“all-chemistry predictive” is not.

## Completion gates

| Gate | Evidence required for closure | Current evidence | Status |
|---|---|---|---|
| One physical input contract | Serialized geometry/materials, `PlasmaBoundaryState`, surface mechanism, units, provenance, and uncertainty replay identically on CPU/GPU | Boundary state exists; high-level API does not yet carry materials, charging, mechanism provenance, or validity | **open** |
| Geometry/AR generality | Same core on open wafer, trench/hole, symmetric/asymmetric pattern, AR ladder, and evolving geometry; no benchmark/AR/region branches; grid and domain convergence | Transport architecture gate spans AR1/4/16; six-step nodal charging descent exists; converged AR/grid/init ladder does not | **partial** |
| Charged-particle/neutral transport | Conservation, no tunnelling, energy error, forward/adjoint reciprocity, collision/re-emission balance, and sample/grid refinement with uncertainty | CPU/CUDA orbit parity is exact on discrete events. A high-sample AR4 audit found state-dependent stepping biased the nonuniform-field adjoint; a common fixed-step ladder removed the observed 8.38-sigma disagreement, and such cases now refuse adaptive stepping. Full step/time-horizon/AR convergence plus reflection/re-emission and collisional transport remain incomplete | **partial** |
| Electrostatics and charging | Compatible deposition/field bases, variable permittivity, conductor equipotential, dielectric storage/conductivity, charge conservation, converged current balance, initialization/grid/sample invariance | Compatible nodal endpoint deposition and checkpoint replay exist, but the former convergence claim was withdrawn because uncertainty-interval overlap had been mistaken for a narrow confidence envelope. The solver now separates update direction from a strict convergence envelope. Physical Q1 charge mode passes capacitance, displacement-continuity, Gauss-law, and restart gates; AR/grid/sample convergence, conductivity, and SEE remain open | **partial** |
| Surface chemistry/state | SiO2 fluorocarbon mechanism with conserved coverages/film/mixed-layer state and sourced energy-angle yields; uniform-flux analytic limits; independent chemistry replay through unchanged core | A new dimensional reduced SiO2/fluorocarbon kernel conserves complex sites and finite polymer inventory, implements the sourced threshold energy-angle law, carries parameter provenance/uncertainty, and refuses unsupported incident species. It explicitly reports omitted crosslinking, resolved complex stoichiometry, F sequencing, redeposition, and mask chemistry; it is not yet coupled to the 3-D engine or experimentally replayed | **partial** |
| 3-D evolving profile | Mass/interface conservation, time-step/grid/ray convergence, topology robustness, charging feedback, and extracted profile metrics | Fast 3-D level-set path and selected parity tests exist; production 3-D charging and complete convergence suite do not | **partial** |
| Notching/bowing/twisting | Charging-off/on causality; absolute profile comparison; symmetric ensemble has zero mean twist with nonzero variance; dense/sparse mean direction and domain convergence | Notch mechanism and frozen pattern sign gates exist; absolute notch energy/depth, stochastic twist, and evolving 3-D coupling remain open | **open** |
| Experimental validation/transfer | Raw-source provenance; preregistered calibration/held-out split; replayable extraction; measurement/digitization/numerical/parameter/model errors separated; held-out profile/AR/recipe predictions | de Boer, Bosch, ALE/cryo/Bosch anchors and literature targets exist at mixed strength; no unified SiO2 held-out transfer result | **open** |
| Differentiability | Finite-difference/complex-step agreement through transport, converged charging fixed point, surface state, and profile loss; stable implicit solve; calibration transfer demonstration | Reduced ALE gradients pass; unified feature and charging gradients do not | **open** |
| Performance | Reproducible warmed CPU/GPU benchmarks; accuracy-matched competitor comparison; scaling with cells/faces/AR/species; device-residency and memory profile; deterministic replay | Exact Warp CUDA parity is gated on an A100; one-off end-to-end tracer speedups were 2.7x-9.8x, but the Python fixed-point loop underutilizes the GPU and no committed benchmark manifest or accuracy-matched competitor result exists | **open** |
| Product behavior | Config-in/result-out API exposes all physical inputs and diagnostics; schema/versioning; checkpoint/restart; invalid-domain refusal; fresh-environment install and examples | `Domain/SF6O2/Process/Result` exists, but lacks the unified physics contract, validity result, UQ, and solver checkpoint schema | **partial** |

## Required validation layers

Every promoted mechanism must pass all applicable layers:

1. **Analytic/manufactured:** normalization, conservation, symmetries, limiting cases, and order of
   convergence.
2. **Numerical cross-check:** independent forward/adjoint or discretization comparison with uncertainty;
   grid, time-step, sample, domain, and initialization ladders.
3. **Mechanism causality:** disabling one physical channel removes its predicted effect without changing
   unrelated physics.
4. **External comparison:** independent code parity where useful, clearly distinguished from experiment.
5. **Experiment:** dimensional profile/rate/state observables under matched conditions, with declared
   calibrated inputs and held-out predictions.
6. **Transfer and validity:** failures outside the supported domain are reported, not silently refitted.

Numerical tolerances must be tighter than the measurement or model uncertainty used in the corresponding
external gate. Each experimental gate owns a preregistered metric and tolerance; simulation benchmarks such
as Hwang-Giapis and ViennaPS are never labeled experimental evidence.

## Completion evidence bundle

The final local release candidate must contain:

- a machine-readable case/config schema and result/validity schema;
- one command that runs the complete verification manifest and emits a signed PASS/FAIL table;
- raw-data provenance and immutable checksums for redistributable data, or reproducible extraction
  instructions where copyright prevents vendoring;
- stored calibration/held-out splits and uncertainty budgets;
- CPU and GPU accuracy/performance reports with hardware/software metadata;
- a clean worktree whose local commits contain every intended change; and
- documentation whose public claims are generated from, or manually reconciled against, the manifest.

The goal is complete only when every table row above is closed by current evidence. Until then the correct
claim is the strongest passing subset, with the remaining rows named explicitly.
