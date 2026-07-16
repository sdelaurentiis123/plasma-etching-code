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
| One physical input contract | Serialized geometry/materials, `PlasmaBoundaryState`, surface mechanism, units, provenance, and uncertainty replay identically on CPU/GPU | `PlasmaBoundaryState` now feeds arbitrary-triangle ballistic and nodal-field 3-D adapters that enforce its SI reference plane and preserve dimensional species/energy/angle events through the surface mechanism. Every built-in continuous boundary density also exposes one deterministic inverse-CDF sampling contract. Serialized materials, reflection, CPU/GPU end-to-end replay, and the high-level product schema remain open | **partial** |
| Geometry/AR generality | Same core on open wafer, trench/hole, symmetric/asymmetric pattern, AR ladder, and evolving geometry; no benchmark/AR/region branches; grid and domain convergence | Transport architecture gate spans AR1/4/16; six-step nodal charging descent exists; converged AR/grid/init ladder does not | **partial** |
| Charged-particle/neutral transport | Conservation, no tunnelling, energy error, forward/adjoint reciprocity, collision/re-emission balance, and sample/grid refinement with uncertainty | CPU/CUDA orbit parity is exact on discrete events. A high-sample AR4 audit found state-dependent stepping biased the nonuniform-field adjoint; a common fixed-step ladder removed the observed 8.38-sigma disagreement, and such cases now refuse adaptive stepping. The 3-D fixed-step Hamiltonian path reproduces ballistic zero field and electrostatic work under refinement. A joint position–velocity Sobol path samples analytic Maxwellian, ion-energy/transverse-Maxwellian, reactor/PIC histogram, and mixture densities in `N` rays rather than an `Nv×Nx` tensor; it resolves the analytic Maxwellian barrier tail to its `1/N` digital-net bound. Replicated scrambled nets propagate signed-current standard error into a strict charging confidence envelope. The charging solve raises its nested Sobol level until estimator uncertainty meets tolerance or the declared maximum is exhausted. Deterministic diffuse transport uses source/target area reciprocity, refuses non-closing form factors, and conserves source = reacted + escaped to about 1e-13. The same operator closes surface-product direct escape, first impact, repeated reaction/re-emission, and final escape. An opt-in same-material redeposition contract now requires bounded/sourced sticking and density, converts captured material to signed interface growth, and refuses cross-material films, pinned growth, and missing launch laws. The charged response cascade closes every landing/escape/unresolved event and merges re-impacts into chemistry. Bounded grazing-specular reflection passes particle/charge/energy/path-reversal gates but remains a sensitivity. Full AR/grid evolution, calibrated reflection, collisional transport, reactive-product branching, and cross-material coatings remain incomplete | **partial** |
| Electrostatics and charging | Compatible deposition/field bases, variable permittivity, conductor equipotential, dielectric storage/conductivity, charge conservation, converged current balance, initialization/grid/sample invariance | Compatible nodal endpoint deposition and checkpoint replay exist, but the former convergence claim was withdrawn because interval overlap had been mistaken for a narrow envelope. Physical Q1 charge passes 2-D/3-D capacitance, displacement-continuity, and Gauss-law gates. Physical time closes charge→field→trajectory→signed current→charge with exact global conservation. Moving-surface remap passes no-op, translation, removal-ledger, first-order refinement, and zero-inventory gates. C3 now wires physical time, response cascades, material-routed chemistry, signed profile motion, state/charge remap, continuation, schema-versioned manifests, and a recovery/error budget. Signed R2 remains the default; the separate two-block profile-relevant contract retains every R2 diagnostic and uses disjoint endpoint ensembles but cannot support an experimental claim. Mixed conductor/dielectric surfaces, conductivity, calibrated electron SEE, co-evolution ladders, and experimental closure remain open | **partial** |
| Surface chemistry/state | SiO2 fluorocarbon mechanism with conserved coverages/film/mixed-layer state and sourced energy-angle yields; uniform-flux analytic limits; independent chemistry replay through unchanged core | A dimensional reduced SiO2/fluorocarbon kernel conserves complex sites and finite polymer inventory, implements the sourced threshold energy-angle law, carries parameter provenance/uncertainty, and refuses unsupported incident species. Its validity result now distinguishes executable manufactured/analytic parameter evidence from evidence that supports prediction in a declared physical domain; the current SiO2 test parameter set is explicitly nonpredictive. Every reduced SiO2 and tabulated Si step now returns an exact material-origin exchange ledger: removed material equals routed outgoing plus chemically unresolved inventory face-by-face, and unresolved products are explicitly ineligible for redeposition. This closes bookkeeping without inventing SiF4/CO/CO2/COF2 branching. A general versioned surface-interaction table supports units-explicit axes, uncertainty, exact serialization replay, leave-one-out interpolation audits, machine-precision endpoint handling, and default extrapolation refusal. The CC-BY OSTI 2589032 archive is checksum-pinned; its byte-exact Si-Cl2-Ar+ sputter, RIE, and ALE product tables replay with MD uncertainty and fixed-condition provenance. A sourced Si-Cl2-Ar+ mechanism runs two evolving 3-D steps through the unchanged transport, species filtering, generic conservative state remap, and level-set engine; at the released 100 eV, normal-incidence, 10:1 flux-ratio node its mean velocity matches the archived yield within 1% and reports prediction-supporting parameter evidence. It refuses unreleased energy, angle, ratio, and species. This demonstrates chemistry extensibility in one narrow second-chemistry domain, not broad Si prediction; SiO2 crosslinking, resolved complex stoichiometry, F sequencing, product branching/transport, and mask chemistry remain explicit omissions | **partial** |
| 3-D evolving profile | Mass/interface conservation, time-step/grid/ray convergence, topology robustness, charging feedback, and extracted profile metrics | The dimensional path carries boundary events through conserved removal/deposition into material-gated signed level-set motion. Material-local remap preserves state to machine precision, fingerprints each handoff, and refuses topology/material/displacement violations. The material router independently advances mask and substrate laws, namespaces their conservative state, merges their ledgers, recomputes material ownership after motion, and round-trips it through checkpoints. Same-material product return can now oppose etch velocity as real growth. Moving-plane and charged planar gates retain their prior analytic behavior. Geometry-native notch/bow and hole centerline/twist measurements pass manufactured refinement gates. Production mask chemistry, cross-material film topology, topology-event physics, and full AR/grid/time/ray convergence remain open | **partial** |
| Notching/bowing/twisting | Charging-off/on causality; absolute profile comparison; symmetric ensemble has zero mean twist with nonzero variance; dense/sparse mean direction and domain convergence | Common C3 coupling exists. Dimensional Poisson arrivals preserve integer count/flux/charge identities; the physical-time ensemble generates independent profiles. Geometry-native left/right notch, bow, centerline, displacement, diameter, and onset-AR metrics pass manufactured gates. C4 packaging enforces checksum/pixel digitization provenance, a maximum-two-parameter committed calibration family, exact held-out coverage, and an uncertainty ladder. C5 packaging executes nested N/2N, paired sample doubling, isotropy, and systematic-direction checks. The actual Nozawa image ingest, held-out solve, N>=30 AR sweep, and stochastic response-branching gate have not run | **open** |
| Experimental validation/transfer | Raw-source provenance; preregistered calibration/held-out split; replayable extraction; measurement/digitization/numerical/parameter/model errors separated; held-out profile/AR/recipe predictions | Krüger calibration/held-out trend facts and both checksum-verified Zenodo Bosch wafer tables are ingested with evidence labels. The 89-position table contributes 7,832 measurements over 88 wafers and preserves 157 original missing values separately from processed thickness; it is reactor/wafer-scale depth/uniformity evidence, not feature-profile validation. Jeon 2022 contributes 54 checksum-verified experimental SiO2 trench depths over six widths. Pixel coordinates replay through stored axis maps, duplicate controls from three independently digitized panels agree within the 35 nm digitization budget, and only the 20% C4F8 continuous-wave curve is preregistered for calibration. Every condition maps to a separately checksum-verified diagnostic-derived neutral/ion control; these are not direct surface-flux measurements. Same-condition width ratios are exposed because coupon exposure cancels algebraically. Pulse/CW ratios are withheld by default because Jeon did not report whether wall time or cumulative RF-on time was shared. A corrected nominal-width one-step diagnostic now scores all six CW widths: non-anchor interval coverage is 4/5 at 20%, 4/5 at the untouched 40% condition, and 2/5 at the untouched 80% condition. These sub-nanometre increments test initial transport/rate shape, not the roughly micrometre-scale evolved experimental profiles, so they do not close validation. Jeong 2023 supplies a checksum-verified 20-minute absolute-time matrix with one anchor and 17 held-outs, but its modeled radical densities still require a reactor-to-reactive-surface flux closure. No full-depth unified-engine held-out prediction has passed | **open** |
| Differentiability | Finite-difference/complex-step agreement through transport, converged charging fixed point, surface state, and profile loss; stable implicit solve; calibration transfer demonstration | Reduced ALE gradients pass; unified feature and charging gradients do not. Deterministic per-cell gather is necessary but not sufficient: estimator selection and adaptive sample levels are now frozen within derivative/root epochs, yet hard hit/escape indicators still move discontinuously with field and geometry for a finite sample set. The latest all-47-column central AR4 map at a 0.05 V radius has condition number about 15,300 (the 0.025 V map was about 72,000), so naive tracer autodiff is not promoted. The gradient gate requires a verified discontinuity treatment—boundary/edge sampling, integral reparameterization, or an equivalent unbiased sensitivity estimator—before implicit differentiation of the fixed point can be claimed | **open** |
| Performance | Reproducible warmed CPU/GPU benchmarks; accuracy-matched competitor comparison; scaling with cells/faces/AR/species; device-residency and memory profile; deterministic replay | Exact Warp CUDA parity is gated on an A100; one-off end-to-end tracer speedups were 2.7x-9.8x, but the Python fixed-point loop underutilizes the GPU and no committed benchmark manifest or accuracy-matched competitor result exists | **open** |
| Product behavior | Config-in/result-out API exposes all physical inputs and diagnostics; schema/versioning; checkpoint/restart; invalid-domain refusal; fresh-environment install and examples | The legacy `Domain/SF6O2/Process/Result` remains labeled. Explicit common-engine APIs now cover uncharged, charged C3, and physical finite-arrival ensembles. Charged continuation preserves geometry/surface-state fingerprint/remapped charge and derives a non-reused seed. Versioned safe NPZ checkpoints support registered conservative surface-state types, reject undeclared arrays/pickle loading, and bind the source manifest by SHA-256. Full config serialization, schema migration, CLI examples, and fresh-environment replay remain open | **partial** |

The Jeon intended-use gate is versioned as `jeon_2022_depth_transfer_v1` and was frozen before the
first unified-engine prediction. It scores five nontrivial calibration ratios, 40 held-out
same-condition width ratios, digitization-only interval coverage, numerical uncertainty, parameter
count, validity, wall time, and accelerated-compute cost separately. Its 24 cross-exposure pulse
ratios and twelve 1 ms reversal directions remain disabled until the exposure basis is independently
resolved; an explicitly labeled development hypothesis may report them but cannot support validation.
Twelve additional electron-density/self-bias controls replay from their source pixels;
independently digitized CW duplicates agree within 4e14 m^-3 and 20 V budgets. A boundary adapter
reproduces the paper's assumed Bohm flux while requiring IEDF and radical composition as explicit
closures. Self-bias is retained as an energy-scale diagnostic, not mislabeled as an IEDF. Shallow
held-out diagnostics have run; no full-depth held-out profile prediction has passed, so the
experimental gate remains open.

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
