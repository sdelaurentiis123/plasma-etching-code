# OSTI / Oehrlein problem map

Research pass: 2026-07-11. This memo separates the papers' actual problem statements from the
claims petch can currently make, then turns each problem into a falsifiable implementation gate.
It is intentionally architecture-facing: feature physics comes first, while the source-plane and
surface-feedback contracts remain compatible with a later reactor-scale model.

## Executive conclusion

The papers describe one connected multiscale problem, not a list of unrelated features:

1. A reactor and sheath produce species-resolved fluxes and joint energy-angle distributions at the
   wafer/source plane.
2. Feature geometry filters those distributions and creates nonlocal electrostatic coupling between
   neighboring features.
3. Each impact changes a surface whose coverage, composition, damage, and mixed-layer depth retain
   history; the altered surface changes later reaction probabilities.
4. Surface consumption, product emission, and wall conditioning feed back to the reactor on a slower
   time scale.

No single mesh should span this entire range. The DOE 2023 semiconductor-plasma report explicitly
recommends hybrid time/spatial slicing and sub-cycling: import reactor fluxes into the feature model
and return effective surface reaction data to the reactor model. That is the right petch boundary.

The first crack is therefore not a reactor rewrite. It is to make the feature solver accept physical
source distributions and MD-derived, stateful surface mechanisms without changing its transport or
geometry code.

## 1. Oehrlein et al. 2024: the field-level problem statement

Source: G. S. Oehrlein et al., *Future of plasma etching for microelectronics: Challenges and
opportunities*, JVST B 42, 041501 (2024), DOI 10.1116/6.0003579, OSTI 2560935.

### Actual problem

Future patterning simultaneously demands atomic-scale precision, extreme aspect ratio, new and
heterogeneous materials, low damage, selectivity, across-wafer uniformity, chamber repeatability,
and lower environmental cost. These requirements interact: chamber state determines incident
fluxes, surface state determines reaction probabilities, and products/seasoning change the chamber.
The review identifies hybrid/coupled digital twins, improved diagnostics, computation, and data as
enablers; it does not claim that a single first-principles solver already closes this loop.

### What petch can legitimately attack

- A common physical input/output contract from reactor controls to source-plane distributions to
  feature profiles.
- Mechanism-based calibration against full profiles rather than per-geometry behavioral fitting.
- Surface state with memory, parameterized by atomistic calculations and experiment.
- Fast uncertainty/sensitivity propagation once gradients are verified.

### What petch must not claim yet

- Predictive reactor chemistry from recipe knobs.
- Broad-material first-principles chemistry.
- Atomic accuracy from a continuum level set alone.
- A digital twin before chamber state, wall state, and uncertainty are represented.

### Gate

One serialized physical case must be capable of carrying: geometry/materials; dimensional species
fluxes; joint energy-angle distributions; waveform/phase metadata; surface mechanism provenance;
and units. Analytic sources and future HPEM/PCMCM/experimental tables must replay through the same
interface.

## 2. OSTI 1802573: pattern-dependent distortion

Source: S. Huang, S. Shim, S. K. Nam, and M. J. Kushner, *Pattern dependent profile distortion
during plasma etching of high aspect ratio features in SiO2*, JVST A 38, 023001 (2020),
DOI 10.1116/1.5132800, OSTI 1802573.

### Actual problem

At aspect ratios above roughly 100, distortion is driven not only by reactor nonuniformity but by
feature-scale statistics and electrostatic interaction. Random differences in the energy, angle, and
sequence of arrivals produce feature-to-feature charging differences. In a symmetric pattern these
produce random tilt with zero ensemble-mean direction. In an asymmetric dense/sparse pattern,
charging produces a systematic horizontal field from the more-positive dense region toward the
less-positive sparse region, causing systematic tilt.

This is a distributional problem. A single smooth mean-charge trajectory cannot validate stochastic
twist, and a single isolated trench cannot validate inter-feature electrostatics.

### Current petch gap

- `charging_backward.py` solves a deterministic mean charging state in a small 2-D edge-array
  geometry.
- It has no ensemble of discrete arrival histories, no variance observable, and no evolving 3-D tilt.
- The finite edge-array excludes some outer faces to avoid an artificial open boundary; that is not
  yet a general dense/sparse pattern model.
- The current `E_defl` comparison does not validate pattern tilt; its AR1--4 correlation against the
  forward reference is presently poor.

### Attack sequence

1. **Mean-field sign gate:** symmetric array has zero mean lateral field; dense/sparse array has the
   paper's field direction. Do this on frozen geometry before profile evolution.
2. **Nonlocality gate:** increase lateral padding and neighbor count until the field/force observable
   is converged; an isolated-cell or boundary artifact must not set the result.
3. **Stochastic gate:** drive the same converged electrostatic solve with discrete, reproducible
   arrival histories. Symmetric ensemble mean tilt must approach zero while its variance remains
   nonzero; asymmetric mean tilt must remain nonzero.
4. **Profile gate:** only then couple lateral force to a 3-D evolving feature and compare tilt/twist
   distributions, not a hand-selected trajectory.

### Kill criteria

- Any systematic symmetric tilt that does not vanish with ensemble size.
- Dense/sparse field direction changes with domain padding or boundary convention.
- Reported twist without a convergence curve in cell size, history count, and pattern extent.

## 3. OSTI 2514378 / 2589032 / 3001885: atomistic data to feature chemistry

Sources:

- A. Kounis-Melas et al., *Deep potential molecular dynamics simulations of low-temperature
  plasma-surface interactions*, JVST A 43, 012603 (2025), OSTI 2514378.
- Associated open model/training/results dataset, DOI 10.34770/rjv6-2w31, OSTI 2589032.
- A. Kounis-Melas et al., *Deep-potential molecular-dynamics simulations of ion-enhanced etching of
  silicon by atomic chlorine*, JVST A 43, 063204 (2025), OSTI 3001885.

### Actual problem

Empirical interatomic potentials do not generalize reliably to the covalent, ionic, and metallic
interactions needed for new plasma-material systems. The papers test whether a DeepMD potential can
reproduce material properties, Ar+ sputter yield, amorphous/crystalline interface depth, spontaneous
Cl etching, ALE behavior, and ion-enhanced Si etch yields. The later work makes the feature-scale
dependencies especially useful: yield versus ion energy, neutral/ion flux ratio, and incidence angle;
steady Cl coverage near 1.25 monolayers; and SiClx mixed-layer thickness.

The open dataset includes the trained potential, training data, example simulations, and tabulated
results. The immediate petch opportunity is to consume the result tables with provenance and
interpolation error. Running DeepMD inside every feature step would be physically unnecessary and
far too slow.

### Current petch gap

- The continuous-etch path still centers on compact empirical yield laws.
- `ale.py` implements the later ALE reduced-order model, but the general evolving feature solver does
  not carry a unified per-surface coverage/mixed-layer/damage state.
- There is no versioned schema for energy-angle-flux-ratio yield and product-branching tables.
- Extrapolation outside a table's validated domain is not centrally prohibited or reported.

### Attack sequence

1. Acquire and checksum the OSTI 2589032 result tables without vendoring the large training corpus.
2. Define a units-explicit `SurfaceInteractionTable` containing material/species, energy, angle,
   flux-ratio/coverage axes, yield, product branching, damage/mixed-layer outputs, uncertainty, and
   source citation.
3. Reproduce the paper tables exactly at their nodes; test interpolation by leave-one-out withholding.
4. Refuse silent extrapolation by default. Permit an explicit policy that records the extrapolated
   fraction of surface impacts.
5. Compare the table path with the current compact law on identical frozen incident distributions.
6. Couple the table to stateful surface chemistry only after the static replay and units gates pass.

### Gates

- Tabulated node replay at floating-point tolerance.
- Interpolation error reported separately from MD/experimental error.
- Rotation/symmetry tests for angle convention.
- Conservation/accounting of incident particles and emitted SiClx products.
- Identical results when a table is loaded from disk or constructed in memory.

## 4. OSTI 2248044 / 2406001 / 1999799 / 2586627: time and depth memory

Sources:

- J. R. Vella et al., *Dynamics of plasma atomic layer etching: Molecular dynamics simulations and
  optical emission spectroscopy*, JVST A 41, 062603 (2023), OSTI 2248044.
- J. R. Vella et al., *A transient site balance model for atomic layer etching*, PSST 33, 075010
  (2024), OSTI 2406001.
- J. R. Vella and D. B. Graves, *Near-surface damage and mixing in Si-Cl2-Ar atomic layer etching
  processes*, JVST A 41, 042601 (2023), OSTI 1999799.
- J. R. Vella and D. B. Graves, *Si-Cl2-Ar+ Atomic Layer Etching Window*, JPCB (2025),
  OSTI 2586627.

### Actual problem

ALE is intrinsically transient. Alternating Cl2 and Ar+ exposures change the surface in time and in
depth. A top monolayer exchanges material with a subsurface mixed layer; ion energy sets damage and
mixing depth; product distributions evolve over the bombardment step. The 2023 experiment/MD paper
also identifies an unresolved discrepancy: experiment implies more Cl mixed into the layer than the
MD procedure predicts. The transient site-balance and ALE-window models are reduced descriptions of
this state, not evidence that memory can be discarded.

### Current petch state and gap

- `ale.py` and `ale_diff.py` already reproduce the published ALE-window reduced model and gradients.
- That is a zero-dimensional process calculation. It is not yet a field of surface states attached to
  an evolving geometry.
- The main feature solver does not advect/remap coverage, subsurface Cl, damage depth, or cumulative
  fluence as interfaces move.
- The published MD/experiment discrepancy must become uncertainty/model discrepancy, not a fitted
  constant hidden inside a yield.

### Attack sequence

1. Add state-conservation tests to the existing zero-dimensional ALE implementation: chlorine/site
   bounds, dose saturation, zero-flux invariance, and limiting continuous-exposure behavior.
2. Define the minimal surface state: top-layer Cl, mixed-layer Cl, mixed-layer/damage depth,
   cumulative ion fluence, and optional product inventory.
3. Advance that state on a static collection of surface elements under prescribed flux histories and
   reproduce the zero-dimensional model when all elements receive identical flux.
4. Specify remapping invariants before coupling to level-set motion: bounded coverage, conserved
   extensive inventory where applicable, and grid-refinement convergence.
5. Couple to one 2-D evolving ALE trench. Only after conservation and resolution gates pass should
   the state be moved into the GPU/3-D path.

### Gates

- Published ALE energy window and etch-per-cycle anchors remain passing.
- Pulse subdivision invariance: splitting an unchanged dose into smaller steps converges to the same
  state and removal.
- Static-uniform surface reproduces the zero-dimensional solver.
- Moving-interface inventory error decreases under refinement.
- Model-discrepancy output exposes the MD-versus-experiment Cl-mixing gap.

## 5. Reactor scale: preserve the coupling boundary now

Primary precedents are MCFPM/HPEM/PCMCM and Baer et al. (SISPAD 2010). Reactor/equipment simulation
provides species fluxes and ion angular/energy characteristics at the wafer; the feature model
computes local surface response. DOE's 2023 report adds the reverse coupling: effective sticking,
product emission, wall condition, and surface composition feed back to reactor chemistry through
hybrid sub-cycling.

The petch boundary should therefore have two objects:

1. `PlasmaBoundaryState`: position/time/phase-dependent species flux and joint velocity/energy-angle
   distribution, with sheath/source-plane metadata and units.
2. `SurfaceFeedbackState`: area-averaged consumption probabilities, emitted-product fluxes and
   distributions, charge/current response, material/coverage fractions, and uncertainty.

Analytic sources are constructors of `PlasmaBoundaryState`, not special branches deep in transport.
Future reactor solvers, reactor surrogates, diagnostics, and proprietary data should all produce the
same object. The feature engine remains independently testable.

## 6. Ranked local program

1. **Charging invariants and convergence** -- current work. Establish analytic boundary gates,
   forward/backward reciprocity, current residuals, and pattern-domain convergence.
2. **Pattern mean-field sign experiment** -- smallest direct attack on OSTI 1802573; frozen geometry,
   symmetric versus dense/sparse.
3. **OSTI surface-table schema and replay** -- start with the released result tables from 2589032;
   no feature coupling until interpolation/provenance gates pass.
4. **Spatial surface-state prototype** -- lift the already-gated ALE model onto static surface
   elements, then an evolving 2-D interface.
5. **Stochastic charging ensemble** -- distributional twist/tilt after the deterministic mean field is
   converged and boundary-independent.
6. **Reactor boundary objects** -- refactor sources behind the common physical contract; later ingest
   HPEM/PCMCM tables or a reactor surrogate without changing feature physics.

This ordering makes each layer independently falsifiable and keeps speed work honest: optimize only
after the physical observable and its convergence gate are defined.

## 7. Commercial wedge: an AI-native etch laboratory

The commercial target is not feature parity with an incumbent process simulator. Closed tools are
strongest when a customer already has a large, structure-specific calibration library. The opening is
where data are sparse, materials/processes are new, charging and history matter, and the model must
transfer beyond the exact structure used for calibration.

### Product claim to earn

**Turn sparse proprietary plasma-surface data and metrology into a transferable physical process
model, with gradients and uncertainty, then predict a held-out geometry or condition.**

This is stronger and more testable than “AI etch simulator.” AI should serve four physical jobs:

1. **Atomistic compression:** ML interatomic potentials and reduced models turn expensive MD/QC into
   energy-angle-state interaction tables.
2. **Physics calibration:** differentiable simulation infers uncertain mechanism parameters from full
   contours and time-series diagnostics while preserving units and constraints.
3. **Experimental design:** sensitivity and uncertainty identify the next experiment that most reduces
   predictive uncertainty, rather than merely adding another SEM to a library.
4. **Reactor acceleration:** a surrogate may map controls to `PlasmaBoundaryState`, but must be trained
   or corrected against high-fidelity reactor physics and carry out-of-distribution diagnostics.

AI must not silently replace charge conservation, transport, electrostatics, surface balances, or
interface kinematics. A learned component needs a physical domain, conservation/bounds tests,
provenance, uncertainty, and a fallback or refusal policy outside its domain.

### First credible buyer benchmark

Use a new or cold-start chemistry/material system with deliberately sparse data:

- fit physical surface/charging parameters on one subset of structures, aspect ratios, or process
  conditions;
- predict held-out structures/conditions without per-structure refitting;
- compare against a behavioral per-structure fit and a derivative-free calibration baseline;
- report profile error, calibration cost, uncertainty calibration, and wall-clock time;
- demonstrate that adding an MD table or one new experiment improves multiple held-out cases.

The benchmark should expose the data-efficiency advantage directly. Matching a heavily calibrated
closed tool on its home dataset would not demonstrate the wedge.

### Data architecture for proprietary chemistry

Keep private data out of solver code. Each mechanism package should contain versioned raw-data hashes,
units, processing scripts, fitted parameters or tables, covariance/model discrepancy, allowed domain,
and citations/ownership metadata. Public physics tests must run without private data; private validation
adds gates rather than altering algorithms. This makes collaboration and eventual deployment possible
without leaking customer or laboratory data.

### Commercial kill criteria

- Calibration improves the fitted structure but not a held-out structure.
- A learned component cannot identify when inputs leave its validated domain.
- A result depends on geometry-indexed corrections that cannot be expressed as physical mechanisms.
- Gradient calibration does not beat a fair derivative-free baseline in data or compute for the chosen
  high-dimensional contour problem.
- The model cannot explain which data or mechanism changed a prediction.

Meeting these criteria is how petch becomes an AI-native counterpart to the leading equipment and
simulation companies: not by imitating their accumulated calibration decks, but by making new physical
knowledge compound faster and transfer farther.
