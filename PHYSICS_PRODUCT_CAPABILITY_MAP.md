# Physics-to-product capability map

Audit date: 2026-07-15. This document reconciles the current two solver paths and fixes the order in
which they may be unified. It is a decision record, not a validation claim.

## Product target

The first product is a SiO2-first feature-scale forward solver. It accepts geometry, material fields,
species-resolved dimensional plasma boundary distributions, and sourced surface-interaction mechanisms.
It returns the evolving 3-D profile, surface state, species-resolved incident observables, charge/field
when enabled, uncertainty and convergence diagnostics, and an explicit validity/refusal result.

The first commercial proof is **calibrate declared uncertain physical inputs on structure N, then predict
held-out structure N+1**. Low error on the calibration structure alone is insufficient. Target-profile to
recipe optimization and a predictive reactor model are later layers.

## Reconciliation: there are currently two engines

1. `src/petch/api.py` exposes `Domain`, `SF6O2`, `Process`, and `Result`, but `Process.run()` calls the
   legacy `src/petch/threed.py::run_etch_3d`. That path contains a monolithic flag-selected collection of
   transport, chemistry, reflection, charging, redeposition, masking, and interface-update models.
2. `src/petch/feature_step_3d.py` is the intended common engine. It has SI boundary states, additive
   material level sets, material-local state, conservative remap, self-consistent dielectric charging,
   explicit validity, and refusal of unsupported topology events. It is not yet the public product path
   and does not yet contain every physical mechanism used by the legacy demonstrations.

Therefore, green legacy tests do not establish the common product, and immediately redirecting the public
API would remove capabilities. Migration must happen mechanism by mechanism through one physical
contract, with the old and new scores reported separately until each capability is re-earned.

## Capability matrix

| Product capability | Governing physics / inputs | Common-engine state | Legacy-only state or defect | External gate needed | Decision |
|---|---|---|---|---|---|
| Plasma boundary | Species flux and joint position-energy-angle-phase density at a physical reference plane | `PlasmaBoundaryState` is dimensional and source-agnostic | Legacy source laws are embedded in flags/parameters | Diagnostic/PIC distribution replay and flux/current invariants | Keep common contract; migrate sources, never benchmark-shaped source weights |
| Ballistic transport | Collisionless characteristics, visibility, conservation, electrostatic work for charged species | Arbitrary-triangle first-hit/field transport and deterministic gather exist | MC and fixed angular paths coexist in the monolith | Open wafer, trench and hole; AR/grid/ray/domain ladders | Keep both forward and deterministic estimators behind one error-controlled operator |
| Neutral re-emission | Surface integral equation with material/species reaction and escape probabilities | Conservative diffuse radiosity closes source = reacted + escaped | Legacy Knudsen/MC closures include AR-shaped parameters | Step coverage/ARDE plus species balance | Promote radiosity; do not port AR-specific attenuation laws as governing physics |
| Ion reflection | Material-, energy-, angle-, roughness-dependent scattering kernel | A material-tagged grazing-specular sensitivity now closes particle/charge/energy ledgers, preserves full lineage through charged re-impact, and contributes every landed event to chemistry in both ordinary and charging/profile modes. Field-free ions reuse the certified zero-field tracer; periodic moving profiles preserve the same hard-visibility topology across wrapped seams. Its three parameters are literature-bounded, not calibrated | ViennaPS-like coned-cosine law and stochastic energy loss | Differential reflection distributions and profile ablation | Keep the bounded common-engine sensitivity; require sourced material tables and held-out morphology before a predictive claim; no universal cone law |
| Product redeposition | Species-resolved production, flight, reaction/sticking, and exact material balance | Surface mechanisms return an exact material-origin ledger; named populations close that ledger; the conservative surface-emission operator closes first escape, repeated impact/re-emission, reaction, and final escape. An opt-in, provenance-bearing v1 converts bounded sticking into signed same-material interface growth and reports emitted/deposited/escaped material exactly. Cross-material capture, pinned-surface growth, missing launch laws, and undeclared target materials refuse. The pinned Si DeepMD yield still refuses transport because its source has no differential emitted energy-angle law. Reactive SiO2 branching remains unresolved | Scalar velocity-derived products, fixed sticking, and a charging/etch-stop foot-band suppression hack | Time-resolved necking/taper and trench-to-hole transfer | Calibrate the explicit v1 law on time-resolved morphology; implement a new material layer before allowing cross-material films; never port the legacy foot-band branch |
| SiO2 fluorocarbon surface state | Conserved sites, fluorination/complexes, polymer inventory, ion-assisted removal, product branching | Reduced dimensional state exists; current parameter set is explicitly nonpredictive | Belen/ViennaPS-shaped calibrated closure produces old SF6/O2 results | Time-resolved contours and held-out process/geometry transfer | Extend state only when a sourced observable requires it; expose unknown rates as calibrated inputs |
| Si SF6/O2 surface state | Coupled F/O site balance, direct F chemical removal, physical sputtering, ion-enhanced removal, neutral loss and target-material conservation | `BelenSiliconSF6O2Mechanism` runs the published pseudo-steady law through the common event/radiosity interface. The same available-site fraction controls neutral re-emission and chemistry; the engine converges that coupling before profile motion. Every parameter has a source and bound, removed Si closes its ledger, and volatile SiFx routing stays explicitly unresolved. Source-correct reflection improved the exposed Figure-9 development RMSE from 3.545 to 2.557 micrometres; a three-level moving-profile gate is conservative and sample-contracting but is not experimental validation | The old generic-SiO2 adapter and legacy Belen monolith remain replay-only | Independent profile/flux-ratio/temperature experiment after the exposed de Boer points; grid/sample/fixed-point refinement; product branching | Promote the common mechanism as the sole development path; do not call it predictive until a new held-out experiment passes and do not add photons without a measured flux/yield |
| Mask evolution | Separate material state, sputter/oxidation/deposition, selectivity, facet motion | `MaterialMechanismRouter3D` now dispatches exposed face IDs to independent mechanisms and state fields, merges signed material/product ledgers, moves mask and substrate level sets with their own velocities, recomputes ownership after motion, and round-trips namespaced state through checkpoints. Manufactured two-material sputter gates pass. A production mask law is still missing | Static mask pinning and region classification | Mask-facet time evolution and high-O2 transfer | Supply material-specific mask evidence through the router; never reuse substrate chemistry silently |
| Charging/notching | Poisson with permittivity/conductors, deposited current, dielectric storage/conductivity, self-consistent current balance | Physical-time charge evolution, full charged surface-response/re-impact, conservative moving-surface charge remap, and C3 charge/profile co-evolution now run through one engine. Signed `CCA-2026-07-13-R2` remains the default. An explicit draft profile-relevant path compares two fresh physical-time blocks using disjoint endpoint ensembles and gates field/current/delivered-flux/profile-increment changes while retaining every R2 diagnostic; it does not authorize an experimental claim | Low-AR HG/notch lineages are separate and include reference-emulation closures | Absolute matched profile; grid/sample/init invariance; charging off/on | Review/sign or reject the separately versioned profile-stationarity contract, then run the held-out C4 campaign; never reopen frozen-map root solvers |
| Bowing/tilt/twist | Coupled transport, reflection/redeposition, pattern electrostatics, stochastic finite-count arrivals | Dimensional Poisson arrivals convert flux×area×time to compact integer counts inside the conservative charge update. Constant-boundary physical-time C3 and the public ensemble evolve independent 3-D profiles. Geometry-native centerline/displacement/onset metrics and an executable nested N/2N, paired sample-doubling, isotropy/systematic-direction scorer now exist. Surface-response branching remains conditional-mean unless a stochastic response law is supplied, and no N>=30 AR campaign has run | No common evolving 3-D validation gate | Symmetric zero-mean/nonzero-variance twist; dense/sparse systematic tilt; time-resolved bow/neck | Run the preregistered C5 ensemble/isotropy/refinement gates; never present one realization as a deterministic forecast |
| Interface evolution | Hamilton-Jacobi motion, additive material fields, conservative state transfer | CR-2 reinitialization, material-local motion/remap and topology refusal exist | Legacy first-order/FSM path pins masks and merges materials | Plane/sphere manufactured gates plus grid/time/profile convergence | Common engine is authoritative |
| Chemistry extensibility | Versioned material mechanism with provenance, uncertainty and refusal | Narrow Si-Cl2-Ar+ table replay proves architecture beyond SiO2 | Bosch/ALE/cryo demonstrations remain separate reduced solvers | Second chemistry through identical end-to-end product schema | Claim chemistry-extensible, not universal chemistry prediction |
| Reactor coupling | Reactor/sheath outputs define feature boundary; surface returns effective reaction/sticking/product fluxes | Inbound state exists; feedback contract does not | None coherent | Held-out wafer-position/process transfer | Add a `SurfaceFeedbackState`; do not claim a multiscale digital twin yet |
| UQ and refusal | Numerical convergence plus parameter, measurement and model-form uncertainty | Step validity exists | Public API omits it | Coverage/refusal on preregistered held-out cases | Validity is mandatory product output, not documentation |
| Differentiable calibration | Sensitivities through transport, state, converged fixed point and interface/profile loss | Only reduced ALE gradient is verified | None end-to-end | Finite-difference agreement and structure-N to N+1 calibration transfer | Do not sell feature-scale backpropagation until boundary terms are handled |
| Public product API | Versioned config/result/checkpoint schemas using the common engine | `PhysicalProcess`, `PhysicalChargingProcess`, and the finite-count ensemble call the common engine explicitly. Uncharged and charged run manifests now bind schema, geometry arrays, material IDs/layers, boundary distributions, mechanism/redeposition provenance, operator controls, device, and numerical settings. C3 itemizes inline replay/horizon/bounce recovery, priced tail absorption, and conservation bounds. Continuation preserves geometry, remapped charge, state and fingerprint; the versioned `allow_pickle=False` checkpoint binds its source manifest by SHA-256 | Friendly `Process` API calls legacy engine and labels that provenance | Fresh-install replay, accuracy-matched CPU/CUDA manifest, and checkpoint migration/version-refusal gates | Keep the explicit APIs alongside the legacy facade; do not silently redirect `Process.run()` |

## External bar set by recent work

- The 2023 DOE semiconductor-plasma report says feature-scale reaction probabilities are still commonly
  selected empirically and calibrated to SEMs, limiting quantitative range. It recommends a hybrid
  reactor/feature exchange through incident fluxes and effective surface feedback rather than one mesh
  spanning sub-nanometres to chamber scale.
- Krueger et al. (JVST A, 2024, DOI `10.1116/6.0003554`) already optimize physical MCFPM reaction
  parameters against a final SEM and report process-trend transfer. Their partial transfer and roughly
  20 CPU-hour single-feature cost define the baseline: petch must demonstrate more data-efficient,
  faster held-out transfer with honest model-form omissions, not merely fit one SEM.
- The 2026 K-SPEED sputtering study (DOI `10.1063/5.0304157`) calibrates against time-resolved trench
  profiles and transfers to a hole, resolving mask facets, entrance necking, taper, and geometry-dependent
  redeposition. This is the stronger near-term validation form.
- Huang et al. (JVST B, 2020, DOI `10.1116/1.5132800`) separates stochastic, zero-ensemble-mean twist in
  symmetric patterns from systematic dense/sparse electrostatic tilt. A single smooth trench cannot
  validate either claim.

## Implementation order and stop rules

1. **Freeze the claim boundary.** Keep legacy demonstrations runnable but label them original-path.
   Remove no legacy capability until the corresponding common-engine evidence exists.
2. **Create the common product schema and adapter.** Geometry, materials, boundary state, mechanisms,
   numerical targets, validity and diagnostics must serialize. Initially expose it alongside the legacy
   API; do not redirect `Process.run()` silently.
3. **Close the mass/material loop.** The generic router and same-material return path are implemented.
   The remaining evidence work is a production mask law, reactive product branching, and an explicit new
   material layer for cross-material coatings. Retire the legacy foot-band suppression rather than
   translating it.
4. **Add sourced energetic reflection.** Use material/energy/angle/roughness tables with interpolation
   audit, uncertainty and refusal outside support.
5. **Run one preregistered time-resolved transfer campaign.** Calibrate only declared unknown physical
   inputs on early trench profiles; predict later times and a held-out geometry/process. Report numerical,
   digitization, parameter and model-form error separately.
6. **Add charging only where its experimental ablation is resolvable.** Require converged current balance
   and charging-off/on profile causality. Do not spend an open-ended campaign matching an internal residual
   that does not change the target profile above its error budget.
7. **Then verify sensitivities and calibration speed.** Gradients are promoted only after boundary-motion
   terms and the converged fixed point pass independent finite-difference gates.

Every investigation gets one external observable, one causal hypothesis, a preregistered numerical error
budget, and a compute ceiling. Non-monotonic progress on the gating observable triggers a method review,
not more samples. No aspect-ratio, benchmark-name, expected-answer, or surface-region branch may enter the
common governing physics.
