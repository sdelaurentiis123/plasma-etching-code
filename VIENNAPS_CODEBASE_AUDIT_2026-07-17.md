# ViennaPS codebase audit and petch superset boundary

Date: 2026-07-17
Status: read-only source audit; implementation proposals remain gated behind the active Krüger verdict

Audited upstream: ViennaPS 4.6.1, commit
`2956ed587984c6dc38be24c6e2390e10c9b2f0a7` (2026-07-02).

Related documents:

- [`COMPETITIVE_DIFFERENTIATION_ROADMAP_2026-07-17.md`](COMPETITIVE_DIFFERENTIATION_ROADMAP_2026-07-17.md)
- [`NUMERICS_CALIBRATION_AMR_PLAN_2026-07-17.md`](NUMERICS_CALIBRATION_AMR_PLAN_2026-07-17.md)
- [`PHYSICS_FIRST_UNIFIED_ENGINE_2026-07-17.md`](PHYSICS_FIRST_UNIFIED_ENGINE_2026-07-17.md)

## 1. Executive verdict

ViennaPS is presently the stronger production topography framework. It has mature C++ and Python
interfaces, sparse level-set storage, several CPU/GPU ray engines, broad process examples, geometry
utilities, serialization, and cross-platform CI.

Petch is already the more ambitious feature-physics architecture: dimensional wafer distributions,
electrostatic charging, charged trajectories, certified response lineage, conservative surface and
product inventories, stochastic twisting, and strict calibration/validation/provenance contracts.
Those are verified operators, not yet evidence of superior untouched wafer prediction.

The objective is a behavioral strict superset, not an indiscriminate source merge:

1. run official ViennaPS unchanged as an independent reference;
2. reproduce overlapping scientific cases through petch's common engine;
3. borrow clean architectural ideas and independently implement equations from primary sources;
4. preserve petch's conservation, uncertainty, charging, and claim contracts;
5. earn superiority through held-out experiments and cost/error measurements.

## 2. Capability catalog

| Area | ViennaPS 4.6.1 | Petch implication |
| --- | --- | --- |
| Core | Header-only C++20 with Python bindings; 2-D/3-D level-set topography and MC flux | Match routine usability and robustness |
| Process loop | Remesh, trace flux, solve coverages, compute velocity, advect, remap | Use as an independent architecture reference |
| Geometry | Nested material level sets, Boolean operations, masks, top layers, HRLE sparse narrow bands | Vienna is currently stronger in geometry ergonomics |
| Time integration | Engquist–Osher default, Euler/RK choices, CFL and adaptive controls | Add a parity/convergence matrix to petch |
| Flux backends | CPU disk, CPU triangle/Embree, GPU disk, GPU triangle/OptiX, GPU line in 2-D | A common petch backend factory is high value |
| Etch models | SF6/O2, HBr/O2, SF6/C4F8, CF4/O2, detailed fluorocarbon, IBE, Faraday cage, generic particles | Crosswalk only validation-demanded mechanisms |
| Deposition | TEOS/PECVD, particle/polymer deposition, ALD | Useful conformance and future breadth targets |
| Cyclic process | Bosch composition and generic pulse/purge atomic-layer strategy | Borrow the process-composition pattern; add evidence-backed ALE separately |
| Other physics | Wet etch, selective epitaxy, neutral transport, oxide regrowth, oxidation with diffusion/flow/mask bending | Vienna is substantially broader outside plasma feature etch |
| Testing | 38 test directories, 41 C++ tests, 4 Python tests, 29 example families, multi-platform CI | Petch needs a similarly routine release league |
| Charging | No built-in feature Poisson charging, charge evolution, SEE, or field-deflected trajectories found | Major petch differentiator if experimentally validated |
| AMR | One fixed `gridDelta` on sparse HRLE storage; no feature-adaptive refinement found | Sparse narrow band and AMR are separate builds |

Primary source locations include:

- process loop: `include/viennaps/process/psFluxProcessStrategy.hpp`;
- domain/material level sets: `include/viennaps/psDomain.hpp`;
- profile integration controls: `include/viennaps/process/psProcessParams.hpp`;
- fluorocarbon deposition/removal: `include/viennaps/models/psFluorocarbonEtching.hpp`;
- energetic reflection/redeposition: `include/viennaps/models/psIonBeamEtching.hpp`;
- examples: `examples/stackEtching/` and `examples/boschProcess/`.

The immutable source is available at
`https://github.com/ViennaTools/ViennaPS/tree/2956ed587984c6dc38be24c6e2390e10c9b2f0a7`.

### 2.1 Source-backed WP1 geometry boundary

The WP1 distinction between sparse bands and true AMR is established directly by upstream source:

- `psDomain.hpp:169-199` inserts nested material level sets into one domain, and `184-188` requires
  every inserted level set to have the domain's single `gridDelta`;
- `psDomain.hpp:420-423` exposes the underlying HRLE grid and that one global spacing;
- `psFluxProcessStrategy.hpp:302-363` expands the narrow band, regenerates the surface used by the
  flux engine, traces, computes velocity, advects, and maps coverage back after motion;
- `psAdvectionHandler.hpp:131-187` moves coverage through the level-set representation during that
  advection;
- `psProcessParams.hpp:52-72` provides spatial/temporal schemes, CFL control, and optional adaptive
  **time** stepping, but no refine/coarsen hierarchy.

Therefore ViennaPS provides sparse HRLE narrow-band evolution on a uniform Cartesian spacing. That is
a valuable WP1 reference and likely the correct first petch optimization, but it is not feature-adaptive
mesh refinement. No octree, block hierarchy, 2:1 balance, coarse/fine flux, or conservative
prolongation/restriction path was found in the audited public source.

The actionable WP1 sequence is consequently:

1. preserve the existing petch state, remap, topology, and conservation contracts behind a geometry
   backend interface;
2. add a uniform-`dx` sparse narrow-band backend and require identical short-burst answers to the dense
   reference at matched `dx`;
3. make surface extraction and transport-cache invalidation explicit when the active band changes;
4. only then add block-structured AMR with the R3 refine/coarsen transfer contract from the numerical
   plan;
5. compare sparse uniform and AMR paths at matched error, not merely at matched nominal cell size.

This is a clean-room architectural lesson. Current ViennaHRLE/ViennaLS implementation source is GPL-3.0
and is not a copy source for a differently licensed core.

## 3. Important corrections to prior petch shorthand

ViennaPS can shrink and close an opening. Its fluorocarbon model explicitly grows polymer once its
coverage threshold is reached, and examples represent the polymer as a separate top material level
set. Therefore “Vienna cannot neck” is false.

The defensible distinction is narrower: petch can make necking depend on a richer dimensional
mask/polymer state, retain material/product inventories, couple charging, and report conservation and
validity. Those advantages matter only after untouched profile prediction is demonstrated.

Another stale claim is that current ViennaPS retraces the source flux four to eight times per profile
step. The audited implementation traces it once during an ordinary step; iterative coverage
initialization occurs at startup. Old projected 4--6x exact-mode speedups based on repeated tracing
must not be repeated.

## 4. Architecture worth learning

High-value patterns to reproduce independently around petch's authoritative state path are:

- a `ProcessModel` composed from surface model, velocity field, optional geometry model, and callbacks;
- a common flux-engine interface with automatic CPU/GPU backend selection;
- a compact `Process.apply()`-style Python workflow;
- typed parameter structures with complete metadata serialization;
- explicit surface-state remapping through level-set motion;
- a generic atomic-layer pulse/purge driver;
- explicitly different fast geometric and physical ray-tracing modes;
- CPU/GPU parity tests for plane, trench, reflection, deposition, and material cases;
- GDS import, geometry factories, and restart/serialization support;
- sparse narrow-band evolution before general AMR.

These patterns must converge on one petch process context. They must not create another legacy engine
beside the common boundary/transport/state/profile path.

## 5. Chemistry crosswalk

Vienna's reduced models are useful references, not automatically authoritative parameter sources:

- **SF6/O2 and HBr/O2:** algebraic etchant/passivation coverages with chemical, ion-enhanced, and
  sputter terms;
- **fluorocarbon:** etchant/polymer coverages, explicit polymer deposition/removal, and
  material-dependent densities/yields;
- **CF4/O2:** F/O/C state variables, with visibly incomplete sections that should not be copied as
  settled chemistry;
- **ion-beam etch:** energy/angle-dependent reflection and redeposition weighting;
- **Bosch:** composition of etch and deposition steps rather than one universal chemistry;
- **atomic-layer process:** sound pulse/purge infrastructure, but its supplied physical model is ALD,
  not a ready-validated ALE chemistry.

Petch should independently source the underlying equations and promote only mechanisms required by
experimental evidence. Every promoted law retains dimensional units, material/product ledgers,
parameter provenance, bounds, and validity.

### 5.1 Public model/backend matrix

GPU availability is model-specific, not a property of every Vienna process:

| Public model family | Main declared state/channel | Audited ray backend status |
| --- | --- | --- |
| SF6/O2 and HBr/O2 | etchant/passivation coverage; ion, neutral, sputter response | CPU and GPU implementations |
| SF6/C4F8 | etch plus separate polymer deposition/removal | CPU and GPU implementations |
| Detailed fluorocarbon | etchant/polymer/polymer-etch coverage; material-specific deposition and etch | CPU implementation found |
| CF4/O2 | F/O/C coverages and energetic-ion response | CPU implementation found |
| Ion-beam etch | energetic reflection plus redeposition weight | CPU and GPU implementations |
| Generic single/multi-particle | configurable neutral/ion flux-to-rate callback | CPU and constrained GPU variants |
| Single-particle ALD | pulse/purge coverage, desorption and optional diffusion | CPU and ballistic GPU variants |
| TEOS PECVD | neutral/ion deposition response | CPU and GPU implementations |

The detailed fluorocarbon model is especially important for necking comparisons: at
`psFluorocarbonEtching.hpp:213-246` it switches among polymer deposition, mask sputtering, polymer
removal, and substrate etch according to material and coverage. Petch conformance must compare that
declared reduced law separately from any richer mask/polymer mechanism; different equations are not a
parity failure.

## 6. Performance evidence

The existing `viennaps-accel` work measured a controlled 17.6x acceleration of the level-set phase
and about 1.58x projected end-to-end speedup from that phase alone. This is valuable but does not show
that Vienna-level accuracy is free:

- current ViennaPS already performs one ordinary source trace per step;
- exact-preserving acceleration is likely nearer 2x unless ray transport also changes;
- the local coarse petch/de Boer match was not refinement-stable;
- `dx=0.25` matched by coincidence while refinement approached the Vienna ballistic result.

The source record is `viennaps-accel/notes/dx-convergence-verdict.md`.

Additional source constraints keep the runtime claim honest:

- upstream documents GPU ray tracing and oxidation diffusion as experimental
  (`README.md:158-160`);
- GPU auto-selection occurs only when the chosen model supplies a GPU implementation
  (`psProcess.hpp:303-317`); otherwise the ordinary automatic path selects CPU disks;
- GPU triangle/line/disk adapters convert geometry to single-precision buffers, so backend parity must
  include grazing/shared-edge and material-ID cases rather than only smooth planes;
- level-set advection remains outside the ordinary ViennaPS GPU ray-tracing path;
- the 17.6x result is a controlled phase microbenchmark of a separate dense GPU
  advect/redistance pipeline, not an integrated upstream end-to-end result.

The local profiling suggests an exact-preserving integrated speedup on the order of 2x, but that is a
planning estimate, not an earned benchmark. WP6 must record matched-error wall time separately for each
model/backend pair; a GPU SF6/O2 timing cannot stand in for the CPU-only detailed fluorocarbon case.

## 7. Strict-superset conformance suite

“Strict superset” is earned only when the following table is executable and versioned:

| Case | Shared truth | Required petch result | Vienna role |
| --- | --- | --- | --- |
| Translating plane | analytic displacement and inventory | convergence plus exact ledger | independent baseline |
| Straight trench conductance | analytic/asymptotic limits | flux/order/grid convergence | compare ray backends |
| Sticking/re-emission sweep | probability conservation | source = reaction + escape | compare surface transport |
| Energetic reflection | number/energy/angular ledger | reciprocity and refinement | compare IBE channel |
| Polymer deposition | mass and opening evolution | material-resolved conservative growth | compare fluorocarbon model |
| Multi-material advection | manufactured interface motion | no owner loss or invented material | compare geometry robustness |
| Bosch composition | declared cycle sequence | same step semantics and state continuity | compare process composition |
| CPU/GPU parity | identical declared operator | tolerance/bitwise gates as appropriate | compare backend maturity |
| Charging trench | petch physical contract | refined Poisson/trajectory/charge result | expected petch-only capability |
| Stochastic twist | ensemble symmetry/statistics | stable distribution, no fake deterministic claim | expected petch-only capability |

Matching Vienna is not sufficient when an analytic or manufactured truth exists. Both codes are
scored against truth first and against one another second.

### 7.1 WP6 execution matrix and gates

Each overlapping conformance case has three deliberately separate scores:

1. **truth score:** analytic or manufactured error, conservation, and refinement order;
2. **Vienna parity score:** difference between version-pinned ViennaPS and petch when the governing law,
   geometry, source distribution, seed mode, and observable are actually shared;
3. **product score:** wall time, peak memory, unattended completion, recovery/refusal counts, and output
   provenance at matched truth error.

The minimum backend matrix is:

| Axis | Required variants |
| --- | --- |
| Dimension | 2-D and one genuinely 3-D geometry |
| Vienna surface representation | CPU disk and CPU triangle; GPU disk/triangle only for models that implement them |
| Sampling | fixed seed replay and independent-seed ensemble statistics |
| Profile numerics | at least three `gridDelta` values; timestep/CFL refinement |
| Material complexity | one material, nested mask/substrate, deposited top layer |
| Motion sign | removal, deposition, and one alternating cycle |

A parity row is invalid rather than failed when the equations differ. Examples include comparing
Petch's charged trajectories to Vienna's uncharged rays, or its inventory-resolved mask chemistry to a
Vienna algebraic coverage law. Those belong in causal or experimental comparisons, not in a numerical
parity metric.

WP6 should begin with the smallest discriminating set: translating plane; straight-trench neutral
conductance at multiple sticking coefficients; energetic reflection with number/energy scoring;
multi-material removal; explicit polymer deposition and reopening; one Bosch cycle; then CPU/GPU
parity. Charging and twisting are petch-only extension rows after overlap closes.

### 7.2 Validation claim boundary

The audited Vienna tree has substantial software verification—38 test directories, 41 C++ test files,
4 Python tests, 29 example families, and cross-platform CTest CI—but no repository-level suite of
digitized, preregistered held-out wafer profiles was found. Its README also warns that the SF6/O2 model
may change and that newer versions need not reproduce the displayed v3.6.0 image
(`README.md:193-201`).

Accordingly:

- ViennaPS is a mature numerical/process-framework baseline;
- its example output is not experimental ground truth;
- agreement with Vienna earns overlap/conformance evidence, not wafer validation;
- disagreement must first be decomposed into equation, sampling, mesh, and implementation differences;
- only untouched experimental scoring can earn the predictive-superiority claim for either engine.

## 8. Independent adapter boundary

The Vienna comparison backend should remain a separate executable/process:

```text
declared geometry + boundary + model parameters
                    |
          +---------+---------+
          |                   |
       petch run          ViennaPS run
          |                   |
          +---------+---------+
                    |
        common geometry/observable scorer
```

Every comparison records Vienna version/SHA, compiler, backend, grid, seed, parameter file, input
checksums, output geometry, wall time, and hardware. No benchmark-specific correction enters either
solver.

## 9. Licensing boundary

Current ViennaPS is GPL-3.0; versions before 4.3 were MIT. A deliberate product/license decision is
required before linking or copying current GPL implementation into a differently licensed core.

Safe default:

- run official ViennaPS unchanged in a separate environment/process;
- exchange declared geometry/configuration and VTP/CSV observables;
- independently implement equations from cited scientific papers;
- preserve notices and exact provenance;
- use old MIT history only after pinning and verifying the specific historical source;
- obtain legal review before distribution if direct code reuse is desired.

## 10. Ordered, bounded follow-up

After the active Krüger verdict:

1. build the black-box Vienna benchmark adapter;
2. land the manufactured conformance suite;
3. add a common petch process context/backend factory;
4. implement periodic, conservative geometric remapping and its refinement suite;
5. benchmark sparse narrow-band evolution;
6. add true AMR only after the uniform-grid and remap contracts pass;
7. crosswalk fluorocarbon first, then one SF6/O2 or HBr/O2 mechanism demanded by validation;
8. publish the accuracy/runtime/reliability matrix without claiming superiority prematurely.

## 11. Claim boundary

Today the accurate statement is:

> ViennaPS has the more mature production topography kernel. Petch has a broader evidence-aware
> feature-physics architecture. Petch becomes a strict, superior superset only when overlapping
> capabilities are at least as reliable and fast, and its additional physics improves untouched
> experimental predictions.
