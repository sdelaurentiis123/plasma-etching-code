# AMR and sparse narrow-band implementation audit

Date: 2026-07-17
Status: implementation design only; no simulation was run and no held-out Krüger outcome was read
Scope: the common 3-D feature engine in `src/petch`, not the legacy 2-D or benchmark-only solvers

Related documents:

- [`NUMERICS_CALIBRATION_AMR_PLAN_2026-07-17.md`](NUMERICS_CALIBRATION_AMR_PLAN_2026-07-17.md)
- [`KRUEGER_2024_R19_RESPONSE_CHECK_REPORT_2026-07-17.md`](KRUEGER_2024_R19_RESPONSE_CHECK_REPORT_2026-07-17.md)
- [`VIENNAPS_CODEBASE_AUDIT_2026-07-17.md`](VIENNAPS_CODEBASE_AUDIT_2026-07-17.md)
- [`VIENNAPS_TOPOLOGY_CONFORMANCE_2026-07-17.md`](VIENNAPS_TOPOLOGY_CONFORMANCE_2026-07-17.md)

## 1. Executive decision

Do **not** begin with a general octree, an adaptive Poisson solve, or a new surface evolution engine.
The minimum safe target is:

1. preserve the current uniform-grid engine as the reference operator;
2. put a narrow geometry-backend seam underneath that same engine and prove exact uniform parity;
3. implement fixed-resolution block-sparse narrow-band storage at one global `dx`;
4. add one block-structured 2:1 refinement level, initially `10 nm -> 5 nm` for the Krüger cell;
5. require every physical interface, material junction, and its full numerical stencil to remain on the
   finest level;
6. coarsen only far-field gas or single-material blocks in v1;
7. certify the result against uniform 5 nm manufactured and short-burst references before one long run.

This is true AMR in the volume, but intentionally conservative AMR at the surface. It avoids the
hardest failure mode—a coarse/fine boundary cutting the zero contour—until the block hierarchy,
ledgers, checkpointing, periodicity, and surface remap have all passed independently.

The first release will reduce level-set memory, velocity-extension work, redistancing work, and
far-field topology work. It will **not** automatically reduce the number of surface triangles or the
particle/radiosity cost, because the physical surface remains fine. Selective coarsening of smooth
surface regions is a later, separately gated milestone.

## 2. Current engine architecture

The current authoritative geometry is a dense, nodal, uniform Cartesian representation:

```text
FeatureGeometry3D
  phi: dense float64[nx,ny,nz]
  material_levelsets[m]: dense float64[nx,ny,nz]
  material_id: dense int[nx,ny,nz]
  one scalar dx
        |
        v
skimage marching cubes (`threed.extract_mesh_3d`)
        |
        +--> triangle hard-visibility transport / radiosity / charging
        |
        +--> material-local surface chemistry and triangle state
        |
        v
triangle normal velocity
        |
nearest-triangle extension to a dense 4 dx band
        |
Godunov advection + CR2/FSM/skfmm redistance on dense arrays
        |
new marching-cubes surface
        |
material-local conservative KNN surface-state remap
```

The architecture is more separable than the size of `feature_step_3d.py` suggests:

| Concern | Current authority | AMR implication |
| --- | --- | --- |
| Material geometry | `FeatureGeometry3D` in `src/petch/feature_step_3d.py` | Dense arrays and scalar `dx` are the main seam to generalize |
| Uniform advection/redistance | `src/petch/threed.py` | Keep permanently as the reference implementation |
| Surface extraction | `extract_mesh_3d` in `src/petch/threed.py` | Add a backend-neutral surface object; do not change transport inputs |
| Face ownership/normals/topology | helpers in `src/petch/feature_step_3d.py` | Must be expressed through geometry sampling/topology methods |
| Surface state | triangle fields returned by mechanisms | Can remain unchanged if remap receives old/new physical surfaces |
| Surface remap | `conservative_remap_surface_state` in `feature_step_3d.py` | Must become shared geometry infrastructure before repeated AMR remeshing |
| Surface charge remap | `src/petch/surface_charge_remap_3d.py` | Reuse the same geometric transfer weights, with its distinct removal closure |
| Profile observables | `src/petch/profile_observables_3d.py` | Replace direct dense slicing with hierarchy line/slice sampling |
| Checkpoints | dense arrays in `charging_checkpoint_3d.py` and campaign scripts | Add a versioned sparse/AMR schema; retain v1 readers |
| Poisson charging | uniform Q1 in `src/petch/charging_poisson_3d.py` | Not AMR-compatible today; keep charging on uniform geometry until its own gate |

Direct dense-grid coupling is concentrated rather than ubiquitous. Within the package, most direct
geometry access is in `feature_step_3d.py`; the other principal consumers are
`profile_observables_3d.py`, `charging_coevolution_3d.py`, `charging_checkpoint_3d.py`,
`nozawa_replay_3d.py`, and `krueger_replay_3d.py`. Transport and chemistry already consume surfaces
and dimensional boundary states, so they should not be rewritten.

## 3. First-principles target representation

### 3.1 Authoritative state

Per-material signed-distance fields remain authoritative. The combined field and owner are derived:

```text
phi(x)      = max_m phi_m(x)
owner(x)    = argmax_m phi_m(x), when phi(x) >= 0
gas         = phi(x) < 0
```

The sparse hierarchy stores fixed-size blocks. A block is one of:

- `far_gas`;
- `far_material(material_id)`;
- `active`, carrying nodal signed distances for every relevant material plus halo data.

Any block containing a gas/solid crossing, a material/material junction, or an interface stencil is
`active`. Therefore a far block has one unambiguous label; topology is not inferred from a single
centroid sample.

### 3.2 Safe first AMR rule

Let `h_f` be the finest spacing. The present engine redistances through `4 h_f`, uses a one-cell
Godunov stencil, and allows an accepted outer move up to a declared fraction of one cell. The v1
fine-band half-width should therefore be at least:

```text
4 h_f redistance + 1 h_f advection/MC halo + accepted displacement + safety
```

Use a conservative initial declaration of **8 fine cells** and prove it by a width sweep. Refine
before an interface can enter the outer buffer. Coarsen only when every material interface is outside
a larger hysteresis radius (for example 12 fine cells). Correctness must not depend on the eventual
optimized values.

For v1:

- base spacing: 10 nm;
- one refined level: 5 nm;
- refinement ratio: exactly 2;
- block size: configurable; `8^3` cells is a performance candidate, not a physical contract;
- boundary blocks may be partial, because the 130 x 20 nm Krüger period is not divisible by an
  arbitrary block size;
- duplicate periodic endpoint planes remain one physical plane and must be canonicalized;
- a sign-changing or material-junction block may never be coarsened.

This exploits the already-proved grid fact that 10 and 5 nm preserve the Krüger domain while nominal
20 nm does not.

### 3.3 Deterministic refinement

The first refinement indicator is geometry-only:

- signed distance to every interface;
- material junction proximity;
- predicted maximum interface displacement before the next rebuild;
- periodic seam ownership.

Do not initially refine from a raw Monte Carlo flux, charge, coverage, or velocity gradient. A noisy
indicator would make the mesh sample-dependent and could turn sampling noise into a persistent
operator bias. Physics-based indicators may be added later only with deterministic filtering,
hysteresis, checksum-bound thresholds, and refinement-schedule stability tests.

## 4. Large work packages and exact subtasks

### WP-AMR0 — Freeze the uniform reference and promotion gates

**Purpose:** AMR needs a truth operator before it needs code.

Existing files:

- `scripts/krueger_2024_multiresolution_audit.py`
- `tests/test_krueger_2024_multiresolution_audit.py`
- `src/petch/threed.py`
- `src/petch/feature_step_3d.py`

Subtasks:

1. Bind every numerical-response or calibration artifact to the geometry/operator epoch.
2. Preserve a current-operator uniform 10 nm and uniform 5 nm initial `0.5 s` reference.
3. Preserve a late-state uniform 5 nm frozen operator and `0.1--0.5 s` short-burst reference.
4. Record stage wall time, peak resident/GPU memory, dense node count, surface triangle count, and
   time split among extraction, transport, chemistry, extension, advection, redistance, remap, and
   observables.
5. Keep all Krüger held-out outcomes outside this campaign; this is base-calibration numerical
   evidence only.

Promotion gate:

- both reference windows close their ordinary ledgers and topology rules;
- exact source/config/checkpoint hashes exist;
- the performance profile shows which cost AMR can actually remove.

Stop rule: no 60 s AMR run is permitted from this work package.

### WP-AMR1 — Introduce a geometry-backend seam with uniform parity

**Purpose:** isolate geometry storage/operations without touching transport or chemistry.

**2026-07-17 increment:** the dependency-neutral state and read-only uniform backend exist. Both
initial and post-advection surface extraction plus face ownership now route through one private
uniform-backend bridge that returns writable C-order copies of the legacy arrays. The existing
active-surface fingerprint is unchanged. Exact array/dtype/order/writeability coverage and the
feature/checkpoint/replay cluster pass `90` tests. Domain sizing, normals, topology, observables,
advection, and redistance are intentionally not routed yet.

New source:

- `src/petch/feature_geometry_backend_3d.py`
- `tests/test_feature_geometry_backend_3d.py`

Existing source to modify:

- `src/petch/feature_step_3d.py`
- `src/petch/profile_observables_3d.py`
- `src/petch/__init__.py`

The backend contract should expose:

- physical domain extent, origin, periodic axes, and finest local spacing;
- exact surface extraction as vertices/faces/centroids/areas;
- sampling of combined and per-material signed distances at points;
- material ownership at points;
- velocity extension, advection, redistance, and topology signature;
- line/slice sampling for observables;
- a deterministic geometry fingerprint and operator-epoch payload;
- an explicitly bounded dense diagnostic view, never an implicit production conversion.

Subtasks:

1. Wrap the existing `FeatureGeometry3D` as `UniformFeatureGeometryBackend3D`.
2. Route `_face_material_ids`, `_surface_gas_normals`, surface extraction, domain sizing, topology,
   and profile observables through the contract.
3. Keep `advect_3d`, `reinit_cr2`, `reinit_fsm`, and `extract_mesh_3d` unchanged as the uniform
   implementation.
4. Preserve the public `FeatureGeometry3D` API while downstream users migrate.

Tests/gates:

- uniform backend versus the pre-refactor engine is bitwise identical for geometry no-op,
  fingerprints, periodic projection, and checkpoint round trip;
- one translating-plane step has identical arrays, mesh, surface state, ledgers, and diagnostics;
- one multi-material junction step has identical owner labels and topology classification;
- all existing `test_feature_step_3d.py` tests remain green.

Kill condition: if uniform parity cannot be made exact, do not begin sparse storage; shrink the seam.

### WP-AMR2 — Shared geometric transfer operator

**Purpose:** prevent repeated remeshing from becoming a hidden diffusion/operator change.

**2026-07-17 standalone increment:** `TriangleSurface3D` now provides immutable/order-sensitive
surface identity, exact closest point/distance, whole-triangle periodic images, strict material-local
queries, deterministic ties, and certified centroid-radius candidate pruning. Its bounded path agrees
with brute force in the manufactured suite (`11` tests). It is not wired into either remap. Seam-
crossing triangles must presently be unwrapped by the caller, and the candidate scan is still `O(N)`;
BVH/AABB acceleration and overlap weights remain this work package's next increments.

New source:

- `src/petch/surface_mesh_3d.py`
- `src/petch/surface_overlap_remap_3d.py`
- `tests/test_surface_overlap_remap_3d.py`

Existing source to modify:

- `src/petch/feature_step_3d.py`
- `src/petch/surface_charge_remap_3d.py`
- `tests/test_feature_step_3d.py`
- `tests/test_surface_charge_remap_3d.py`

Subtasks:

1. Factor the exact point-to-triangle certification and periodic images into one reusable
   `TriangleSurface3D` object (**standalone correctness increment complete; BVH still open**).
2. Define one sparse old-face to new-face geometric transfer matrix.
3. First retain the current KNN reconstruction behind that interface and prove no behavioral change.
4. Then add local tangent-chart triangle-overlap weights for coincident/retriangulated and small-motion
   patches.
5. Apply field-specific closures on top of the same weights:
   - conservative extensive surface inventories;
   - bounded intensive coverages with no new extrema;
   - signed charge with positive/negative ledgers and etched-charge removal.
6. Canonicalize surface ordering so far-field block refinement cannot change a physical surface
   fingerprint solely through block enumeration order.

Manufactured tests:

- unchanged mesh/no motion: bitwise state;
- same plane with unrelated triangulation;
- periodic seam translation;
- translating plane and curved cylinder/sphere refinement order;
- material junction with zero cross-material borrowing;
- discontinuous surface field with no new extrema and L1 convergence;
- 100 repeated refine/coarsen/remesh round trips with reported accumulated diffusion;
- uniform recession and growth with exact charge/inventory removal or creation closure;
- old/new overlap rows and columns close their declared area ledgers.

Promotion gate: local field error converges and every extensive ledger closes to roundoff. Global
conservation alone is insufficient.

### WP-AMR3 — Fixed-resolution block-sparse narrow band

**Purpose:** prove sparse storage and kernels before introducing more than one `dx`.

New source:

- `src/petch/block_levelset_3d.py`
- `src/petch/block_surface_extract_3d.py`
- `src/petch/block_topology_3d.py`
- `tests/test_block_levelset_3d.py`
- `scripts/block_levelset_manufactured_audit_3d.py`

Existing source to modify only after standalone gates pass:

- `src/petch/feature_geometry_backend_3d.py`
- `src/petch/feature_step_3d.py`

Subtasks:

1. Implement active/far block storage at one global spacing, with partial boundary blocks and
   periodic neighbor maps.
2. Pack active blocks plus deterministic halos for CPU reference operations.
3. Implement blockwise velocity extension, Godunov advection, and monotone redistance using the same
   equations as the uniform kernels.
4. Extract marching surfaces only from cells owned by active blocks; canonicalize/deduplicate vertices
   and faces across block boundaries.
5. Implement component/topology classification over active cells and uniform far blocks. A dense
   boolean diagnostic view is acceptable during development but not as the production topology path.
6. Rebuild the band deterministically with hysteresis and record active block/node fractions.

Manufactured tests:

- static and translating plane crossing several block boundaries;
- expanding/receding sphere and cylinder with convergence order;
- periodic plane and trench seam equivalence;
- multi-material junction motion;
- keyhole open -> sealed -> reopened under the existing event policy;
- exact hard-visibility triangle closure across every block seam;
- dense/sparse signed-distance and surface Chamfer/Hausdorff comparison;
- dense/sparse topology and material-volume comparison;
- deterministic checkpoint/restart and block-order permutation;
- active-band width sweep proving the selected width is on a result plateau.

Acceptance gates:

- same topology classification as the uniform grid at every checkpoint;
- no surface cracks, duplicate-area contribution, or visibility leak;
- uniform-spacing depth/opening and area-weighted velocity within the measured uniform discretization
  uncertainty;
- material and surface ledgers close;
- at least 3x geometry-memory reduction and 2x level-set-stage speedup on a representative deep
  feature, or a documented no-go before further optimization.

### WP-AMR4 — One-level 2:1 block AMR with a finest interface

**Purpose:** obtain 5 nm interface accuracy without a full-domain 5 nm volume grid.

New source:

- `src/petch/amr_levelset_3d.py`
- `tests/test_amr_levelset_3d.py`
- `scripts/amr_manufactured_audit_3d.py`

Existing source:

- `src/petch/feature_geometry_backend_3d.py`
- `src/petch/block_levelset_3d.py`
- `src/petch/block_surface_extract_3d.py`

Subtasks:

1. Add a 10 nm base level and one 5 nm refined level with 2:1 balance.
2. Refine every sign-changing/material-junction block plus the declared fine buffer.
3. Refine ahead of the projected accepted interface move; refuse a step if the interface could leave
   the fine buffer.
4. Forbid coarsening of any block that intersects a physical surface, material junction, or live
   triangle state.
5. Prolong per-material signed distance into newly activated fine blocks, then redistance while
   freezing the represented zero set.
6. Restrict only far, uniform-sign blocks in v1. This avoids a refine/coarsen surface-volume
   correction in the first authority path.
7. Record refinement/coarsening events, cell counts, work units, and hysteresis decisions in the
   operator epoch and checkpoint.

Manufactured tests/gates:

- plane, sphere, and junction cross block/refinement boundaries without a change in speed or topology;
- 100 refine/coarsen cycles away from the interface leave the physical surface and surface state
  bitwise unchanged;
- periodic refinement decisions and duplicate endpoint planes match;
- block-order and restart decisions are deterministic;
- AMR surface, material volume, and observables converge to uniform 5 nm;
- no result depends on a coarse/fine surface transition, because such a transition is forbidden in v1.

Kill condition: if keeping the whole interface fine does not reduce matched-error runtime or memory
enough to justify the hierarchy, stop before multi-level surface AMR. The next performance target
would be transport/surface-mesh reduction, not more volume-grid machinery.

### WP-AMR5 — Common-engine integration and durable checkpoints

**Purpose:** make AMR a backend of the same engine, not another benchmark solver.

Existing source to modify:

- `src/petch/feature_step_3d.py`
- `src/petch/profile_observables_3d.py`
- `src/petch/charging_checkpoint_3d.py`
- `src/petch/physical_api.py`
- `src/petch/charging_coevolution_3d.py` (AMR refusal/dispatch only at this stage)
- `src/petch/__init__.py`

New tests:

- `tests/test_amr_feature_step_3d.py`
- `tests/test_amr_checkpoint_3d.py`

Subtasks:

1. Select `uniform`, `sparse_uniform`, or `amr_2to1` through one declared geometry backend input.
2. Keep boundary construction, hard-visibility transport, radiosity, mechanisms, material exchange,
   surface products, adaptive outer time stepping, and validity reporting identical.
3. Use finest local `dx` for CFL and displacement gates.
4. Add checkpoint schema v2 for block payloads, refinement masks, periodic maps, hysteresis state, and
   backend/operator epoch; retain safe v1 uniform readers.
5. Make observables sample the hierarchy in physical coordinates instead of assuming dense slices.
6. Report active blocks/cells, finest-surface fraction, remap events, mesh fingerprints, and per-stage
   wall time on every step.
7. Refuse charging on AMR geometry until WP-AMR6 passes; do not silently project to a different field
   operator.

Gates:

- uniform backend still passes the complete existing suite;
- AMR and uniform use the same physics manifest apart from declared numerics;
- interrupted/resumed AMR is identical to uninterrupted execution;
- a backend change creates a new operator epoch and cannot reuse a response/calibration artifact from
  another epoch.

### WP-AMR6 — Electrostatics/charging compatibility

**Purpose:** restore the full unified physics set without pretending the current uniform Q1 system is
already adaptive.

Current blockers:

- `NodalPoissonSystem3D` assembles one uniform rectilinear Q1 grid;
- triangle charge projection assumes each marching-cubes triangle lies in one potential cell;
- charged trajectory interpolation consumes one uniform nodal spacing;
- charging checkpoints and stationarity gates assume one nodal shape.

The lowest-risk bridge is an explicitly **independent uniform electrostatic grid** coupled to the AMR
surface, not an immediate octree Poisson solver.

Existing source to modify:

- `src/petch/charging_poisson_3d.py`
- `src/petch/boundary_transport_3d.py`
- `src/petch/charging_coevolution_3d.py`
- `src/petch/surface_charge_remap_3d.py`
- `src/petch/charging_checkpoint_3d.py`

New tests:

- `tests/test_amr_charging_3d.py`

Subtasks:

1. Generalize triangle-to-Q1 load projection so triangles crossing field-cell boundaries are split or
   integrated cell-by-cell without charge loss.
2. Keep potential interpolation on the declared uniform field grid while profile geometry is sparse
   or refined.
3. Remap signed surface charge through the shared overlap operator and preserve the removal ledger.
4. Refine the independent electrostatic grid until potential, deflected flux, B1/B2 diagnostics, and
   profile velocity are invariant.
5. Only then consider adaptive Q1 assembly with hanging-node constraints as a separate project.

Gates:

- exact surface-to-node charge conservation;
- uniform-profile-grid versus AMR-profile-grid charging agreement at the same field grid;
- field-grid refinement convergence;
- exact hard-visibility final audit and existing charging stationarity contract;
- no AMR Nozawa/notching claim before these gates pass.

### WP-AMR7 — GPU packing, Krüger short certification, and promotion

**Purpose:** earn performance only after CPU/reference correctness.

New/modified source:

- packed-block Warp kernels in `src/petch/block_levelset_3d.py` and
  `src/petch/amr_levelset_3d.py`;
- `scripts/krueger_2024_amr_short_burst_audit.py`;
- `tests/test_amr_cpu_gpu_parity.py`.

Subtasks:

1. Pack active blocks into stable arrays plus integer neighbor maps; do not expose Python dictionary
   order to GPU reductions.
2. Port advection/redistance/extension only after the CPU block operator passes.
3. Verify CPU/GPU geometry, surface, topology, and ledger parity.
4. Run the same base-only initial 0.5 s and late 0.1--0.5 s windows at uniform 5 nm and AMR 10/5 nm.
5. Compare depth/opening increments, surface velocity area-L1/RMS/worst-feature errors, surface-state
   integrals, topology, wall time, memory, and work units.
6. Promote AMR to a calibration authority only after matched-error evidence; bind any new response
   model and final base confirmation to the AMR operator epoch.

Promotion gates:

- AMR lies inside uniform 5 nm numerical uncertainty for every declared short-burst observable;
- no experimental parameter changes between the paired runs;
- topology and ledgers match;
- end-to-end wall time improves by at least a predeclared useful amount (recommended initial target:
  25%) and memory improves materially;
- only after all of the above may one clean long base confirmation run begin.

## 5. Test matrix before any long AMR run

| Layer | Minimum manufactured evidence | Failure meaning |
| --- | --- | --- |
| Backend seam | uniform no-op and one-step bitwise parity | refactor changed the operator |
| Shared remap | retriangulation, seam, junction, round trips, recession/growth | AMR would diffuse or misplace state |
| Sparse band | plane/sphere/junction/keyhole plus width sweep | band or block seams change geometry/topology |
| 2:1 hierarchy | interface crossing block boundaries, deterministic restart | refine/coarsen is not path-independent |
| Surface mesh | no cracks, area duplication, or hard-visibility leaks | transport operator is invalid |
| Observables | hierarchy sampling agrees with uniform physical coordinates | validation scorer changed with backend |
| Charging | independent-grid charge/field/profile refinement | AMR cannot yet support notching claims |
| Krüger | initial and late short paired bursts | no authority for a 60 s run |

## 6. Principal risks and controls

| Risk | Why it matters here | Binding control |
| --- | --- | --- |
| Remap-induced operator drift | The R17/R19 history already shows path feedback can amplify small remap changes | shared overlap operator; round-trip diffusion gate; operator epoch binding |
| Coarse/fine surface cracks | A triangle seam can become a hard-visibility particle leak | no coarse/fine interface crossing in v1; canonical surface seam tests |
| Lost cavity/topology information | A sparse far field can erase a sealed component | uniform-label far blocks plus sparse connectivity; keyhole open/sealed/reopen test |
| Material mixing at junctions | `max/argmax` ownership is sensitive near equal level sets | junction always fine; material-local prolongation/remap; no cross-material weights |
| Periodic endpoint duplication | This previously caused real remap failure | canonical physical endpoint, periodic block-neighbor tests, checksum-bound policy |
| Stochastic mesh decisions | Repeated sample error could steer geometry as frozen samples steered charging | geometry-only v1 indicator; deterministic hysteresis; no raw flux indicator |
| Checkpoint non-reproducibility | Adaptive layouts can depend on iteration/dictionary order | sorted block keys, serialized hysteresis, deterministic reductions and restart test |
| False performance promise | Transport may dominate even after volume work shrinks | stage profile before build; matched-error end-to-end gate; stop if gain is immaterial |
| Charging overclaim | Current Q1 and triangle coupling are uniform-grid contracts | explicit AMR charging refusal until WP-AMR6 |
| Licensing contamination | Vienna's HRLE implementation is GPL-3.0 | behavioral reference only; independent implementation or explicit product/license decision |

## 7. What AMR will and will not solve

AMR can:

- make a 5 nm (and later local 2.5 nm) interface affordable without fine cells throughout the full
  volume;
- resolve narrow openings, necks, corners, and material junctions without calibrating away coarse
  geometry error;
- reduce redistance, extension, topology, checkpoint, and memory cost;
- make paired fine confirmations and multi-fidelity calibration practical.

AMR cannot by itself:

- repair an incorrect chemistry law or reactor boundary;
- make a stale response model valid across an operator change;
- reduce hard-visibility/radiosity cost if every surface triangle remains fine;
- replace robust surface-state/charge remapping;
- validate a held-out experiment;
- make charging adaptive before the electrostatic coupling is refined and verified.

## 8. Recommended execution order

```text
WP-AMR0  uniform current-operator references and cost profile
    |
WP-AMR1  exact uniform backend seam
    |
WP-AMR2  shared certified/overlap surface transfer
    |
WP-AMR3  fixed-dx sparse narrow band
    |
WP-AMR4  one 10 nm -> 5 nm 2:1 hierarchy, interface always fine
    |
WP-AMR5  common-engine integration + durable checkpoints
    |
    +--------------------------+
    |                          |
WP-AMR7 Krüger base-only        WP-AMR6 independent-grid charging
short certification             compatibility
    |                          |
one authoritative base run     Nozawa/notching AMR eligibility
```

The immediate implementation milestone is **WP-AMR1 plus the standalone half of WP-AMR2**. That is
bounded, reversible engineering and directly addresses the operator-drift risk. Sparse storage should
not enter `advance_feature_step_3d` until the uniform backend is exactly identical and the shared
surface transfer passes its manufactured matrix.

## 9. Definition of done for the first AMR release

The first AMR release is complete only when all of the following are true:

1. the uniform backend remains the bitwise reference and all existing tests pass;
2. one common engine selects uniform, sparse-uniform, or 10/5-nm AMR numerics by a declared input;
3. the physical surface and all its stencils remain at 5 nm in AMR v1;
4. every remap/refine/coarsen event closes material, surface inventory, and charge ledgers applicable
   to the enabled physics;
5. periodic, topology, hard-visibility, checkpoint/restart, and CPU/GPU gates pass;
6. AMR reproduces uniform 5 nm manufactured and base-only short-burst results at lower matched-error
   cost;
7. the operator epoch prevents reuse of incompatible response/calibration artifacts;
8. charging is either independently certified through WP-AMR6 or explicitly refused;
9. no held-out experimental datum has influenced numerical refinement, parameters, or promotion;
10. only then is one long authoritative base run justified.
