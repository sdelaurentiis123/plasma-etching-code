# ViennaPS topology conformance: pinch-off, enclosed cavities, and reopening

Date: 2026-07-17
Status: read-only upstream source audit plus bounded local ViennaLS manufactured run; no petch engine
files modified by this audit

Related documents:

- [`VALIDATION_FIRST_SUPERSET_CAMPAIGN_2026-07-17.md`](VALIDATION_FIRST_SUPERSET_CAMPAIGN_2026-07-17.md)
- [`VIENNAPS_CODEBASE_AUDIT_2026-07-17.md`](VIENNAPS_CODEBASE_AUDIT_2026-07-17.md)
- [`NUMERICS_CALIBRATION_AMR_PLAN_2026-07-17.md`](NUMERICS_CALIBRATION_AMR_PLAN_2026-07-17.md)

## 1. Executive verdict

ViennaPS does not need a special topology solver to pass through ordinary pinch-off and reopening.
Its material interfaces are implicit level sets, so connectivity may change during the same CFL-limited
Hamilton--Jacobi evolution used for every other step. The ordinary process loop does not stop merely
because an open gas feature becomes an enclosed cavity.

That is the production behavior petch should match. It does **not** imply that Vienna supplies the
entire petch contract:

- Vienna does not classify enclosure/reopening as a first-class process event or emit a conservation
  receipt for it.
- Vienna's ordinary coverage propagation copies values through sparse-grid point lineage. It is not a
  conservative old/new triangle-overlap remap and does not conserve a finite surface inventory by
  construction.
- An enclosed surface receives no direct external first-hit flux when a watertight triangle backend is
  used, but what happens to that surface under intrinsic chemistry, desorption, diffusion, or an
  analytic velocity remains a model choice.
- `ignoreVoids=true` can freeze enclosed interfaces; the ViennaPS default is `false`.

The strict-superset target is therefore:

```text
Vienna's natural implicit topology continuation
                         +
petch's explicit event classification, accessibility gate,
material-local conservative state transfer, and refusal ledger
```

## 2. Immutable source provenance

The campaign plan names ViennaPS commit
[`2956ed587984c6dc38be24c6e2390e10c9b2f0a7`](https://github.com/ViennaTools/ViennaPS/tree/2956ed587984c6dc38be24c6e2390e10c9b2f0a7),
dated 2026-07-02. Its CMake project version is `4.6.1`, and it pins ViennaLS `v5.8.3`.

There is a provenance nuance: the immutable Git tag `v4.6.1` currently resolves to the earlier commit
`1dd15f43896bb053836e5c9351a992c962dfab12` (2026-06-25), not to `2956ed...`. The topology-relevant
ViennaPS controls and the pinned ViennaLS implementation are unchanged between these snapshots. The
later commit changes ordinary flux/desorption plumbing, not the level-set topology mechanism. Future
black-box runs should record `2956ed...` explicitly and label it `4.6.1+2956ed`, rather than relying on
the moving meaning of a version string.

| Component | Immutable revision | Role in this audit |
| --- | --- | --- |
| ViennaPS campaign snapshot | `2956ed587984c6dc38be24c6e2390e10c9b2f0a7` | process loop and parameters |
| ViennaPS `v4.6.1` tag | `1dd15f43896bb053836e5c9351a992c962dfab12` | release-tag cross-check |
| ViennaLS `v5.8.3` | `c3a1a2ad5bb0b05f75cd4beb71b540f6c544a287` | implicit interface and void detection |
| ViennaCore `v2.2.1` | `640411618bef44717487f25ec7417b4f6daaabf7` | point-data translation |
| ViennaHRLE `v1.1.2` | `8aa4aa7ac4a72d101bee8a140dbf64e20963e305` | sparse level-set storage |
| ViennaRay `benchmark` branch at audit time | `6cd664ceddebfa2ef67bda944e68886d880f889e` | source-plane first-hit transport |

ViennaPS 4.3 and newer is GPL-3.0. This report records behavior and an independent conformance design;
it does not copy Vienna implementation into petch.

## 3. What Vienna actually does

### 3.1 Ordinary process loop

The flux-process strategy repeatedly performs:

```text
extract current surface
        ↓
trace source and optional surface flux
        ↓
update coverages and calculate velocity
        ↓
copy coverages onto the top sparse level set
        ↓
advect one CFL-limited level-set step
        ↓
regenerate surface and copy coverages back
```

The loop condition is only process time, and the advection handler returns success after a normal
step. There is no connectivity test that turns a gas-cavity enclosure into early termination:

- [ordinary loop and step order](https://github.com/ViennaTools/ViennaPS/blob/2956ed587984c6dc38be24c6e2390e10c9b2f0a7/include/viennaps/process/psFluxProcessStrategy.hpp#L249-L380);
- [single-step ViennaLS advection wiring](https://github.com/ViennaTools/ViennaPS/blob/2956ed587984c6dc38be24c6e2390e10c9b2f0a7/include/viennaps/process/psAdvectionHandler.hpp#L41-L128).

ViennaLS stores a signed interface in a sparse narrow band. Its advection rebuild finds new sign
crossings from neighboring grid values; that representation naturally permits components to merge or
separate without following triangle identities. This is why pinch-off is not a triangle-correspondence
problem in Vienna's geometry operator:

- [advection rebuild and point lineage](https://github.com/ViennaTools/ViennaLS/blob/c3a1a2ad5bb0b05f75cd4beb71b540f6c544a287/include/viennals/lsAdvect.hpp#L238-L428);
- [CFL integration and material transitions](https://github.com/ViennaTools/ViennaLS/blob/c3a1a2ad5bb0b05f75cd4beb71b540f6c544a287/include/viennals/lsAdvect.hpp#L431-L611).

### 3.2 Enclosed-interface policy

ViennaLS can mark surface points not belonging to the externally connected positive component. When
`ignoreVoids=true`, those points receive zero advection rate. When false, the selected model velocity
is evaluated normally:

- [connected-component void marking](https://github.com/ViennaTools/ViennaLS/blob/c3a1a2ad5bb0b05f75cd4beb71b540f6c544a287/include/viennals/lsMarkVoidPoints.hpp);
- [rate suppression only when `ignoreVoids` is enabled](https://github.com/ViennaTools/ViennaLS/blob/c3a1a2ad5bb0b05f75cd4beb71b540f6c544a287/include/viennals/lsAdvect.hpp#L431-L530);
- [public `setIgnoreVoids` semantics](https://github.com/ViennaTools/ViennaLS/blob/c3a1a2ad5bb0b05f75cd4beb71b540f6c544a287/include/viennals/lsAdvect.hpp#L893-L899).

ViennaPS exposes that switch and defaults it to `false`:

- [ViennaPS `AdvectionParameters::ignoreVoids = false`](https://github.com/ViennaTools/ViennaPS/blob/2956ed587984c6dc38be24c6e2390e10c9b2f0a7/include/viennaps/process/psProcessParams.hpp#L52-L94).

The correct petch analogue is not a universal `freeze_cavities` flag. It is channel-specific:

- direct external ion/neutral/electron source on a sealed surface: zero by hard visibility;
- local stored-film reaction or surface diffusion: evolve if the declared model permits it;
- internal desorption/re-emission: evolve and transport inside the cavity if that channel is enabled;
- purely geometric manufactured velocity: evolve unless the test explicitly requests the frozen-void
  variant.

### 3.3 Why direct external flux is zero

Vienna's default random source launches particles from the source-side bounding plane. The triangle
kernel asks Embree/OptiX for the nearest intersection and processes that surface hit before any
reflection. A watertight cap therefore shields every surface inside a sealed cavity:

- [source origins on the source-side bounding plane](https://github.com/ViennaTools/ViennaRay/blob/6cd664ceddebfa2ef67bda944e68886d880f889e/include/viennaray/raySourceRandom.hpp#L25-L85);
- [nearest-intersection and surface-hit processing](https://github.com/ViennaTools/ViennaRay/blob/6cd664ceddebfa2ef67bda944e68886d880f889e/include/viennaray/rayTraceKernel.hpp#L153-L335);
- [ViennaPS triangle surface regeneration](https://github.com/ViennaTools/ViennaPS/blob/2956ed587984c6dc38be24c6e2390e10c9b2f0a7/include/viennaps/process/psCPUTriangleEngine.hpp#L89-L155).

This is a structural zero for exact first-hit triangle visibility, subject to the mesh being closed and
the intersection implementation being numerically sound. A disk backend is an approximate surface
representation and should be scored with an estimator tolerance rather than assumed exact.

### 3.4 `RemoveStrayPoints` is not the ordinary cavity path

ViennaLS also provides an explicit cleanup that keeps one selected connected surface and removes the
others. The ViennaLS test constructs a substrate with an internal spherical cavity and demonstrates
that direct use of `RemoveStrayPoints` removes that internal interface:

- [ViennaLS cleanup implementation](https://github.com/ViennaTools/ViennaLS/blob/c3a1a2ad5bb0b05f75cd4beb71b540f6c544a287/include/viennals/lsRemoveStrayPoints.hpp);
- [manufactured internal-cavity cleanup test](https://github.com/ViennaTools/ViennaLS/blob/c3a1a2ad5bb0b05f75cd4beb71b540f6c544a287/tests/RemoveStrayPoints/RemoveStrayPoints.cpp).

The higher-level ViennaPS `Domain::removeStrayPoints()` first inverts the top surface, applies that
cleanup, and inverts back. It is aimed at disconnected solid fragments, and it is an explicit user
operation rather than part of every `Process.apply()` step:

- [ViennaPS domain cleanup wrapper](https://github.com/ViennaTools/ViennaPS/blob/2956ed587984c6dc38be24c6e2390e10c9b2f0a7/include/viennaps/psDomain.hpp#L319-L338).

Petch must not fill or delete a resolved cavity under a generic cleanup label. Unresolved subcell
bubbles may be treated only under the already declared resolution/ledger rule; resolved topology is a
physical event.

## 4. Bounded manufactured run performed locally

### 4.1 Fixture

An original 2-D ViennaLS-only probe was compiled against the exact dependency revisions above. It did
not require VTK, ViennaRay, a GPU, or a long process model.

| Item | Value |
| --- | --- |
| Grid | uniform `dx = 0.25`, domain `[-10,10] x [-10,10]` |
| Initial solid | half-plane below `y=0` |
| Gas feature | narrow neck of width `2.0` joined to a circular chamber of radius `3.0` |
| Deposit | constant normal velocity `+1.0` for `1.5` time units |
| Reverse etch | constant normal velocity `-1.0` for `3.5` time units |
| Scheme | first-order Engquist--Osher, CFL ratio `0.25` |
| Event detector | ViennaLS `MarkVoidPoints`, recomputed after checkpoints |
| Local execution | approximately `0.01 s` after compilation on the current Mac |

The fixture is deliberately keyhole-shaped. Uniform deposition closes the narrow neck while leaving
a resolved chamber; reverse etching removes the bridge and reconnects the chamber to the external gas.

### 4.2 Result

| State | Phase time | Connected-component count | Marked enclosed points | Interpretation |
| --- | ---: | ---: | ---: | --- |
| Initial | 0.0000 | 2 | 0 | chamber connected through open neck |
| Last recorded open deposit state | 0.8881 | 2 | 0 | neck narrowing |
| First recorded sealed state | 1.0486 | 3 | 84 | gas cavity formed naturally |
| End of deposition | 1.5000 | 3 | 68 | cavity retained; process did not stop |
| Reverse etch, default void policy | 0.5000 | 3 | 84 | enclosed interface also advected |
| Reverse etch, `ignoreVoids=true` | 0.5000 | 3 | 68 | enclosed interface remained frozen |
| First recorded reopened state | 1.8048 | 2 | 0 | external gas reconnected naturally |
| End of reverse etch | 3.5000 | 2 | 0 | open topology continued normally |

The important evidence is qualitative and structural, not the coarse-grid event time:

1. closure occurred without a special mesh operation;
2. the level-set integrator continued after closure;
3. the same interface reopened under reversed motion;
4. the void-policy switch changed only whether the enclosed interface itself moved.

The event times must be refined before being used as a numerical benchmark. At `dx={0.25,0.125,0.0625}`
the enclosure/reopening times should converge and the final surfaces should be compared by Chamfer and
critical-dimension errors.

### 4.3 Local petch policy test status

The current targeted petch topology-policy tests pass:

```text
4 passed, 63 deselected in 1.01 s
```

They verify explicit refusal versus `continue_gas_cavity`, event labels, exact zero direct source on
manufactured sealed faces, conservative state diagnostics, and continued refusal of solid/material
component changes. The open-to-sealed-to-open geometry in that test is injected at the advection
boundary. It is therefore a valid policy/transport/remap test, but it is **not yet** the paired
end-to-end proof that petch's own geometric advection spontaneously pinches and reopens the fixture.
That end-to-end case remains required below.

## 5. Surface-state and coverage semantics

Vienna's geometry and coverage paths should not be conflated.

### Geometry

The sparse signed field is authoritative. Triangles/disks are regenerated views for transport, so
triangle identities may change freely and topology remains robust.

### Coverage

Before advection, ViennaPS maps each surface coverage value to a sparse level-set point through a
translator. During `Reduce`/rebuild, ViennaLS records one old source point ID for each new sparse point
and copies its scalar/vector data. After advection, ViennaPS regenerates the surface and reads the level-
set values back through the new translator:

- [coverage to/from level-set mapping](https://github.com/ViennaTools/ViennaPS/blob/2956ed587984c6dc38be24c6e2390e10c9b2f0a7/include/viennaps/process/psAdvectionHandler.hpp#L131-L187);
- [old-point selection during rebuild](https://github.com/ViennaTools/ViennaLS/blob/c3a1a2ad5bb0b05f75cd4beb71b540f6c544a287/include/viennals/lsAdvect.hpp#L258-L423);
- [literal indexed data copy](https://github.com/ViennaTools/ViennaCore/blob/640411618bef44717487f25ec7417b4f6daaabf7/include/viennacore/vcPointData.hpp#L312-L355).

This is efficient and avoids matching arbitrary triangle centroids, but it has different semantics
from petch's finite material inventories:

| Property | Vienna point lineage | Required petch authoritative remap |
| --- | --- | --- |
| Intensive coverage stays bounded | usually, by copied values/model update | yes, explicitly gated |
| Integrated finite inventory conserved | not guaranteed | yes, material-local ledger |
| Retriangulation invariance | geometry yes; surface sample values depend on grid lineage | yes within order |
| Merge/split ownership | one selected source point | overlap/declared new-surface closure |
| Removed material itemized | no generic surface ledger | yes |
| AMR prolong/restrict | not supplied by fixed-grid lineage | required |

The lesson is not “copy Vienna's remapper.” It is:

1. keep state attached to an Eulerian interface representation rather than ephemeral triangle IDs;
2. use a spatial acceleration structure for correspondence;
3. add conservative overlap/area semantics where the state represents a finite inventory;
4. make topology and AMR transfers explicit events with receipts.

## 6. Small paired petch/Vienna conformance suite

All cases below are manufactured and should run locally. No Krüger held-out data and no long GPU run
is needed.

### T1 — Geometry-only keyhole close/reopen

Run the exact fixture from Section 4 in ViennaLS and petch at three grid levels.

Record:

- first enclosure and first reopening times;
- gas, solid, and material component counts;
- neck width and chamber area/volume versus time;
- surface Chamfer distance at matched physical times;
- timestep count, wall time, and any recovery/refusal.

Gate: both converge toward the same morphological normal-motion solution. A topology transition is an
accepted event, not a terminal success or an unclassified exception.

### T2 — Sealed-cavity hard visibility

Freeze matched open, just-sealed, and reopened geometries. Launch a collimated external source using
Vienna's CPU triangle backend and petch's exact hard-visibility backend.

Record external first-hit flux separately on cap, outer wall, and enclosed surface.

Gate:

- enclosed direct flux is exactly zero for the deterministic triangle fixture, or statistically
  consistent with zero for a randomized estimator;
- cap plus escape closes the source ledger;
- reopening restores nonzero interior access;
- disk-backend leakage, if measured, is reported rather than used as triangle truth.

### T3 — Intensive stripe through topology

Initialize a bounded coverage stripe on one chamber wall, close and reopen geometrically, and compare
the transported field.

Gate:

- no values leave declared bounds;
- translation/retriangulation away from the event converges;
- petch conserves every extensive inventory separately;
- Vienna's integral drift is measured as comparator behavior, not silently granted as truth.

### T4 — Disappearing/created surface closure

Use a bridge whose old surface area disappears at merger and new surface appears at separation.

Gate: petch assigns disappeared inventory to the declared removal ledger and initializes genuinely new
surface through the declared physical closure. No nearest-neighbor value may teleport across the
cavity.

### T5 — Multi-material pinch

Repeat with mask and substrate level sets.

Gate: solid and per-material component counts remain unchanged across the allowed gas event; owner IDs
never average across the junction. Any material merger/fragmentation remains a structured refusal.

### T6 — Periodic translation

Translate the same neck across a lateral periodic seam and repeat T1--T3.

Gate: event time, geometry, accessibility, and ledgers are invariant under the translation. This test
specifically detects a nonperiodic nearest-surface/remap search.

### T7 — Cleanup controls

Apply direct ViennaLS `RemoveStrayPoints` and ViennaPS `Domain::removeStrayPoints()` as separate,
explicit controls. Never include either in the ordinary T1 path.

Gate: the report distinguishes an intentionally removed disconnected interface/component from natural
level-set continuation.

## 7. Runtime expectations

| Run | Expected budget | Hardware |
| --- | --- | --- |
| ViennaLS T1 geometry at one grid | `<1 s` | local CPU |
| Three-level ViennaLS/petch T1 | seconds to a few minutes | local CPU |
| T2 small triangle-flux fixture | seconds to `<1 min` per backend | local CPU |
| T3--T7 | minutes total when kept 2-D/small extruded 3-D | local CPU |
| Full regression | existing suite scale | local CPU |

No external GPU is justified for this topology-conformance block. GPU comparison begins only after the
CPU reference cases close and uses the same fixed geometries and observables.

## 8. What Vienna teaches us, and what it does not

### Adopt

- implicit level-set topology should continue without human intervention;
- triangles are transport views, not persistent material elements;
- sparse narrow-band evolution is the first scaling target before general AMR;
- an explicit void-advection control is useful for models whose internal-interface physics differs;
- ordinary geometry operations should finish in seconds on manufactured cases.

### Retain from petch

- distinguish gas enclosure/opening from solid or material fragmentation;
- hard-visibility accessibility certification;
- material-local extensive and intensive remap modes;
- signed material/product/charge ledgers;
- periodic topology and remap invariance;
- structured refusal for unsupported physics, not for routine topology;
- exact checkpoint/provenance at the event.

### Do not copy as authoritative truth

- point-ID copying for a finite surface inventory;
- a disk flux result as exact cavity visibility;
- removal of a resolved interface under generic cleanup;
- a version-only comparator without an immutable SHA;
- a successful completion as evidence of conservation or experimental accuracy.

## 9. Ordered implementation/conformance decision

1. Keep the current petch `continue_gas_cavity` policy and its stronger refusal taxonomy.
2. Add/run the non-injected T1 end-to-end petch geometry case; this is the immediate missing proof.
3. Run T2 on frozen surfaces before resuming a long experimental checkpoint.
4. Complete T3/T4 before calling the current nearest-surface correction production-grade.
5. Use T1--T4 as mandatory gates for conservative triangle-overlap remapping.
6. Use the same suite again for sparse narrow bands and later AMR; AMR must reproduce the uniform-grid
   result rather than redefine it.
7. Only then resume the archived Krüger state under a bounded wall budget.

The practical conclusion is favorable: the cavity problem is not an impossible etch-physics problem.
It is a bounded topology/state-transfer maturity gap. Vienna shows that geometry continuation itself is
routine; petch's opportunity is to make that routine continuation conservative, auditable, periodic,
and compatible with richer physics.
