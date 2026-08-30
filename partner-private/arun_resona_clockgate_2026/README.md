# Arun / Resona clock-gate geometry handoff

Private partner data, intentionally co-located with the engine in the private
multiphysics repository. Do not mirror this directory or its raw geometry to a
public/shared remote.

## Received artifact

- `raw/CLKGATE X1 25x Reduced 3D.STL`
- SHA-256: `5df028e30436ca2313fd58a79a3240b4db7775ccf62301b87866612109f901fa`
- Binary STL; 19,334 bytes

## Verified geometry result

The received mesh contains 208 vertices and 385 triangular facets. Nine facets
are exact zero-area export artifacts. Removing only those facets produces
`derived/clkgate_x1_topology_repaired.stl`, a closed, manifold, consistently
oriented mesh with 376 facets and ten disconnected solid components. The
relative signed-volume change is `2.13e-15`, so this is a topology-preserving
repair, not a geometric redesign.

The file is a planar footprint extruded along its first coordinate axis. Its
overall file-unit extents are approximately `0.0300 x 0.0988 x 0.0628`.

## Recovered source, units, and scale

The scale is no longer an open guess. The binary STL header identifies
`CLKGATE_X1`, layer 11. The public FreePDK45/Nangate source GDS contains that
exact cell, and layer 11 is Metal 1. All ten polygons and all 104 unique planar
vertices match the STL after one isotropic affine transform. The maximum
bidirectional mismatch is 0.19 nm after physical scaling, which is binary-STL
float32 roundoff rather than a geometric discrepancy.

The recovered workflow is:

- the source Metal-1 cell is `2.47 x 1.57 um` on a `5 nm` GDS grid;
- the STL was authored in millimetres and must receive the normal Nanoscribe
  `1000x` millimetre-to-micrometre import conversion;
- the actual geometry enlargement from the source GDS is exactly `40x`;
- the resulting printed geometry is `98.8 x 62.8 um` in plan and `30 um` tall;
- the enlarged layout grid is `0.2 um`;
- the narrowest source inter-polygon gap is `65 nm`, hence `2.6 um` after
  enlargement; the common `70 nm` Metal-1 tracks become `2.8 um` wide.

The filename's `25x` therefore cannot be a geometry scale. It is consistent
with a Nanoscribe 25x objective/process preset. STL itself has no unit field,
so the millimetre label remains an evidence-backed workflow inference; the
physical dimensions are independently pinned by the exact 40x source match.

Projection through the unique extrusion axis gives a solid plan area of
`0.00275824` file-units squared inside a `0.00620464` file-units-squared
bounding box: a scale-independent fill fraction of `0.444545`. The ten
component plan areas are preserved separately in the audit so physical areas
follow immediately once one length scale is supplied.

The complete machine-readable receipt is `results/geometry_audit.json`; the
human-checkable rendering is `results/geometry_preview.png`.

## Working pattern-transfer interpretation

The default transfer polarity is that the printed STL solids protect silicon,
so the target is a raised silicon replica of the Metal-1 wiring. Invert that
polarity only if Arun wants Metal-1-shaped trenches. The 30 um printed polymer
is a tall entrance mask: at the minimum 2.6 um opening, its entrance aspect
ratio is about 11.5 before any silicon is etched. Feature transport and mask
sidewall interactions therefore matter even though the lateral dimensions are
large compared with the published nanograting reference.

## Selected first SF6/O2 transfer recipe

The evidence-selected first shot is the undercut-minimizing Miao et al. Oxford
Plasmalab 100 cryogenic ICP recipe:

- substrate: silicon;
- mask: 30 um printed polymer footprint, optionally with a thin angled Cr cap;
- chuck temperature: `-110 +/- 2 C`;
- pressure: `8 mTorr`;
- ICP/source power: `1000 W`;
- RF platen power: `10 W`;
- flows: `52 sccm SF6`, `8 sccm O2`;
- first-shot duration: `2 min`.

Miao et al. measured `3.5--3.9 um` depths in two minutes across polymer, Cr,
SiO2, and Cr-on-polymer masks on 400 nm-pitch gratings. The polymer selectivity
was about 15, while Cr exceeded 500; the same recipe reached 10.6 um after ten
minutes with Cr-on-polymer. The direct transfer expectation for the first shot
is therefore a roughly `3.5--3.9 um` silicon depth, with only
`0.23--0.26 um` polymer loss if Arun's printed polymer shares that selectivity.
That mask-loss estimate is a conditional feasibility bound, not a material
identity claim.

This recipe wins over the faster 92/8 sccm literature branch because the latter
reported a bottle-shaped profile and 40--60 nm sidewall undulation at its
high-rate point. The fidelity challenge should first minimize lateral error;
throughput is a tie-breaker, not the objective.

The machine-readable selection and the minimal nine-run derivative board are
in `results/sf6_o2_recipe_board.json`.

The exact ten-polygon footprint has also been converted into the common
periodic layered 3-D engine geometry without boxes or a hand-redrawn mask. At
the production `0.2 um` mesh the minimum opening has 13 cells and the domain is
`495 x 315 x 193` nodes (30,093,525 voxels). A `0.4 um` pilot build is finite,
periodic, and contains gas, silicon, mask, and base materials. The receipt is
`results/feature_geometry_audit.json`.

## Remaining confirmations and measurements

1. Confirm raised-silicon versus trench polarity and the desired silicon depth.
2. Confirm the actual silicon etcher can hold `-110 C` and identify its model.
3. Record exact-run DC self-bias or platen voltage, helium backside pressure,
   and sample position if the tool exposes them.
4. Measure pre/post polymer height and one blanket or wide-open Si depth on the
   same run. Those two numbers identify tool-specific selectivity and absolute
   rate without fitting the patterned result.
5. Preserve a pre-etch print metrology image and a post-etch cross-section for
   the blind profile score.

No computational ambiguity remains in the STL scale or mask footprint. A
unique tool-specific absolute profile still depends on the etcher boundary and
the actual printed-polymer response; until those are measured, simulation must
propagate them as declared physical intervals rather than silently tune them.

## Reproduction

The generic audit implementation and tests live in the private working branch
of `plasma-etching-code`:

```text
python scripts/audit_stl_geometry.py \
  "partner-private/arun_resona_clockgate_2026/raw/CLKGATE X1 25x Reduced 3D.STL" \
  --audit-output "partner-private/arun_resona_clockgate_2026/results/geometry_audit.json" \
  --repaired-output "partner-private/arun_resona_clockgate_2026/derived/clkgate_x1_topology_repaired.stl"
```

Reproduce the source/scale match after obtaining the public Nangate GDS:

```text
python scripts/audit_source_scale.py \
  --reference-gds /path/to/nangate-stdcells.gds --check

python scripts/audit_feature_geometry.py \
  --reference-gds /path/to/nangate-stdcells.gds --check
```

## Primary sources

- FreePDK45/Nangate GDS: https://github.com/mflowgen/freepdk-45nm
- FreePDK45 layer map: https://github.com/The-OpenROAD-Project/OpenROAD-flow-scripts/blob/master/flow/platforms/nangate45/FreePDK45.lyt
- Nanoscribe CAD unit convention: https://www.nanoscribe.com/en/contact-support/support/cad-model-creation/
- Miao et al. 2016: https://doi.org/10.1109/JMEMS.2016.2593339
- Wu et al. 2011: https://doi.org/10.1016/j.mee.2010.11.055
- Zhang et al. 2023: https://doi.org/10.3390/mi14040846
