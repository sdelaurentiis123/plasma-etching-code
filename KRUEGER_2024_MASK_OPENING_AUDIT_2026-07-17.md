# Krueger 2024 mask-opening observable audit

Date: 2026-07-17
Scope: base-development artifacts only; no held-out profile was opened; no transport or profile
evolution was run.

## Result

The final reopened R1.8 checkpoint has a definition-correct development value

`w_m = 39.46998 nm`.

That value is not yet grid-authoritative. A clean, same-operator 10/5 nm endpoint comparison is
still required before this observable can freeze or steer held-out validation.

One historical reporting defect and one reference-frame hazard were found:

1. The paper assigns `w_m = 0` while the mask is clogged. The R1.8 history instead retained the
   nonzero width of the sealed internal gas pocket at all 53 enclosed-state checkpoints. This does
   not affect the final reopened scalar, but those 53 transient values are not the paper's
   observable.
2. The R1.8 mesh is translated: the substrate top is 1.8 um. Calling the old metric with its generic
   1.0 um default measures the wrong vertical band and produces 38.89991 nm. With the declared
   1.8 um reference, the executed union-field value is 39.46637 nm. The certification tool now
   infers the benchmark geometry from the audit and refuses when it is absent; it never silently
   assumes 1.0 um.

## Primary definition

Primary source: [Krueger et al., JVSTA 42, 043008 (2024)](https://cpseg.eecs.umich.edu/pub/articles/JVSTA_42_043008_2024.pdf),
doi:10.1116/6.0003554.

- Section V and Fig. 7 define `w_m` as the narrowest horizontal width through the mask, including
  deposited polymer. Its physical role is necking/clogging.
- Table IV declares the base target as 45 nm.
- Section VIII states that a fully clogged feature has `w_m = 0`, irrespective of the remaining
  internal pocket width.
- Section IV says the reference MCFPM used 1 nm cubic voxels. Section VII reports final optimizer
  error at roughly one voxel or less. Grid-scale metric motion is therefore part of the reference
  method too; it must be priced, not hidden.

The authoritative geometric construction is consequently:

1. Reconstruct the two mask sidewalls from the **mask material level set**, not the union field.
2. Measure their horizontal separation at common physical height.
3. Minimize that separation through the mask thickness.
4. Average equivalent periodic-extrusion rows for this nominal line trench and report their spread.
5. Qualify the result by exterior-to-mask-bottom gas connectivity. The lower witness is the first
   resolved plane inside the mask, not a gas node below the substrate; the latter would falsely
   close a shallow subcell etch. If disconnected, report exactly zero while retaining the
   enclosed-pocket width as a separate diagnostic.

No persistence filter, rolling average, quantile, morphological opening, or soft geometry is
allowed. Those operations could erase a real short neck. The dense cross-section scan evaluates
the same piecewise-bilinear level-set reconstruction; it changes no geometry.

## R1.8 5 nm evidence

Inputs:

- `results/krueger_2024_base_calibration_r18/mixed_operator_topology_continuation/audit.json`
- `results/krueger_2024_base_calibration_r18/mixed_operator_topology_continuation/checkpoint.npz`
- generated diagnostic:
  `results/krueger_2024_base_calibration_r18/mixed_operator_topology_continuation/opening_certification.json`

Final checkpoint:

| Quantity | Value |
| --- | ---: |
| Audit-reported endpoint | 39.46637 nm |
| Mask-level-set minimum | 39.46998 nm |
| Dense reconstruction minimum | 39.46998 nm |
| Minimum location | z = 2.325 um |
| Widths across three interior extrusion rows | 37.04055, 42.29175, 39.07763 nm |
| Cross-extrusion span | 5.25119 nm = 1.050 cells |
| Combined-union contrast value at correct 1.8 um reference | 39.46637 nm at z = 2.325 um |
| Combined-union minus mask-layer value | -0.00361 nm |
| Final exterior-to-mask-bottom connectivity | open |

The exact nodal and dense reconstructed minima coincide. The remaining cross-extrusion spread is
not averaged out of the evidence; its half-span is 2.626 nm and must accompany the development
scalar until the refinement audit closes.

Topology interval:

- enclosure accepted at 56.939835049 s;
- reopening accepted at 58.663673663 s;
- 53 history checkpoints lie in the enclosed interval;
- all 53 retained a nonzero pocket width, but the paper-defined `w_m` is zero there.

The final state is reopened, so its nonzero 39.46998 nm endpoint is conceptually valid.

## The apparent flicker

The highlighted sequence at steps 1927--1929 is

`40.45280 -> 38.89000 -> 40.47013 nm`.

After removing the local time-linear trend, the middle-row excursion is 1.57270 nm, or 0.315 of a
5 nm cell. It occurs 0.7594 s after the reopening event, not on the event itself. It is the
second-largest one-step excursion in the whole open-state trajectory. The largest is 1.79923 nm,
or 0.360 cell, at 47.33072 s, well before either topology event.

Therefore the evidence does **not** support a topology-specific failure. It supports intermittent
selection of a grid-scale throat. We retain the true minimum and price the excursion; we do not
smooth it away.

## Existing 10 nm contrast

The completed R1.7 10 nm checkpoint at the same parameter pair gives:

| Quantity | 10 nm | R1.8 5 nm development |
| --- | ---: | ---: |
| Definition-correct final `w_m` | 42.84680 nm | 39.46998 nm |
| Largest one-step flicker | 5.82118 nm (0.582 cell) | 1.79923 nm (0.360 cell) |
| Cross-extrusion span at final throat | 0 nm (one interior row) | 5.25119 nm |

The final scalar difference is 3.37682 nm, less than one 5 nm cell, and absolute flicker decreases
with refinement. This is encouraging development evidence, not a refinement proof: the 5 nm path
contains documented operator transitions and a gas-cavity continuation interval, whereas the
10 nm path does not.

The completed initial 0.5 s 10/5 nm audit also supplies the shallow-opening regression that the
original below-substrate connectivity test lacked:

| Grid | Etch depth | Paper-qualified opening | Lowest mask-plane witness |
| --- | ---: | ---: | ---: |
| 5 nm | 6.83362 nm | 87.81703 nm | z = 1.805 um, open |
| 10 nm | 6.82863 nm | 87.78575 nm | z = 1.810 um, open |

At 10 nm the etch depth is smaller than one cell, yet the mask is plainly open. The new throat
criterion classifies both checkpoints correctly. Recomputing the opening from their saved
checkpoints changes the recorded values by less than `3e-8 nm`; no profile rerun is needed.

## Production integration

The pilot now:

- measures mask widths from `material_levelsets[2]`;
- continues to measure feature widths/depth from the combined `geometry.phi`;
- reports paper-qualified `mask_opening_nm`, which is exactly zero only when the exterior gas
  component fails to reach the first resolved interior mask plane inside the declared opening ROI;
- retains `mask_pocket_width_nm` even while sealed;
- reports throat height, per-extrusion-row widths, cross-extrusion span, mask-exit height, and the
  connectivity Boolean on every checkpoint.

This is an observable-only change. It does not alter transport, chemistry, remap, level-set motion,
or topology acceptance.

## Binding next gate

Before a **final authority** endpoint may use opening:

1. Use the integrated mask-layer, throat-qualified metric for every new checkpoint. The existing
   initial 10/5 geometry is recertified offline; it need not be rerun.
2. A future clean 5 nm authority endpoint and its 10 nm comparison must use the same operator and
   observable implementation. The base freeze still requires the 5 nm result within its one-cell
   target gate. No result here opens a held-out profile or relaxes that contract.

The reusable diagnostic is `scripts/krueger_2024_opening_certification.py`; its manufactured open,
sealed-pocket, and flicker contracts are covered by
`tests/test_krueger_2024_opening_certification.py`.
