# Krüger 2024 surface-normal failure closure

Date: 2026-07-18  
Status: engine defect closed at both preserved real checkpoints; base authority remains incomplete

## Finding

The late 5 nm and 10 nm runs did not fail because diffuse particles were launched inside solid,
because the periodic horizon was too short, or because the etch physics had become undefined. The
particle origins had negative level-set values (gas), and exact float64 rays reached genuine
gas-to-solid crossings. Five target triangles in each mesh had their gas normals reversed.

The old orientation rule sampled the level set one quarter-cell on either side of the triangle.
That is nonlocal at a thin fold: a probe can cross a neighboring interface, place both samples in
the same phase, and choose the wrong side from their small value difference.

The corrected rule keeps the marching-cubes geometric normal but orients it from the exact local
derivative of the trilinear nodal field:

```text
phi > 0 in solid
grad(phi) points locally into solid
gas normal = oriented triangle normal aligned with -grad(phi)
```

An exactly ambiguous critical point refuses instead of inventing an orientation.

## Evidence

| Checkpoint | Old failure | Backward faces | Corrected full visibility | Next physical step |
| --- | ---: | ---: | --- | --- |
| 5 nm, step 203, 12.443006 s | 12 / 36,064 solid-facing rays | 5 / 4,508 | pass; 0 replayed rays | accepted step 204; exact ledger; 0.0556-cell motion |
| 10 nm, step 315, 43.606900 s | 20 / 12,928 solid-facing rays | 5 / 1,616 | pass; one completed long-wrap replay | accepted step 316; exact ledger; 0.2933-cell motion |

The corrected visibility audits took 0.436 s and 2.047 s on CPU. The complete one-step engine
continuations took 483.297 s at 5 nm and 92.351 s at 10 nm. They deliberately wrote no replacement
authority artifacts and read no held-out profile data.

The manufactured trilinear-fold regression passes, and the repository suite is green with
974 passed and one expected skip. Machine-readable evidence is in
`results/curated/krueger_2024_normal_orientation_closure_2026-07-18.json`.

## Consequence

The engine can advance through both formerly fatal geometries without deleting triangles,
softening visibility, dropping particles, changing chemistry, or loosening conservation. This
closes the transport-correctness blocker only. It does **not** turn either partial trajectory into
a calibration endpoint and does not authorize a held-out prediction by itself.

The next expensive action remains one completed base endpoint under one checksum-bound operator
epoch. Its base-only depth/opening result decides whether the fixed R1.9 pair freezes or earns the
single allowed current-epoch correction. No held-out oxygen or power profile should be opened first.
