# Krueger 5 nm authority run: topology refusal

## Status

This run is **not** an authority result and must not be resumed. It used engine revision
`a15634f` and stopped safely after step 168 at 11.376815 s physical etch time. The last recorded
geometry had 155.697 nm depth and 76.709 nm opening. Charge and material ledgers remained closed,
and no incomplete trajectory was accepted.

## Refusal

The hard-visibility certification found one of 36,856 rays reaching a solid-facing intersection
after exhausting the 1,024-wrap exact replay budget. Exact checkpoint replay showed that the ray
did not suffer a floating-point or horizon failure: it passed through an unmatched interior edge
in the marching-cubes surface, travelled inside solid, and hit the far side of the feature.

## Root cause

`extract_mesh_3d` removed all triangles below an arbitrary float32-scaled area floor. At this
checkpoint, 29 real positive-area sliver triangles fell below that threshold. Removing them
converted a watertight 4,636-face scalar-field contour into a 4,607-face mesh with 56 unmatched
interior edges. Keeping all positive-area faces restored zero unmatched interior edges. The
failed ray then hit the correct near-side face immediately with positive incoming cosine and no
replay.

An independent duplicate-facet issue at exact multi-material grid crossings was also exposed by
the new topology certification. Exact opposite-winding duplicates now cancel as an oriented
surface pair; same-winding duplicates reduce deterministically to one face.

## Repair and evidence

Revision `418fb28`:

- retains every representable positive-area triangle;
- disables degenerate marching-cubes output at extraction;
- cancels indexed duplicate facets by orientation;
- refuses unmatched interior or non-manifold edges while permitting declared domain boundaries;
- applies the same topology contract to CPU and GPU extraction.

Verification at the repair revision: 53 focused tests passed, followed by the full suite with
968 passed and 1 skipped. The repaired production ray hits the near-side surface with zero wraps.

## Restart rule

Because the surface operator changed, the authority campaign restarts from zero on `418fb28` or a
descendant containing it. The checkpoint in this directory is retained only for forensic replay.
