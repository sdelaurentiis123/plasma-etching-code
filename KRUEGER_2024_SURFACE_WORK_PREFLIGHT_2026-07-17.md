# Krüger 2024 surface-work preflight

Date: 2026-07-17  
Status: bounded performance decision; no calibration, profile evolution, or held-out reveal

## Question

The fixed-`dx` sparse-volume audit showed that Krüger is surface-work dominated. Before building a
new reduced engine, this preflight asked two narrower questions:

1. Does the already-implemented exact line-extruded diffuse operator accelerate the real Krüger
   surface without changing its physical equations?
2. Can the existing full-3-D operator remove repeated immutable event work without changing any
   flux, yield, material ledger, or checkpoint semantics?

## Exact extrusion result

The zero-motion 10 nm base geometry certified as an exact extrusion: 724 triangles reduced to 181
cross-section segments across two identical periodic strips. Its maximum diffuse balance residual
was `5.73e-15`, versus `1.03e-13` for the existing eight-ray 3-D estimator. The exact construction
was scientifically clean, but its single paired CPU timing was `30.262 s`, versus `28.085 s` for the
current path: a `7.75%` penalty.

This is not evidence that symmetry is invalid. It is evidence that the current Python crossed-string
construction does not earn a production runtime switch on this mesh. The sampled full-3-D Krüger
operator remains unchanged. Full 3-D also remains mandatory for holes, twisting, asymmetric masks,
array edges, and along-line fields.

## Immutable event-work result

The CUDA profile showed that a large face-resolved ion measure was repeatedly reconstructed,
revalidated, and evaluated under the same frozen yield law as it passed through active-surface and
material routing. The accepted repair changes no mathematical operator:

- a face selection validates its complete index map, then reuses the already-certified immutable
  energy, cosine, flux, position, and direction payloads without recomputing physical validity;
- built-in frozen yield laws memoize their immutable per-face rate on that event measure;
- cache state is deliberately excluded from dataclass fields, equality, and serialization.

On a one-million-event manufactured measure, the new subset was bitwise identical and `1.50x`
faster; four identical yield requests were bitwise identical and `4.03x` faster. These are hot-path
microbenchmarks, not whole-engine speedup claims. The repository-wide gate passes: 955 tests passed,
one expected capability test skipped.

## Decision

Keep the unified full-3-D production operator. Do not promote the existing exact diffuse extrusion
backend for runtime. Retain the operator as a deterministic reference and eligibility check. Accept
the immutable event-work optimization, then return to the original validation sequence:

1. paired 5 nm common-refinement versus indexed-remap confirmation;
2. one clean current-epoch fine base anchor at the fixed R1.9 parameter pair;
3. freeze or earn one safeguarded response direction;
4. reveal held-outs once and score experimental/Vienna overlap.

Machine-readable evidence is in `results/krueger_2024_surface_work_preflight/audit.json`.
