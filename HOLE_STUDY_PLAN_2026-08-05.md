# 200:1 HAR hole study — assembly plan (2026-08-05)

A profile study of an extreme-aspect-ratio contact/memory hole etched from an
imported CAD geometry, assembled from pieces already built and gated. This
document is the run spec and the honesty contract; it introduces no new physics.

## What already exists (and its gate)

| piece | module | gate |
|---|---|---|
| STL import → signed distance field | `src/petch/stl_import.py` | SDF exact to the file's own tessellation (≤ 1e-12 beyond faceting); ASCII/binary bitwise agreement; non-watertight rejected (13 gates, `02e4a62`) |
| axisymmetry detection → profile | `stl_import.extract_axisymmetric_profile` | cylinder 6e-4 out-of-roundness accepted, prism 0.16 rejected; frustum r(z) to the faceting bound |
| body-of-revolution transport | axisymmetric operator | Clausing benchmark **1.295 % @ 100:1, 0.656 % @ 200:1** |
| two-component ion beam | `src/petch/iadf_two_component.py` | 13 gates; measured core/tail widths reproduced (Kim 2025), Krüger Fig-4 round trip inside the digitisation band; analytic cone acceptance to 1e-10 (`b79e968`) |
| surface chemistry | `chemistry_deck.py` + mixed layer | Deck 1 bitwise-equivalent to the validated construction; ledgers < 1e-9 |
| metrics | `scripts/regrade_neck_metrics.py` | top CD / neck CD / neck depth, synthetic recovery to grid tolerance |

## The physics that makes this study worth running

At AR 200 the acceptance half-angle of the hole is `arctan(1/200) = 0.286°`.
The measured ion angular distribution is two-component — a thermal core plus a
sheath-collision tail — and the split decides the answer
(`src/petch/iadf_two_component.py`, table reproduced from its gates):

| beam model | cone acceptance @ AR 200 | sidewall share |
|---|---|---|
| core only (T⊥ = 0.044 eV, measured core) | 0.860 | 0.140 |
| two-component, tail fraction 0.50 | 0.501 | 0.499 |
| two-component, tail fraction 0.65 | 0.393 | 0.607 |

Adding the measured tail multiplies AR-200 sidewall delivery by **3.6–4.3×**.
A study that ignores it is wrong at the bottom of the hole by that factor. The
tail *fraction* is the one declared open parameter (`[VERIFY]`; deriving it from
a sheath collision operator is the S2 rung of
`RESEARCH_IADF_SUBDEGREE_AND_REACTOR_2026-07-29.md`).

## Run spec

**Geometry.** Import the supplied STL; run `diagnose_mesh` (watertight,
consistently oriented) and `extract_axisymmetric_profile`. Route on the measured
deviation: ≤ 1 % out-of-roundness → axisymmetric operator (the Clausing-gated
path); otherwise full 3-D. Record the deviation in the deliverable — the routing
decision is a receipt, not an assumption.

**Beam.** `build_two_component_boundary` at the declared reference set. Sweep
the tail fraction over {0.0, 0.35, 0.50, 0.65} as the *headline sensitivity* —
the deliverable reports a band, not a point, because that is the honest state of
the parameter.

**Chemistry.** Deck 1 (fluorocarbon/SiO₂) — the arm carrying the validation
record. Report the validated envelope explicitly (below).

**Series.**
1. Transport-only reference: Clausing transmission vs AR at 50/100/150/200 —
   reproduces the published benchmark inside the study itself.
2. Static delivery vs depth: ion, hot-neutral, and thermal flux to the etch
   front vs AR at the imported geometry, with and without the tail.
3. Profile evolution to the target depth, per tail fraction.

**Reported observables.** Bottom CD, top CD, neck CD and its depth, taper angle,
depth vs time, and the delivery curves from series 2. Metrics via
`regrade_neck_metrics` so the definitions match the community convention
(top/neck/z_neck), not a bespoke aperture number.

## Honesty appendix (ships with the study)

- **Validated:** transport geometry at 200:1 against the analytic Clausing
  solution (0.656 %); beam width against a 0.1°-resolution measurement; surface
  chemistry against a published trench experiment — mask survival exact, mouth
  equilibrium at 50.9 nm vs 45 nm target, late-time closure ratio inside the
  reference band.
- **Extrapolated:** surface chemistry beyond AR ≈ 17, where no matched
  profile data exists; the depth channel currently runs +29 % on the validation
  case (`RESULTS_CASCADE_FUNNELLING_2026-08-05.md` locates it in the ion-energy
  channel magnitude, which the source mechanism's author flags in his own
  thesis).
- **Not modelled:** in-feature charging in this pipeline. For a conventional
  plasma tool this bounds deep-AR claims; the preregistered gate ladder is
  `RESEARCH_CHARGING_DEEP_AR_VALIDATION_2026-08-04.md` (hard ceiling: in-feature
  potential never exceeds the maximum ion energy and is nearly AR-independent
  above AR ≈ 17). A neutral-beam source is largely exempt from this term, which
  is the regime where the validated subset covers the dominant mechanisms.
- **Declared open:** tail fraction (swept, not fitted); crosslink probability
  (unpublished); Chang class-2 roll-off interpolation (our form within 0.065
  absolute of the digitised curve).

## Cost

Series 1–2 are frozen-geometry: minutes each, CPU. Series 3 is the only box
item — one run per tail fraction. Preflight locally with
`scripts/preflight_krueger.sh` first; deploy by clean `git archive` (never patch
a live tree — three false relaunches on 2026-08-04).
