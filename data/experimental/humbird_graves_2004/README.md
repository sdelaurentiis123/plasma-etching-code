# Humbird--Graves fluorocarbon/silicon MD surface-state board

This board recovers transient surface inventories and silicon etch yields
shown by David Humbird and David Graves in the UC Berkeley FLCC research
seminar *Molecular Dynamics Simulations of Plasma-Surface Interactions and
Etching* (23 February 2004).  It is primary-author evidence, but it is not a
peer-reviewed archival figure.  The related peer-reviewed mechanism papers
are DOI
[10.1063/1.1644338](https://doi.org/10.1063/1.1644338),
[10.1063/1.1736321](https://doi.org/10.1063/1.1736321), and
[10.1063/1.1769602](https://doi.org/10.1063/1.1769602).

The board is useful because the longer seminar trajectories expose states
that an absolute-depth closure must reproduce:

- at 20 eV, CF2/Ar+ = 9/1 produces continuing C/F deposition and almost no
  silicon removal;
- at 200 eV, the same boundary accumulates a thick carbon-rich inventory,
  while the silicon yield decays strongly with fluence;
- adding one or two thermal F atoms per Ar+ makes the layer more permeable
  and sustains silicon removal;
- the surface contains a finite-residence Si--C transport layer and a deeper
  Si--F reaction front.  Silicon is not created and volatilized in one
  bookkeeping step.

These are classical molecular-dynamics constraints, not atomic truth:
potential error, finite-cell effects, accelerated fluence, and missing
ensemble uncertainty all remain.  The curves grade topology and transient
surface response.  They do not supply reactor fluxes or validate oxide
chemistry.

## Reproduce the digitization

Download the two source rasters without renaming or committing them:

```bash
curl -L \
  'https://image1.slideserve.com/3215203/thermal-cf-2-ar-9-1-si-impacts-l.jpg' \
  -o /private/tmp/humbird_graves_slide20.jpg
curl -L \
  'https://image1.slideserve.com/3215203/thermal-f-cf-2-200-ev-ar-l.jpg' \
  -o /private/tmp/humbird_graves_slide21.jpg
```

Then replay source hashes, dimensions, axis transforms, every selected
full-resolution trace pixel, the committed CSV, and the JSON manifest:

```bash
python scripts/digitize_humbird_graves_2004_seminar.py \
  --slide-20 /private/tmp/humbird_graves_slide20.jpg \
  --slide-21 /private/tmp/humbird_graves_slide21.jpg \
  --check
```

For visual QA, add
`--overlay-dir tmp/humbird_graves_2004_digitization`.  The overlays are
deliberately untracked; the copyrighted source rasters are not redistributed.

The CSV retains pixel centers and conservative digitization allowances.  They
cover raster reading, not the source MD model discrepancy.
