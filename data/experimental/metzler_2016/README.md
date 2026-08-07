# Metzler 2016 cyclic C4F8/Ar-ion Si and SiO2 surface board

Primary source: D. Metzler, *High Precision Plasma Etch for Pattern
Transfer: Towards Fluorocarbon Based Atomic Layer Etching*, PhD thesis,
University of Maryland (2016).  The official repository record is
[DRUM b01bb413](https://drum.lib.umd.edu/items/b01bb413-6152-492a-98ec-0abe752ea240).

This is a **retrospective mechanism-validation board**, not a value-blind
target.  The figures were read before this board was frozen.  It can reject
an incorrect surface closure and can grade a closure whose parameters were
fixed independently; it cannot support a claim that these particular values
were prospectively held out.

## What is measured

The experiment separates precursor supply from ion activation:

1. deposit a C4F8-derived fluorocarbon film;
2. deplete the precursor for 12 s;
3. expose the film/substrate stack to low-energy ions from a steady Ar
   plasma.

The reported reactor condition was 200 W ICP power, 10 mTorr, 50 sccm Ar,
and a 10 °C sample.  Biases of approximately -5 to -15 V produced nominal
maximum ion energies of 20--30 eV.

`figures6_5_6_6_cyclic_depth.csv` contains all 42 measured markers from:

- Figure 6.5: Si and SiO2 removal per cycle versus deposited C4F8-film
  thickness at 20, 25, and 30 eV for a 40 s ion step;
- Figure 6.6: Si and SiO2 removal per cycle versus 20, 40, and 60 s ion
  steps at 25 and 30 eV for a nominal 5 Å deposited film.

The source raster has two open Si markers at the same 40 s/30 eV condition
in Figure 6.6.  Both are retained as plotted replicates.  No experimental
error bars or statistical semantics are reported; the CSV's uncertainty is
digitization placement only.

`figure6_9_cycle_averaged_yield.csv` contains all 25 markers from Metzler's
own ion-normalized surface-response plot:

- panel (a): Si and SiO2 removal per incident Ar ion versus fluorine in the
  deposited film per incident Ar ion at 25 eV and a 40 s ion step;
- panel (b): the same ratios at 25 and 30 eV for a nominal 5 Å film while the
  ion step changes from 20 to 60 s.

The direct response is especially useful because both axes share the same
incident-ion denominator.  It constrains the fluorine-supply versus
ion-capacity closure without requiring us to invent an absolute current for
this reactor.  The source does not report error bars or enough details of its
film-to-F conversion to reverse the plotted ratios into an independently
measured wafer flux or film density.  The board therefore preserves the
author-derived ratios as plotted and forbids that reverse inference.

`figures6_14_6_15_xps_cycle_state.csv` contains all 32 markers from the XPS
trajectories at 0, 5, 15, and 40 s during the 25 eV ion step:

- the film F/C ratio derived from its C 1s components;
- `delta F/C`, the source's difference between F/C obtained from the C 1s
  and F 1s/C 1s ratios and its relative proxy for F in the substrate;
- the CFx C 1s and F 1s intensities.

The contrast between 5 and 11 Å films is a model-form gate.  The thin film
rapidly loses C-bound fluorine while `delta F/C` rises above one; the thick
film defluorinates but its substrate proxy remains below roughly 0.06.  A
single well-mixed film reservoir with no finite transfer/mixing depth cannot
reproduce both trajectories.

## Atomic-accounting boundary

The deposited-film quantity is ellipsometric optical thickness.  The thesis
does not publish a time-resolved density and C/F composition that would
convert every plotted Å into absolute C and F areal inventories.
`delta F/C` is an XPS-derived relative proxy, not an absolute number of
substrate F atoms.  These quantities constrain state evolution but must not
be placed directly into an atom ledger.

Figure 6.9 supplies a total-incident-ion-normalized response, but the board
still does not publish species-resolved neutral/ion wafer fluxes or an ion
energy-angle distribution.  It therefore validates a surface closure under
declared cyclic boundary conditions; it does not identify a complete
reactor-to-wafer boundary for Krüger or Yoshie.

## Pixel replay and visual audit

Download the official thesis without renaming it:

```bash
mkdir -p tmp/sources/metzler_2016
curl -L -o tmp/sources/metzler_2016/thesis.pdf \
  https://api.drum.lib.umd.edu/server/api/core/bitstreams/227fdb28-6ea7-4a8d-92ae-5c521f2d1e0b/content
```

Then verify the PDF checksum, render pages 141, 142, 147, 157, and 158 at 240 dpi,
verify the raster checksums and dimensions, reproduce all three CSVs and the
manifest, and optionally create full-resolution PIL overlays:

```bash
python scripts/digitize_metzler_2016_fc_ale.py \
  --overlay-dir tmp/metzler_2016_overlays
```

The replay metadata and strict claim boundary are frozen in
`digitization_manifest.json`.  Source rasters and overlays are not committed.
