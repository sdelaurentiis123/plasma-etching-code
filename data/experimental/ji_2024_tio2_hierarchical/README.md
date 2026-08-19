# Ji 2024 same-gas TiO2/Cr morphology-response boards

This directory digitizes Figures 3(b1-b3) and 4(a-f) of Ji et al.,
*Micromachines* **15**,
1160 (2024), DOI `10.3390/mi15091160`. The source is CC BY 4.0. Source pixels
are not committed; the numerical table retains every 600-dpi marker center,
axis calibration, PDF/render hash, and a conservative digitization bound.

With the checksum-matched source PDF, reproduce the render and overlay:

```bash
pdftoppm -f 6 -l 6 -singlefile -png -r 600 \
  /path/to/ji_2024_tio2_hierarchical.pdf /tmp/ji_page6_600
python scripts/digitize_ji_2024_figure3.py \
  --source-pdf /path/to/ji_2024_tio2_hierarchical.pdf \
  --render /tmp/ji_page6_600.png \
  --overlay /tmp/ji_figure3_overlay.png

pdftoppm -f 7 -l 7 -singlefile -png -r 600 \
  /path/to/ji_2024_tio2_hierarchical.pdf /tmp/ji_page7_600
python scripts/digitize_ji_2024_figure4.py \
  --source-pdf /path/to/ji_2024_tio2_hierarchical.pdf \
  --render /tmp/ji_page7_600.png \
  --overlay /tmp/ji_figure4_overlay.png
```

The original-resolution overlay was visually inspected. At fixed `350 W` ICP,
`10 mTorr`, `40 C`, and `40/10/5 sccm SF6/CHF3/O2`, increasing RF power from
`90` to `210 W` changes the upper-triangle height from `143.7` through a
`347.8 nm` maximum at `180 W` to `336.1 nm`; reduces the tip radius from
`99.9` to `9.8 nm`; and narrows the interfeature gap from `96.0` to `18.0 nm`.
The last response is especially discriminating: a removal-only solid cannot
make the gap strictly narrower. A model of this experiment needs an evolving
Cr mask and a retained/deposited passivating surface-volume channel in
addition to ion-assisted removal.

At fixed source/RF power, Figure 4 provides an independent pattern-loading
axis. The printed SEM annotations show that the `350--750 nm` gap cases form a
compact morphology cluster, whereas the `100` and `70 nm` points have shorter
upper regions and steeper upper angles. The paper's prose places the boundary
at `100 nm`, but the digitized `100 nm` point has already shifted; the data
therefore bracket a transition between `100` and `350 nm` and do not identify
a sharp threshold. For the conditional Freddie board, only the inferred
`320 nm`-wide/`80 nm`-gap member is at or below the changed Ji point; the
`280 nm`-wide/`120 nm`-gap member lies in the unsampled interval. Even those
flags are sensitivities, not transferable Oxford thresholds.

This is a topology and response validation board for the `800 nm` EBD-TiO2,
`60 nm` Cr, `200 nm` CD, `300 nm` pitch Ji experiment. The paper does not state
the etch time for the RF-power board and does not publish species-resolved
fluxes, yields, or passivation thickness. These values therefore do not set an
Oxford NPG80 coefficient, spacing threshold, or Freddie profile.
