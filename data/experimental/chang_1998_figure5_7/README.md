# Chang 1998 Figure 5.7 Cl2+ poly-Si yield

This directory recovers the missing intercept of the Cl2+ poly-Si yield law
plotted in Chang's MIT thesis, Figure 5.7 (PDF and print page 113). Chang
prints the fitted slope, `0.22 Si/ion/sqrt(eV)`, but not the threshold. A
checksum-bound 600-dpi Poppler render and an original-pixel PIL/NumPy replay
recover the x-intercept at `sqrt(E/eV) = 5.099`, or `E_th = 25.999 eV`.

The resulting source law is

```text
Y_Cl2+ = 0.22 * max(sqrt(E_eV) - sqrt(25.999), 0) Si/ion.
```

The plotted Cl2+ points are attributed by Chang to Balooch et al. and were
measured with an approximately `1e-4 Torr` chlorine background. Chang
explicitly interprets the elevated yield as a highly chlorinated surface.
Accordingly this law supports Cl2+ removal on highly chlorinated poly-Si; it
does not support bare-Si sputtering, incidence-angle response, oxygenated or
passivated surfaces, or extrapolation outside the plotted 26-625 eV support.

Reproduce and checksum-check the digitization with:

```bash
python scripts/digitize_chang_1998_figure5_7.py --check \
  --overlay /tmp/chang_figure5_7_overlay.png
```

The source thesis PDF is already version-controlled at
`research_sources/chang_thesis.pdf`; the derived raster is not duplicated.
