# Wise 1996 direct Boltzmann/electron-pressure closure audit

This is a direct check of the bulk closure using three independently measured Figure-3 fields. Potential is aligned at the axis because its gauge is arbitrary; no slope, temperature, reactor, or feature parameter is fitted.

- measured axis-to-outer potential drop: `-6.383 V`
- local-Te pressure-gradient reconstruction: `-7.680 V`
- local-Te unweighted MAPE, excluding the gauge point: `5.07%`
- constant-axis-Te sensitivity MAPE: `3.82%`
- digitized density FWHM: `7.53 cm`
- independent Miller typical FWHM interval: `7.0--9.0 cm` (inside: `True`)

The closure is supported at the few-percent unweighted level, but the paper prints no marker error bars, so this is not labeled a formal uncertainty-weighted pass. It validates neither a knobs-to-state solve nor wafer flux or etch depth.
