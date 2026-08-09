# Mahorowala 1998 photo-assisted depth identifiability

## Verdict

The residual can be written as a 105-nm-equivalent photon target, but broadband amplitude overlap is not spectral closure. A subsequent line-resolved Cl I audit finds negligible atomic emission in Du's 104.82--106.67 nm response band and strong unmeasured 109--120 nm lines. This scalar audit therefore grants no photo-assisted depth closure.

The conserved no-PAE board has 15.52% MAPE. Expressing its positive residuals in units of Du's independently measured 90--244 Si/photon response gives a 106-nm-equivalent flux spanning 4.23e+11--7.75e+13 cm^-2 s^-1. That is a target-variable identity, not a fitted validation score or a spectrally supported source.

| run | source W | bias W | flow sccm | observed nm | baseline nm | missing nm | required 106 nm flux, 244 yield | required 106 nm flux, 90 yield | max beta |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 400 | 80 | 100 | 393.8 | 291.0 | 102.8 | 2.81e+13 | 7.61e+13 | 0.0061 |
| 2 | 550 | 80 | 175 | 459.4 | 354.7 | 104.7 | 2.86e+13 | 7.75e+13 | 0.0044 |
| 3 | 250 | 80 | 25 | 250.0 | 208.3 | 41.7 | 1.14e+13 | 3.09e+13 | 0.0041 |
| 4 | 400 | 20 | 175 | 175.0 | 117.5 | 57.5 | 1.57e+13 | 4.26e+13 | 0.0032 |
| 5 | 550 | 20 | 100 | 162.5 | 124.4 | 38.1 | 1.04e+13 | 2.82e+13 | 0.0015 |
| 6 | 550 | 80 | 25 | 296.9 | 289.1 | 7.8 | 2.12e+12 | 5.75e+12 | 0.0004 |
| 7 | 250 | 80 | 175 | 265.6 | 264.1 | 1.5 | 4.23e+11 | 1.15e+12 | 0.0001 |
| 9 | 400 | 140 | 25 | 362.5 | 341.4 | 21.1 | 5.77e+12 | 1.56e+13 | 0.0014 |
| 10 | 400 | 20 | 25 | 112.5 | 102.4 | 10.1 | 2.77e+12 | 7.51e+12 | 0.0006 |
| 11 | 250 | 20 | 100 | 134.4 | 98.6 | 35.7 | 9.76e+12 | 2.65e+13 | 0.0031 |
| 13 | 250 | 140 | 100 | 362.5 | 347.4 | 15.1 | 4.13e+12 | 1.12e+13 | 0.0014 |

## Guardrail

Kemaneci reaction 18 produces a 139-nm excitation proxy. Du's absolute surface yield is restricted to 104.82--106.67 nm, so the code does not apply that yield at 139 nm. The reported proxy fractions combine an unknown shortwave branch, radiative escape, and spectral surface response; they are target requirements, not calibrated reactor constants. The line-resolved follow-up also forbids treating all sub-120-nm photons as though they carried Du's 105-nm response.

The independent wavelength result is also carried as a negative bound: using 10 Si/photon as a nonpredictive 130--139 nm analog requires as much as 4.67x the transparent one-pass Kemaneci reaction-18 proxy. Thus 139-nm light alone is not an amplitude closure for all runs under that response.

## Exact experiment that closes the board

Measure the absolute sub-120-nm spectrum at the wafer synchronously with the species-resolved ion flux/IED under each Table-2.2 condition, then measure the RF-waveform-resolved Si photo-etch response with ions optically blocked/unblocked. Those independent boundaries make the depth board out-of-sample. No feature depth may be used to choose their normalization.
