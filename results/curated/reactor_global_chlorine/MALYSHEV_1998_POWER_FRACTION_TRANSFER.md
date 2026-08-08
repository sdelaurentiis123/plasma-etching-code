# Malyshev 1998 constant power-fraction transfer

## Verdict

One effective source-to-plasma fraction, `0.364647`,
was inferred from the 300 W volume-average electron density only. With that
fraction frozen, the untouched 500 W density error is
`+1.76%`. This closes the power-scaling
residual numerically, but it is **not a direct absorbed-power validation**:
Malyshev reports forward TCP power into the matching network and no density
uncertainty in this article. The independently predicted Eq.-11 Cl2 density is
`-12.28%` relative to the measured
value and **passes**
the source's about +/-25% absolute-density accuracy band.

| role | source W | absorbed W | ne error | energy-proxy error | Eq.-11 Cl2 error | Cl2 accuracy gate | axial ion flux m-2 s-1 |
|---|---:|---:|---:|---:|---:|:---:|---:|
| training | 300 | 109.39 | +0.00% | +14.16% | n/a | n/a | 3.390e+19 |
| held out | 500 | 182.32 | +1.76% | +4.98% | -8.96 pp | PASS | 5.779e+19 |

## Interpretation boundary

- Neither temperature, dissociation, the held-out 500 W condition, nor any
  feature depth selected the fraction.
- The held-out density transfer is descriptive because the source provides no
  uncertainty for that observable here.
- The held-out Eq.-11 Cl2 observable is graded against the source's about
  +/-25% absolute-density accuracy statement; this is not relabeled as sigma.
- The bounded wall ratio/temperature transfer remains sensitivity evidence, so
  this diagnostic pass does not by itself certify a predictive wafer flux.
- The axial ion flux is still volume-model output, not a wafer boundary; no
  species-resolved IED/IAD exists yet.
- Every prediction/depth support flag remains false.
