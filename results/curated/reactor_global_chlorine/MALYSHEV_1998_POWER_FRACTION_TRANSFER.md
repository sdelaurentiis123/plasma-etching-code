# Malyshev 1998 constant power-fraction transfer

## Verdict

One effective source-to-plasma fraction, `0.357165`,
was inferred from the 300 W volume-average electron density only. With that
fraction frozen, the untouched 500 W density error is
`-1.01%`. This closes the power-scaling
residual numerically, but it is **not a direct absorbed-power validation**:
Malyshev reports forward TCP power into the matching network and no density
uncertainty in this article.

| role | source W | absorbed W | ne error | energy-proxy error | Cl2-proxy error | axial ion flux m-2 s-1 |
|---|---:|---:|---:|---:|---:|---:|
| training | 300 | 107.15 | +0.00% | +16.50% | n/a | 3.407e+19 |
| held out | 500 | 178.58 | -1.01% | +7.70% | -18.71 pp | 5.708e+19 |

## Interpretation boundary

- Neither temperature, dissociation, the held-out 500 W condition, nor any
  feature depth selected the fraction.
- The held-out density transfer is descriptive because the source provides no
  uncertainty for that observable here.
- The remaining temperature and Cl2 residuals show that a power fraction alone
  does not close the reactor physics.
- The axial ion flux is still volume-model output, not a wafer boundary; no
  species-resolved IED/IAD exists yet.
- Every prediction/depth support flag remains false.
