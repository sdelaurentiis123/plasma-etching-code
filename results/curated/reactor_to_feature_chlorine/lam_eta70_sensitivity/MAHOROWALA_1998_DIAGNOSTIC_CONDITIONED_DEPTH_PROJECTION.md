# Mahorowala 1998 diagnostic-conditioned depth projection

This is a reactor-to-surface-plane failure-localization run. No etch rate or depth conditioned the reactor, wafer normalization, or surface law. It is not yet a formal feature-depth prediction.

- Cl+ only surface-plane MAPE: `52.65%`
- measured species-resolved Cl+/Cl2+ surface-plane MAPE: `22.43%`
- power-closed RF-sheath species-resolved surface-plane MAPE: `20.70%`
- one-step joined Cl+/Cl2+ surface-plane MAPE: `20.70%`
- SiCl2 feedback, reflective-wall limit MAPE: `20.90%`
- SiCl2 feedback, reactive-wall limit MAPE: `21.65%`
- all ions treated as Cl+ sensitivity MAPE: `46.79%`
- points inside the measured 35--100 eV surface domain: `6/11`
- formal held-out feature-depth passes: `0`

| run | source W | bias W | flow sccm | observed nm | RF species nm | reflective-product nm | reactive-product nm | Cl+ IEAD mean+/-sd eV | Cl2+ IEAD mean+/-sd eV | Cl/ion |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 400 | 80 | 100 | 393.8 | 286.9 | 278.2 | 278.2 | 85.5+/-17.6 | 93.4+/-13.0 | 100.0 |
| 2 | 550 | 80 | 175 | 459.4 | 332.0 | 324.7 | 319.7 | 76.8+/-16.7 | 70.1+/-12.4 | 81.4 |
| 3 | 250 | 80 | 25 | 250.0 | 246.1 | 227.1 | 240.5 | 134.3+/-19.1 | 129.1+/-13.3 | 131.4 |
| 4 | 400 | 20 | 175 | 175.0 | 98.0 | 96.5 | 95.5 | 40.1+/-5.7 | 38.5+/-4.0 | 97.4 |
| 5 | 550 | 20 | 100 | 162.5 | 108.5 | 105.1 | 105.0 | 34.3+/-5.3 | 33.5+/-3.5 | 83.0 |
| 6 | 550 | 80 | 25 | 296.9 | 344.3 | 304.1 | 330.7 | 75.6+/-16.8 | 68.2+/-12.2 | 84.6 |
| 7 | 250 | 80 | 175 | 265.6 | 238.6 | 235.2 | 233.3 | 135.0+/-18.7 | 131.9+/-13.3 | 121.7 |
| 8 | 550 | 140 | 100 | n/a | 487.5 | 469.1 | 468.5 | 108.9+/-25.3 | 109.1+/-18.7 | 83.0 |
| 9 | 400 | 140 | 25 | 362.5 | 427.7 | 385.3 | 413.9 | 137.2+/-27.5 | 144.1+/-19.1 | 102.7 |
| 10 | 400 | 20 | 25 | 112.5 | 101.9 | 92.7 | 99.1 | 39.7+/-5.7 | 38.0+/-4.0 | 102.7 |
| 11 | 250 | 20 | 100 | 134.4 | 85.6 | 83.8 | 83.9 | 47.3+/-6.1 | 49.3+/-4.4 | 126.5 |
| 12 | 400 | 140 | 175 | n/a | 416.9 | 409.1 | 404.2 | 143.5+/-27.4 | 145.9+/-19.0 | 97.4 |
| 13 | 250 | 140 | 100 | 362.5 | 346.7 | 338.4 | 338.8 | 210.0+/-30.2 | 215.5+/-20.7 | 126.5 |

## Exact blockers to a formal depth grade

- measured per-run or validated bias-power-to-IEAD transfer
- measured Cl2+ incidence-angle response for feature sidewalls
- deterministic in-feature ion/radical transport and evolving geometry
- uncertainties on the center ion current, radical flux, and plasma potential
- same-tool reactor geometry and absorbed-power diagnostics
- self-consistent SiClx collision-power and chlorine-source feedback
- measured chamber-wall Si/Cl coverage state between Lee limits
- surface reaction probabilities for SiClx+ product-ion redeposition
