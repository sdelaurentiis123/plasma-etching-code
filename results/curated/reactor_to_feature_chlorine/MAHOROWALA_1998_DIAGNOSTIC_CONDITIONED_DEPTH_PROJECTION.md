# Mahorowala 1998 diagnostic-conditioned depth projection

This is a reactor-to-surface-plane failure-localization run. No etch rate or depth conditioned the reactor, wafer normalization, or surface law. It is not yet a formal feature-depth prediction.

- Cl+ only surface-plane MAPE: `40.58%`
- measured species-resolved Cl+/Cl2+ surface-plane MAPE: `21.81%`
- power-closed RF-sheath species-resolved surface-plane MAPE: `21.72%`
- one-step joined Cl+/Cl2+ surface-plane MAPE: `21.72%`
- SiCl2 feedback, reflective-wall limit MAPE: `18.17%`
- SiCl2 feedback, reactive-wall limit MAPE: `21.06%`
- full Table-IV SiClx network, reflective-wall limit MAPE: `15.02%`
- full Table-IV SiClx network, reactive-wall limit MAPE: `15.39%`
- all ions treated as Cl+ sensitivity MAPE: `83.06%`
- points inside the measured 35--100 eV surface domain: `6/11`
- formal held-out feature-depth passes: `0`

| run | source W | bias W | flow sccm | observed nm | RF species nm | reflective-product nm | reactive-product nm | Cl+ IEAD mean+/-sd eV | Cl2+ IEAD mean+/-sd eV | Cl/ion |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 400 | 80 | 100 | 393.8 | 365.9 | 356.9 | 356.9 | 113.9+/-23.7 | 124.0+/-16.8 | 100.0 |
| 2 | 550 | 80 | 175 | 459.4 | 421.1 | 413.7 | 408.5 | 96.3+/-21.6 | 93.0+/-16.2 | 81.4 |
| 3 | 250 | 80 | 25 | 250.0 | 308.8 | 289.1 | 303.1 | 173.9+/-23.5 | 178.1+/-17.8 | 131.4 |
| 4 | 400 | 20 | 175 | 175.0 | 133.4 | 131.7 | 130.5 | 47.8+/-7.6 | 46.1+/-5.6 | 97.4 |
| 5 | 550 | 20 | 100 | 162.5 | 149.7 | 145.8 | 145.6 | 40.3+/-7.2 | 38.6+/-4.7 | 83.0 |
| 6 | 550 | 80 | 25 | 296.9 | 437.1 | 394.4 | 423.2 | 95.6+/-21.7 | 90.0+/-16.3 | 84.6 |
| 7 | 250 | 80 | 175 | 265.6 | 300.0 | 296.6 | 294.7 | 176.0+/-23.8 | 181.5+/-17.4 | 121.7 |
| 8 | 550 | 140 | 100 | n/a | 607.9 | 589.4 | 588.8 | 137.0+/-34.2 | 151.1+/-24.7 | 83.0 |
| 9 | 400 | 140 | 25 | 362.5 | 547.6 | 501.8 | 533.3 | 204.8+/-36.3 | 189.6+/-26.7 | 102.7 |
| 10 | 400 | 20 | 25 | 112.5 | 138.6 | 128.0 | 135.4 | 47.3+/-7.6 | 45.3+/-5.6 | 102.7 |
| 11 | 250 | 20 | 100 | 134.4 | 115.7 | 113.7 | 113.7 | 59.8+/-8.4 | 60.9+/-5.7 | 126.5 |
| 12 | 400 | 140 | 175 | n/a | 532.6 | 524.6 | 519.5 | 212.6+/-36.9 | 193.2+/-27.3 | 97.4 |
| 13 | 250 | 140 | 100 | 362.5 | 435.3 | 426.9 | 427.2 | 296.7+/-38.5 | 297.0+/-27.6 | 126.5 |

## Exact blockers to a formal depth grade

- measured per-run or validated bias-power-to-IEAD transfer
- measured Cl2+ incidence-angle response for feature sidewalls
- deterministic in-feature ion/radical transport and evolving geometry
- uncertainties on the center ion current, radical flux, and plasma potential
- same-tool reactor geometry and absorbed-power diagnostics
- self-consistent SiClx collision-power and chlorine-source feedback
- measured chamber-wall Si/Cl coverage state between Lee limits
- surface reaction probabilities for SiClx+ product-ion redeposition
