# Malyshev 1998 deterministic Cl2 EEPF source replay

## Verdict

This board is a **physical source replay and sensitivity, not a validated
knobs-to-wafer model**. Its collision variant is: **legacy siglo hamilton plus comsol nist atomic cl plus isotropic ee**. It solves
the non-Maxwellian two-term EEPF together with particle and power balances and
never selects a coefficient from feature depth. Forward TCP power is not
measured absorbed plasma power, so 30%, 50%, and 70% are all reported rather
than optimized.

| absorbed fraction | source W | E/N Td | ne error | 2/3 mean-E proxy error vs OES | relative-Cl2 proxy error | axial positive-ion flux m-2 s-1 | max closure |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.30 | 300 | 1116.8 | -18.4% | +16.3% | n/a | 2.857e+19 | 6.0e-11 |
| 0.30 | 500 | 1077.9 | -18.5% | +7.2% | -17.1 pp | 4.786e+19 | 1.1e-08 |
| 0.50 | 300 | 1077.9 | +46.7% | +17.2% | n/a | 4.786e+19 | 2.4e-09 |
| 0.50 | 500 | 1047.1 | +42.9% | +8.8% | -21.6 pp | 8.023e+19 | 3.9e-10 |
| 0.70 | 300 | 1056.9 | +112.9% | +18.2% | n/a | 6.725e+19 | 2.1e-09 |
| 0.70 | 500 | 1030.4 | +104.8% | +10.0% | -24.1 pp | 1.128e+20 | 3.0e-10 |

## Use boundary

- The 300 W condition is only a future reactor-diagnostic training candidate;
  500 W is reserved as a held-out reactor diagnostic. No fraction is selected
  by this run.
- The 13.56 MHz electron-heating term uses the Hagelaar--Pitchford
  high-frequency two-term equation and the modern RMS-field convention. At
  this intermediate RF frequency it is a declared quasi-stationary
  local-field sensitivity, not a time-periodic or spatially nonlocal ICP
  heating solution.
- `2/3 <E>` is not the exact OES forward observable. Its error diagnoses EEPF
  shape/chemistry but is not an apples-to-apples temperature validation.
- The Stafford wall regression is extrapolated in Cl/Cl2 from the direct
  `(0.10561, 0.779646)` marker interval to
  `(1e-05, 1.5)`. Every other Stafford domain remains strict.
- Atomic-Cl ionization is included, but electron detachment from Cl- and tracked excited-state kinetics remain absent.
- Isotropic electron-electron Coulomb drift/diffusion is included with the classical Debye logarithm. Electron-ion and anisotropic electron-electron momentum terms remain absent; this is a density-coupling sensitivity, not a complete Coulomb model.
- The global axial flux is not yet a local wafer flux,
  and it carries no species-resolved IED/IAD.
- Raw collision bytes are not committed. Replay identity is the SHA-256 in the
  JSON receipt.

Consequently this board cannot support feature depth. It exists to locate the
reactor-state failure before any feature coupling is attempted.
