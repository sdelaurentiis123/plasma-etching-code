# Malyshev 1998 deterministic Cl2 EEPF source replay

## Verdict

This board is a **physical source replay and sensitivity, not a validated
knobs-to-wafer model**. Its collision variant is: **legacy SIGLO molecular Cl2 plus official COMSOL atomic-Cl momentum and NIST/Hayes measured atomic-Cl ionization**. It solves
the non-Maxwellian two-term EEPF together with particle and power balances and
never selects a coefficient from feature depth. Forward TCP power is not
measured absorbed plasma power, so 30%, 50%, and 70% are all reported rather
than optimized.

| absorbed fraction | source W | E/N Td | ne error | 2/3 mean-E proxy error vs OES | relative-Cl2 proxy error | axial positive-ion flux m-2 s-1 | max closure |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.30 | 300 | 1200.9 | -26.1% | +15.6% | n/a | 2.613e+19 | 1.6e-12 |
| 0.30 | 500 | 1160.9 | -24.8% | +4.6% | -18.4 pp | 4.411e+19 | 8.6e-13 |
| 0.50 | 300 | 1160.9 | +35.4% | +14.2% | n/a | 4.411e+19 | 8.1e-14 |
| 0.50 | 500 | 1130.0 | +35.0% | +3.7% | -22.7 pp | 7.466e+19 | 7.5e-13 |
| 0.70 | 300 | 1139.7 | +99.4% | +13.5% | n/a | 6.237e+19 | 2.1e-12 |
| 0.70 | 500 | 1113.6 | +96.5% | +3.2% | -25.1 pp | 1.057e+20 | 1.8e-12 |

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
- The global axial flux is not yet a local wafer flux,
  and it carries no species-resolved IED/IAD.
- Raw collision bytes are not committed. Replay identity is the SHA-256 in the
  JSON receipt.

Consequently this board cannot support feature depth. It exists to locate the
reactor-state failure before any feature coupling is attempted.
