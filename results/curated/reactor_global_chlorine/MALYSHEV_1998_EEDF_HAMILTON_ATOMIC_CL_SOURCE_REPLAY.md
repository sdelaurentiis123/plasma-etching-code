# Malyshev 1998 deterministic Cl2 EEPF source replay

## Verdict

This board is a **physical source replay and sensitivity, not a validated
knobs-to-wafer model**. Its collision variant is: **legacy siglo hamilton plus comsol nist atomic cl**. It solves
the non-Maxwellian two-term EEPF together with particle and power balances and
never selects a coefficient from feature depth. Forward TCP power is not
measured absorbed plasma power, so 30%, 50%, and 70% are all reported rather
than optimized.

| absorbed fraction | source W | E/N Td | ne error | 2/3 mean-E proxy error vs OES | relative-Cl2 proxy error | axial positive-ion flux m-2 s-1 | max closure |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.30 | 300 | 1127.9 | -16.6% | +12.1% | n/a | 2.894e+19 | 5.6e-14 |
| 0.30 | 500 | 1091.2 | -15.1% | +1.3% | -16.9 pp | 4.880e+19 | 9.4e-13 |
| 0.50 | 300 | 1091.2 | +52.7% | +10.7% | n/a | 4.880e+19 | 1.0e-13 |
| 0.50 | 500 | 1062.9 | +52.1% | +0.5% | -21.4 pp | 8.249e+19 | 1.2e-12 |
| 0.70 | 300 | 1071.8 | +124.8% | +10.0% | n/a | 6.894e+19 | 6.7e-13 |
| 0.70 | 500 | 1048.0 | +121.2% | +0.1% | -23.9 pp | 1.167e+20 | 2.4e-13 |

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
