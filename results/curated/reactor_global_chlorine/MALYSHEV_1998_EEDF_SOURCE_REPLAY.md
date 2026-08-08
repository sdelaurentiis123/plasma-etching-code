# Malyshev 1998 deterministic Cl2 EEPF source replay

## Verdict

This board is a **physical source replay and sensitivity, not a validated
knobs-to-wafer model**. Its collision variant is: **legacy SIGLO molecular Cl2 only**. It solves
the non-Maxwellian two-term EEPF together with particle and power balances and
never selects a coefficient from feature depth. Forward TCP power is not
measured absorbed plasma power, so 30%, 50%, and 70% are all reported rather
than optimized.

| absorbed fraction | source W | E/N Td | ne error | 2/3 mean-E proxy error vs OES | relative-Cl2 proxy error | axial positive-ion flux m-2 s-1 | max closure |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.30 | 300 | 1448.9 | -28.0% | +23.7% | n/a | 2.464e+19 | 4.8e-12 |
| 0.30 | 500 | 1451.9 | -27.3% | +13.4% | -19.2 pp | 4.117e+19 | 2.3e-12 |
| 0.50 | 300 | 1451.9 | +30.8% | +23.9% | n/a | 4.117e+19 | 1.5e-12 |
| 0.50 | 500 | 1459.3 | +29.3% | +13.9% | -23.4 pp | 6.909e+19 | 1.6e-12 |
| 0.70 | 300 | 1456.4 | +91.5% | +24.2% | n/a | 5.788e+19 | 1.0e-12 |
| 0.70 | 500 | 1466.0 | +87.3% | +14.3% | -25.8 pp | 9.730e+19 | 1.6e-13 |

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
- Atomic-Cl ionization, electron detachment from Cl-, and tracked excited-state kinetics are absent.
- The global axial flux is not yet a local wafer flux,
  and it carries no species-resolved IED/IAD.
- Raw collision bytes are not committed. Replay identity is the SHA-256 in the
  JSON receipt.

Consequently this board cannot support feature depth. It exists to locate the
reactor-state failure before any feature coupling is attempted.
