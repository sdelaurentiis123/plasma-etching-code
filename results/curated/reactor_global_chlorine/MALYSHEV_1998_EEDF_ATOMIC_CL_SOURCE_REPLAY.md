# Malyshev 1998 deterministic Cl2 EEPF source replay

## Verdict

This board is a **physical source replay and sensitivity, not a validated
knobs-to-wafer model**. Its collision variant is: **legacy siglo plus comsol nist atomic cl**. It solves
the non-Maxwellian two-term EEPF together with particle and power balances and
never selects a coefficient from feature depth. Forward TCP power is not
measured absorbed plasma power, so 30%, 50%, and 70% are all reported rather
than optimized.

| absorbed fraction | source W | E/N Td | ne error | 2/3 mean-E proxy error vs OES | Eq.-11 Cl2 error | within reported Cl2 accuracy | axial positive-ion flux m-2 s-1 | max closure |
|---:|---:|---:|---:|---:|---:|:---:|---:|---:|
| 0.30 | 300 | 1035.5 | -28.6% | +13.7% | n/a | n/a | 2.538e+19 | 6.9e-13 |
| 0.30 | 500 | 905.8 | -24.8% | +1.8% | -7.8 pp | PASS | 4.361e+19 | 1.3e-12 |
| 0.50 | 300 | 905.8 | +35.4% | +11.2% | n/a | n/a | 4.361e+19 | 1.3e-12 |
| 0.50 | 500 | 772.4 | +43.9% | -0.4% | -19.1 pp | MISS | 7.754e+19 | 2.9e-13 |
| 0.70 | 300 | 818.0 | +107.2% | +9.6% | n/a | n/a | 6.344e+19 | 5.0e-13 |
| 0.70 | 500 | 686.8 | +122.0% | -1.9% | -27.7 pp | MISS | 1.156e+20 | 1.7e-12 |

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
- Malyshev's about +/-25% Cl2 absolute-density accuracy is used as a reported
  accuracy band, not a statistical sigma; digitization error remains separate.
- The bounded Stafford coverage fit is trained only on the direct Cl/Cl2
  interval `(0.10561, 0.779646)`. Its transfer to ratio domain
  `(1e-05, 30.0)` and from 300 to 333 K is sensitivity evidence;
  pressure, power, material, and the direct-marker provenance remain explicit.
- Atomic-Cl ionization is included, but electron detachment from Cl- and tracked excited-state kinetics remain absent.
- Electron-electron and electron-ion Coulomb collisions are absent.
- The global axial flux is not yet a local wafer flux,
  and it carries no species-resolved IED/IAD.
- Raw collision bytes are not committed. Replay identity is the SHA-256 in the
  JSON receipt.

Consequently this board cannot support feature depth. It exists to locate the
reactor-state failure before any feature coupling is attempted.
