# Malyshev 1998 deterministic Cl2 EEPF source replay

## Verdict

This board is a **physical source replay and sensitivity, not a validated
knobs-to-wafer model**. Its collision variant is: **legacy siglo hamilton plus comsol nist atomic cl**. It solves
the non-Maxwellian two-term EEPF together with particle and power balances and
never selects a coefficient from feature depth. Forward TCP power is not
measured absorbed plasma power, so 30%, 50%, and 70% are all reported rather
than optimized.

| absorbed fraction | source W | E/N Td | ne error | 2/3 mean-E proxy error vs OES | Eq.-11 Cl2 error | within reported Cl2 accuracy | axial positive-ion flux m-2 s-1 | max closure |
|---:|---:|---:|---:|---:|---:|:---:|---:|---:|
| 0.30 | 300 | 995.2 | -19.3% | +10.3% | n/a | n/a | 2.816e+19 | 2.5e-14 |
| 0.30 | 500 | 880.5 | -15.9% | -1.3% | -4.8 pp | PASS | 4.792e+19 | 2.9e-13 |
| 0.50 | 300 | 880.5 | +51.3% | +7.9% | n/a | n/a | 4.792e+19 | 4.7e-13 |
| 0.50 | 500 | 762.8 | +57.9% | -3.2% | -15.1 pp | PASS | 8.376e+19 | 5.6e-13 |
| 0.70 | 300 | 803.1 | +129.1% | +6.4% | n/a | n/a | 6.900e+19 | 7.3e-13 |
| 0.70 | 500 | 686.3 | +139.9% | -4.4% | -23.1 pp | MISS | 1.229e+20 | 2.7e-13 |

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
