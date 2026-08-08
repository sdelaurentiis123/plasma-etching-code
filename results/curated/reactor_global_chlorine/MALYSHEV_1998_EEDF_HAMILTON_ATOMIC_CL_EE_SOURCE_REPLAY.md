# Malyshev 1998 deterministic Cl2 EEPF source replay

## Verdict

This board is a **physical source replay and sensitivity, not a validated
knobs-to-wafer model**. Its collision variant is: **legacy siglo hamilton plus comsol nist atomic cl plus isotropic ee**. It solves
the non-Maxwellian two-term EEPF together with particle and power balances and
never selects a coefficient from feature depth. Forward TCP power is not
measured absorbed plasma power, so 30%, 50%, and 70% are all reported rather
than optimized.

| absorbed fraction | source W | E/N Td | ne error | 2/3 mean-E proxy error vs OES | Eq.-11 Cl2 error | within reported Cl2 accuracy | axial positive-ion flux m-2 s-1 | max closure |
|---:|---:|---:|---:|---:|---:|:---:|---:|---:|
| 0.30 | 300 | 983.3 | -20.8% | +14.3% | n/a | n/a | 2.781e+19 | 1.6e-10 |
| 0.30 | 500 | 865.7 | -19.2% | +4.6% | -5.2 pp | PASS | 4.701e+19 | 1.3e-09 |
| 0.50 | 300 | 865.7 | +45.5% | +14.3% | n/a | n/a | 4.701e+19 | 3.9e-10 |
| 0.50 | 500 | 743.2 | +47.6% | +5.9% | -15.9 pp | PASS | 8.141e+19 | 6.7e-11 |
| 0.70 | 300 | 785.4 | +116.3% | +15.1% | n/a | n/a | 6.728e+19 | 1.4e-10 |
| 0.70 | 500 | 662.2 | +119.2% | +7.5% | -24.2 pp | MISS | 1.188e+20 | 1.2e-10 |

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
- Isotropic electron-electron Coulomb drift/diffusion is included with the classical Debye logarithm. Electron-ion and anisotropic electron-electron momentum terms remain absent; this is a density-coupling sensitivity, not a complete Coulomb model.
- The global axial flux is not yet a local wafer flux,
  and it carries no species-resolved IED/IAD.
- Raw collision bytes are not committed. Replay identity is the SHA-256 in the
  JSON receipt.

Consequently this board cannot support feature depth. It exists to locate the
reactor-state failure before any feature coupling is attempted.
