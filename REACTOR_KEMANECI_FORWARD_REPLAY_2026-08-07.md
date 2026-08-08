# Kemaneci 2014 chlorine forward-reaction replay - 2026-08-07

## Result

`build_kemaneci_2014_forward_chlorine_network()` is a conservation-audited
native replay of the 36 non-elastic forward volume reactions printed in
Kemaneci et al. Table 4.

The replay contains ten unique heavy states plus explicit electrons:

- ground `Cl2` (the paper's `v=0`) and `Cl2(v=1--3)`;
- ground `Cl(2P3/2)` represented as `Cl`, plus `Cl(2P1/2)` and `Cl(1P5/2)`;
- `Cl2+`, `Cl+`, and `Cl-`;
- `e`.

Every reaction is exactly atom and charge conserving. Independent tests
evaluate every printed electron-fit row at a fixed temperature, every
gas-temperature row at a fixed heavy-particle temperature, and reject electron
temperatures outside the paper's declared `0.5--10 eV` fit domain.

## Source anomaly resolved explicitly

Kemaneci Table 4 labels the charge-exchange group `(28)--(32)` but prints
`Cl2(v=0--3) + Cl+`, which expands to four reactions. Table 2 also prints both
unqualified `Cl2` and `Cl2(v=0,1,2,3)` in a way that can be misread as duplicate
ground states.

The official COMSOL reproduction resolves the intended topology:

- unqualified `Cl2` is the ground vibrational state;
- there are four charge-exchange reactions, for ground `Cl2` and `v=1--3`;
- 38 forward volume entries include those four plus two elastic channels;
- six reverse excitation reactions obtained by detailed balance bring the
  complete implementation to 44.

The native replay therefore retains source labels 28--31, leaves label 32
absent, and records the printed range as an off-by-one source defect. It does
not invent a fifth molecular state or reaction.

## Deliberate boundary

This rung omits:

1. the two elastic electron channels, because Table 4 points to cross-section
   data rather than printing it;
2. the six reverse vibrational/atomic excitation reactions, because the paper
   says they follow detailed balance but does not print their coefficients;
3. electron-event energies, because rate-fit exponents are not event energies
   and the paper's Figure-10 level convention is not thermochemically complete;
4. wall, flow, charged transport, gas heating, and absorbed-power closures.

Accordingly, the network is a source-transcription and conservation-verification
rung. `has_complete_electron_energy_ledger` is false and an electron-power
request fails closed. It does not support a reactor prediction, wafer flux, or
feature depth.

## Next gates

1. Recover the exact elastic cross-section files used by the official
   implementation and compare them with the primary Griffin/Gregorio sources.
2. Implement detailed-balance reverse rates from explicit degeneracies and
   level energies, then reproduce the 44-reaction source model.
3. Add a separate physical/evaluated energy ledger; do not reuse Table-4 fit
   exponents as thresholds.
4. Preregister published density/temperature boards before solving or grading
   the source reproduction.
