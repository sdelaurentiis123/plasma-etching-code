# Native chlorine electron-energy ledger - 2026-08-07

## Result

The reactor network now distinguishes two physically different electron-loss
objects:

1. a fixed energy transferred per reaction event, multiplied by the event
   rate; and
2. the incident-electron kinetic-energy moment for a particle-removing
   collision, evaluated as `<sigma v E>` from the **same** cross-section table
   that supplies `<sigma v>`.

This is an infrastructure gate, not a completed chlorine power solve. The
NIST attachment table remains support-limited over every measured Malyshev
electron temperature, so the native Lam model still fails closed before
absorbed power can determine `T_e`.

## Why attachment is a separate moment

For a Maxwellian electron population and an electron-target collision,

`R = n_e n_t <sigma v>`

is the particle event rate. If the collision removes the incident free
electron, its kinetic-energy removal rate is

`L_E = n_e n_t <sigma v E>`.

The ratio `<sigma v E>/<sigma v>` is the mean incident energy conditioned on
that reaction occurring. It is not generally `3/2 T_e`: the collision rate is
weighted by both speed and the energy-dependent cross section. For a constant
cross section on complete Maxwellian support, the exact ratio is `2 T_e`.

This term also cannot be replaced by an Arrhenius fit exponent. A fit exponent
describes the temperature dependence of an integrated rate over its fit
domain; it is not an event energy or an energy moment.

The free-electron energy density `3/2 n_e T_e` is an independent Boltzmann
moment. Evolving it alongside electron number does not make collision-
selected kinetic-energy removal appear automatically. A temperature-form
equation can move number-source terms through the chain rule, but it still
requires a declared closure for the energy carried by the created or removed
electrons.

## Source audit

`kemaneci-2014-chlorine-global` supplies a useful detailed source-model
topology:

- chemical, charged-wall, and elastic electron losses are separate;
- the chemistry includes vibrationally excited Cl2 and two atomic Cl excited
  states;
- all printed electron-rate fits are bounded to `0.5--10 eV`.

It does not supply a fundamental thermochemical ledger. Its 500-dpi-audited
Figure 10 places both ground-state Cl2 and ground-state Cl at zero, so the
species-level difference in its Equation 17 cannot reproduce the NIST
`D0(Cl2)=2.4793 eV`. Its ground-state dissociation rate fit contains an
`8.84 eV` exponential parameter, independently demonstrating why fit
parameters and physical event energies must remain separate.

The implementation therefore preserves Kemaneci as a future
source-reproduction tier and uses NIST/evaluated thresholds and collision
moments for the physical tier.

## Executable guarantees

`ElectronMaxwellianCrossSectionRateCoefficient` now evaluates both the rate
and incident-energy moments analytically over each piecewise-linear cross-
section segment. It applies independent high-energy support gates because the
energy kernel `E^2 exp(-E/T_e)` has a materially heavier tail than the rate
kernel `E exp(-E/T_e)`.

`Reaction` accepts the incident-energy moment only when all of the following
hold:

- the reaction removes exactly one incident electron;
- the electron kinetic order is exactly one;
- the rate is driven by the same Maxwellian tabulated-cross-section object;
- no competing fixed energy-per-event value is present.

The network electron-power sum fails closed if any reaction lacks either a
fixed event energy or an allowed energy moment.

The audit also repaired a threshold-node defect: a cross section printed at
exactly the physical threshold is now retained, while a zero node is inserted
only when the threshold lies in an unsampled interval. Tightening the analytic
tests exposed that the previous ~`1e-14` rate assertion had been made vacuous
by pytest's default `1e-12` absolute tolerance. All small-magnitude moment and
rate assertions now set `abs=0` and compare by declared relative error.

## Remaining closure gates

1. Extend or replace NIST Table 16 with primary attachment cross-section data
   that close both rate and energy tails over the Lam `T_e` board, with an
   evidence-backed uncertainty model.
2. Complete evaluated energy losses for molecular/atomic ionization,
   dissociation, excitation, detachment, elastic transfer, and charged-particle
   wall losses without importing fit exponents as thresholds.
3. Add a heavy-particle energy balance or measured gas-temperature boundary.
4. Solve `T_e` from measured/validated absorbed power, then preregister and
   grade reactor observables before emitting wafer fluxes.
5. Only after the reactor boundary passes independently may it be coupled to
   feature depth. Feature depth may not select any reactor parameter.
