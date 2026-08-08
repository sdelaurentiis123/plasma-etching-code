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

It does not supply a fundamental thermochemical ledger. The 600-dpi-audited
Figure 10 places ground-state Cl2 at zero and ground-state Cl at `1.25 eV` per
atom, encoding an approximate `2.50 eV` molecular asymptote. The official
COMSOL reproduction nevertheless uses `ediss=4 eV`, and its ground-state
dissociation rate fit contains an `8.84 eV` exponential parameter. The
evaluated physical value is `D0(Cl2)=2.4793 eV`. These three distinct numbers
demonstrate why a level coordinate, implementation event input, fit parameter,
and physical threshold must remain separate.

The raw COMSOL model also uses Figure 10's absolute `1.35/10.17 eV` atomic
coordinates directly as excitation gaps from ground Cl, despite that ground
state being at `1.25 eV`. Wang's primary fine-structure calculation gives the
physical first gap as `0.109 eV`. The exact COMSOL replay therefore remains a
separate code-verification mode and cannot supply the evaluated energy ledger.

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

Elastic momentum transfer now has two explicit, non-interchangeable modes:

- `STATIONARY_TARGET_ELASTIC_ENERGY_MOMENT` evaluates
  `2 me M/(me+M)^2 <sigma_m v E>` from the same momentum-transfer table. This
  is collision-exact for a stationary heavy target and applies the stricter
  energy-kernel support gate.
- `KEMANECI_ELASTIC_ENERGY_APPROXIMATION` reproduces Equation 18 as
  `3 Te (me/M) <sigma_m v>` and applies the rate-kernel support gate.

For a constant cross section and heavy Cl2 target, the Kemaneci form is about
25 percent below the collision-exact stationary-target moment because
collision events are speed-weighted. Neither mode upgrades an uncertain
cross-section table; the evaluated Cl/Cl2 momentum-transfer evidence remains
open and is documented in `REACTOR_CHLORINE_ELASTIC_GATE_2026-08-07.md`.

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
