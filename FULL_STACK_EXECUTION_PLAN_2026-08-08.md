# Knobs-to-flux-to-depth execution plan - 2026-08-08

## Objective and claim rule

Build one evidence-carrying chain:

`machine knobs -> absorbed plasma power and gas state -> EEDF and reaction rates
-> volume species -> sheath-edge species -> wafer flux/IED/IAD -> evolving
surface state -> feature profile and depth`.

The output is a prediction only when no feature depth selected any upstream
reactor, transport, sheath, or surface parameter. A facility-conditioned mode
may use an independently measured blanket observable, but it must be reported
separately from a first-principles/published-knob mode and tested on untouched
features.

## Current truth about depth and HAR

| capability | status on 2026-08-08 |
|---|---|
| analytic neutral/energetic transport | certified through aspect ratio 200:1 |
| coupled rate-ARDE | executable; 60% decline at AR50 lies inside the broad 43--80% literature scale band, not an absolute validation |
| full evolving HAR profile | open; the straight-wall operator correctly refuses after about 4 nm once taper develops |
| SF5+ surface depth per measured dose | retrospective no-yield-fit transfer, four points, 5.88% MAPE |
| Cl2/Ar+ silicon ALE absolute depth | retrospective cross-source transfer, three points, 12.88% maximum relative error |
| Krueger Ar/C4F6/O2 feature depth | MISS: 346.833 nm versus 825 nm under an underidentified aggregate boundary |
| formal value-blind held-out feature-depth passes | zero |
| ready untouched targets | 49 Yoshie cyclic SF6/C4F8 feature depths; phase boundary and persistent surface memory still open |

Thus the feature mechanics and deep transport remain useful, but “all other
depths are working” is not an allowed summary.

## Gate 1 - collision decks and swarm physics

### Landed now

`electron_collision_deck.py` provides rights-safe BOLSIG/LXCat-format
ingestion. It reads a user-supplied local deck, verifies an optional SHA-256,
preserves process type, products, mass ratio, energy loss, statistical weights,
comments, and every cross-section node, and packages no third-party database
bytes. Structural completeness deliberately does not imply swarm, reactor,
flux, or depth support.

The parser passes manufactured corruption/evidence tests and parsed all 16
legacy SIGLO Cl2 processes in a local-only audit.

The BOLSIG/LXCat deck contains a momentum-transfer row, not the higher
Legendre moments of an elastic differential cross section. Multi-term use
therefore fails closed until an angular closure is declared. The first
Kawaguchi source-reproduction mode may use the paper's explicit isotropic-
scattering assumption; that mode must remain separate from a later provider
backed by angular differential data.

`electron_kinetics.py` now lands the deterministic fixed-grid foundation and
the first physical two-term solve. It carries exact finite-volume EEPF
normalization/mean-energy weights, exact piecewise-linear cross-section rate
and incident-energy moments, support/tail refusal, analytic moment and
normalization JVP/VJP contracts, conservative Scharfetter--Gummel energy
fluxes, excitation scattering-out/in, equal-sharing ionization, attachment,
elastic gas exchange, field heating, and temporal-growth normalization. The
zero-field elastic limit recovers the gas Maxwellian and every returned state
keeps reactor/wafer/depth evidence false.

The independent manufactured-deck receipt at
`results/curated/reactor_global_kinetics/two_term_bolos_oracle_v1.json`
converges monotonically toward local LGPL BOLOS 0.2. At 2400 cells the mean
energy residual is `0.458%`, excitation-rate residual `0.634%`, and weighted
EEPF L1 residual `0.122%`. This passes the numerical two-term operator gate;
it does not pass a physical collision set or any measurement-defined swarm
observable. The current piecewise-constant inelastic reconstruction remains
first order; a piecewise-exponential upgrade is the performance/accuracy path
before large physical decks are graded.

### Required solver physics

A simple local-field EEDF and scalar mobility are insufficient for the direct
2018 board. The compared observables have different transport definitions:

- `Wm` is pulsed-Townsend mean-arrival-time drift;
- `NDL` is pulsed-Townsend spatiotemporal longitudinal diffusion;
- `alpha/N`, `eta/N`, and `(alpha-eta)/N` require the corresponding
  nonconservative steady-state Townsend treatment.

The native solver must therefore use a deterministic multi-term spherical-
harmonic Boltzmann expansion with a density-gradient/spatiotemporal hierarchy
and a separate steady-state-Townsend eigenproblem. Flux mobility may not be
substituted for `Wm`, scalar diffusion may not be substituted for `NDL`, and a
single generic EEDF may not be substituted for the experiment-specific EVDFs.

### Deterministic differentiable architecture

Use fixed-topology, conservative energy finite volumes and a configurable
Legendre order. Each field/composition state becomes a sparse block-banded
bordered solve for the angular harmonics, growth eigenvalue, and normalization.
The density-gradient hierarchy reuses the same factorization for successive
transport orders. Positivity, normalization, particle growth, and energy
closure are residuals of the operator, not post-processing repairs.

Differentiability is analytic and implicit rather than a finite-difference
wrapper. For a linear solve `A(theta)x=b(theta)`, JVPs/VJPs reuse solves with
`A` and `A.T`. Eigenproblems and nonlinear reactor roots add their
normalization/closure equations to a bordered Jacobian and use the same
implicit-function adjoint. This gives deterministic gradients without
backpropagating through iterations.

The public API is batch-first: energy grids and process topology are static;
cross-section interpolation is positivity-preserving on fixed knots; field,
mixture, pressure, power, flow, bias, and geometry conditions carry a leading
batch axis. Conditions are independent and therefore embarrassingly parallel.
A NumPy/SciPy reference backend and compiled GPU backend must implement the
same residual/JVP contract and pass bitwise-or-tolerance parity gates.

Convergence must be demonstrated in energy resolution, upper energy, density-
gradient order, and Legendre order. With the isotropic replay provider,
orders above one test numerical convergence under the source assumption; they
do not upgrade the angular collision evidence. A DCS-backed provider is a
separate later gate.

Acceptance requires:

1. normalization, positivity, particle/growth, and energy ledgers;
2. energy-grid and maximum-energy convergence;
3. manufactured elastic, excitation, attachment, and ionization limits;
4. independent agreement with a standard deterministic solver on
   redistributable manufactured decks;
5. a preregistered no-fit grade against all 52 direct pure-Cl2 measurements;
6. separate residuals for drift, effective ionization, and longitudinal
   diffusion, with the source uncertainty ranges retained.

BOLSIG+ may be used only as a local deterministic scientific oracle: its
official terms forbid commercial use and third-party redistribution. BOLOS is
LGPL and useful as an independent implementation cross-check, but its basic
flux mobility and scalar diffusion outputs do not by themselves reproduce the
measurement definitions above. LXCat bytes remain user-supplied until
contributor permission or a compliant interface exists.

## Gate 2 - EEDF-coupled chlorine volume reactor

Replace supplied `Te` with the converged EEDF and integrate every accepted
collision channel over that distribution. Complete electron energy losses,
including channel thresholds, attachment-selected kinetic energy, elastic
transfer, ionization partition, and charged-wall losses. Then solve particle
and power balance simultaneously with pressure control, molecular feed,
equipment conductance, neutral wall return, and gas temperature.

Acceptance requires exact atom/charge/current/power closure and refusal when a
cross-section tail, channel, wall state, absorbed-power boundary, or equipment
conductance is outside evidence support.

## Gate 3 - knobs and equipment power

Treat generator setpoint, forward/reflected power, and absorbed plasma power as
different objects. For a Lam-equipment-class provider, carry chamber/active
volume, coil/window or electrode losses, pressure control, flow, frequency,
gap, and bias explicitly. Infer a loss model only from electrical or plasma
diagnostics, never from etch depth.

The first target should be the public Lam Alliance chlorine board because the
repository already contains 62 electron-temperature, 27 electron-density, and
38 dissociation markers across power, pressure, and gap. This earns a public
equipment-class model; it does not claim proprietary fidelity to every Lam
platform.

## Gate 4 - charged transport, sheath, and wafer boundary

Replace private-communication ion mobilities with public measured/evaluated
species-resolved transport where possible. Solve electronegative edge losses,
RF sheath dynamics, ion-neutral collisions, and phase-resolved transit with
deterministic characteristics or discrete ordinates in phase/energy/angle.
Emit for every positive ion and reactive neutral:

- absolute flux and covariance;
- energy and angle distribution, not one effective ion;
- sheath/collisionality regime and applicability receipt;
- neutral/radical temperature and angular state;
- source-plane-to-wafer particle and energy closure.

Grade total and `Cl+`/`Cl2+` flux against the existing 24-point GEC-ICP board
before any feature result is opened.

## Gate 5 - feature coupling and first absolute-depth prediction

Pass the immutable species/energy/angle boundary into the common feature
engine. The first end-to-end feature campaign should be chlorine on silicon,
not Krueger, because the public Lam/chlorine reactor boards, Cl2 swarm data,
controlled Cl2/Ar+ surface evidence, and Mahorowala fixed-75-s profile/rate
board overlap in chemistry.

Use controlled-beam surface observables for the surface law and hold all 11
Mahorowala fixed-time rates and 13 profiles out from reactor/surface fitting.
Only after this passes should the chain expand to cyclic SF6/C4F8 and then
C4F6/oxide. Krueger becomes gradeable only if its species-resolved ion and
stable-parent boundary is independently recovered or predicted from a
validated C4F6 reactor deck.

## Gate 6 - full HAR evolution

Replace the straight-cylinder restriction with the planned tapered-profile
self-pair/form-factor kernel, retain exact enclosure and per-step conservation,
add mask evolution, and run grid/form-factor/time-step convergence. The
existing 200:1 transport result remains a valid transport benchmark; a full
200:1 etched profile is a later claim that must pass this gate.

## Performance target

- collision interpolation, sparse symbolic structure, and factorizations
  cached by deck hash, grid, mixture, and gas state;
- one 0-D steady reactor condition in seconds or less after cache warm-up;
- batched knob sweeps vectorized across conditions and sharded without shared
  mutable state;
- feature transport uses deterministic quadrature/radiosity operators; CPU
  reference and compiled GPU paths share conservation and gradient gates;
- gradients come from operator JVP/VJP and implicit adjoints, never stochastic
  estimators or finite differences in the released path;
- surrogates train only on validated governing-solver output and every released
  result remains checkable against that solver.

The shortest defensible route to a first knobs-to-depth result is therefore:

`user-supplied Cl2 deck -> transport-definition-correct swarm solver -> Lam
state validation -> mass-resolved chlorine wafer boundary -> held-out
Mahorowala depth/profile board`.
