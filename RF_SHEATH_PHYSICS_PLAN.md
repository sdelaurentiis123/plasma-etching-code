# RF sheath to feature-boundary physics

Date: 2026-07-11. The production feature source remains unchanged while this path is gated.

## Why the instantaneous arcsine source is not the final engine

Sampling `E = E_Bohm + Vdc + Vrf sin(phi)` assumes an ion instantaneously acquires the sheath voltage
at one RF phase. Its time-averaged density has integrable arcsine singularities at both energy horns.
That is a limiting model. Real ions integrate the electric field during a finite sheath transit; the ratio
of transit time to RF period controls peak separation and phase mixing. Collisions, charge exchange,
waveform harmonics, and a moving sheath edge further change the joint energy-angle-phase distribution.

Primary references reviewed:

- Kawamura, Vahedi, Lieberman & Birdsall, *Ion Energy Distributions in RF Sheaths: Review, Analysis
  and Simulation* (UCB/ERL M98/62, 1998).
- Charles et al., *Absolute measurements and modeling of radio frequency electric fields using a
  retarding field energy analyzer*, Physics of Plasmas (2000): intermediate transit-time phase mixing.
- Brinkmann et al., *Kinetic simulation of sheath dynamics in the intermediate RF regime*,
  arXiv:1305.6786: ion inertia, hysteresis, and agreement between kinetic approaches.
- Huybrechs & Cools, *Generalized Gaussian Quadrature Rules for Singular and Nearly Singular
  Integrals*, SIAM J. Numer. Anal. 47 (2009): endpoint singularities require matched quadrature.

## Implemented first-principles reduced model

`petch.sheath.CollisionlessRFSheath` integrates Bohm-entering ions through

`Phi(x,t) = Vs(t) (x/s)^(4/3)`

using measurable `Vdc`, `Vrf`, frequency, electron temperature, ion mass, and either density or sheath
thickness. When density is supplied, `s` is derived from Child–Langmuir current with Bohm influx.
Velocity-Verlet resolves both RF period and transit time. Current gates establish static energy gain,
Child thickness scaling, and high-frequency phase mixing.

## Unified-engine sequence

1. Define one `PlasmaBoundaryState` containing weighted joint species phase-space distributions and
   provenance. Analytic, sheath-ODE, reactor, diagnostic, and surrogate sources are constructors of the
   same object.
2. Make every forward and adjoint transport engine consume that object. Remove embedded source laws.
3. Add arbitrary voltage waveforms and self-consistent moving sheath edge.
4. Add ion-neutral collisions/charge exchange using cross sections and gas density, not fitted IEDF horns.
5. Gate energy/particle conservation, frequency and collision limits, and measured RFEA distributions.
6. Couple reactor outputs to the same state; a learned surrogate may accelerate the mapping but cannot
   redefine conservation or source normalization.

The instantaneous arcsine constructor remains useful as a named low-transit-time numerical limit. It must
not silently serve as universal production physics.
