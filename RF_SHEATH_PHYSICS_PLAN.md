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

## Boundary-state implementation and current research implications

`petch.boundary_state` now implements the common contract as immutable weighted joint velocity-energy
measures with absolute species flux, charge, mass, RF phase, optional wafer position, reference plane,
and provenance. Both the instantaneous limit and `CollisionlessRFSheath` construct it.

Recent and foundational work sharpens what later constructors must preserve:

- HPEM/MCFPM reactor-to-feature coupling passes spatially resolved species fluxes and energy/angle
  distributions, not recipe knobs directly (Hoekstra, Kushner & Sukharev, OSTI 323619).
- Hybrid bulk/sheath models produce coupled IEDF and IADF in electronegative biased ICPs; energy and angle
  cannot generally be factored (OSTI 1399615).
- Measured sheath waveforms and density can drive fast virtual IED sensors including ion-neutral
  collisions (Bogdanova et al., arXiv:2012.14882).
- Tailored and chirped waveforms deliberately reshape IEDFs, so arbitrary measured voltage waveforms must
  be first-class inputs rather than reduced to one `Vrf` (Lanham & Kushner, DOI 10.1063/1.4993785;
  Giesekus et al., arXiv:2509.01171).
- Charge exchange creates secondary peaks and couples pressure to IEDF/IADF (Georgieva et al.,
  Phys. Rev. E 69, 026406).
- Reactor wall seasoning changes radical and ion flux composition wafer-to-wafer, requiring boundary
  state provenance and time/wafer indexing (Agarwal & Kushner, DOI 10.1116/1.2909966).

Accordingly, no downstream feature module may reconstruct missing energy-angle correlations from scalar
means. Surrogates may predict a boundary state, but conservation, normalization, support, and uncertainty
remain explicit fields and gates.

`petch.boundary_transport` is the first real consumer. It tensors any species measure with spatial
quadrature, preserves normalized probability, absolute flux, and three-dimensional kinetic energy, and
feeds the nodal particle tracer without reconstructing a source law. The identical adapter passes open-
wafer gates for `Ar+` and neutral `CF2`. Migration is partial: existing production forward and adjoint
charging paths still contain embedded analytic sources.
