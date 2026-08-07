# Conserved global-reactor kernel and argon verification gate

Frozen on 2026-08-07 after the Lee--Lieberman equations, Table 3, Figure 3,
and Figure 8 were visually inspected, but before any reactor code was written,
any curve was digitized, or any reactor-model result was generated.

This is the first independent in-house reactor slice. It is intentionally
argon-only at the validation boundary. The implementation must remain
chemistry-agnostic so later O2, Cl2, and C4F6 decks use the same conservation
and unit machinery.

## Claim boundary

Passing this gate means:

1. the reaction-network and cylindrical-loss numerics preserve declared atoms
   and charge;
2. an unfit implementation can reproduce the published Lee--Lieberman argon
   global-model trend; and
3. the kernel is eligible for a separately preregistered comparison with
   independent reactor measurements.

It does **not** mean that a 0-D model is first-principles in its rate
coefficients, that the Lee--Lieberman paper is an independent validation of
itself, or that the Krüger C4F6 boundary has been recovered.

## Frozen source authority

- C. Lee and M. A. Lieberman, *Global Model of Ar, O2, Cl2 and Ar/O2 High
  Density Plasma Discharges*, UCB/ERL M94/49, revised 28 November 1994,
  `lee-lieberman-1994-global`.
- Source PDF SHA-256:
  `f2049e7041984d658d23688e8e8112a8d8e8a524172a8d2e335be8fde7fc2e23`.
- Reference reactor: cylinder length 7.5 cm, radius 15.25 cm, neutral
  temperature 600 K.
- Published Figure 3 condition: absorbed power 1000 W, flow 35 SCCM,
  recombination coefficient zero.

Every imported coefficient must retain its original unit, converted SI value,
source location, and evidence class. No coefficient may be selected from a
target electron temperature, density, flux, etch rate, or feature depth.

## Stage A: structural gates

For the chemistry-agnostic kernel:

- species names are unique and every reaction participant exists;
- volume reactions and closed wall-return events conserve every declared
  element and total electric charge, with electrons represented explicitly;
- any open feed, pump, deposition, or wall-retention term is represented as a
  declared boundary exchange and excluded from a closed-system claim;
- rate-law kinetic orders are explicit and are not inferred from reaction
  stoichiometry;
- centimetre-cubed-per-second to SI conversion is exact to floating-point
  roundoff (`1 cm3/s = 1e-6 m3/s`);
- reaction rates, source vectors, and power terms reject nonfinite or negative
  densities and nonpositive electron temperature;
- atom and charge source ledgers from a closed reaction network have normalized
  residual at most `1e-12`;
- an analytic first-order decay and an analytic bimolecular invariant agree
  with the kernel to relative error at most `1e-12`;
- cylinder volume, physical area, effective loss area, and diffusion length
  agree with their analytic definitions to relative error at most `1e-14`;
- all tests are deterministic.

## Stage B: Lee--Lieberman argon reproduction

The argon deck must use the visually audited Table 3 rate laws without fitting:

- ground-state ionization:
  `1.23e-7 exp(-18.68/Te) cm3/s`;
- metastable excitation:
  `3.71e-8 exp(-15.06/Te) cm3/s`;
- metastable step ionization:
  `2.05e-7 exp(-4.95/Te) cm3/s`;
- electron superelastic quenching:
  `2.0e-7 cm3/s`;
- metastable pooling:
  `6.2e-10 cm3/s`;
- Bohm wall loss and metastable diffusive wall loss use the paper's declared
  forms, with every additional transport input sourced before use.

The Figure 3 argon curve will be digitized from the archived source rendering
with a pixel-to-data calibration manifest and native-resolution visual
inspection. It is withheld from coefficient selection. At common pressures
inside the visible source support:

- electron temperature must decrease monotonically with pressure;
- mean absolute percentage error must be at most 10%;
- maximum absolute percentage error must be at most 20%;
- the solver's normalized particle and power residuals must each be at most
  `1e-8`;
- every returned density and flux must be positive and finite.

Failure is recorded as a failed reproduction, not repaired by changing a
published rate. Any missing transport input is resolved from a primary source
and recorded as a new evidence item before the run.

## Stage C: independent reactor validation

Figure 8 combines three experiments with different geometries, powers, flows,
and observables. It may not be graded as one anonymous curve. Before using any
of those points, the corresponding Ra, Mahoney, or Oomori operating conditions
must be recovered from its primary source and a condition-specific gate must be
frozen in an addendum.

Only a condition not used to choose a rate, wall coefficient, absorbed-power
fraction, or transport closure may be called held out.

## Link to feature depth

No reactor result from this board may be used to fit Krüger's 825 nm endpoint.
The path to a predictive depth claim is:

`recipe -> validated species densities/fluxes -> validated sheath/IEAD ->
unchanged feature physics -> depth/profile`.

Until the C4F6 molecular deck and the relevant reactor/sheath boundary are
validated on non-Krüger measurements, Krüger remains a published-boundary
sensitivity forecast rather than a reactor-derived absolute-depth prediction.
