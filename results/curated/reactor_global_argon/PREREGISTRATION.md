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

## Execution record

Stage A and Stage B were executed on 2026-08-07. Stage B passed at both
endpoints of the source's published ion-wall-energy range:

- `5 Te`: MAPE `9.206%`, maximum APE `14.669%`;
- `8 Te`: MAPE `9.210%`, maximum APE `14.669%`;
- strictly decreasing electron temperature with pressure in both cases;
- maximum normalized particle/power residual `1.44e-14`.

The digitization manifest, per-point results, machine-readable grade, and
claim-boundary report are adjacent to this preregistration. Stage C remains
unexecuted and requires condition-specific primary experimental inputs.

## Stage C1 addendum: Mahoney 1994 independent argon ICP

Frozen on 2026-08-07 after the primary Mahoney paper and its Table I/Figure
11 were visually audited, and before any Mahoney-condition reactor result was
generated.

### Withheld source data and condition

- Primary source: `mahoney-1994-planar-icp`, DOI `10.1063/1.357672`.
- Source PDF SHA-256:
  `acd59d5def6373a81f1ec73248608b7c67a8769b010f2287fa05b04fd8cc61b7`.
- Grounded cylindrical plasma boundary: radius `0.114 m`, length `0.137 m`.
- Observable location: peak electron density on axis at `z = 0.050 m`.
- Withheld rows: every 100 W row printed in Table I: `10 mTorr` cryogenic
  pump; both `20 mTorr` cryogenic/mechanical repeats; `50` and `100 mTorr`
  mechanical-pump rows.
- Withheld observables: peak electron density and bulk electron temperature.
- The experiment reports net generator power, not calorimetric absorbed
  plasma power. The first board therefore uses `100 W` as a declared
  all-net-power-absorbed **upper-bound scenario**; it does not fit an
  absorption fraction.
- Neutral-gas temperature is unpublished. Both `300 K` and `600 K` endpoints
  are run as a declared sensitivity bracket. The bracket may not be narrowed
  from the target data.
- Both endpoints of Lee and Lieberman's published `5–8 Te` ion-wall-energy
  range are run. No endpoint may be selected after seeing the comparison.

### Frozen gates

Every combination of gas-temperature endpoint and ion-wall-energy endpoint
must:

1. converge with maximum normalized particle/power residual at most `1e-8`
   and return positive finite densities, temperatures, and fluxes;
2. predict electron temperature decreasing strictly along the
   `10, 20-mechanical, 50, 100 mTorr` sequence;
3. achieve bulk-electron-temperature MAPE at most `30%` and maximum APE at
   most `50%` across all five printed rows;
4. predict nondecreasing peak density along that same pressure sequence;
5. place `model peak density / measured peak electron density` within
   `[1, 5]` for every printed row. This interval is taken directly from the
   paper's stated electron-versus-ion diagnostic discrepancy, not from a
   reactor result; and
6. reproduce the normalized density-pressure shape with log RMSE at most
   `ln(2)` after both model and experiment are normalized to the 10 mTorr
   row. This normalization is a grade only and does not alter a model input.

The duplicate 20 mTorr rows are both scored for temperature and absolute
density. The mechanical-pump repeat is used in the ordered pressure-shape
sequence so no duplicate pressure is silently averaged.

Passing means the unchanged argon closure survives one independent
condition-specific plasma-state board within the source's explicit diagnostic
and missing-boundary limitations. It does **not** validate net-to-absorbed
power transfer, promote the transport closure to `validated_model`, supply a
C4F6 deck, or close Krüger depth.

### Stage C1 execution record

Executed on 2026-08-07 without changing the frozen gates or source-backed
argon coefficients. **Verdict: FAIL.**

- Every numerical solve converged with maximum normalized balance residual
  below `4.9e-15`.
- Electron-temperature and density-pressure trend gates passed at every
  neutral-temperature/wall-energy corner.
- Normalized density-shape log RMSE was `0.350–0.387`, inside `ln(2)`.
- At `600 K`, the electron-temperature MAPE was `29.86%`, inside the frozen
  `30%` limit; at `300 K` it was `37.78%`.
- The absolute center/model to measured-peak electron-density ratio was
  `4.05–19.09` across all corners, violating the frozen `[1,5]` interval.

The failure is preserved. Treating Mahoney's net `100 W` as plasma absorption
is an upper-bound scenario, and the paper explicitly leaves match/coil
dissipation and neutral-gas temperature unmeasured. Pressure-response shape
and conservation do not authorize fitting those missing boundaries from the
target density.

The implementation response is an explicit RF power-boundary contract:
measured plasma absorption can support prediction; forward-minus-reflected or
DirectDrive output power requires an independently measured downstream
hardware-loss closure. No reaction rate was changed to rescue this board.

## Stage C1 diagnostic addendum: target-inverted absorbed power

Frozen on 2026-08-07 after the Stage C1 FAIL and before running the
inversion. This is explicitly **not a validation gate**: it uses the withheld
Mahoney density to diagnose whether one constant RF transfer fraction could
explain the failed upper-bound board.

For each of the same four neutral-temperature/wall-energy corners and all
five Table I rows:

1. leave every chemistry, geometry, transport, and wall coefficient
   unchanged;
2. solve in log absorbed power for the powers at which model center density
   equals `2x` and `5x` the measured peak electron density, matching the
   paper's stated electron-versus-ion diagnostic interval;
3. divide each inferred absorbed power by the reported `100 W` net RF power;
4. intersect all five per-row `[P(2x), P(5x)]` intervals with the physical
   `0–100 W` net-power ceiling to test whether a single condition-independent
   transfer fraction exists; and
5. compare that target-inverted intersection to Hopwood's separately
   measured `70–90%` planar-ICP coupling range only as a nonportable external
   context band.

The inversion may not select a coefficient, relabel the Mahoney board PASS,
or become an absorbed-power provider. A narrow intersection means a hardware
measurement could close the board; no intersection means the missing
boundary cannot be represented by one constant transfer fraction under that
sensitivity corner.

### Stage C1 diagnostic execution record

Executed on 2026-08-07. Every sensitivity corner admitted a single constant
transfer-fraction interval satisfying all five rows:

- `300 K`, `5 Te`: `21.9–26.2%`;
- `300 K`, `8 Te`: `26.3–30.4%`;
- `600 K`, `5 Te`: `41.4–50.4%`; and
- `600 K`, `8 Te`: `49.7–59.0%`.

None overlaps Hopwood's separately measured `70–90%` context band. This does
not select the smallest or largest interval. It shows that one independently
measured constant hardware-transfer fraction could reconcile the pressure
series inside each declared corner, while the unresolved neutral temperature
and ion-wall energy move the inferred range by more than a factor of two.

The Stage C1 validation verdict remains **FAIL**. The inverted intervals are
target-informed diagnostics and are prohibited as production boundary inputs.
