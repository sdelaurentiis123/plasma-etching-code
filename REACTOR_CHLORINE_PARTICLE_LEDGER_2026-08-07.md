# Native chlorine particle ledger - 2026-08-07

## Result

`src/petch/reactor_global/chlorine_particle_model.py` is the first open,
pressure-controlled chlorine reactor solve in the native stack. It is a
particle-balance rung, not yet a knobs-to-flux predictor.

For a supplied electron temperature and a versioned charged-wall transport
provider, the model solves seven positive unknowns:

1. `n_Cl2`, `n_Cl`, `n_Cl2+`, `n_Cl+`, `n_Cl-`, and `n_e`;
2. the common volumetric exhaust-loss frequency required by the pressure
   controller.

The seven independent residuals are the five heavy-species balances,
quasineutrality, and the neutral pressure constraint. The implementation uses
the conserved Lee-Lieberman volume deck, exact cylindrical neutral Robin loss,
explicit `2 Cl -> Cl2` wall return, species-resolved positive-ion wall loss and
neutralization return, molecular feed, and neutral exhaust.

The fixed test condition closes every normalized solve residual below
`3.1e-15`. Independent integrated audits close:

- inlet chlorine atoms versus pumped neutral chlorine atoms;
- electron volume production versus positive-ion wall current;
- `n_Cl2 + n_Cl = p/(k_B T_g)`;
- each positive-ion axial wafer flux as its supplied axial flux velocity times
  its solved volume density.

These are structural verification checks, not comparison with experiment.

## Why pressure and flow are not a residence-time input

The model takes molecular flow and controlled pressure separately. It solves
the exhaust frequency after dissociation and wall return are known. In the
fixed regression condition, the solved exhaust frequency is more than `1.5x`
the pre-chemistry shortcut `Q_particles/(n_neutral V)`. The latter silently
assumes one feed molecule remains one exhaust particle and therefore fails as
soon as `Cl2 -> 2 Cl` is appreciable.

This resolves the bookkeeping issue identified on the Malyshev Lam board:
published flow and controlled pressure are valid open-system constraints, but
they are not permission to invent a pre-solved residence time.

The present controller uses a common well-mixed volumetric exhaust frequency
for `Cl` and `Cl2`. Species-dependent molecular-flow conductance remains a
future equipment closure and is one reason the model is not predictive.

## Evidence and physics boundaries

- Lee and Lieberman, `lee-lieberman-1994-global`, supplies the particle deck,
  quasineutrality structure, positive-ion wall return, and the source-model
  assumption that negative ions vanish at the sheath edge.
- Chantry, `chantry-1987-wall-diffusion`, supplies the partially reactive wall
  transport structure; petch evaluates the exact cylindrical Robin roots.
- Lymberopoulos/Economou and Ramamurthi/Economou,
  `economou-1995-2002-cl-transport`, supply the chlorine reduced diffusivity
  and confirm that neutral pumping and ion neutralization return are separate
  ledger terms.
- Stafford/Guha wall measurements supply the condition-scoped recombination
  boundary and expose the incident-velocity assumption.

Every scalar condition input carries units, source, evidence kind, and an
optional uncertainty. SCCM conversion requires explicit standard pressure and
temperature. Charged transport is recomputed from the current density state
inside every nonlinear ledger evaluation. A fixed provider remains available
only for declared source-reproduction and sensitivity checks.

The first dynamic provider reproduces the exact Lymberopoulos/Economou ion
mobilities and composes them with the Lee–Lieberman electronegative edge
closure. It exposes mobility, collision frequency, ion-speed convention,
momentum mean free path, ambipolar diffusivity, Bohm speed, electronegativity,
and axial/radial edge factors separately for `Cl2+` and `Cl+`. It is not
predictive: the papers' constant collision cross sections came from a private
communication, carry no uncertainty, and the identical mobilities were used
at conflicting ion temperatures (`0.12 eV` in 1995 and `300 K` in 2002).

An evaluated-dissociation deck now preserves all non-dissociation
Lee–Lieberman particle rows but replaces its one lumped `e + Cl2 -> e + 2Cl`
fit with the eight Hamilton state-resolved cross-section integrals. The
replacement is exact over the frozen `0.3--5 eV` Maxwellian table and carries
each state's excitation energy. The legacy deck remains unchanged for exact
source reproduction. This is still not a complete predictive chemistry deck:
the remaining attachment, detachment, neutralization, and species-resolved
ionization rows retain their separate evidence boundaries.

A separate detailed source-reproduction rung now replays Kemaneci's 36
printed non-elastic forward reactions with ten heavy states plus electrons.
It resolves the paper's charge-exchange numbering anomaly on the record and
enforces the printed `0.5--10 eV` domain. It is not yet the complete 46-feature
source model: two elastic channels and eight detailed-balance reverses remain
open, and no fit exponent has been imported as an electron-event energy.

Neutral wall transport is also now a versioned state-dependent provider. On
every nonlinear ledger evaluation it obtains the current `n_Cl/n_Cl2`,
evaluates the condition-scoped wall-response provider, and re-solves the exact
cylindrical Robin roots. The Stafford Figure-8 empirical providers reproduce
the committed direct-data regressions and report fit/leave-one-out residuals,
but remain non-predictive because the source has no absolute measurement
uncertainty and omits individual marker powers. They fail outside the measured
surface, ratio, pressure, power, and 300 K domains; in particular they cannot
be silently transplanted to Malyshev's 333 K Lam wall.

## Deliberate fail-closed verdict

`ChlorineParticleSolution.supports_prediction` is always `False` at this rung.
The solution reports the missing closures explicitly:

1. electron power balance (so `T_e` follows absorbed power rather than being
   supplied);
2. condition-input evidence where measurements are absent;
3. predictive neutral diffusivity, incident-velocity, and conditioned-wall
   evidence with uncertainty;
4. predictive species-resolved charged transport (the current provider is an
   exact source-model replay, not measured/evaluated transport);
5. equipment-specific exhaust/conductance when the common volumetric throttle
   model is inadequate.

The neutral transport operator is exact for homogeneous radical decay with a
Robin wall. Ion neutralization is presently an integrated wall-return source
in the zero-dimensional ledger; a spatially inhomogeneous boundary-source
solution is a later sophistication rung. Nothing in this implementation is
fitted to the Malyshev dissociation board or to any etched depth.

## Next gates

1. Replace the private-communication ion mobilities with public
   measured/evaluated momentum-transfer data and uncertainty when available;
   until then preserve this source-reproduction boundary.
2. Close physical electron energy losses before solving `T_e` from absorbed
   power. The first infrastructure rung now separates fixed event energies
   from the same-cross-section `<sigma v E>` moment required by electron-
   removing collisions; see `REACTOR_CHLORINE_ENERGY_LEDGER_2026-08-07.md`.
3. Only then grade the Lam Alliance dissociation board as held-out reactor
   evidence and emit a feature-plane boundary.
