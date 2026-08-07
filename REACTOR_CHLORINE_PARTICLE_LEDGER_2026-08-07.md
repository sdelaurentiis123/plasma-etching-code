# Native chlorine particle ledger - 2026-08-07

## Result

`src/petch/reactor_global/chlorine_particle_model.py` is the first open,
pressure-controlled chlorine reactor solve in the native stack. It is a
particle-balance rung, not yet a knobs-to-flux predictor.

For a supplied electron temperature and species-resolved charged-wall
transport, the model solves seven positive unknowns:

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
temperature. Charged transport is supplied as separate axial/radial flux
velocities for `Cl2+` and `Cl+`; the solver does not invent them from pressure.

## Deliberate fail-closed verdict

`ChlorineParticleSolution.supports_prediction` is always `False` at this rung.
The solution reports the missing closures explicitly:

1. electron power balance (so `T_e` follows absorbed power rather than being
   supplied);
2. condition-input evidence where measurements are absent;
3. predictive neutral diffusivity and incident-velocity evidence;
4. predictive species-resolved charged transport;
5. equipment-specific exhaust/conductance when the common volumetric throttle
   model is inadequate.

The neutral transport operator is exact for homogeneous radical decay with a
Robin wall. Ion neutralization is presently an integrated wall-return source
in the zero-dimensional ledger; a spatially inhomogeneous boundary-source
solution is a later sophistication rung. Nothing in this implementation is
fitted to the Malyshev dissociation board or to any etched depth.

## Next gates

1. Build an evidence-bearing Lee charged-transport provider from sourced ion
   momentum transfer and ambipolar diffusion, with the common-edge-factor
   assumption visible.
2. Replace the legacy neutral dissociation channel with the already frozen
   Hamilton state-resolved rates while preserving a source-reproduction mode.
3. Close physical electron energy losses before solving `T_e` from absorbed
   power.
4. Only then grade the Lam Alliance dissociation board as held-out reactor
   evidence and emit a feature-plane boundary.
