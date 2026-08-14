# Current-driven moving-RF reactor stack — 2026-08-13

## Outcome

The deterministic pure-Ar stack now has a physically explicit path

`absorbed bulk condition -> global plasma -> axial Bohm flux -> moving RF sheath -> phase-resolved collisional wafer boundary`.

The previous collisional rung solved all ion collision orders in a static
phase-conditioned Child field.  The new operator advances ions through a
time-dependent Poisson field with a moving electron front.  RF phase is part
of the kinetic state through every elastic and resonant charge-exchange event.
The collision series is summed by one sparse absorbing solve, not particle
Monte Carlo.

This is a reactor-boundary breakthrough, not a claimed Krüger depth match.
The code still refuses feature depth because the remaining evidence gaps are
physical rather than numerical.

## What landed

1. `PeriodicCurrentDensity` represents an arbitrary zero-mean Fourier sheath
   current with explicit evidence: measured, validated circuit, or assumed.
2. `TurnerChabertCurrentDrivenSheath` implements Turner and Chabert equations
   1--19: charge excursion, waveform-dependent `xi`, maximum width and
   voltage, moving electron front, instantaneous potential and field, and
   deterministic ion trajectories.
3. The exact common-current-scale JVP is analytic: `s_max ~ J^3` and
   `V_max ~ J^4`.
4. `DeterministicMovingCollisionalRFSheath` solves on fixed
   `phase x direction x position x energy x transverse fraction` ordinates.
   A sparse `(I-Q.T)` solve closes every ion elastic/CX collision order and
   differentiates exactly with respect to neutral density.
5. `DeterministicCurrentDrivenArgonReactorToWaferModel` connects the global
   plasma to the sheath without using the volume-average electron density at
   the wrong plane: `n_s = Gamma_i,axial / u_B`.
6. The common `PlasmaBoundaryState` export carries the arriving Ar+ IEAD/IAD
   and the explicitly lower-bound fast-neutral branch.

The nonlinear sheath closure is cached once per condition.  This reduced the
13 focused sheath/kinetic gates from roughly 141 s to 11 s without changing
the equations or results.

## Frozen numerical audit

The machine-readable receipt is
`results/curated/current_driven_argon_reactor_stack/audit.json`; regenerate it
with:

```bash
python scripts/audit_current_driven_argon_reactor_stack.py
```

The manufactured audit uses 500 W absorbed bulk power, 2 mTorr Ar, 500 K,
15 cm radius, 7.5 cm length, and a declared two-harmonic 2 MHz current.  No
depth or IEAD target is used.

| Ledger | result |
|---|---:|
| global particle/power maximum residual | 3.394e-15 |
| global-to-sheath Bohm-flux seam | 0 |
| Child-current residual | 3.808e-16 |
| charge-voltage residual | 1.328e-16 |
| sparse kinetic-solve residual | 5.421e-19 |
| probability-ledger residual | 0 |
| maximum two-body collision-energy residual | 7.718e-16 |

At this declared condition the stack produces a 453.165 V maximum sheath,
0.8386 mm maximum width, 4.659e20 m^-2 s^-1 source Ar+ flux, and a resolved
234.458 eV mean wafer ion energy with 1.151 degree RMS impact angle.  Those are
manufactured-condition outputs,
not equipment validation data.

The 48-to-96-step temporal-characteristic gate passes: mean energy changes
0.00538%, RMS angle 0.00240%, expected collision count 0.0789%, and arrival
probability 0.000189%.  The 12-to-24 RF-phase ordinate gate also passes:
mean energy changes 0.684%, RMS angle 1.104%, expected collision count 0.249%,
and arrival probability 0.000639%.  A production claim must additionally
carry condition-specific position/energy/transverse-ordinate receipts and an
independent measured IEAD grade.

## Physics and evidence gates that remain open

1. **Low-energy Ar+--Ar experimental grade.** The deterministic operator now
   uses the Phelps/LXCat isotropic/backscatter decomposition below 400 eV,
   enforcing `Qm = Qi + 2 Qb` and using the physical event cross section
   `Qi + Qb`.  This closes the previous zero-angle fallback as a model.  It
   does not replace an independent measured IEAD grade on the target tool.
   The audit routes 0.01505 collisions per source ion through this model.
2. **Repeated fast neutrals.** The first born-neutral flight is resolved;
   subsequent neutral--neutral collisions remain a ledger.  The audit reports
   2.317e-4 unresolved collisions per source ion.
3. **Generator and matching network.** The model consumes current at the
   sheath reference plane.  Generator forward watts cannot determine it
   without the matching-network, stray, electrode, and plasma impedances.
4. **Equipment validation.** The frozen audit deliberately uses assumed bulk
   transport and current, so the resolved-ion predictive flag is false.
5. **Molecular chemistry and surface kinetics.** Pure-Ar ion transport cannot
   identify Krüger's CFx+ mixture, neutral fluorocarbon supply, or the surface
   removal/polymerization law.  Consequently `supports_feature_depth` stays
   false and the certified Krüger published-boundary result remains
   346.83264081620524 nm versus 825 nm measured.

## Experimental package needed for knobs -> depth

For one machine and one recipe family, measure or export:

- absorbed source power (or a calibrated generator-to-absorbed transfer),
  pressure, flow, gas temperature, and chamber dimensions;
- de-embedded electrode/sheath current waveform or the full RF network and
  plasma-load diagnostics needed to calculate it;
- mass-resolved wafer ion flux plus an IEAD/IAD condition for validation;
- neutral radical flux/density and gas composition at the wafer;
- blanket etch/deposition rates on the actual material and mask; and
- cross-sectional SEMs with exact time, geometry, and recipe provenance.

That package makes the boundary identifiable.  It does not justify fitting a
surface yield to one feature depth; the blanket observables set the surface
law, and held-out feature profiles grade depth.
