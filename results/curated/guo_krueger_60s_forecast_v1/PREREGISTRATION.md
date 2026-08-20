# Preregistration — 60 s Guo/Kwon to Krueger no-fit feature forecast

Frozen before launching any 60 s trajectory from the current implementation.
The experimental endpoint may be opened only by a separate scorer after the
runtime audit is sealed.

## Scientific question

What depth, mask opening, and terminal profile are produced when the
source-fixed finite-fluence Guo/Kwon Si/O/C/F/V surface mechanism is advanced
through the unchanged common feature engine under Krueger's published wafer
boundary, without an ion-flux multiplier, oxide-yield multiplier, or
feature-depth fit?

## Numerical authority

The production member selected by the committed deterministic-prefix gates is:

- 10 nm uniform grid;
- 0.015625 s nominal time step, 3840 nominal steps over 60 s;
- deterministic extruded two-dimensional neutral exchange with analytic
  occlusion;
- 16 source positions, three face quadrature points, 8 x 16 neutral
  directions, and 16 ion azimuths;
- compressed joint IEAD at 250 eV / 0.25 degree bins;
- common-refinement surface-state remap;
- adaptive displacement controller with 0.35-cell target and 0.75-cell hard
  bound;
- explicit continuation only for gas-cavity topology changes;
- exact published aggregate ion flux, normalization 1.0;
- source-fixed 2.5 nm Guo translating-layer capacity;
- `blind_execution=true` so target observables are absent from runtime audits
  and plots.

The 0.015625 s step is the coarse member of the passing
0.015625/0.0078125 s time pair and the production step used by the passing
10/5 nm spatial pair.  It was selected without the 825 nm endpoint.

## Physical closures

Run these declared aggregate-ion identity endpoints independently:

1. `nominal_unresolved`: Krueger's aggregate ion row deposits energy but no
   guessed C/F atoms;
2. `all_cf2`: every aggregate ion is declared CF2+;
3. `all_cf3`: every aggregate ion is declared CF3+.

The second and third cases are envelope endpoints, not candidate reactor
mixtures.  No mixture fraction may be inferred from the experimental depth.
The nominal case is the primary published-boundary sensitivity.

## Terminal outcomes

A run is complete only if it reaches 60 s or the common engine returns a
declared `terminal_feature_clogged` event.  Wall-budget checkpoints are not
scientific endpoints and must resume from their exact conservative state.
Every runtime artifact must report `experimental_outcomes_read=false`.

After completion, freeze hashes for `audit.json`, `checkpoint.npz`, and the
execution receipt.  Only then may a separate scorer compare the prediction
with Krueger's depth/profile observables.

## Claim boundary

Even a numerical match would remain a transfer sensitivity until the missing
species-resolved positive-ion composition, stable C4F6 flux, high-energy
surface response, and mask coefficients are independently identified.  The
composition endpoints determine whether the unpublished ion identity can span
the measured depth; they do not identify which endpoint the reactor produced.
