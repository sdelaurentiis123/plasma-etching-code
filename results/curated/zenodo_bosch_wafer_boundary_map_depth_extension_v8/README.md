# Bosch v8 conservative wafer-ion boundary map

## Scientific question

Bosch v7 transferred mean depth across whole held-out calibration lots but
lost to a trivial mean-map baseline on radial point and normalized-shape
accuracy. A calibration-only residual decomposition localized 86.2% of the
squared spatial error to one repeatable machine-coordinate map and a smaller
C4F8 platen-Vpp edge mode. V8 tests whether those terms can be expressed as a
physical wafer-boundary operator rather than an output-depth correction.

## Frozen operator

The preregistered law is a positive complete-real-Zernike field applied only to
positive-ion wafer transmission:

```text
u(r,phi,Vpp) = exp(static Zernike field + z(Vpp)*dynamic field)
T = u * integral(Jion_base dA) / integral(Jion_base*u dA)
```

The finite-volume normalization conserves total positive-ion current exactly.
The law does not change ion energy, neutrals, v7 wall memory, or the frozen
Belen and La Magna/Garozzo surface mechanisms. No depth correction is applied.

Twenty candidates were evaluated: static maximum orders 1-10, each with no
dynamic term or a frozen order-2 dynamic term driven only by measured C4F8
platen Vpp RMS. Selection was the lowest coefficient count passing every
absolute, whole-lot, identifiability, stability, and slope gate.

## Selected candidate

- static maximum order: 9;
- dynamic maximum order: 2;
- coefficient count: 59;
- exact full-fit Jacobian rank: 59/59;
- full-fit condition number: 9.954;
- maximum parameter correlation: 0.714;
- no coefficient bound contact;
- maximum raw log field: 0.14744, below the frozen `ln(2)` ceiling.

This is a smooth machine-specific calibration for the measured SPTS tool and
coordinate frame. The 59 coefficients are not claimed as universal plasma
constants.

## Exact replay result

The selected candidate was replayed through the exact reactor and unchanged
surface recurrence rather than graded from the accelerated tensor fit.

Full calibration:

- silicon mean MAE: 0.223668 um;
- silicon mean MAPE: 0.511315%;
- silicon point RMSE: 0.387056 um;
- normalized shape RMSE: 0.620074%;
- oxide mean MAE: 0.039242 um;
- selectivity MAPE: 6.015456%.

Whole-lot leave-one-out:

- silicon mean MAE: 0.270179 um versus 0.338486 um global-depth baseline;
- silicon point RMSE: 0.431199 um versus 0.486585 um mean-map baseline;
- normalized shape RMSE: 0.635267% versus 0.636619% mean-map baseline;
- within-lot slope MAE: 0.077810 um/wafer versus 0.082903 v7 gate.

All exact replay gates pass. The shape advantage is real but narrow, which is
why the exact and refinement gates are binding rather than optional.

The independent 144-point wall/ion tensor-midpoint audit reaches at most
0.036988 of any frozen gate. At the selected parameters, the exact-versus-
accelerated differences reach at most 0.008727 of a full-calibration gate and
0.009672 of a whole-lot gate.

Increasing the cylindrical grid from 24x24x16 to 32x32x24 changes no observable
by more than 0.043559 of its frozen gate; the denser Zernike certification grid
also preserves the field bound. The refinement gate passes.

## Target firewall and remaining step

The chronological heldout outcome remains unopened:

```text
heldout_outcomes_read = false
heldout_prediction_written = false
eligible_for_prediction_seal = false
```

The next step is to construct and hash the chronological heldout prediction
from process inputs only, commit it with all code/data/parameter/fold hashes,
and only then allow a separate scorer to parse the heldout numeric outcome.
If that heldout score loses the frozen baselines, v8 fails despite its
calibration and whole-lot results.

## Artifacts

- `preregistration.json` and `PREREGISTRATION.md`: pre-code scientific freeze;
- `exact_wall_ion_response_table.npz`: 13x13x75x89 exact response tensor;
- `full_calibration_capacity.json`: all 20 full-calibration candidates;
- `calibration_fit.json`: whole-lot fits, stability, and selection;
- `interpolation_validation.json`: independent tensor-midpoint audit;
- `exact_replay.json`: exact full/whole-lot replay and refinement receipt.
