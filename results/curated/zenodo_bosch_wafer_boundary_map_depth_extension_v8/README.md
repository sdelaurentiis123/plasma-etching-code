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

## Chronological heldout result

The 20-record chronological prediction was generated from process inputs only,
hashed, committed, and pushed before the heldout numeric outcomes were parsed.
The sealed prediction SHA-256 is
`56ed2429832fe77280762fbca86cb6ffa4de3fd9687aa84f3b5cfd4ca99a3b1a`.
A separate post-seal extractor then found measured 89-point wafer maps for 13
of the 20 heldout process records.  Seven process records have no measurement
in the source asset and are reported as missing rather than imputed.

On the 13 measured heldout wafers:

- silicon wafer-mean MAE: 0.238906 um;
- silicon wafer-mean MAPE: 0.541245%;
- silicon 89-point RMSE: 0.329751 um;
- normalized radial-shape RMSE: 0.261855%;
- oxide wafer-mean MAE: 0.023908 um;
- selectivity MAPE: 3.161040%;
- silicon wafer-mean correlation: 0.873410.

The frozen absolute gates all pass.  The physics prediction also beats the
calibration-global mean depth baseline (0.400599 um mean MAE), the frozen
calibration mean-map point baseline (0.477852 um point RMSE), and its radial
shape baseline (0.294344%).  At the wafer-bootstrap level, the depth and point
advantages remain positive across the 95% intervals.  The shape advantage is
smaller: its bootstrap interval crosses zero, so it is a point-estimate win,
not evidence of a statistically resolved shape improvement.  The predicted
within-lot depth drift also beats a zero-slope baseline.

The target firewall records that the prediction did not change after reveal.
This is a successful heldout validation of absolute reactor/wafer depth,
selectivity, radial transfer, and drift for one Bosch tool.  It does not by
itself validate feature charging, sidewall angle, scallop geometry, or ARDE.

## Artifacts

- `preregistration.json` and `PREREGISTRATION.md`: pre-code scientific freeze;
- `exact_wall_ion_response_table.npz`: 13x13x75x89 exact response tensor;
- `full_calibration_capacity.json`: all 20 full-calibration candidates;
- `calibration_fit.json`: whole-lot fits, stability, and selection;
- `interpolation_validation.json`: independent tensor-midpoint audit;
- `exact_replay.json`: exact full/whole-lot replay and refinement receipt;
- `heldout_prediction.json` and `heldout_prediction_seal.json`: pre-reveal
  prediction and cryptographic seal;
- `heldout_score.json`: post-seal heldout score, baselines, bootstrap, drift,
  missing-outcome ledger, and target-firewall receipt.
