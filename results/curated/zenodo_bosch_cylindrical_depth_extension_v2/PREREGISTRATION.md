# Bosch cylindrical depth extension v2

This is a calibration-visible model-form extension of the frozen v1 Bosch
gate.  The chronological heldout outcomes remain sealed.

The extension is necessary, not discretionary.  Across the 75 calibration
wafers, even an oracle shared radial curve has 2.128943% normalized Si-map
RMSE, above the frozen 2% heldout gate.  A separate oracle radial curve for
every wafer still has median 2.053148% RMSE.  A shared unconstrained `(x,y)`
map reaches 0.621753%, demonstrating a stable azimuthal signature.  The exact
replayable calculation is `axisymmetry_model_form_audit.json`.

V2 changes only the equipment-transfer model:

- the measured-waveform 0-D inventories are lifted through a positive,
  conservative cylindrical `(r,phi,z)` finite-volume operator;
- each effective species may carry a positive exponential Fourier source with
  complete orders selected sequentially up to order four;
- one shared wall-conditioning law may depend on the declared lot preparation,
  but independent lot/wafer offsets are forbidden;
- all equipment parameters have explicit bounds and are selected by
  leave-one-lot-out calibration;
- the Belen and La Magna/Garozzo surface mechanisms, v1 split, baselines, and
  acceptance gates remain unchanged.

No cylindrical fit was run before this file was frozen.  The mixed official
outcome CSV is forbidden to fitting code.  A hash-sealed heldout prediction is
required before any heldout numeric field is parsed.
