# Preregistration - SPTS Bosch reactor-to-wafer absolute-depth gate

Source: Sayyed et al., Zenodo 17122442, CC BY 4.0.

This gate uses the official 5 Hz machine record as the input and reserves the
last two experiment days, 2024-08-21 and 2024-08-22, for chronological
reactor-to-wafer transfer.  The split is frozen from process identifiers only;
the preregistration script never opens either measurement CSV.

This is an execution-held-out test, not a pre-exposure blind test.  The wafer
measurement files existed in the repository before this protocol was written.
That limitation is carried in the machine-readable manifest and may not be
silently upgraded to a blind-prediction claim.

## Frozen path

The accepted implementation path is:

1. decode the measured gas, pressure, source, platen, and thermal waveforms;
2. advance a deterministic, phase-resolved zero-dimensional reactor state;
3. condition the sheath and ion-energy boundary on the measured platen
   waveform;
4. propagate the boundary through an axisymmetric wafer-transfer operator;
5. advance the unchanged Belen silicon surface law with an explicit conserved
   C4F8-derived film-memory state.

Direct depth regression from date, lot, wafer number, experiment identifier, or
the held-out target is forbidden.  One set of equipment-transfer parameters is
shared across all wafers.  Model selection uses leave-one-lot-out validation on
the first eight experiment days only.

## Held-out board

The source process file contains 76 calibration records and 20 later process
records.  Available 89-point measurements will be joined only after a model
receipt is hash-sealed; missing source measurements remain missing rather than
being manufactured.

The held-out scorer grades wafer-mean silicon depth, the 89-point silicon depth
map, oxide-mask loss, selectivity, and within-lot drift.  Its absolute limits and
baseline requirements are frozen in `preregistration.json`.  Statistical
intervals bootstrap wafers, not the many spatial points on one wafer.

## Claim boundary

A pass validates absolute wafer depth, selectivity, radial transfer, and drift
for this SPTS Bosch process.  It does not validate feature charging, sidewall
angle, scallop shape, or aspect-ratio-dependent etching.  Those remain separate
feature-profile gates.
