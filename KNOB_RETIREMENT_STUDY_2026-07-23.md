# Knob-retirement study K24-DEKNOB-1 (declared 2026-07-23, before any study runs)

Open-book post-reveal study (held-out observations already unsealed by R5); its
discipline is the DECLARED SUCCESS CRITERIA below, written before execution.

## Configuration under test (vs the frozen R5 authority)

1. `oxide_etch_yield_scale` knob (0.5868…) REMOVED. Substrate/complex/polymer ion yield
   laws use `energy_model="deposited_in_layer"` (Sigmund form: yield proportional to
   ZBL nuclear energy deposited within the reactive layer; `ion_energy_deposition.py`;
   layer depth 1.5 nm, literature band 1-3 nm, ratio-insensitive <3 percent). The single
   absolute-rate constant per channel (Sigmund Lambda, entering as the reference-anchored
   yield) is calibrated ON THE BASE CONDITION ONLY. The fitted knee is superseded.
2. `oxygen_half_saturation_flux_m2_s = 4.0e20` (m^-2 s^-1; from the figure-16a oxygen
   flux scale, declared here, not tuned per-case).
3. Mask crosslinked-growth fraction recalibrated on base (unchanged role).
4. Everything else identical to the frozen R5 authority operator.

## Declared success criteria (all simultaneously)

- Base: opening and depth within +/-5 nm of (45, 825) after recalibrating exactly two
  constants (mask fraction, Lambda-scale) by the receipted Jacobian/preview procedure.
- Power sweep (the parameter-free prediction): depth ratios r(4/6) in [0.84, 0.94] and
  r(8/6) in [0.97, 1.06] (experiment/MCFPM: 0.888, 1.007; R5 blind values 0.85, 1.21 --
  the 8 kW ratio is the decisive number).
- Oxygen sweep: clog at 0.5 preserved; necking-absent at 2.5 preserved; depth rank
  maximum at 1.5 OR 1.5-to-2.5 increase < 10 nm (the saturation claims).
- No refusals; conservation ledgers exact as always.

## Failure disposition

Any criterion missed is recorded as-is; the yield-scale knob is NOT reintroduced --
the miss becomes the mixed-layer v1 target (RESEARCH_MIXED_LAYER_DESIGN_2026-07-23.md).
Blind credit for any of this requires a fresh preregistered campaign on new sealed data.
