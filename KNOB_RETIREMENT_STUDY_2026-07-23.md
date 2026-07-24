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

## Execution log

### Stage 1 — Lambda ladder (box 45653143, epoch d23190f, complete 2026-07-24 ~06:30 UTC)

All runs at fraction 0.8765 (R5 frozen value), declared operator, 10 nm grid:

| Lambda | opening (nm) | depth (nm) |
|--------|--------------|------------|
| 1.0    | 62.30        | 233.40     |
| 2.5    | 69.33        | 522.08     |
| 3.0    | 74.01        | 612.63     |
| 3.5    | 74.35        | 700.55     |

Depth is linear in Lambda (local slope 175.8 nm/unit at 3.0-3.5). Opening rises then
plateaus at ~74.3 nm. The base opening target is 45 +/- 5 nm; the fraction constant can
move opening by only ~-2.4 nm at its lower bound (probe column ~67 nm/unit, bound
0.84 vs 0.8765). **The base opening gate is therefore unreachable under the derived
law with the two permitted constants — recorded as a criterion miss per the failure
disposition (mixed-layer v1 target). The yield-scale knob stays retired.**

### Stage 2 — declared calibration (recorded 2026-07-24 ~06:50 UTC, BEFORE any stage-2 run completed)

Receipted-Jacobian solve for the depth gate with fraction at its opening-minimizing
bound: **fraction = 0.84, Lambda = 4.272** (predicted base ~(72.5, 825.0)). Queue:
dk-final (base), power 4/8 kW, oxygen 0.5/1.5/2.5 — all at this single calibration.
The power depth-ratios r(4/6) and r(8/6) are computed against dk-final and are the
decisive parameter-free prediction; the oxygen claims are scored per the declared
criteria above. Results to be appended verbatim when the queue drains.
