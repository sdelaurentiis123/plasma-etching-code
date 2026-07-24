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

### Stage 2 — results and scorecard (queue drained 2026-07-24 ~10:50 UTC)

All runs at the declared pair (0.84, 4.272); audits archived in
`results/curated/knob_retirement_stage2/`:

| run | opening (nm) | depth (nm) |
|-----|--------------|------------|
| dk-final (base, 6 kW, O2 1.0) | 69.83 | 828.61 |
| power 4 kW | 56.45 | 769.22 |
| power 8 kW | 68.88 | 844.67 |
| oxygen 0.5 | 0.000 | 241.25 |
| oxygen 1.5 | 64.09 | 764.64 |
| oxygen 2.5 | 69.58 | 785.80 |

Scorecard against the declared criteria:

1. **Base depth: PASS** — 828.61 vs 825 +/- 5 (first Jacobian shot, no iteration).
   **Base opening: MISS** — 69.83 vs 45 +/- 5, exactly as declared unreachable in the
   stage-1 record (written before any stage-2 run). Mixed-layer v1 target.
2. **Power sweep: PASS (both, the decisive prediction)** — r(4/6) = 0.928 in
   [0.84, 0.94] (experiment 0.888; R5 blind 0.85); r(8/6) = **1.019** in [0.97, 1.06]
   (experiment/MCFPM 1.007; R5 blind 1.21). The 6->8 kW saturation that the frozen R5
   physics missed by 20 percent is now *emergent* from the unscaled deposited-energy
   yield law plus the neutral-supply ceiling — zero fitted parameters in the ratio.
3. **Oxygen sweep: clog at 0.5 PASS** (opening exactly 0.000); **necking-absent at
   2.5 PASS** (opening 69.6, no pinch); **saturation shape MISS** — depth rank
   maximum sits at ratio 1.0 (828.6) not 1.5, and the 1.5->2.5 change is +21.2 nm
   (criterion < 10). The declared half-saturation constant (4e20) over-suppresses the
   1.0->1.5 limb: the real curve rises to 1.5 then flattens; ours dips then re-rises.
4. **No refusals: PASS** — all six runs completed under the declared operator.

Disposition (per the preregistered failure clause): the yield-scale knob stays
retired — the power-saturation physics it was masking is now derived. The two
misses (base opening, oxygen-saturation shape) transfer to the mixed-layer v1
program (`src/petch/mixed_layer.py`, committed standalone at 7d2466b), where the
oxygen channel is an element ledger rather than a Langmuir constant: the
`oxygen_half_saturation_flux_m2_s` parameter becomes the next retirement target.
Blind credit for any of this still requires a fresh preregistered campaign on new
sealed data (Sawin/Yin beam database is the declared candidate).
