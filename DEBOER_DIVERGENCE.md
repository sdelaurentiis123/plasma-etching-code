# de Boer AR>20 divergence: characterization + charging-throttle test (2026-07-07)

Mission: characterize the de Boer deep-trench floor divergence precisely, then test ONE
hypothesis (a charging-derived floor ion-flux throttle). Scripts:
`scripts/deboer_divergence_arm.py` (runner) + `scripts/deboer_divergence_plot.py` (curves+figure).
Data: `deboer_divergence.npz`. Figure: `viz/deboer_divergence.png`.
Experiment gate: de Boer/Blauw normalized floor rate 1.0 / 0.43 / 0.29 / 0.20 @ AR 0/10/20/40.

## Headline: the banked divergence NARRATIVE was inverted by two artifacts

The working memory said "the EXPERIMENT's floor collapses above AR~20 and our simulation's does
not." Measured honestly, it is the OPPOSITE: **the simulation's floor collapses; the experiment
sustains ~0.20 to AR40.** Two artifacts had hidden this:

1. **Time-limited runs read as a stall.** The banked KNEE curve (deboer_final.npz, t_end=10)
   ends at AR19 -- that is end-of-run, not a stall. Run at t_end=45 the same config reaches
   AR31 (62 um of the 64 um substrate); on the tall stack (Lz=100, GPU) at t_end=90 it etches
   the FULL 92 um substrate to AR46+. Nothing in the model stops; it just gets slow.
2. **r0 normalization bias.** Long runs record depth coarsely (first record already at AR~3), so
   the AR<2 reference rate r0 is understated and the whole normalized curve inflates (this made
   the banked KNEE curve look like it matched at 0.43 @ AR10). Fix: r0 from a fine-cadence
   (record-every-step) short run of the same config, median smoothed rate over 0.5<AR<2.

## STEP 1 -- the corrected divergence (KNEE = the matched de Boer process config)

Tall stack (W=2 um, Lz=100, sub_top=94 -> AR ceiling 46), dx=0.25, 40k ions/neutrals,
2 seeds + fine-cadence r0 run, RTX 4090 (~80-90 s per 360-step run; identical physics to CPU).

| AR | experiment | petch KNEE | delta | sim/exp |
|---:|---:|---:|---:|---:|
| 5  | 0.715 | 0.758 | +0.043 | 1.06 |
| 8  | 0.544 | 0.545 | +0.001 | 1.00 |
| 10 | 0.430 | 0.309 | -0.121 | 0.72 |
| 15 | 0.360 | 0.223 | -0.137 | 0.62 |
| 20 | 0.290 | 0.148 | -0.142 | 0.51 |
| 25 | 0.267 | 0.116 | -0.193 | 0.43 |
| 30 | 0.245 | 0.066 | -0.179 | 0.27 |
| 40 | 0.200 | 0.051 | -0.149 | 0.25 |
| 44 | 0.200 | 0.056 | -0.144 | 0.28 |

- **Divergence onset: AR ~ 10** (matches within noise at AR<=8, -0.12 at AR10).
- **Magnitude: normalization-independent.** From AR8 to AR40 the experiment decays 2.7x
  (0.544->0.20); the simulation decays ~11x (0.545->0.05). No choice of r0 changes that ratio.
- DEFAULT (ViennaPS-regime, etchant-rich) shows the same deep collapse (0.073 @ AR30,
  RMSE 0.145 vs 0.177 for DEFAULT over AR5-44) -- the collapse is not a process-param artifact.
- Honest re-gate: KNEE RMSE(AR5-44) = 0.145 (the banked "RMSE 0.041 to AR20" does not
  survive the r0 correction + longer runs).

## STEP 2 -- hypothesis test: charging floor throttle Q(AR). VERDICT: REFUTED (FAIL)

Implemented opt-in (`Flags.floor_charge_throttle`, default OFF; 26-test suite green):
multiply the floor ion flux by Q(AR) from `charging_general.floor_charge_throttle_profile`
(first-principles table AR1-4 + deep-AR survivor values Q~0.25/0.60/0.50 at AR 8/15/25 from the
C12 charging runs, flat 0.50 extrapolation beyond AR25). Physical basis: cryo SF6/O2 grows a
dielectric SiOxFy passivation film that can trap charge even over grounded Si. Wired into
`mc_flux_3d_coupled` (threed.py `_apply_floor_charge_throttle`) before the coverage fixed point.

Result (same tall-stack protocol, before vs after):

| metric | KNEE | KNEE+throttle |
|---|---:|---:|
| nr @ AR20 | 0.148 | 0.168 |
| nr @ AR30 | 0.066 | 0.065 |
| nr @ AR40 | 0.051 | 0.066 |
| nr @ AR5 / AR8 | 0.758 / 0.545 | 0.618 / 0.412 |
| RMSE(AR5-44) | 0.145 | 0.142 |
| max AR reached | 46.8 | 45.7 |

- **No collapse "emerges"** -- the collapse was already in the baseline, with the opposite sign
  to what a flux throttle can address (sim UNDER-etches the deep floor).
- **The throttle is inert at deep AR** (differences within MC noise for AR>=15): the
  etchant-starved (cal_F=1.5) floor is NEUTRAL-limited, so halving the floor ion flux barely
  moves the rate. It only bites at mid-AR where the etch is ion-limited -- and there it makes
  the match WORSE (0.545->0.412 @ AR8 vs experiment 0.544).
- Charging (as an ion-flux reducer) is therefore NOT the missing physics for the de Boer
  AR>20 gap. The gap needs a mechanism that SUSTAINS the deep floor -- consistent with the
  earlier deboer_floor.py finding: the Belen rate has no large coverage-independent etch channel
  (spontaneous/thermal SF6 etch of Si at cryo + ion-assist recovery, or a direct F+Si channel
  that survives coverage starvation) rather than any further transport/charging suppression.

## Status of the opt-in code

- `src/petch/params.py`: `Flags.floor_charge_throttle` (default False -- default behavior
  unchanged; pytest 26/26 green).
- `src/petch/threed.py`: `_apply_floor_charge_throttle` + hook in `mc_flux_3d_coupled`.
- `src/petch/charging_general.py`: `floor_charge_throttle_profile` + `_THROTTLE_*` table.
The flag stays available as the "dielectric-passivated floor" experiment; it is NOT part of the
de Boer configuration.
