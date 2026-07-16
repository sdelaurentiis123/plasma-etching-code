# de Boer AR>20 divergence: characterization + charging-throttle test (2026-07-07)

> **Evidence correction (2026-07-15):** this document treated the normalized
> `1/.43/.29/.20` Blauw/Clausing calculated curve as measured de Boer data. That identification and
> all experimental-pass language are withdrawn. The numerical mechanism study is retained as legacy
> development history only; see `ARDE_LEGACY_UNIFIED_RECONCILIATION_2026-07-15.md` for the corrected
> evidence map and common-engine migration result.

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

## STEP 3 -- hypothesis test: PASSIVATION-LINKED wall loss + the transport-artifact finding (2026-07-07)

Hypothesis (mission): the uniform Knudsen wall-loss scale fits the knee but starves the deep floor;
a FRONT-CONCENTRATED wall loss (full on the freshly-exposed contested band near the etch front,
reduced 5-20% on the passivated SiOxFy sidewall column above it) should do BOTH. Implemented opt-in
(`Flags.knudsen_front_loss`, default OFF; pytest 26/26 green): a per-slice effective loss that relaxes
from the full base scale (shallow / all-fresh walls) toward `knudsen_passive_frac*base` (deep /
mostly-passivated) over an AR band -- see `knudsen._front_loss_scale`.

**First, a premise correction that dominates the result.** The banked AR46 divergence curves
(STEP 1/2) run `neutral_transport="mc"`; the Knudsen path with `knudsen_wall_loss_scale` was NEVER
exercised in them (it lives in the notching/evolving calibration). Testing the hypothesis therefore
means running the de Boer geometry under `neutral_transport="knudsen"`. Doing so exposes the real
lever: **the deep-floor collapse is largely a Monte-Carlo deep-floor UNDER-SAMPLING artifact, not a
missing wall-loss term.** Deterministic Knudsen transport already sustains the floor.

Tall stack (Lz=100, sub_top=94 -> AR46), 2 seeds + fine-cadence r0, KNEE process params, CPU
(~270 s/run; Knudsen is deterministic + fast). Same r0/normalization protocol as STEP 1
(analyze_tall). r0 is identical across bases (15.0 um/t): at AR<2 the open floor sees full flux
regardless of wall loss.

| AR | experiment | MC banked (STEP1) | Knudsen front-OFF (wls 2.9) | Knudsen front-ON base 2.9 | Knudsen front-ON base 2.6 |
|---:|---:|---:|---:|---:|---:|
| 5  | 0.715 | 0.758 | 0.667 | 0.667 | 0.667 |
| 8  | 0.544 | 0.545 | 0.500 | 0.500 | 0.583 |
| 10 | 0.430 | 0.309 | 0.333 | 0.500 | 0.500 |
| 15 | 0.360 | 0.223 | 0.333 | 0.333 | 0.500 |
| 20 | 0.290 | 0.148 | 0.217 | 0.233 | 0.239 |
| 30 | 0.245 | 0.066 | 0.200 | 0.233 | 0.233 |
| 40 | 0.200 | 0.051 | 0.200 | 0.233 | 0.233 |
| 44 | 0.200 | 0.056 | 0.200 | 0.225 | 0.233 |
| **RMSE(AR5-44)** | -- | **0.145** | **0.069** | **0.066** | **0.063** |

**Two-sided gate:**
- (a) PRESERVE the knee (AR<=8): the front-loss band is full through AR<=band_W, so at fixed base
  the knee is byte-identical with the flag OFF (AR8 0.500 both, base 2.9). vs the MC banked 0.545
  the Knudsen knee is 0.500 (delta -0.044) -- slightly softer than MC's near-perfect knee but well
  inside the honest agreement band; co-tuning base to 2.6 restores AR8 0.583 (sim/exp 1.07). **PASS.**
- (b) SUSTAIN the deep floor (target sim/exp @ AR30-40 from ~0.25 toward >0.6): Knudsen sim/exp is
  0.82 @ AR30 / 1.00 @ AR40 (front-OFF) and 0.95 / 1.17 (front-ON). **PASS (far exceeds >0.6).**

**Verdict -- both gate sides PASS, but the credit is split and honest:**
1. The DOMINANT fix is the transport integrator: MC -> deterministic Knudsen drops RMSE(AR5-44)
   from **0.145 to 0.069** and sustains the floor to ~0.20 (matching the experiment). The banked
   "sim floor collapses 11x while the experiment sustains" narrative was substantially an MC
   deep-floor sampling artifact -- consistent with [[reconcile-craig-into-petch]] ("petch-MC
   under-samples deep floor in static eval"). The sustained tail is real (at AR40 the front is at
   z=14.5 um, far from the z=0 domain floor) and is carried by the ion-sputter floor (Ysp) plus the
   deterministic neutral tail, not a domain-bottom artifact.
2. The passivation-linked front-loss is a REAL but SUBDOMINANT, knee-safe improvement: cleanly
   ablated at fixed base 2.9 it lifts the mid+deep tail toward the experiment (AR10 0.333->0.500,
   AR30 0.200->0.233) and improves RMSE 0.069 -> 0.066 without touching the knee; co-tuning the base
   reaches 0.063. It does NOT make-or-break the gate (Knudsen alone passes) and it cannot fix the
   one residual shape error -- the AR10-20 slope is slightly too steep vs the experiment's gentle
   roll, and because the front-loss reduces loss monotonically with depth it lifts the deep tail
   (already on target) as much as the mid, introducing a small AR10-15 bump when pushed. So the
   hypothesis is directionally CORRECT (inert passivated walls -> flatter tail) but is a refinement,
   not the missing physics; the missing physics was Monte-Carlo variance, not chemistry.

Figure: `viz/deboer_passivation_wall_loss.png` (experiment vs MC-banked vs Knudsen OFF/ON/co-tuned).

## Status of the opt-in code

- `src/petch/params.py`: `Flags.floor_charge_throttle` and `Flags.knudsen_front_loss` (both default
  False -- default behavior unchanged; pytest 26/26 green). Front-loss params: `knudsen_front_band_W`
  (=8.0 feature-widths fresh band), `knudsen_passive_frac` (=0.5), `knudsen_front_ar_pass` (=15.0).
- `src/petch/threed.py`: `_apply_floor_charge_throttle` + hook in `mc_flux_3d_coupled`;
  `knudsen_front_loss` wired into `mc_flux_3d_knudsen` (per-slice loss array to `knudsen_face_flux`).
- `src/petch/knudsen.py`: `_front_loss_scale` (per-slice passivation-linked scale) +
  `conductance_profile` accepts a scalar OR per-slice `wall_loss_scale`.
- `src/petch/charging_general.py`: `floor_charge_throttle_profile` + `_THROTTLE_*` table.
Both flags stay available as opt-in physics experiments; neither is on in the default de Boer config.
The actionable de Boer takeaway is the transport choice: `neutral_transport="knudsen"` (deterministic)
should be preferred over `"mc"` for HARC deep-floor ARDE, where MC variance collapses the floor.
