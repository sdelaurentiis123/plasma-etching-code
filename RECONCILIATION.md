# Reconciling petch with Craig's plasma_sim (2026-06-29)

Absorbed Craig Xu Chen's `plasma_sim` (Apple-Metal DDA solver) physics into petch, CUDA-first,
then stress-tested where the two engines diverge and scored both against the de Boer wafer.

## What was ported into petch (from `plasma_sim`)

All as `neutral_transport` options, coupled to petch's Belen coverage fixed point:

- **`dda`** — deterministic discrete-ordinates grid-march (`src/petch/dda.py`). The clean,
  noise-free deep-AR transport. **Validated:** direct sky-view rolloff matches the analytic
  cylinder view factor (0.05 vs 0.059 at AR2), and the full ARDE reproduces ViennaPS (below).
- **`knudsen`** — 1-D molecular-flow conductance tail (`src/petch/knudsen.py`). Runs, but the
  floor reaction sink is under-weighted at small W → too flat. **Needs calibration.**
- **`radiosity` + `radiosity_solver='gmres'`** — matrix-free GMRES on petch's form-factor
  radiosity (Craig's better-conditioned solve). Radiosity still over-corrects (flat) — known.

The Belen SF6/O2 chemistry is the same model in both; petch's `belen.py` already subsumes it
(petch adds faithful MC ion reflection). petch's earlier-removed `cal_f` knob was NOT re-added
(intentional reversal on 2026-06-18).

## Numerics cross-validation (petch-MC vs Craig-DDA, matched chemistry)

`scripts/cross_validate_dda.py` — only the transport method differs. Headline: at a cheap 8k-ray
test budget petch-MC **starves the deep floor** (nr→0 by AR~5.5–6), while Craig's deterministic
DDA gives a clean rolloff. Max divergence: 0.60 (hole, AR4) / 0.76 (trench, AR6) — dominated by
MC under-sampling, not a model disagreement.

## Static ARDE — the DDA fixes it (`scripts/dda_vs_mc_arde.py`, `viz/dda_vs_mc_arde.png`)

Clean W=0.5 µm trench, instantaneous nr(AR)=V_floor/V_field:

| AR | petch-MC (200k rays, static) | **petch-DDA** | ViennaPS ref | de Boer wafer |
|----|------|------|------|------|
| 2  | 0.54 | 0.99 | — | — |
| 6  | 0.08 | 0.86 | — | — |
| ~8.6 | ~0.04 | **0.727** | **~0.73** | — |
| 10 | 0.03 | 0.65 | — | 0.43 (W=2 µm) |

- **petch-DDA ≈ ViennaPS** (0.727 vs 0.73 @ AR8.6) — the port reproduces the SOTA ballistic rolloff.
- **petch-MC's static deep-AR collapse is under-sampling** — even 200k rays starve the floor in a
  single static eval (the deep-floor solid angle ~1/AR). The DDA, being deterministic, doesn't.
  (Evolving etches accumulate over steps so MC is milder there — consistent with petch-MC tracking
  ViennaPS within ~0.1 in the evolving sweeps.)
- **Both ballistic engines sit ABOVE the wafer** (petch-DDA 0.65 / ViennaPS ~0.73 vs wafer 0.43 at
  AR10) — the structural gas-conductance/charging gap, not a petch bug.

## de Boer wafer scorecard (Craig's harness, Metal)

`bench/industrial_validation.py` (Craig's calibrated defaults, W=2 µm trench, AR [0,10,20,40]):
- DDA (`gmres`): RMSE **0.072**, open rate 1.301 µm/min (= Gomez 1.3, tuned)
- Knudsen: RMSE **0.048** ✅ (passes the 0.05 wafer gate)
- source: RMSE 0.072

Craig's engine matches the wafer because its defaults are calibrated to it (open rate tuned to Gomez,
knee within RMSE 0.05–0.07). petch (DDA or MC) sits in the *ballistic* regime near/above the wafer
knee; the de Boer "process" knobs do NOT bring it down to the wafer's deep-AR floor (see below).

## Speed (RTX 3090, CUDA) — real depth-matched test

petch is **GPU-resident and CPU-independent (~1 s)**; ViennaPS's wall is CPU-bound (its level-set
advection is ~40% CPU), so the *ratio* depends entirely on the ViennaPS host CPU:

| ViennaPS host | ViennaPS wall | petch wall | speedup |
|---|---|---|---|
| 32-core @ 4.4 GHz (fast, **depth-matched**) | 7–12 s | 0.9–1.2 s | **7.3–10.4× (median 8.3×)** |
| mid/typical CPU (documented) | 19–25 s | ~1.3–1.8 s | ~14× |
| 6 weak vCPUs (CPU-starved) | 52–64 s | ~1 s | 46–66× (**artifact**) |

So the honest floor is **~8×** (ViennaPS on a fast many-core CPU, depth-matched to 3–7%); the
advantage only grows on weaker hosts. The 46–66× from the weak-vCPU box was the inflation trap and is
NOT a real gain. `scripts/vps_sweep.py` (rate range widened so it actually depth-matches).

## petch-DDA at the de Boer width (W=2 µm, CUDA)

With the Warp gather, DDA is **fast on CUDA (~0.1 s per static eval)**. petch-DDA static ARDE at W=2 µm:

| AR | petch-DDA nr | de Boer wafer |
|----|------|------|
| 10 | 0.58 | 0.43 |
| 20 | 0.23 | 0.29 |

It straddles the wafer (above at AR10, below at AR20), RMSE ~0.11 with default ballistic params
(Craig's *calibrated* engine: 0.072). **Tested 2026-06-30: the de Boer "process" knobs do NOT close
the gap.** A narrow (0.8°) IADF funnels ions to the floor → *gentler* knee (0.77 @ AR10, wrong way);
3× etchant starvation barely moves it (0.69). The wafer's flat high-AR tail (0.29 @ AR20, 0.20 @ AR40)
is the **structural frontier** — it needs missing physics (surface charging / a coverage-independent
etch channel), not a transport, IADF, or flux knob. So petch-DDA reproduces the clean *ballistic* ARDE
(matching ViennaPS); it does not, with current knobs, reach the real wafer's deep-AR floor. (The figure
is the right panel of `viz/dda_vs_mc_arde.png` / `dda_vs_mc_arde.npz`.)

## Knudsen tail: bug found, fixed, calibrated to the wafer (2026-07-02)

The flat-profile bug in petch's Knudsen port is fixed: mesh-face centroids sit ON the zero contour, so
floor faces' slice index rounded into the gas-free slice below the last valid one — the reaction sink
landed on an invalid slice and never entered the tridiagonal solve. Fix: snap every face to its nearest
gas-carrying slice before accumulating (`knudsen.py`; plasma_sim's band-cell sampling is immune, which is
why Craig's worked). Verified A/B vs plasma_sim's `_neutral_flux_knudsen`: petch now decays properly.

**Benchmarked to the real wafer** (static nr(AR) at the de Boer points, W=2 µm, petch chemistry,
sweeping the one knob `knudsen_wall_loss_scale`):

| wls | nr @ AR [0,10,20,40] | RMSE vs wafer |
|---|---|---|
| 1.0 | 1.0 / 0.54 / 0.29 / 0.14 | 0.062 |
| **1.3 (new petch default)** | 1.0 / 0.47 / 0.25 / 0.11 | **0.053** |
| 1.85 (plasma_sim default) | 1.0 / 0.38 / 0.19 / 0.08 | 0.082 |

Wafer = 1.0 / 0.43 / 0.29 / 0.20; gate ≤ 0.05. **RMSE 0.053 is a near-pass and ~2.4× better than any
ballistic config (~0.13)** — the Knudsen conductance channel is the first petch transport that bends
toward the real wafer. Honest labels: (a) this is a 1-knob *calibration*, not a prediction (same status
as plasma_sim's 1.85 on its own chemistry); (b) the residual is *structural in shape* — the single
floor-sink model decays linearly with depth while the measured tail flattens (AR40: 0.11 vs 0.20), so
no value of the knob passes the gate. The identified next physics: a self-limiting (starvation-coupled)
floor sink or distributed sidewall loss, which is what produces a flattening tail.

## Follow-ups

- ✅ Warp-ify the DDA neutral gather — done (`_dda_gather_kernel`, ~14× over numpy on CPU, 0.1 s/eval on CUDA).
- ✅ Tested petch-DDA with de Boer process knobs (narrow IADF, etchant starvation) at W=2 — they do NOT
  close the wafer gap (narrow IADF goes the wrong way; the deep-AR floor is structural).
- ✅ Knudsen floor sink: bug fixed (contour-snap), calibrated to the wafer (wls=1.3, RMSE 0.053 — near-pass,
  2.4× better than ballistic). Remaining shape residual is structural (linear decay vs flattening tail).
- Diagnose petch radiosity's over-correction.
- ✅ Tested the self-limiting-sink hypothesis (`knudsen_sink='field'`: clamp the sink at the field
  coverage so it stops growing as the floor starves). Result: RMSE 0.050 at wls=2.0 — marginally better
  than local (0.053) but the tail *shape* is unchanged (~2× decay per AR-doubling vs the wafer's ~1.45×).
  **Falsified**: the steepener is not the neutral-sink feedback. The concentration profile is ~1/(1+kAR);
  nr decays faster because of the multiplicative ion-side decay (m_i drops with AR; starved-regime rate
  ∝ θ_F ∝ conc × ion term).
- The real frontier, re-sharpened: the wafer's flat tail needs something that keeps the *deep-floor etch
  efficient* at high AR — the coverage-independent etch channel from FINDINGS (a reaction-limited floor
  term that doesn't die with θ_F), NOT a transport/sink knob, and NOT charging (charging cuts deep ion
  flux and would steepen the tail further). (The old `cal_F` knee-fix was removed; etchant starvation
  alone is insufficient given ion-enhanced floor etching.)
