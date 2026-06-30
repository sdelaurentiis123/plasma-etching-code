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

Craig's engine matches the wafer because its defaults are calibrated to it. petch matches the wafer
with its de Boer *process* params (narrow IADF / etchant-starved) — documented separately, not
re-derived here (the W differs, and petch-DDA's numpy gather is too slow for W=2 µm deep ARs).

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

It straddles the wafer (above at AR10, below at AR20), RMSE ~0.09 with default ViennaPS-regime params
(Craig's *calibrated* engine: 0.072). petch matching the wafer needs its de Boer process params.

## Follow-ups

- ✅ Warp-ify the DDA neutral gather — done (`_dda_gather_kernel`, ~14× over numpy on CPU, 0.1 s/eval on CUDA).
- Calibrate the Knudsen floor sink; diagnose petch radiosity's over-correction.
- Re-run petch with de Boer *process* params (narrow IADF / etchant-starved) through the W=2 DDA scorecard
  to show petch can match the wafer (not just the ViennaPS-regime defaults).
