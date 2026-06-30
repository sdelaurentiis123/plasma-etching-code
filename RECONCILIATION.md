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

## Speed (RTX 3090, CUDA)

petch sub-second (0.98–1.30 s) vs ViennaPS-OptiX. The benchmark box reported 46–66× but that is
**inflated** — the box had only 6 weak vCPUs, so ViennaPS's CPU-side level-set (~40% of its time)
was starved (52–64 s vs the validated 19–25 s), and the run wasn't depth-matched. **Honest,
comparable speedup remains ~14×** (balanced box, depth-matched, as documented in the README).

## Follow-ups

- Warp-ify the DDA neutral gather (currently numpy) for CUDA speed — then `dda` is fast too.
- Calibrate the Knudsen floor sink; diagnose petch radiosity's over-correction.
- petch-DDA static ARDE at W=2 µm for a direct petch-vs-wafer scorecard line.
