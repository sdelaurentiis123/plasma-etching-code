# Night log — autonomous session

Summary of what changed this session (newest first). Full detail in `FINDINGS.md`; explainers in `docs/`.

## Speed: now FASTER than ViennaPS-GPU (was 2x slower)
- Same-engine d=6 hole ~9um: **ours 10.18s vs ViennaPS-GPU 11.6s** (was 22.6s -> **2.2x self-speedup**).
- The earlier "23x faster" was a GPU-vs-CPU artifact -- corrected.
- SOTA research found the root cause: **dense-vs-sparse, not CPU-vs-GPU.** ViennaPS's level set is CPU+sparse
  (HRLE) and does ZERO global reinit per step; our dense skfmm-every-step reinit (42%) was self-inflicted.
- Shipped: **narrow-band reinit** (skfmm `narrow=`, ~3-5x, EXACT in band, zero accuracy cost) +
  **GPU advection** (Warp kernel, phi on-device across all CFL substeps, matches numpy to 1e-6).
- The speedups are pure -- verified accuracy-neutral. Profile now flux(48%)+reinit(40%) bound.

## Accuracy: all of ViennaPS's physics, 3D ARDE rmse ~0.05-0.08
- Replicated EVERY ViennaPS mechanism: Belen coverages, exact Russian-roulette weighted transport,
  coverage-dependent sticking, 1-neighbor flux smoothing. Flux smoothing was the key small-hole fix.
- 2D matched (rmse 0.016); 3D brackets ViennaPS (rmse ~0.05-0.08). Residual = smoothing-neighborhood
  calibration (edge-mesh vs disk-radius), a tuning detail, not missing physics. alpha knob added for it.

## Physics BEYOND ViennaPS (it omits these)
- **Full ion-energy-distribution (IED) yield integration**: mean / Gaussian(=ViennaPS) / **bimodal** sheath
  (real low-freq-bias distribution). Jensen bias <Y> < Y(<E>) verified (~0.4% at 100eV, 4.4% near threshold).
- **Etch-product redeposition**: product emitted from faces, re-sticks on lower sidewalls -> passivation ->
  **taper** (verified: taper 0.13 -> 0.26). flags.redeposition (default off).

## Research (all documented as HTML in docs/, both repos)
- Physics-constants dig, SOTA-speed dig, AND a REAL-EXPERIMENT dig. Three new explainer pages:
  `docs/experimental-validation.html`, `docs/performance.html`, `docs/physics-grounding.html`.
- Honest reality: ViennaPS validates against ONE dataset (Belen SEMs) by fitting fluxes per condition.
  REAL-wafer benchmark targets identified: Gomez/Belen 2004, de Boer 2002 cryo RIE-lag (Blauw Knudsen
  S_F~0.47, ER(AR)=1.0/0.43/0.29/0.20), Hoekstra-Kushner microtrench. Flamm/Donnelly F-Si chemistry.

## The honest frontier (to beat REAL wafers, not just ViennaPS)
1. **Surface charging** -- notching, AR-dependent ion deflection. Biggest gap. ViennaPS also lacks it.
2. Real bimodal IEDF tail sensitivity (we now have the IED machinery).
3. Flux-dependent reaction probability; fluxes are fit not predicted (feature-scale sim, not plasma-to-wafer).

## State
All committed + pushed (both repos). No GPU boxes left running. Next levers ready: charging (the
real-wafer differentiator), GPU iFIM reinit + radiosity neutral flux (push further past ViennaPS),
and running the de Boer ARDE real-wafer validation on a box.

## CAPSTONE: first REAL-wafer validation (honest gap found)
Ran our model vs the de Boer/Blauw cryo experiment (NOT ViennaPS). Our trench STALLS at AR~7.75 and
over-starves (~7x slower floor than the wafer at AR~8); RMSE 0.22 vs experiment. So: we match ViennaPS,
but do NOT yet match real wafers -- our ballistic transport over-predicts ARDE vs real Knudsen molecular
flow. Path (quantified): true Knudsen transport + condition calibration (lower betaE) + charging.
This is the honest answer to "real physics accuracy": not there yet against wafers; clear path documented.
