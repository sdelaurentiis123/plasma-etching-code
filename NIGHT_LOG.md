# Night log — autonomous session

Summary of what changed this session (newest first). Full detail in `FINDINGS.md`; explainers in `docs/`.

## GPU flux smoothing: profiled the flux internals, modest win (honest)
- After FSM+warm-start the loop is flux-bound (~70%). Profiled the flux INTERNALS
  (flux_internal_profile_3d, M1 -- host ops are device-independent): on a deep 54k-face mesh the
  cost is HOST, not GPU traversal -- **smooth_flux 16.7ms/call** (np.add.at scatter over 80k edge
  pairs) x3/step, source-gen 1.5ms, coverage 1.1ms. So the lever was host-op offload, NOT a faster
  BVH (cuBQL would not have helped).
- Shipped `smooth_flux_gpu` (Warp atomic-add scatter `_smooth_scatter`), matches numpy to 3e-7
  (validated on the Warp CPU backend), auto-on under CUDA. **HONEST: only ~5% overall** (same-box A/B:
  GPU 3.07s vs CPU 3.20s, flux 1.92 vs 2.13, depth identical) -- the per-call host->device transfers
  (pairs/weights/src/num/den, 80k each, rebuilt every call) eat most of the savings, and averaged over
  the run the mesh is smaller than the 54k deep-mesh profile. Kept (clean accuracy-neutral win).
- Real remaining flux lever = the per-step host orchestration itself (array allocs, transfers,
  .numpy() syncs spread thin) -> CUDA graph capture, which needs a device-resident loop (no mid-step
  .numpy()). That's a refactor, not a one-liner. Next target after flux: mesh (marching-cubes CPU, ~18%).
- Could ~2x the smoothing win by caching pairs/weights on-device ONCE per step (shared across the 3
  smooth calls) instead of per-call -- deferred (moderate plumbing, ~5% more).

## Warm-start coverage: accuracy WITHOUT speed loss (the tradeoff was an artifact) -- VALIDATED
- Q: "no way to get the accuracy without losing speed?" A: yes -- the 8 neutral MC traces/step were
  8 COLD restarts of the same coverage fixed point. The geometry (transport visibility) is fixed
  within a step; only sticking changes. Seed each step from the PREVIOUS step's coverage (front moves
  <1 cell/step) -> same converged fixed point in 1-2 iters not 4. `warm_start_coverage` flag.
- **ACCURACY (clean 30k-ray/dx0.25 box, non-bottoming):** WARM n_fp=1 (40 neutral launches) matches the
  cold n_fp=8 converged truth EXACTLY (mean/max 0.000um). Same-budget COLD n_fp=2 drifts 1.57um.
  Warm n_fp=1 is MORE converged than the old cold n_fp=4 (which still under-etches). Accuracy-neutral
  by construction. (First CPU check at 8k rays was noise-dominated/inconclusive -- needed the clean regime.)
- **SPEED (combined w/ FSM reinit):** NEW (warm n_fp=1 + fsm GPU) vs OLD default (cold n_fp=4 + skfmm CPU):
  **7.15s -> 3.41s = 2.09x**, depth identical 20.00um. Loop now flux-THROUGHPUT-bound (72%); reinit 3.6%.
- Next lever: flux throughput (cuBQL BVH 1-line, wavefront, or move smoothing/coverage host-ops to GPU).
  Hardened: warm-start skips its KDTree seed if centroids non-finite (CFL-blowup safe).

## GPU campaign Wave 2: FSM reinit kills the reinit bottleneck (RTX 3090, validated + box killed)
- **GPU Jacobi Godunov-Eikonal reinit (`reinit_fsm`, `reinit_method='fsm'`)** replaces CPU skfmm AND
  the biased PDE `reinit_gpu`. The Godunov solve enforces |grad phi|=1 EXACTLY at the fixed point;
  the monotone min-update is unconditionally stable (no forward-Euler blowup).
- **Correctness (validated on the Warp CPU backend, zero box cost):** |grad phi| mean 0.988 / p95 1.049
  (the old PDE reinit was 1.32 -- masked-front bias GONE); distance vs skfmm max 0.61 cells; contour
  preserved 100%; full-etch depth delta 0.000um.
- **Speed (RTX 3090, d=6 hole, 30 steps, 30k rays):** reinit line **4.03s -> 0.07s = 57.9x**; overall
  **16.18s -> 12.30s = 1.32x**, depth identical (18.00um, accuracy-neutral). Reinit went 25% -> 0.6%
  of the loop. The ~40% reinit bottleneck is eliminated; loop is now flux-dominated (~62%).
- **HONEST: radiosity-as-default REFUTED for speed.** Profiler measured radiosity 0.48x (SLOWER:
  34.6s vs 16.5s MC) at this moderate AR -- the form-factor build (F x 64 rays) + host scipy sparse
  solve (160 matvecs) costs more than 8 GPU neutral traces. Radiosity's speed win is HARC-only (where
  MC bounce counts explode); it stays the accuracy/deep-floor tool, not the default. Depth identical.
- **n_fp exposed via `par['n_fp']`** (coverage fixed-point iters; each = 2 neutral launches). Sweep:
  n_fp=3 borderline-safe (depth delta 0.25um, 1.18x), n_fp<=2 drifts (1.75-4.25um). Default stays 4.
- Next lever (own wave, accuracy-sensitive): cut the flux work -- warm-start coverage across steps /
  Anderson acceleration to drop n_fp, the now-dominant 62%. Box destroyed; no instances running.
- Key research finding: Warp `mesh_query_ray` is a SOFTWARE BVH (not OptiX/RT-cores) -> headroom
  inside Warp (cuBQL, wavefront) without losing autodiff; OptiX+SER is last-resort (breaks autodiff).
  Full 7-wave plan in `GPU_CAMPAIGN.md`.

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

## RADIOSITY built — the de Boer gap was MC under-sampling; we now BRACKET the real wafer
- de Boer betaE sweep: lowering betaE does NOT fix the over-starvation (RMSE 0.19-0.25, ~0.13 at AR10 vs
  wafer 0.43). Ray-count test: 30k->120k rays lifted AR10 rate 0.129->0.205 -> it's MC UNDER-SAMPLING of
  the deep floor, not a transport-model error.
- Built **deterministic radiosity neutral flux** (`mc_flux_3d_radiosity`, flags.neutral_transport='radiosity'):
  form-factor matrix by single-bounce ray casting + exact multi-bounce linear solve -> the deep floor gets
  its true conductance flux with NO under-sampling. Same build the SOTA research flagged as the speed lever.
- Result vs de Boer wafer: radiosity reaches AR=33 (vs MC's 10) and lands in the gentle-ARDE regime, BUT
  over-corrects (1.0 at AR10 vs wafer 0.43; 0.47 at AR20 vs 0.29). So **MC too steep, radiosity too flat,
  the real wafer is between** -- bracket-then-tune, same pattern that closed the ViennaPS gap.
- NEXT: tune/debug the radiosity over-recirculation (or blend MC<->radiosity) to nail the de Boer ARDE.
  Mechanism is right (deterministic conductance, differentiable); calibration isn't yet.

## RADIOSITY VALIDATED correct (the de Boer over-feed was a finite-trench artifact)
- Diagnosed the radiosity over-feed: face-emission rays used unoriented normals AND the de Boer TRENCH is
  open-ended in y (Ly=5um) -> rays escape the trench ends -> no lateral confinement -> no conductance.
- Fixed normals (_gas_normals, into-gas). On a CLOSED hole, radiosity m_F now DECAYS correctly with depth:
  0.656 (top) -> 0.256 (AR3) -> 0.031 (AR8) -> 0.003 (floor AR12) = it captures the Knudsen conductance
  EXACTLY (vs MC's under-sampling). **The radiosity build is CORRECT.**
- So: the de Boer (real ~infinite trench) needs PERIODIC-y BCs; radiosity is validated for closed geometry
  (holes) and should improve the ViennaPS-3D HOLE ARDE (deep floor no longer under-sampled).

## Clear next steps (documented, ready to execute)
1. Periodic-y BCs for trench radiosity -> proper de Boer trench comparison.
2. Re-run ViennaPS-3D hole ARDE with radiosity (does exact conductance match/diverge from ViennaPS?).
3. Re-validate redeposition (it shared the normal bug, now fixed).
4. Surface charging -- the biggest missing real-wafer physics.
