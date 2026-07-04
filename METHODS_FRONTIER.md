# Frontier methods for the GPU differentiable charging engine (research 2026-07-04)

One-line finding: **no public differentiable feature-scale plasma-etch simulator exists** (commercial
TCAD is finite-difference/DoE-calibrated, not differentiable) -> greenfield. The transferable methods
are mature in three adjacent fields: **differentiable rendering** (our MC-charging flux IS a transport
integral over a moving level-set surface), **learned Krylov preconditioning** (our Poisson inner loop),
and **implicit/fixed-point differentiation** (our steady-state charging equilibrium). All 15 anchor
citations verified.

## Top 3 to adopt first

### #1 Path-Replay Backprop + boundary-term gradients for the MC-charging flux
Our MC particle-in-field charging over a moving level-set is mathematically identical to
**differentiable rendering of an SDF** -- a solved, production problem in graphics.
- Path-Replay Backpropagation (Vicini/Speierer/Jakob, SIGGRAPH 2021, ships in Mitsuba3/Dr.Jit):
  differentiates an MC transport integral in LINEAR TIME, CONSTANT MEMORY by replaying the random walk
  with the same RNG seed backward instead of taping every event. https://rgl.epfl.ch/publications/Vicini2021PathReplay
- Reparameterized boundary term (Loubet 2019; Bangaru 2020; path-space implicit-surface SIGGRAPH 2024):
  naive autodiff SILENTLY DROPS the surface-motion gradient d(charge map)/d(surface position) at the
  moving level-set wall; the change-of-variables recovers it correctly + low-variance.
- Grazing-incidence variance control (diff. radiation transport 2026, arXiv:2605.06779): pathwise AD
  explodes as 1/cos(beta) at grazing walls (our AR>20 floor-collapse regime); primal-preserving
  stop-gradient when |v.n|<0.2 kills outliers ~6 orders, forward pass unchanged.
Why: constant-memory autodiff = the difference between OOM and runs-on-one-GPU; the boundary term is the
biggest correctness upgrade for the charging gradient. Ships today, maps 1:1 onto our level-set + MC.

### #2 Implicit / deep-equilibrium differentiation of the steady-state charging
Steady-state charging (Poisson + surface-charge balance -> self-consistent potential) IS a fixed point.
Differentiate via the implicit function theorem -> exact gradients in O(1) MEMORY from ONE linear solve,
vs unrolling hundreds of relaxation sub-iterations through the tape (the memory blocker).
- FNO-DEQ (arXiv:2312.00234); discrete-adjoint linear-solve pattern (JAX-BTE npj 2025; Ceviche custom-vjp).
Wrap the SOR relaxation as an implicit layer; composes with #1 and Warp's tape.

### #3 Neural-preconditioned CG to replace red-black SOR Poisson (HYBRID, not surrogate)
The Poisson solve is the per-timestep bottleneck. A ~25k-param matrix-free CNN preconditioner cuts CG to
**~21 iters vs ~616 for plain CG (~30x iteration reduction, 3-13x wall-clock)** on EVOLVING MIXED-BC
geometry -- exactly our setting. arXiv:2310.00177 (Lan et al.); provable-convergence framing
arXiv:1906.01200; learned-SPAI GPU fallback (SpMV-only, Warp-friendly) arXiv:2510.27517.
CRITICAL: keep the outer CG driving true residual to tolerance -> speed WITHOUT the near-exponential
rollout error that kills bare FNO/DeepONet solvers in a time loop. Preconditioner, never surrogate.

## GPU engineering patterns to copy (WarpX, arXiv:2101.12149, verified)
- Cell-sort particles every ~4 steps -> **7.5x** on charge deposition (cache locality). Warp radix sort / HashGrid.
- Direct global `wp.atomic_add` into the charge grid (skip private buffers) -- fast once sorted.
- Fuse field-gather + push into one kernel -> ~1.6x memory + ~25% speed.
- Everything device-resident; never `.numpy()` in the hot path (managed-mem overhead <0.2%).
- ~1e6-1e8 particles is all charging needs -> memory-comfortable on one GPU; store state in fp16/bf16 SoA.
- Warp tile programming (1.5+): shared-mem cooperative deposit + on-device cuFFTDx -> FFT-Poisson without
  leaving the kernel (likely simpler/faster than multigrid on our tiny micron grids).

## Differentiable-through-discrete-events (the hard part, for collisions/surface reactions)
- Adjoint DSMC (Caflisch/Yang, arXiv:2207.11579, 2026 arXiv:2603.20946): the ONLY mature method for
  gradients through discrete stochastic collision/scatter/neutralize events -- hybrid pathwise +
  score-function estimator. Exact recipe for gradients w.r.t. cross-section/sheath/yield params.
- Estimator taxonomy (Mohamed JMLR 2020, arXiv:1906.10652): reparameterize continuous scattering angles;
  score-function for discrete branches; control variates to tame variance.
- JAX-in-Cell (arXiv:2512.12160): fully differentiable PIC via plain autodiff -- but COLLISIONLESS
  (confirms field-solve+pusher are "free"; collisions are the value-add).

## Chemistry / multiscale coupling (the Graves moat, method side)
- Neural Master Equation (Nath, Vella, GRAVES, Mesbah, npj 2025, 10.1038/s41524-025-01677-4): dP/dt=W(theta)P,
  W a NN from atomistic data, probability-conservation structurally enforced. Demonstrated Si ALE + RIE.
  Replace hand-tuned per-site kinetics with an NME block; gradients flow from final-profile objective to
  kinetic params. No code -> reimplement in Warp/JAX. THIS is the differentiable atomistic->feature bridge.
- beta-VAE for the EAD interface (arXiv:2109.01406): compress 1800-dim energy-angle dists to 2-D latent,
  0.3-0.9 ms/eval, replaces TRIDYN at runtime -- differentiable BC knob.
- DP-GEN active learning (arXiv:2203.00393): ensemble-variance trigger for on-the-fly surrogate refinement.
- CAUTION: MD differentiability is fragile (chaotic gradient blow-up); make the SURROGATE the
  differentiable bridge, don't backprop through long reactive-MD trajectories.

## Calibration / inverse design (the payoff of differentiability)
- Warp IS our autodiff engine (source-to-source reverse-mode over GPU kernels via wp.Tape). Record
  level-set advect + MC deposit + flux kernels on one tape; tape.backward() -> d(loss)/d(cal_F, IADF, sheath).
- FDTDX (arXiv:2412.12360): closest architectural sibling -- GPU JAX FDTD + AD + inverse design with a
  TIME-REVERSIBLE gradient (recompute states backward, don't store every step -> the memory fix for long
  etch marches). Verified 10x faster than Meep, 415x faster than Ceviche at 288M cells. Our "GPU+diff
  beats CPU incumbents" precedent.
- Differentiable LArTPC detector sim (arXiv:2309.04639): proof that entangled stochastic-sim params can be
  co-calibrated by gradient, not one-at-a-time DoE. Directly = co-calibrating our charging+transport+rate params.
- MMD/kernel-score calibration for inexact stochastic sims (arXiv:2411.05315): fit DISTRIBUTIONS of
  profiles to CD-SEM with UQ that accounts for model misspecification (our high-AR frontier).

## Genuinely new (2024-2026) worth watching
- Indirect Neural Corrector (NeurIPS 2025, arXiv:2511.12764): inject learned correction INTO the governing
  equation (source/velocity term), not onto the output state -> error amplification bounded (up to +158% R2).
  The how-to for a stable+differentiable hybrid: feed any neural corrector into the level-set velocity/flux.
- Path-space differentiable rendering of level-sets (SIGGRAPH 2024): d(flux)/d(surface) correct+low-variance.
- PDE-Refiner / PCFM / "predict change not states": stable long rollouts, hard-constraint projection
  (charge neutrality, etched volume), dphi/dt learning -- for if/when we add a neural etch-front stepper.

Hype flags (don't build on): generic FNO plasma surrogates, image-to-image profile CNNs, chamber-level
etch-rate regressions -- black boxes with no transport gradients / wrong altitude.
