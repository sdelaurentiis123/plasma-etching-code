# GPU speed campaign — scoped, iterative run plan

Goal: "lightning fast" — beat ViennaPS-GPU wall by a wide margin while staying as-accurate
(3D ARDE rmse ~0.05-0.08, depth within MC noise) AND differentiable. Every wave:
**implement -> correctness/parity run -> speed run (profiler) -> decision gate -> commit+push**.
Run `profile_steps_3d.py` after every change; it is the feedback loop.

Ruler (fixed): `head_to_head_3d.py` (ours vs ViennaPS-GPU, matched d=6 hole). Baseline to beat:
ours ~10.18s vs ViennaPS ~11.6s. Accuracy ruler: ViennaPS-3D hole ARDE (rmse ~0.05-0.08).

Box discipline: spin ONE box, run as many waves as possible that session, kill at end
(`vastai destroy instance <id>`). Don't idle a box.

## KEY RESEARCH FINDING (reshapes the strategy)

**Warp `wp.mesh_query_ray` is a SOFTWARE stack-based BVH on CUDA cores — NOT OptiX, NOT RT-cores**
(verified in Warp native `mesh.h`/`bvh.h`). Consequences:
- The "RT-core-class 200M rays/s" was a software walker. Flux is NOT hardware-raytrace-limited ->
  real headroom inside Warp (cuBQL backend, wavefront restructuring) WITHOUT leaving autodiff.
- OptiX SER does NOT apply to current code. True RT-core traversal = an OptiX rewrite that
  **sacrifices Warp differentiability** -> last resort only (Wave 7), defer.
- Strategic read: items that make the loop GPU-resident + take reinit off the critical path
  (kill `.numpy()` syncs, GPU fast-sweep reinit, extension-velocity lazy reinit) are the cleanest
  path to decisively beating ViennaPS-GPU, and all preserve the differentiable edge.

Ranked action list (from the dig), mapped to waves below:
1 kill per-step `.numpy()` readbacks  | 2 GPU fast-sweep reinit | 3 extension-vel + lazy reinit |
4 reuse form-factor matrix across coverage | 5 radiosity = neutral default | 6 Anderson on coverage |
7 CUDA graph capture | 8 wavefront neutral loop | 9 cuBQL BVH (1-line) | 10 F/O on streams |
11 mixed-precision traversal | 12 OptiX+SER (defer, breaks autodiff).

---

## Wave 0 — Baseline + instrument (run FIRST; no code change)
Runs: `profile_steps_3d.py` (per-step breakdown, MC vs radiosity neutral, depth delta) +
`head_to_head_3d.py` (fix the ruler this session).
Gates: **G0a** radiosity faster AND depth delta <0.5um? -> radiosity becomes neutral default (Wave 1
is then already done; items 4+5). **G0b** confirm bottleneck ranking (expected reinit>flux>mesh>host)
to order the rest.

## Wave 1 — Radiosity as neutral default  [items 4,5 — biggest EASY flux win, already coded]
The coupled-MC path re-traces full neutral MC 8x/step (n_fp=4 x {F,O}). `mc_flux_3d_radiosity`
already builds the form-factor matrix ONCE then does cheap sparse re-solves (item 4 for free) AND is
noise-free at the deep floor (fixes de-Boer accuracy). Promote it to the default neutral path for
speed configs.
Runs: covered by `profile_steps_3d.py` G0a (faster?) + an ARDE run (accuracy-neutral?).
Gate **G1:** radiosity faster AND ARDE rmse held -> set default. Else keep MC, do item 6 (Anderson)
in Wave 5 instead.

## Wave 2 — Kill the host syncs + GPU reinit  [items 1,2,3 — the 40% reinit + loop-resident win]
This is the structural win. Three coupled changes:
- (1) replace per-step `.numpy()` readbacks (`fl.numpy()`, `fi.numpy()`, CFL `np.max(V)`, coverage
  residual) with on-device reductions into 1-element arrays. Prereq for graphs.
- (2) GPU parallel Fast Sweeping reinit in Warp (Detrixhe-Gibou-Min diagonal/red-black Godunov,
  fixed 8 sweeps, enforces |grad|=1) replacing CPU skfmm. FIXES the masked-front |grad|=1.32 bias
  that made the old `reinit_gpu` unshippable.
- (3) extension velocity `grad F . grad phi = 0` (Adalsteinsson-Sethian) via Warp Jacobi sweeps
  replacing nearest-face `extend_velocity_gpu` -> |grad phi|=1 preserved -> lazy reinit (every K
  steps) becomes SAFE -> reinit 40% -> ~8%.
Runs:
- NEW `reinit_correctness_3d.py` — fixed deep-hole phi: |grad phi| error in band vs skfmm-narrow +
  full-etch depth match. MANDATORY before trusting GPU reinit.
- NEW `extension_lazy_3d.py` — sweep K (reinit_every) with the proper extension; depth drift vs K=1.
- `profile_steps_3d.py` — reinit % and host % must drop.
Gate **G2:** depth matches skfmm-narrow within MC noise AND reinit+host wall drops. Else revert to
skfmm-narrow (already SOTA). Don't ship faster-but-wrong reinit.

## Wave 3 — CUDA graph capture the step  [item 7 — only pays AFTER Wave 2 removes syncs]
`wp.ScopedCapture` / `wp.capture_launch` the per-step kernel sequence; static-unroll the fixed
coverage iters + CFL substeps inside the capture; `wp.capture_while` for data-dependent counts (no
host sync). Composes with `wp.Tape`.
Runs: `profile_steps_3d.py` (launch overhead gone) + accuracy-neutral depth check.
Gate **G3:** wall drops, depth unchanged. (Blocked until Wave 2 makes the loop device-resident —
any skfmm/scipy call breaks capture.)

## Wave 4 — Flux throughput inside Warp  [items 9,8 — no OptiX, keeps autodiff]
- (9) cuBQL BVH backend: `wp.Mesh(..., bvh_constructor="cubql")` — 1-line try. ~3x on large meshes;
  our meshes are ~2500 faces so VERIFY the win at our scale + check gradient support before committing.
- (8) if flux still dominates: wavefront restructure of the neutral RR loop (per-bounce kernels +
  SoA path state + compaction) — HPG2013 measured 36-221% on incoherent multi-bounce. Bigger lift.
Runs: NEW `bvh_backend_3d.py` (cubql vs default: flux wall + per-face parity + gradient check) +
`profile_steps_3d.py`.
Gate **G4:** faster AND per-face flux parity AND gradients intact.

## Wave 5 — Coverage iteration count  [item 6 — Anderson acceleration]
Anderson acceleration (window m=2-5, thin-QR LS) on the coverage fixed point: ~8 -> ~3 iters.
Differentiate the converged point via IFT (gradients unchanged). Complements Wave 1's per-iter cost
cut with an iteration-count cut.
Runs: NEW `coverage_anderson_3d.py` — iters-to-converge + depth/ARDE vs Picard. `profile_steps_3d.py`.
Gate **G5:** fewer iters, depth+ARDE held.

## Wave 6 — Re-validate accuracy (prove "as-accurate AND faster")
Runs: ViennaPS-3D hole ARDE re-run (rmse stays ~0.05-0.08) + `head_to_head_3d.py` (NEW wall vs 11.6s
= headline) + `validate_deBoer_radiosity.py` with periodic-y BCs for the trench.
Gate **G6:** rmse held + wall beats ViennaPS-GPU by target margin -> SHIP + update docs/memory.

## Wave 7 — OptiX + SER (DEFER; last resort, breaks autodiff)
Only if Waves 1-6 don't close the gap. Rewrites flux against OptiX 8 `optixReorderThread` for true
RT-core traversal + SER (up to 2x on incoherent paths). COST: leaves `wp.mesh_query_ray`, loses
Warp's automatic adjoint -> must hand-write the differentiable intersection/tally adjoint. Stage last.

---

## Parallel accuracy track (fold into a box session; not speed)
periodic-y BCs for trench radiosity (proper de Boer) | re-validate redeposition (shared the fixed
into-gas normal bug) | surface charging (biggest missing real-wafer physics; needs wafer data).

## Dependency graph
```
Wave 0 (profiler + ruler) ─ G0 ranks ─┬─> Wave 1 radiosity default (items 4,5) ─ G1 ─┐
                                       │                                              │
                                       └─> Wave 2 kill syncs + GPU reinit (1,2,3) ─ G2 ┼─> Wave 3 CUDA graphs (7) ─ G3 ─┐
                                                                                       │                                │
   Wave 4 cuBQL/wavefront (9,8) ─ G4 ───────────────────────────────────────────────┘                                │
   Wave 5 Anderson (6) ─ G5 ─────────────────────────────────────────────────────────────────────────────────────────┼─> Wave 6 re-validate ─ G6 ─> SHIP
                                                                                       Wave 7 OptiX+SER (12) — defer ───┘
```
One box session = Wave 0 + 1 + (2 if research lands the FSM kernel) + cheap tries (4 cuBQL), then
Wave 6 re-validate, then kill. Graphs (3) + Anderson (5) + wavefront (8) are follow-on sessions.
