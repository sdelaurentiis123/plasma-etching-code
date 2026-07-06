# petch frontier loop — autonomous research/build charter

Standing autonomous loop (authorized 2026-07-06). Each cycle pushes one frontier of the open,
GPU, differentiable, feature-scale etch simulator, benchmarks it honestly against published data,
visualizes it, smoke-tests it, documents it, and commits. Runs across wake-ups; survives compaction.

## First principles (the whole thing in one line)
Etch is electrons + bonds evolving forward under ion/neutral/photon flux. We model that transport +
surface kinetics FAST (GPU), DIFFERENTIABLE (autodiff calibration/inverse design), at FEATURE SCALE,
and we GATE every claim against published data. The moat = open + GPU + differentiable + charging +
reactor↔atom coupling, which no open tool has.

## Cycle structure (every iteration)
1. **Pick target** from the queue (below) or a newly discovered gap.
2. **Research (literature-first).** Depth search, primary sources, working URLs. Refute before build.
3. **Implement / prototype.** CPU first for mechanism proof; box (CUDA) for fine-grid/long runs.
4. **Benchmark + GATE.** Every claim = script + published-data comparison + explicit PASS/fail.
   No claim ships without a number and a gate. Negative results are results — document them.
5. **Data viz.** A figure per cycle (into `viz/`, and `docs/*.html` when it's a headline).
6. **Smoke test.** `pytest` the touched paths; confirm default behavior unchanged.
7. **Document + commit.** Update the cycle log here + the relevant finding doc; commit LOCAL.

## Guardrails (hard rules)
- **Boxes:** spin up vast.ai only when a cycle needs CUDA; **kill when done + VERIFY the account is
  empty**; switch boxes FAST if slow to boot; pick many/fast vCPUs for parallel benches.
- **Honesty:** never over-claim "closed." Per-metric PASS/fail. If a fix is refuted, say so.
- **Do NOT touch** `viennaps-accel/` or `plasma_sim/` (vendored/reference dirs).
- **Repos:** `plasma-etching-code` is public/solo — commit freely, **do not push unless asked**.
  The `~/chip-etch` docs repo is shared/PRIVATE (fetch+rebase before any push; only `docs/` files).
  NEVER commit Resona details to either repo.
- **Commit trailer:** `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>` +
  `Claude-Session:` line. Commit each cycle's work; keep local.

## Target queue (reprioritize as findings land)
- [active] **C1 — reconcile the two charging solvers** (charging_general vs charging2d edge_array):
  quantify the floor-flux gap (0.16 vs 0.22) and attribute it to the old solver's sheath launch
  plane (`boundary_um=3.7`) + adaptive survivor-gated integrator + rf_bursts source. Output = the
  spec for the kinetic engine.
- [next]  **C2 — fine-grid GPU kinetic charging** (device-resident, cell-sort, CUDA-graph, phase-
  resolved Boltzmann-gated source) to produce HG's spatially-selective anti-shadowing → close the
  floor over-charge. CUDA box. Gate: AR4 floorV→33, floorFlux→0.22 WHILE neighbor stays ~39.
- [queued] **C3 — new chemistry frontier** (research agent running): ALE synergy curves
  (Kanarik/Gottscho), cryo chemistry, HAR SiO2. Pick one with a clean published benchmark; build a
  fast + differentiable reduced kinetic module. Gate vs the paper's EPC/coverage numbers.
- [queued] **C4 — differentiable calibration demo**: wp.Tape through charge→field→trajectory→profile;
  invert Te/V_rf/flux against a target. The payoff of differentiability.
- [backlog] task #42 DDA flux-conservative; task #43 notch wired into etch vs Nozawa/Fujiwara.

## Cycle log
(newest first; each entry: target, what was done, gate result PASS/FAIL, artifacts, commit)

- **C3 (2026-07-06, in progress):** ALE chemistry module (`src/petch/ale.py`, `CHEMISTRY_FRONTIER_ALE.md`).
  Research scout picked directional plasma ALE (Si/Cl₂/Ar⁺), Vella–Graves 2025 ROM — CPU, differentiable,
  unoccupied. CAUGHT a bad equation: the scout's `J_Si = J_Ar·Y_Si·(1−θ₂)` gives ZERO etch in the
  15–20 eV window (Y_Si=0 there) — that's the sputter-leak term only; real ALE removal is the SiCl/SiCl₂
  channels. Pulled the real paper (OSTI 2586627), agent extracting the exact Si-removal flux + ODEs +
  cycle schedule. GATE PENDING: EPC 0.67/0.75/1.0 Å/cyc at 15/17.5/20 eV, rise at 22.5/30 eV.

- **C1 (2026-07-06, DONE):** reconciled charging_general vs charging2d edge_array. GATE = explain the
  floor over-charge. RESULT: **the electron source model brackets HG's floor** —
  | convergence | old (sheath launch plane) | new (z=1 launch) | HG |
  |---|---|---|---|
  | nit 120 | 21.5 V / 0.31 | 38.4 V / 0.19 | 33 V / 0.22 |
  | nit 400 | 26.4 V / 0.28 | 44.2 V / 0.18 | 33 V / 0.22 |
  | nit 700 | 25.6 V / 0.32 | 45.2 V / 0.17 | 33 V / 0.22 |
  Old solver's `boundary_um=3.7` sheath launch plane + rf_bursts OVER-delivers floor electrons
  (plateaus at 25.6 V, below HG); the general engine's z=1 launch UNDER-delivers (climbs to 45 V, past
  HG). HG's 33 V / 0.22 sits between → the faithful phase-resolved source IS the lever (C2 spec earned).
  Also: at short convergence e_traced==e_geom (0.124); as the floor charges, traced drops below
  geometric — a positive feedback (more charge → stronger entrance barrier → fewer e⁻ → more charge).
  Artifacts: `scripts/charging_solver_reconcile.py`, `viz/charging_solver_reconcile.png`,
  `FLOOR_OVERCHARGE_FINDING.md`. Commit: local.
