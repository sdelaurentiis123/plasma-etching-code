# petch frontier loop — autonomous research/build charter

Standing autonomous loop (authorized 2026-07-06). Each cycle pushes one frontier of the open,
GPU, differentiable, feature-scale etch simulator, benchmarks it honestly against published data,
visualizes it, smoke-tests it, documents it, and commits. Runs across wake-ups; survives compaction.

## Superseding operating contract (2026-07-12)

The external intended-use demonstration, rather than an internal residual in isolation, defines the
active objective. The first demonstration is Jeon 2022 SiO2 **dimensionless depth transfer**: calibrate
only the 20% C4F8 continuous-wave width curve, then predict the untouched widths, gas fractions, and
pulse conditions, including all twelve 1 ms pulse-response directions. The versioned scorecard is
`src/petch/validation_demo.py`. It was frozen before the first unified-engine prediction run.

Every proposed physics or numerical change now follows this bounded loop:

1. derive the change or cite a primary source, and state which error term it should reduce;
2. pass the nearest analytic/conservation/limit invariant;
3. immediately widen one ring: neighboring widths, aspect ratios, geometries, and the external score;
4. compare with the last accepted checkpoint, reporting improvements and regressions separately;
5. promote only if the broader evidence improves, or document an explicit, useful tradeoff;
6. after two attempts that do not improve the external objective, stop and change method or research
   the missing physics rather than extending the same campaign.

Numerical approximations are allowed only with a measured convergence error below the relevant
experimental/model uncertainty. Unknown surface chemistry remains an explicit closure input with
provenance and validity limits; it must not be disguised as derived physics. Geometry-, benchmark-,
aspect-ratio-, or condition-specific correction branches are prohibited.

Charging is no longer an open-ended `confidence max < 0.15` campaign. Its accepted fixed-map
checkpoint is retained. For the Jeon demonstration, charged and uncharged predictions are first
compared on the narrowest/deepest structure. Charging work resumes only if that profile sensitivity
is material relative to the declared error budget, and then converges against profile error rather
than a detached balance statistic.

No remote accelerator is launched until a reproducible local bottleneck is named. Each remote run
needs a wall-time and dollar ceiling in advance, must preserve an accepted checkpoint, and ends by
verifying that no instance remains. Nothing is pushed without an explicit request.

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
- [active] **D1 — Jeon SiO2 depth-transfer baseline**: run the unified engine without using held-out
  observations, expose numerical/model/parameter-evidence errors separately, and identify the one
  dominant product-facing miss. Gate with `jeon_2022_depth_transfer_v1`.
- [next] **D2 — bounded SiO2 closure calibration**: source or derive fixed yield/threshold/angular
  forms; calibrate no more than three global surface parameters on the five nontrivial 20% C4F8 CW
  ratios; predict the 64 nontrivial held-out targets without condition-specific corrections.
- [next] **D3 — charging relevance ablation**: charged versus uncharged narrow/deep Jeon structure;
  resume charging convergence only when the predicted profile effect exceeds the error budget.
- [historical] **C1 — reconcile the two charging solvers** (charging_general vs charging2d edge_array):
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

- **C11 (2026-07-07, DONE — mechanism PASS under the fully-derived table; AR4 trend = coupling frontier):**
  petch-computed charging table (Q/Vf/E_defl vs AR from the theorem-correct solver, zero knobs) wired
  into the notch via surface_charging="petch". GATE A (charging-specific mechanism) PASS: off=0 at all
  AR, resolved notches AR>=2 (0.207/0.265 um at AR2/3 — larger than the HG-closure mode). GATES B/C
  fail at AR4: the notch forms then is ERASED in late overetch (0.076@100%OE -> 0.004 final) — the
  table's documented +15-20% AR4 over-charge (Vf 39.7 vs 33; Q 0.40 vs 0.22) amplifies through the
  overetch redep/erosion balance. Default stays "hg" (validated); "petch" is the first-principles path
  with the AR4 sensitivity as the named frontier. E_defl(AR4)=28.0 eV lands exactly on HG's 28.

- **C9 (2026-07-07, CONCLUDED — the charging module is pure-physics; joint residual ~20-35% attributed):**
  the full first-principles rebuild: corrected integrator (C6) + mirror-image tracer BCs + log
  current-balance dynamics + physical -10Te bound + DERIVED ion sheath source (exact at 400 kHz:
  bathtub+Bohm+anticorrelation, P(E<33)=0.434=analytic) + electron derivation pushed to its honest
  ceiling (projection bug found+fixed: cos^0.35→cos^1.0; isotropic-top→cos^0.88; HG's Vlasov gives
  0.6 — the residual is thick-sheath time-resolved dynamics, agent hunting Ootera's theory).
  Field reversal REJECTED as explanation by mass-ratio flux balance (collapse thermal flux covers
  ions 2.6x). Final config table (AR4, all crutches OFF):
    pure_lambert: floor 32.4✓/0.318/11.3/27.7 | pure_cos06: 39.7/0.259/10.3/35.3 |
    HG-faithful: 44.8/0.256/13.4/39.9✓ | HG: 33/0.22/7/39.
  Every observable individually reachable with faithful physics; foot peak 61V≈HG 60; all structure
  (dipole, splits, focusing, deflection feedback) emergent. Remaining joint gap = 2D-projection
  conventions + Vlasov EADF breadth + rectangular mask idealization — named modules, not knobs.
  Retrofit queued: wire petch-computed charging into the notch arc (task #43).

- **C10 (2026-07-07, DONE — both configs + kill-test PASS):** Bosch DRIE scalloping, SEM-gated
  (`src/petch/bosch.py`, BOSCH_BENCHMARK_SPEC.md). The first open, SEM-gated, forward-predictive Bosch
  model (VLSet-AE 2026 = inverse/measurement, cited+distinguished; ViennaPS v3.6.0 example ungated;
  arXiv 2606.11247 names the niche). Config R (Ayon 1999): depth 28.6 (28.2±2.8), pitch 440 (434±43),
  scallop 140 (140±35 — exact), undercut 220 (250±50). Config S (Tillocher 2021): pitch 60.0 (60.8±6),
  scallop 15 (≤30), D→60 µm ✓. Cross-config s_R/s_S = 9.3 (≥4). The DATA forced the physics: sequential
  punch-then-iso mechanics (analytic s = r−sqrt(r²−(p/2)²) = 140 predicted, grid confirmed); the
  simultaneous model provably gives 32 nm. Deferred: ARDE 0.82 sub-gate (needs the transport tier).
  3/3 tests; suite 26/26. Artifacts: viz/bosch_scallops.png, tests/test_bosch.py.

- **C5 (2026-07-06, DONE):** CRYO etch chemistry (`src/petch/cryo.py`) — the hot 2023-2026 chemistry
  (3D-NAND/DRAM), unoccupied in open source. Temperature-dependent physisorption: a condensed etchant
  layer builds as T drops (residence time ~exp(E_ads/kT)) and multiplies the etch rate. Langmuir isotherm
  θ(T)=Kp/(1+Kp), Kp=A·exp(E_ads/kT); rate ER(T)=R_base·(1+gain·θ). E_ads=0.4 eV FIXED from independent
  measurement (Antoun Sci.Rep.11,357: 0.406 eV), not free-floated. GATE PASS vs CF4/H2 pseudo-wet (Small
  Methods 2024, doi 10.1002/smtd.202400090): ER(+20°C)=2.32 (target 2.3 plateau), ER(−60°C)=3.76 (target
  3.76), ratio 1.62 (target 1.6×). Honesty: the folk "rate doubles" anchor is really 1.6× and from CF4/H2
  not CHF3 — gated on the true number. Cross-checks: HF cryo-ALE ~3.2× (diff system), C4F8 residence-time
  cliff −120 vs −110°C reproduced. Smooth/differentiable (drops into the ale_diff pattern). 5/5 cryo tests,
  14/14 full chemistry suite. Artifacts: `scripts/cryo_window.py`, `viz/cryo_window.png`, `tests/test_cryo.py`.

- **C2 (2026-07-06, DONE — premise REFUTED with box data, box killed+verified-empty):** the original
  premise ("a fine-grid GPU kinetic engine closes the floor over-charge") is WRONG, shown on a rented
  RTX 4090 (contract 44044631, ~1 GPU-hr, destroyed, account verified 0 instances):
  1. **GPU barely helps here:** the tracer round-trips host↔device each iteration, so at W16 the 4090 is
     only ~35% faster than CPU (100 iters/31 s). A real GPU win needs the full device-resident rewrite
     (persistent arrays + cell-sort + CUDA-graph) — not worth it for this problem.
  2. **Fine grid is NOT the fix:** old sheath-source solver floorV = 26.7 (nit800) → 27.4 (nit1600) at
     W16 — a plateau; convergence and (per memory) W32 don't move it.
  3. **The electron-delivery knob feeds the WALLS, not the floor:** raising `open_wall_boost` 1.0→3.6 made
     the floor WORSE (44→47 V) and drove the edge line negative (over-fed) — wrong lever.
  Honest conclusion: the floor over-charge is set by the electron **launch geometry / EAD** delivering
  the right floor flux. Old solver (boundary_um=3.7 sheath launch) over-delivers → floor 27 (best base,
  ~18% low); general engine (z=1 launch) under-delivers → floor 45. The precise fix is a source-LAUNCH
  calibration (which launch plane + EAD delivers HG's 0.22 floor flux) — small-grid, CPU-cheap,
  differentiable. **NOT** a GPU build. Requeued as C6 (below). Net: spent ~1 GPU-hr to kill a multi-hour
  wrong path — good trade.

- **C8 (2026-07-06, DIAGNOSED + infra built; stable Poisson solver is the open build — CHARGING_POISSON_PLAN.md):**
  ROOT CAUSE of the floor over-charge fully pinned by 2 primary-source research passes + a line-referenced
  code audit: it is the **LOCAL vs GLOBAL charge->surface-potential map**, NOT the interior PDE and NOT SEE.
  HG solve Laplace-in-gas too (same as us); the difference is they map deposited charge to surface potential
  GLOBALLY and ε-aware (Coulomb superposition / full variable-ε Poisson), while we map it LOCALLY
  (`Vs+=net` per cell). Local map -> flat lateral potential -> no inward-bending field -> electrons hit the
  floor by pure geometry -> floor climbs to 37 V. Confirmed: at fine grid e_traced -> 0.124 = geometric exactly.
  My earlier "long-range fringing field" hypothesis was WRONG (the focusing is a LOCAL in-trench well effect).
  BUILT (committed, correct, backward-compatible): `GROUND` material + `add_grounded_substrate()` (oxide
  dielectric stack on grounded Si) + Poisson substrate BCs. NOT WORKING YET: `field_model="poisson"`+substrate
  is UNSTABLE (walls run to -1000s V, no focusing) — unbounded `rho` accumulation with no physical scaling.
  This is genuine MCFPM-level numerics. Concrete path in CHARGING_POISSON_PLAN.md: physical-unit charge
  (ρ·h²/ε₀), σ-sheet on interface cell, capacitance-matched substrate, stable damped charge dynamics, uneven
  conductor equipotential -> then remove the band-aid knobs. SEE ruled out for the floor (recaptured).
- [next] **C9 — implement the stable first-principles Poisson** per CHARGING_POISSON_PLAN.md (physical units +
  interface-σ + capacitance matching + stable dynamics). Gate: floorV->33, e_traced>geometric, all knobs OFF.

- **C6 (2026-07-06, DONE — real numerics fix, first-principles, no fudge; GPU-ready):** the floor
  over-charge is (partly) an INTEGRATOR bug, found by thinking physically. Tell: traced floor electron
  flux (0.115-0.130) was ≤ pure geometric shadowing (0.124) — physically impossible if the integrator
  is right, since electrostatic focusing must give MORE than geometry. The coarse leapfrog (dt_cap=0.45
  cell/step) was losing focused floor-bound electrons. Parametrized the integrator resolution in BOTH
  the numba CPU tracer AND the Warp-CUDA kernel (trace_dt/trace_dt_field/trace_steps; GPU-accelerated).
  Result (AR4, no focus knob):
    dt=0.45: floorV 46.0, e_traced 0.130 (≈geom)
    dt=0.15: floorV 41.6, e_traced 0.154 (ABOVE geom 0.124 — the anti-shadowing now EMERGES from the trace)
    dt=0.08: floorV 42.1, e_traced 0.152 (converged)
  So refining the integrator makes HG's electrostatic focusing appear from first principles (no fitted
  term), and drops the floor 46→42 V. Backward-compatible (default still 0.45; 4/4 charging tests pass).
  Also ruled out SEE as the floor fix by physics: a 1-5 eV secondary from a +33 V floor is recaptured
  (can't escape the well) → zero net electrons; SEE only redistributes from low-V regions.
  RESIDUAL (42→33): the genuine reduced-model limit — at W16 the coarse FIELD grid gives too weak a
  focusing lever. The honest close is fine grid + the corrected integrator on GPU (C2 tested fine grid
  with the BROKEN integrator; re-test warranted). The earlier `insulator_e_focus` knob is kept opt-in
  (default 0) as a documented phenomenological stand-in for this residual, NOT presented as the fix.
- [next] **C7 — fine-grid GPU floor close**: re-run the AR sweep at W32/48 with the corrected integrator
  (dt=0.15) on a CUDA box; test whether the now-emergent focusing + finer field lever lands floorV→33.

- **C4 (2026-07-06, DONE):** DIFFERENTIABLE ALE (`src/petch/ale_diff.py`, torch) — the moat payoff.
  Reverse-mode autograd through the whole cyclic site-balance chemistry: `dEPC/dE` and `dEPC/dparams`.
  Speed: integrate the coverage transient finely (~45 s) then add the constant bare-Si sputter-leak tail
  analytically → ~40× fewer autograd steps, forward parity with numpy ale.py <2%. GATE PASS: autograd
  `dEPC/dE` matches central finite-diff to 4 decimals (0.061/0.389/0.393 at 17.5/22.5/27.5 eV);
  gradient-based INVERSE DESIGN recovers the ion energy from a target EPC (self-consistent <0.1 eV;
  solved E=20.39 eV for target EPC=1.0 Å/cyc). 4/4 tests pass. No open feature-scale etch tool exposes
  chemistry gradients (ViennaPS: no ALE; Kushner: closed/CPU/non-diff). Artifacts:
  `scripts/ale_calibrate.py`, `viz/ale_calibrate.png`, `tests/test_ale_diff.py`. (Done live, no idle wait.)

- **C3 (2026-07-06, DONE):** ALE chemistry module (`src/petch/ale.py`, `CHEMISTRY_FRONTIER_ALE.md`).
  Directional plasma ALE (Si/Cl₂/Ar⁺), Vella–Graves 2025 ROM — CPU, differentiable, unoccupied in open
  source. Process: scout picked it → CAUGHT a bad equation (scout's `J_Si=J_Ar·Y_Si·(1−θ₂)` gives ZERO
  window etch, Y_Si=0 there; that's only the sputter-leak term) → pulled the real paper (OSTI 2586627),
  agent transcribed Eqs 6–17 verbatim → fixed: Si leaves in THREE channels (SiCl+SiCl₂ chemical-sputter
  via `(θ₁+θ₂)` + bare-Si sputter), real step times t_dose 0.112 s / t_barr 113.5 s (fluence/flux).
  GATE PASS (nothing tuned): EPC 17.5→0.76 (ROM 0.7), 20→0.86 (0.9), 30→4.74 (4.8); the 15–20 eV
  **window** + **self-limitation loss** above 20 eV + **synergy 100%→20%** all emerge from the yields
  alone. 15 eV floor 0.36 vs 0.6 runs low (near-threshold sensitivity) — the one soft point. 5/5 ALE
  tests pass. Artifacts: `scripts/ale_window.py`, `viz/ale_window.png`, `tests/test_ale.py`.

  Next chemistry (same skeleton, queued as C5): temperature-dependent physisorption term → cryo
  SiO₂/Si₃N₄ (CHF₃ "rate doubles +20→−60 °C", JAP 133,113306); thermal Al₂O₃ ALE Arrhenius unit-test.

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
