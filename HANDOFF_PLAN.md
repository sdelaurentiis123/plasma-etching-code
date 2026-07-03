# Handoff plan: close the HG charging gates, then production charging

Written 2026-07-03 for the next agent taking over. Everything below is measured or committed;
no aspirational numbers. Read `ROBUST_PHYSICS_MODEL_PLAN.md` (architecture + literature),
`CHARGING_DEEP_ISSUES.md` (what is falsified), and this file (execution order) before coding.

## 0. Ground truth: repo + working-tree state

- `plasma-etching-code` (this repo, solo, PUBLIC, remote `github-personal:sdelaurentiis123/plasma-etching-code`):
  - Everything through `e6c5eac` is committed. Commits through `19c40e7` are pushed; do NOT push
    again without the user asking. Never force-push.
  - UNCOMMITTED (working tree, in `src/petch/charging2d.py`): the two W0 numerics fixes below.
    Commit them once the reduced gate re-run is read and recorded here.
- `~/chip-etch` (shared with Craig, PRIVATE): clean, synced. Only ever add `docs/` files there;
  never touch `viennaps-accel/` or `plasma_sim/`. Fetch+rebase before any push.
- Open tasks: #42 (DDA flux-conservative march), #43 (notch-depth gate vs Nozawa/Fujiwara),
  #44 (this W0 stabilization, in progress).
- Operational rules that keep biting: no pushes unless asked; kill any vast.ai box when done and
  verify the account is empty; literature-first when stuck (user directive); every claim needs a
  script + published-data gate + PASS/fail line.

## 1. What was just fixed (uncommitted) and why

Codex's audit (`ROBUST_PHYSICS_MODEL_PLAN.md` section W0) found the edge-array solver
"residual-limited": current residual 0.24-1.1 vs gate 0.08, tail/final floor flux disagreeing.
Two root causes were found and fixed today, both numerics, not physics:

1. **Laplace was masked-Jacobi, not red-black GS.** `solve_edge_array_charging.laplace()`
   computed the neighbor average once per sweep and then applied it to both colors, which is
   damped Jacobi: it diverges for omega > 1 (verified: omega=1.88 overflows) and was therefore
   left at omega=1.0 — badly under-converged on the wide nonperiodic domain (nx=214 at the
   reduced setting). Fix: recompute the neighbor average inside the color loop (true red-black),
   exactly like the proven periodic `solve_trench_charging.laplace()`, and default omega=1.88.
2. **The residual gate was measuring shot noise.** The per-surface residual came from a single
   4x-n_per_iter snapshot; at n=1200 the trench sees ~360 particles, so 1-sigma shot noise is
   ~0.05 — the same size as the 0.08 gate. Measured proof: AR4 seed-104 run gave residual 0.067
   while the 8-seed saved run gave 0.10-0.33 scatter. Fix: residual is now the tail-averaged net
   current over the last k iterations (k = n_iter//3, ~40x the samples); the old snapshot is kept
   as `diag["residual_snapshot"]`.

Verification state: `tests/test_charging_edge_open.py` passes (4/4).

Reduced gate re-run WITH the fixes (`python scripts/charging_reduced_gate.py`, defaults:
edge_array, line_of_sight, analytic source, W16/mouth80/n1200/it100), 2026-07-03:
floor RMSE `0.079` fail (model ABOVE HG by ~+0.08 at every AR); survivor `0.0000` pass;
tail-avg residual max `0.220` fail; foot-E max rel err `0.624` fail (now overshoots at LOW AR:
24.4 eV at AR1 vs HG 15, while AR4 is 26.4 vs 27.5 — the deep end is close, the shallow end is
hot); foot-flux ratio `1.28` pass; Vpoly rising 3.3->25.0 vs HG 6->39, max rel err `0.452` fail;
Matsui pass (0.619); 0-D closure pass. Codex's pre-fix baseline (same config): RMSE 0.077,
residual 0.329, foot-E err 0.592, Vpoly err 0.504 — i.e. the numerics fix re-baselined the
numbers slightly; the old ones came from unconverged fields and are not a regression reference.

## 2. Measured diagnosis (do this fix first): PR insulator floating balance

Per-surface signed tail residuals at AR4 (n1200/it100) isolate the imbalance:

| config | floor | pr | poly_edge | poly_neighbor |
|---|---|---|---|---|
| default (`insul_vmin_Te=1`), seed 104 | +0.029 | **-0.142** | +0.021 | +0.048 |
| default, seed 7 | +0.052 | **-0.146** | +0.006 | +0.057 |
| `insul_vmin_Te=5`, seed 104 | -0.009 | **+0.134** | -0.004 | +0.012 |
| default, `n_iter=300`, seed 104 | +0.001 | **-0.325** | +0.016 | +0.001 |

Reading:

1. **Floor and conductor residuals converge to ~0 with iterations** — the RB-SOR fix made the
   rest of the solve honest. Conductors are just SLOW: Vneighbor 23.5 (it100) -> 27.1 (it300),
   still climbing toward HG's 39 V. So part of the Vpoly "miss" is unconverged conductor
   charging, not physics. Increase n_iter (or the conductor relax step, which is divided by
   `neighbor_exposed_area`) until Vedge/Vneighbor plateau; re-measure Vpoly gate only then.
2. **The PR clip is the real charge-conservation violation and it is NOT a knob.** Pinned at
   -1*Te the resist drains electrons forever (-0.14, growing to -0.33 by it300 as the trench
   charges up and pulls more electrons in); bounded at -5*Te it overshoots the other way
   (+0.13). The fix is a per-cell (or per-segment) floating insulator balance: PR potential must
   settle where local ion current equals local electron current (the erfc transverse-energy
   rejection already provides the electron cutoff mechanism the balance needs). Note commit
   94662b7 introduced the -Te bound as "physical"; it is not — it was masking drift. Replace it
   with the self-consistent balance, keep a wide safety bound (say -10*Te) only as a numerical
   guard, and re-run the residual table above. Kill criterion: if pr residual still exceeds 0.08
   tail-averaged after the floating balance, stop and re-derive the PR electron acceptance.
3. **Expected knock-ons once PR conserves charge:** the ~+0.08 uniform floor-flux excess should
   drop (those drained electrons belong in the trench balance), Vpoly should rise further, and
   the low-AR foot-E overshoot should be re-measured before touching anything else.
4. The `edge_extra_e` line-of-sight top-up residual is small (+/-0.02) — leave it as diagnostic;
   its removal is W2 as planned.

## 3. Physics order after W0 (from ROBUST_PHYSICS_MODEL_PLAN, sharpened)

Work these strictly in order; each has a gate and a kill criterion. Reduced setting first
(W16/mouth80/n1200), one high-stat confirm (W32/mouth237/n8000, seeds 0-7) only after reduced
passes — high-stat is ~25 min/config on this M1, use a vast.ai box only if iterating.

1. **W2 — explicit pattern electrostatics** (the known dominant miss; CHARGING_DEEP_ISSUES
   proves scalar boundary currents are insufficient and the imposed line-to-line bias moves foot
   energy eV-scale). Generalize `_build_edge_array_geometry` to include the real open field on
   the left (already partially there), remove the `line_of_sight` top-up from the production
   path, and let the open-area electron supply arrive as traced particles. Gates: HG floor RMSE
   <= 0.05; neighbor poly potential rising 6->39 V within 30%; foot energy rising 15->27.5 eV
   within 30%; Kushner sanity (symmetric pattern: no net lateral field; open/dense: dense-to-open
   field sign).
2. **W1 — table source interface** (pure refactor gate: replay == analytic within seed noise).
   Do it when touching the source anyway; it unblocks HPEM/PCMCM tables later.
3. **W3 — SEE in the pattern model** (sign test first; PR-only PMMA yields exist in
   `_pmma_see_yields`). Gate: AR4 fixed geometry, SEE lowers positive potential, total charge
   conserved. Note `solve_edge_array_charging` currently raises on see_model != none — port the
   trench SEE branch.
4. **W5 — one shared 3-D charging interface** in `threed.py` (deposit current -> field solve ->
   tracer query), replacing the per-mode HG table hook. The HG table stays as fallback closure.
5. **W4 — Cl2/ALE surface-state chemistry** (independent of charging; see
   CHEMISTRY_EXPANSION_PLAN.md; anchors already digitized in the robust plan).

Then the standing tasks: #43 notch-depth vs Nozawa (JJAP 34,2107)/Fujiwara (JJAP 34,2095) using
the multi-material etch-stop that already passed the shape gate (corr 0.92,
`notching_depth_result.npz`), and #42 DDA flux-conservative march (quarantine lift).

## 4. Reporting discipline

- Every gate run: exact command + env in the output, seeds recorded, PASS/fail per gate line.
- Docs (`docs/*.html`, one-pager, Experiments page) currently carry the honest baseline numbers
  (closure 0.039 config and conductor 0.060 trade). Do not update docs from reduced runs; only
  from the high-stat config with saved npz.
- When a claim is retracted or a number re-baselined, say so in RECONCILIATION.md — that file is
  the research log and has carried every retraction so far.
