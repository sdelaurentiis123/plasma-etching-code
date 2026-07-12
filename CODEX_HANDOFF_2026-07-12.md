# Handoff to Codex / GPT-5.6 — petch ARDE + differentiability + charging

Date: 2026-07-12. From: a Claude Code session. Repo:
`/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code`. Branch:
`codex/unified-engine-root-fixes`. **Local commits only — do not push.** Suite: `pytest -q` = 305
passed, 1 skipped. Read `AUTONOMOUS_PROGRESS.md` (running log) and `ARDE_PHYSICS_REFERENCE.md`
(sourced analytic targets) for full detail; the approved plan is `.claude/plans/wiggly-launching-fern.md`
(also summarized below).

This continues the Codex reconciliation in `CLAUDE_CODE_HANDOFF_2026-07-12.md` (three lineages; the
common 3-D engine `feature_step_3d.py` / `feature-3d-v1` is the only target for new physics; every
historical benchmark is de-earned until re-proven through it). All discipline there still holds.

## What this session established (three pillars, all demonstrated + gated)

1. **ARDE neutral transport is validated + converged + deterministic + GATED.** The common engine
   reproduces free-molecular floor flux three independent ways (adjoint gather, forward+QMC, and an
   independent pure-numpy analytic ray-trace `scripts/arde_mc_reference.py`), converging in grid
   (>=~25 cells/opening) and angle (AMR), bit-identical on re-run. First ARDE tests in the suite:
   `tests/test_arde_transport.py`. See `ARDE_PHYSICS_REFERENCE.md`.

2. **de Boer SF6/O2 ARDE from REAL coupled chemistry through the engine** (`scripts/deboer_feature3d.py
   --deboer`). Reparameterizes the validated coupled kernel `ReducedSiO2FluorocarbonMechanism` for Si-F
   (complex=SiF_x built by F / removed by ion; polymer=O-passivation gating F access). Calibrating ONE
   physical knob (F sticking s_F=0.08) on the knee (AR10,20; RMSE 0.033) PREDICTS the held-out AR40
   floor = 0.166 vs exp 0.20. Genuinely ion-assisted (ion_flux=0 -> rate exactly 0). AR40 runs in ~10 s
   (the feared "mesh wall" does not exist). Retires the additive `scripts/deboer_two_channel.py`.

3. **Differentiable calibration (the moat) + self-consistent charging both demonstrated.**
   - Calibration gradient `d(floor flux)/d(sticking)` is EXACT via implicit diff of the radiosity fixed
     point (one adjoint solve), matching FD to ~1e-7 (`scripts/diff_calibration_gradient.py`, gated).
     Multi-param recovery shows the data-efficiency wedge: adjoint solve-count flat in K, FD grows O(K)
     (`scripts/diff_calibration_demo.py`). Chemistry-parameter gradients are solved; only GEOMETRY/shape
     gradients hit the discontinuous ray-hit boundary (the real "differentiability open").
   - Self-consistent charging loop converges (I+ = I- exactly, self-consistent potential) on a grounded
     dielectric: `scripts/charging_selfconsistent_demo.py`; gate
     `tests/test_feature_step_3d.py::test_feature_step_solves_charge...`.

## Reproduce (all CPU, deterministic)

```bash
pytest -q                                        # 305 passed, 1 skipped
pytest -q tests/test_arde_transport.py           # ARDE transport + reactive + calibration-gradient gates
python scripts/deboer_feature3d.py --deboer      # coupled de Boer, calibrate knee predict AR40 (~1 min)
python scripts/diff_calibration_gradient.py      # exact sticking gradient vs FD (~1e-7)
python scripts/diff_calibration_demo.py          # adjoint-vs-FD calibration wedge
python scripts/charging_selfconsistent_demo.py   # self-consistent charging converges (I+=I-)
```

## Remaining plan, in order (pick up at item 1)

1. **Periodic-trench charging fix (small, HIGH value), then charging on features (Phase 2b).**
   `scripts/charging_selfconsistent_demo.py` runs on a flat plane because a periodic-cell TRENCH trips
   `src/petch/charging_poisson_3d.py::lump_triangle_sheet_charge_3d` line ~200: mesh verts sit exactly
   on the cell boundary and float32 rounding exceeds the `1e-10` tolerance ("triangle vertices lie
   outside the nodal grid"). Fix = widen that boundary tolerance to a float32-appropriate value (or
   clamp `normalized_vertices` to `[0, cell_shape]`); verify the 305 suite stays green. Then enable
   charging on a real trench/notch (Fujiwara/Nozawa poly-on-oxide): build `NodalPoissonSystem3D` from
   the geometry (pattern `tests/test_feature_step_3d.py::_plane_poisson_system` + the demo's
   `poisson_from_geometry`), pass `charging_poisson_system` + `potential_origin/spacing` +
   `trajectory_fixed_dt` (multi-step needs a `charging_system_builder`). Gate: converged current
   balance, charging-off/on causality (the electron-shading dipole: upper walls negative, floor
   positive), then the deep-AR de Boer residual (0.166 vs 0.20) + sub-degree ion IAD.

2. **Full-chain gradient** — extend the exact adjoint from `d(flux)/d(param)` through
   `surface_kinetics.advance` to `d(etch rate)/d(sticking, yields, O-passivation)`, so the actual de
   Boer RATE is calibrated with gradients (not just transport flux). The chemistry `advance` is pure
   numpy (Strang-split exact sub-operators) and cheap (no ray tracing), so a hybrid works: analytic
   radiosity adjoint for the expensive transport part + cheap chemistry differentiation (FD reuse of
   the same fluxes, or a torch/jax re-expression like `ale_diff.py`). Then run calibrate-N/predict-N+1
   vs a derivative-free baseline (Krueger 2024) — the commercial-wedge proof.

3. **Geometry/shape gradients** — the genuinely hard part (discontinuous hit/escape). Approach:
   edge/boundary reparameterization (Li et al. `10.1145/3272127.3275109`; Loubet et al.
   `10.1145/3355089.3356510`) on the form factors, FD-gated. Spike + go/no-go before a full build.

4. **AMR** (`src/petch` has none — uniform `dx` grid always) — de-prioritized: the mesh wall was
   refuted (AR40 in 10 s). Needed for 3D/many-feature/gradient scale, not the current regime.

5. **GPU + accuracy-matched speed benchmark** vs ViennaPS (`device=` already threads to Warp CUDA;
   needs an A100).

6. **Phase 3** — Jeon SiO2 held-out transfer (`scripts/jeon_unified_baseline.py`), 3D holes (not
   slots), second chemistry through the unchanged contracts.

## Key files

- `src/petch/feature_step_3d.py` — the common engine (`solve_feature_3d`, `advance_feature_step_3d`).
- `src/petch/surface_kinetics.py` — the coupled ion+neutral mechanism template (`_substrate_step` is
  the synergy: ion removes complex the neutrals built; `advance` -> etch velocity).
- `src/petch/neutral_radiosity_3d.py` — the linear radiosity fixed point (sparse GMRES; the calibration
  adjoint differentiates this).
- `src/petch/charging_coupled_3d.py` — `solve_dielectric_charging_steady_3d` (the self-consistent loop).
- `src/petch/charging_poisson_3d.py` — Q1 Poisson + `lump_triangle_sheet_charge_3d` (item-1 fix here).
- `scripts/deboer_feature3d.py`, `scripts/diff_calibration_gradient.py`,
  `scripts/diff_calibration_demo.py`, `scripts/charging_selfconsistent_demo.py`,
  `scripts/arde_mc_reference.py`, `scripts/deboer_arde_static.py`.
- `tests/test_arde_transport.py`, `tests/test_feature_step_3d.py` (charging gates),
  `tests/test_ale_diff.py` (the torch FD-gate autodiff template).

## Discipline (unchanged)

Local commits only, no push. Every step: one causal observable + independent judge + numerical error
budget + compute ceiling; validated / converged / deterministic. No AR / benchmark / region branch in
governing physics. Unknown surface rates are DECLARED calibrated inputs with provenance/uncertainty
(the honest structure: calibrate on N, predict held-out N+1 — never fit an AR-shaped formula or a
per-AR parameter). Keep `scripts/arde_mc_reference.py` (particle MC) as the non-differentiable
ground-truth reference. Do not touch `viennaps-accel/` or `plasma_sim/`. Commit trailers:
`Co-Authored-By:` and `Claude-Session:` (or your own equivalent).
