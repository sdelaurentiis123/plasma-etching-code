# Autonomous progress log (single writer)

Started 2026-07-12, continued from the Codex reconciliation handoff. Branch:
`codex/unified-engine-root-fixes`. **Local commits only, no push** (per user). One active writer.

Goal restated by the user: reach the OSTI/Oehrlein roadmap targets with a self-consistent engine
that **matches correct aspect ratios in a correct first-principles physics way** — no shortcuts, no
AR-shaped fudge; adaptive mesh/phase-space refinement (AMR); fast + GPU-accelerated at the end.

## Verified state (independently audited, not trusted from docs)

- **Branch is sound — KEEP it.** Forensic audit (Fable, read-only): 128 commits, one linear 24h
  chain, zero deletions, no concurrent-writer scars, docs *under*-claim. Rollback only re-introduces a
  known geometry bug and deletes real tested work. `pytest -q` = **296 passed, 1 skipped**.
- **The green tests are real conservation/analytic gates but contain ZERO aspect-ratio validation.**
  The one "aspect ratio ladder" test (`test_boundary_transport.py:89`) fires a *vertical monodirectional*
  beam into a straight trench and asserts floor flux = 1.0 at AR 1/4/16 — true by geometry, cannot show
  ARDE. All de Boer/ARDE work is loose scripts against the *legacy* monolith and was de-earned in the
  reconciliation. **"Matches aspect ratios first-principles" is unbuilt in the common engine.**

## Direction

1. **de Boer SF6/O2 ARDE through `feature-3d`** (radical-transport-limited → isolates AR transport from
   chemistry). First: a static neutral-flux-vs-AR gate that shows the floor-flux collapse and proves it
   CONVERGES under refinement. Then the rate closure with radical sticking as a *declared* input;
   calibrate low-AR, predict held-out AR40.
2. Jeon SiO2 depth-transfer (harness already built). 3. Charging only if it moves the profile above the
   error budget. (OSTI_OEHRLEIN_PROBLEM_MAP.md §6 ordering, adjusted so the pure-transport gate leads.)

## Work in progress

- **`scripts/deboer_arde_static.py`** (new): static ARDE floor-flux gate through the common engine —
  half-Maxwellian *flux* neutral source (analytic cosine law, no fitted angular closure) → ballistic
  face-gather (`gather_boundary_state_ballistic_3d`) → diffuse molecular-flow radiosity
  (`solve_diffuse_neutral_radiosity_3d`, one physical sticking `s`). Conservation exact (balance ~1e-13).

- **CONFIRMED: a FIXED angular quadrature flatlines the floor flux at high AR; the AMR requirement is
  real.** At s=1 (pure line-of-sight shadowing) a fixed 5-node quadrature flatlines transmission ~0.53
  for AR≥2 (~17x too high at AR16). Root cause, literature-confirmed (Coburn-Winters; JVST A 35
  05C301): the floor-reaching acceptance cone is ~arctan(1/A), so a fixed N-node quadrature aliases it
  once A≳N; angular samples must scale ∝A.

- **GPU-ready high-AR path found: QMC-sampled source + forward first-hit tracer.** The forward tracer
  batches every angular atom x source position into ONE Warp kernel; sampling the flux density with
  N=2^L Sobol points (`thermal_neutral_qmc_boundary_state`) concentrates rays by the physical cosine
  measure so the acceptance cone is resolved by raising N. It is **N-converged and fast** (AR16 in ~1s,
  stable across N=2^14..2^18): T(AR8)=0.065, T(AR16)=0.0355.

- **VALIDATED (three independent ways).** The common engine's ballistic neutral transport reproduces
  first-principles geometric shadowing (s=1). The two open items from the earlier correction are both
  closed:
  1. **Method reconciliation:** the adjoint gather CLIMBS to the forward+QMC value as it resolves
     (AR1: nt 24->48->72 gives 0.302->0.342->0.348 vs forward+QMC 0.345, ~1%). The earlier "8-14%
     disagreement" was the adjoint being angular-under-resolved (nt=24); properly refined, the two
     independent estimators agree.
  2. **Exact reference:** an INDEPENDENT pure-numpy analytic ray-trace of the same box
     (`reference_floor_transmission`, no engine code) gives 0.303/0.193/0.110/0.059 at AR1/2/4/8 —
     matching the opposed-strip view factor. The engine CONVERGES to it under GRID refinement
     (AR1: dx 0.02->0.0125->0.008 gives 0.345->0.295->0.310 vs ref 0.303). The coarse-grid dx=0.02
     over-prediction (~14% at AR1, ~3% at AR8; worse at low AR because 5 cells/opening staircases the
     walls) is grid error, not an engine bias.
  So getting a correct ARDE number requires BOTH adequate grid (>=~12 cells/opening) AND angular AMR
  (QMC-refined or nt ∝ A; a fixed quadrature flatlines/over-predicts the deep cone). Both are
  error-driven, no shortcuts. The earlier committed "premature 2% at AR1" (405e555) was the coincidence
  of nt=24 + dx=0.02 landing near the approximate analytic; corrected in 9eafb1f and closed here.

- **Committed gate:** `tests/test_arde_transport.py` — the first ARDE physics test in the suite:
  monotone ARDE collapse, engine == independent reference within grid tolerance, and the AMR-necessity
  regression (a fixed coarse quadrature over-predicts the deep flux; QMC does not).
  Analytic targets + citations: `ARDE_PHYSICS_REFERENCE.md`.

## Reactive s<1 family (qualitative demonstration, not yet a validated gate)

Engine forward+QMC, dx=0.02, floor transmission vs AR for sticking s (Coburn-Winters regime):

```
  AR  s=0.03  s=0.1   s=0.5   s=1
   1  0.928   0.797   0.466   0.345
   2  0.871   0.660   0.300   0.212
   4  0.717   0.436   0.160   0.119
   8  0.487   0.220   0.077   0.062
```

Correct behavior: monotone-decreasing in AR for every s; higher s -> lower floor flux (more wall
consumption); gentlest collapse at s->0 (T(AR8)/T(AR1)=0.53, approaching the Clausing conductance
limit) steepening as s rises. s=0.5 slightly steeper than s=1 is the expected exp(-alpha*s*A^2)-vs-1/A
crossover, not a defect. This is a DEMONSTRATION that the reactive re-emission machinery
(diffuse radiosity with one physical sticking coefficient) produces the right physics; a validated
reactive gate needs an independent multi-bounce reference (next step).

## Roadmap to the de Boer product (remaining)

1. Reactive gate: independent multi-bounce numpy reference (diffuse re-emission with sticking) ->
   promote the s-family to a validated pytest gate like the s=1 one.
2. de Boer SF6/O2 rate curve: rate model = f(ion flux, radical floor flux, sticking); sticking is a
   DECLARED calibrated input with provenance/uncertainty. Calibrate on low-AR points, PREDICT held-out
   AR40. Report grid/ray/digitization/model error separately.
3. GPU: run the forward+QMC path with device="cuda" (already threaded); accuracy-matched speed report.
4. Then Jeon SiO2 depth-transfer, then charging (only if it moves the profile above the error budget).

## Guardrails honored

Single writer. Local only. No benchmark/AR/region branch in governing physics. Sticking is the one
declared physical surface input with provenance. Refinement is error-driven. Device-agnostic transport
(Warp CUDA path ready for the GPU deployment).
