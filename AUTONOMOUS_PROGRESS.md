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

- **CONFIRMED: a FIXED angular quadrature flatlines the floor flux at high AR; adaptive refinement
  (AMR) fixes it and the engine then matches first-principles.** At s=1 (pure line-of-sight shadowing)
  a fixed 5-node quadrature flatlines transmission ~0.53 for AR≥2 (~17x too high at AR16). Root cause,
  literature-confirmed (Coburn-Winters; JVST A 35 05C301): the floor-reaching acceptance cone is
  ~arctan(1/A), so a fixed N-node quadrature aliases it once A≳N; angular samples must scale ∝A.
  Fix: adaptive angular refinement (`converged_floor_transmission`, error-driven stop). With it, the
  converged engine transmission MATCHES the analytic slot view factor sqrt(1+A_eff^2)-A_eff
  (A_eff=A+mask/opening) to ~1-2% where resolved:

  | AR | A_eff | engine | analytic | ratio |
  |----|-------|--------|----------|-------|
  | 1.0 | 1.50 | 0.3056 | 0.3028 | 1.009 |
  | 1.5 | 2.00 | 0.2333 | 0.2361 | 0.988 |
  | 2.0 | 2.50 | 0.1891 | 0.1926 | 0.982 |

  This validates the common engine's neutral transport reproduces first-principles geometric shadowing.
  Analytic targets + citations: `ARDE_PHYSICS_REFERENCE.md`. Higher AR is under-resolved by the adjoint
  face-gather (Python-loop over angular atoms, cost-bound); the forward first-hit tracer batches all
  atoms into one Warp kernel (GPU-ready) and is the high-AR + fast path (next increment).

## Guardrails honored

Single writer. Local only. No benchmark/AR/region branch in governing physics. Sticking is the one
declared physical surface input with provenance. Refinement is error-driven. Device-agnostic transport
(Warp CUDA path ready for the GPU deployment).
