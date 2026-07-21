# Engineering plan: noise-free diffuse transport (2026-07-21)

Goal: eliminate form-factor sampling noise at its source so Krueger calibration observables
become well-posed, per RESEARCH_CHAOS_SOLUTIONS_2026-07-21.md. Two paths, phased.

## Phase 1 — Extruded-authority mode (the operator already exists; days -> hours)

Discovery: `deterministic_exchange_2d` (Hottel crossed-strings: EXACT unobstructed exchange
lengths, deterministic partial-visibility refinement, reciprocity by construction, escape =
unassigned row remainder) + `extruded_exchange_3d` (strict bridge to extruded triangle
surfaces, fail-closed extrusion certification, reciprocity/partition tests green) are already
in the tree as an authority backend (`form_factor_backend=deterministic_extruded_2d`,
requires periodic_lateral=True — the Krueger cell qualifies).

Plan:
1. Pilot flag `--radiosity-backend` (tooling change; exposes the existing option).
2. Local coarse smoke (dx=0.02 failed before on quadrature grounds — retry; else dx=0.01
   short run on CPU): confirm the full evolution LOOP preserves extrusion within the
   certification tolerance once the QMC y-noise source is removed (ballistic + chemistry are
   y-symmetric in exact arithmetic; certification is fail-closed if drift accumulates).
3. Full 10 nm Krueger base at the R1.9 pair, extruded backend: DETERMINISTIC, ZERO sampling
   noise, no rays parameter at all. Endpoint = the mean-field profile — exactly the
   self-averaged observable the literature prescribes (and what 2-D Kushner-type sims
   implicitly compute). Grid ladder (10/5 nm) is then the only convergence axis.
4. R5 protocol on this authority: calibrate the pair on extruded observables (grid-converged),
   freeze, run held-out conditions extruded; the existing narrow-cell 3-D ensemble data
   quantifies the fluctuation band around the mean-field claim. Categorical held-out gates
   evaluated on the extruded profiles.

Validation gates: extrusion certification every step (fail-closed); reciprocity + row-sum
(escape) ledgers; crossed-strings vs QMC-ensemble-mean agreement on a static mid-etch
checkpoint (the ensemble mean of our 6 chaos runs is the reference band); 10-vs-5 nm ladder.

Risks: (a) evolution drifts out of extrusion tolerance -> add a declared per-step y-projection
(mean-field closure, documented) or fall back to Phase 2; (b) crossed-strings partial
visibility at extreme neck -> covered by its deterministic refinement; verify at t~57 s state.

## Phase 2 — General-3-D analytic-unoccluded x visibility-ratio estimator (the moat build)

For wide cells, LER studies, holes/200:1 where extrusion does not hold: Baum-style analytic
unoccluded factors (contour/Stokes form) x stratified-ray visibility fraction (ratio
estimator; 1-2 orders variance reduction), hierarchical refinement with error bounds later.
Implementation-numerics research in flight (adjacent-edge singularities, periodic images,
culling); plan finalizes when it reports. Validated against: crossed-strings on extruded
geometry (exact cross-check!), existing ARDE analytic gates, QMC at high rays.

## Phase 3 — Statistics layer (from the solutions synthesis)

hetGP/replication upgrade of the BO harness; Mullins-scaled regularization option (length
below LER correlation length, never tuned to hide instability); wide-cell (>=5 xi) runs as
physical-fidelity confirmation and LER-study geometry; shadowing/EKI reserved for the
differentiable roadmap.

Sequencing: Phase 1 now (target: extruded base running today); Phase 2 after its research
lands and Phase 1 verdict known; Phase 3 alongside partner work.
