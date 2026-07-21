# NEXT STEPS (pinned 2026-07-21, post quadrature finding)

Standing context: the 2026-07-20/21 session proved the diffuse form-factor quadrature
(rays_per_face=8) was unconverged at the mask neck — 2.75 nm endpoint scatter from the seed
alone, and refinement 8->32 moved the opening onto the experimental target with no knob
changes (noise experiment receipts: results noise-arm0/1/2, protocol R4). Every priority
below inherits that lesson: sampling depth is a physics-grade parameter and must carry a
convergence receipt.

## 1. Adaptive rays (top numerics priority after the R4 reveal)

Nobody in feature-scale etch simulation does adaptive transport sampling; graphics and
radiative heat transfer solved the pattern decades ago (hierarchical radiosity, Hanrahan
1991; variance-driven adaptive sampling; adaptive angular quadrature). Build order:

1. **Progressive per-face QMC with a split-sample variance gate.** The Sobol block sampler
   already supports (index_start, index_stop) continuation: shoot 8 rays/face, compare two
   interleaved half-sets' form-factor rows, double the face's rays (16, 32, 64 cap) until the
   declared row tolerance passes. Deterministic, receipted (per-face counts into the audit),
   necks refine themselves, open walls stay cheap: rays-32 accuracy at ~12 average rays.
2. **Aperture importance sampling** for neck/throat-dominated faces (next-event-estimation
   analog; changes the estimator -> own validation ladder).
3. **Hierarchical/clustered radiosity** O(N log N) — required for 200:1 face counts.
4. **Adjoint-weighted sample allocation** (research tier; dovetails with differentiability).

## 2. Active pipeline (R4 -> reveal)

rays=32 calibration bracketing (cand2 yield 0.574 in flight) + rays=64 convergence check ->
freeze (tooling rebound to K24-PETCH-R4, committed) -> sealed oxygen/power transfer runs ->
reveal + scorecard -> 5 nm post-hoc confirmation later. If cand2 misses: one more secant
interpolation inside the 6-run budget.

## 3. Partner deliverables (unblocked, this week)

STL->level-set importer (~1 day, first thing that breaks otherwise); 200:1 static ladder
(transport starvation curve anchored to the industry 1.3%-at-100:1 figure + charging
potential vs depth; needs NO calibration); time-to-depth manufacturability curves;
axisymmetric hole mode (biggest lever for 200:1 dynamics); LER propagation study (needs
wider 3-D cells + local-shape work).

## 4. Transport hygiene package

Seam stuck-loop fix (dispositioned, not fixed; 200:1 corridors multiply it); rays
convergence policy folded into every future operator declaration; polymer-ledger clip path
manufactured-test sweep.

## 5. Calibration/UQ harness v2

Noise-aware gating (replicates when observables near-critical), identifiability
diagnostics (R3 showed opening barely identifiable — detect after 3 points, not 7),
multi-fidelity (10 nm explore / 5 nm confirm), adjoint inner loop later (Tier 2).

## 6. After the reveal

Speed stack port (transport reuse + GPU-resident, then hierarchical radiosity + angular/
spatial AMR); Nozawa sealed-notching protocol draft (next validation flag, plays to the
charging differentiator); SEE for >20:1 charging; deck/site claim refresh against the
converged operator (the 0.19% claim must cite the frozen R4 numbers, not the pre-repair
lineage); methods-paper skeleton: certified transport + graduated refusals + the
quadrature case study.
