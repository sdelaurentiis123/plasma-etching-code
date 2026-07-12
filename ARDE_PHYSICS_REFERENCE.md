# ARDE neutral-transport: analytic validation targets

The reference physics the common-engine ARDE gate (`scripts/deboer_arde_static.py`) must reproduce.
These are first-principles identities / established models, not fits. Regime: free-molecular
(Knudsen, Kn >> 1) interior, diffuse (cosine) re-emission on non-reacting walls — the regime
Coburn & Winters argue holds inside etch features.

Notation: aspect ratio `A = L/w` (slot width `w`, depth `L`); `s` = wall/floor reaction (sticking)
probability; `Gamma0` = open-field incident neutral flux.

## 1. Non-sticking transmission (Clausing factor), s -> 0

- Cylindrical tube (radius r, length L): long-tube asymptotic `W -> (8/3)(r/L) = (4/3)(D/L)`;
  all-length interpolation `W ~= 1 / (1 + (3/4)A)`, `A = L/D`.
- Long 2D slot (width w, depth L, infinite third dim): `W ~= (w/L)[ln(L/w) + C]`, `C = O(1) ~= 3/2`.
  Slot decays as `ln A / A` (slower than the tube's `1/A`) — the 2D-vs-round signature.

## 2. Reactive-surface ARDE bottom flux, walls+floor react with probability s

Normalized bottom flux `Gamma_b/Gamma0` (approx normalized etch rate):

- `s -> 0`: follows the passive Clausing transmission (power-law, gentle). Fitting form
  `Gamma_b/Gamma0 ~= 1 / (1 + (3/4) s A)`.
- `0 < s < 1`: transverse random walk gives `~A^2` wall collisions; survival `~ (1-s)^{N}, N~A^2`
  -> `Gamma_b/Gamma0 ~= exp(-alpha s A^2)`. Larger `s` => steeper collapse.
- `s -> 1` (pure line-of-sight shadowing): only direct mouth->floor trajectories survive; bottom
  flux = geometric view factor of the floor from the cosine mouth. Acceptance half-angle
  `theta_max = arctan(w/L) ~= 1/A`; `Gamma_b/Gamma0 ~ 1/A`.

The family is MONOTONE decreasing in `A` and STEEPENS with `s`, interpolating from the Clausing
conductance curve (`s->0`, gentlest) to pure geometric shadowing (`s->1`). This monotone, s-ordered
family is the correct convergence target for the transport solver.

## 2a. Exact geometric target used in the gate (s = 1)

For an infinite 2D slot, the Hottel crossed-strings view factor between the mouth strip and the
directly opposed floor strip is exact:

    T_geom(A) = sqrt(1 + A^2) - A          (asymptote ~ 1/(2A))

In the feature the flux-limiting slot runs from the mask top to the floor, so the effective aspect
ratio includes the mask: `A_eff = A + mask/opening`, and the target is `T_geom(A_eff)`. The gate's
`--validate-geometric` mode checks the converged engine transmission against this.

Verified (2026-07-12, common engine + adaptive angular refinement, dx=0.02um, opening=0.10um,
mask=0.05um):

| AR | A_eff | engine | sqrt(1+A_eff^2)-A_eff | ratio |
|----|-------|--------|-----------------------|-------|
| 1.0 | 1.50 | 0.3056 | 0.3028 | 1.009 |
| 1.5 | 2.00 | 0.2333 | 0.2361 | 0.988 |

Agreement to ~1-2% where the angular quadrature resolves the acceptance cone.

## 3. Numerics: adaptive angular refinement (AMR) is REQUIRED

The floor-reaching acceptance cone has half-angle `~arctan(1/A) -> 1/A`. A fixed uniform angular
quadrature of `N` directions has resolution `~pi/N`; once `1/A < pi/N` (i.e. `A >~ N`) the entire
floor-reaching cone falls BETWEEN quadrature nodes -> the deterministic solver misses the direct flux
(artificial zero) or aliases onto one node (biased flux). The observed symptom in a fixed 5-node
quadrature: floor transmission FLATLINES (`~0.53` at s=1 for `A>=2`, ~17x too high at A=16) instead
of following `sqrt(1+A_eff^2)-A_eff`.

Required fix: angular samples must scale `∝ A`, or use importance/adaptive angular sampling
concentrated near vertical, or Monte Carlo with ray count raised so `O(10^2-10^3)` rays fall inside
the shrinking cone. The forward first-hit tracer batches all angular atoms into one Warp kernel
(GPU-ready), so high angular resolution is affordable there; the adjoint face-gather loops over atoms
in Python and is angular-resolution-cost-bound. High-AR + speed => forward tracer / importance sampling
on GPU. Deviation from monotonicity or from the `s->0` Clausing and `s->1` view-factor limits at high
AR is a numerics (angular-resolution / ray-starvation) artifact, not physics.

## 4. de Boer / Blauw experimental context

- de Boer et al., J. Microelectromech. Syst. 11, 385 (2002): SF6/O2 cryo Si; etch rate falls
  monotonically with AR; used to set DRIE guidelines.
- Blauw et al., JVST B (2001-2002): SF6 Si DRIE rate-vs-AR fit well by a Knudsen neutral-transport
  model -> ARDE is neutral-transport (conductance) limited, not ion-limited, over the measured range.
- Range ~ AR 20-40; high-AR normalized rate ~0.3-0.5 and below.

## Sources

- Coburn & Winters, Appl. Phys. Lett. 55, 2730 (1989).
- Gottscho, Jurgensen & Vitkavage, J. Vac. Sci. Technol. B 10, 2133 (1992).
- Berman, J. Appl. Phys. 36, 3356 (1965); Steckelmacher, Rep. Prog. Phys. 49, 1083 (1986).
- de Boer et al., J. Microelectromech. Syst. 11, 385 (2002); Blauw et al., JVST B (2001-2002).
- "Role of neutral transport in aspect ratio dependent plasma etching...," JVST A 35, 05C301 (2017).

Confidence: the tube 8/3 constant, s->0/s->1 limits, monotone-and-steepens-with-s, arctan(1/A)
acceptance, and the N ∝ A angular-resolution requirement are firmly established. The slot additive
constant C (~3/2) and the exact alpha in exp(-alpha s A^2) are geometry/model-dependent prefactors:
treat the SCALING (slot ~ln A/A; reactive ~exp(-alpha s A^2)) as load-bearing, not the O(1) constants.
The precise Blauw volume/page and exact high-AR rate numbers were not opened from the primary PDFs;
verify before quoting in print.
