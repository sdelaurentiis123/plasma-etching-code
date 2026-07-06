# C10 spec — Bosch DRIE scalloping, SEM-gated (open+gated+differentiable = unclaimed)

Research verdict (2026-07-07): no published open, quantitatively SEM-gated, feature-scale Bosch
scallop benchmark exists. ViennaPS's `boschProcess` example is demonstrative (arbitrary units, no
gate); Ertl & Selberherr 2010 (the academic 3D reference) is a parameter study with NO experimental
validation; the strongest experimental gate in the literature is Volland & Rangelow JVST B 20,3111
(2002) (profile agreement, closed code). An open differentiable Bosch model gated on published
scallop numbers is unclaimed space.

## The cycle model (standard emulation, same as Zhou'04/Ertl'10/ViennaPS)
per cycle: (1) conformal polymer deposition thickness t over the profile; (2) directional ion
punch-through removes t on up-facing surfaces (floor clears, sidewalls stay coated); (3) etch of
exposed Si: isotropic bite R (F neutrals) + directional advance; (4) repeat N cycles.

## Closed-form arc geometry (derivation; the reduced gate abstraction)
- pitch p = vertical etch per cycle = D_total/N
- scallop depth s = R − sqrt(R² − p²/4) ≈ p²/(8R) (small-scallop limit)
- undercut U ≈ R of the FIRST cycle (unpassivated lateral bite under the mask)
- scallop ∝ t_etch (both p and R scale with etch-step time) — the observed linear scaling
- SELF-CONSISTENCY on the anchor dataset: Ayon s=140nm, p=434nm → fitted R=238nm vs measured
  U=250nm — closes to 5%. Passing s+p+U simultaneously is a physics test, not curve fitting.

## Gate Config R — long-cycle "rough" (Ayon 1999 JES 146,339 via McVittie NNIN deck, firm slide text)
STS-HRM shallow-trench recipe: 3.5 s SF6/O2 etch cycles, 65 cycles (6 min), 2 µm trench:
  1. total depth D = 28.2 µm ± 10%       (4.7 µm/min × 6 min)
  2. scallop pitch p = 434 nm ± 10%      (28.2 µm / 65)
  3. scallop depth s = 140 nm ± 35 nm    (slide text "140 nm scallops")
  4. mask-edge undercut U = 250 nm ± 50  (slide text)
  (secondary: wall angle 90.2° ± 0.5°; PR selectivity 76:1)
Calibrate per-cycle iso bite R≈240 nm + vertical clear≈434 nm from the published endpoints, then
score 1-4 with NO further tuning.

## Gate Config S — short-cycle "smooth" (Tillocher 2021, Micromachines 12,1143, open access, firm)
Ultrafast switching: 500 ms SF6 etch / 50 ms passivation, 1000 cycles, 10 µm trench:
  1. D = 60.8 µm ± 10%                   (firm text, Fig 4)
  2. p = 60.8 nm ± 10%                   (D/N)
  3. s ≤ 30 nm                           (upper bound, "residual roughness")
  4. ARDE: depth(4µm)/depth(10µm) = 49.8/60.8 = 0.82 ± 0.05 (both firm)

## Cross-config assertion (kills accidentally-right models)
s(3.5 s etch) / s(0.5 s etch) ≥ 4 within one model family.

## Bonus/optional
- Code-to-code: Ertl & Selberherr Fig 8/9 (2.5 µm hole, 20 cycles, their Tables 1-2 params) — beat
  their 2-day/8-core runtime. Open PDF: https://www.iue.tuwien.ac.at/pdf/ib_2009/hashed_links/ep4PPErIJjnr4Y_us.pdf
- Qualitative: ARDSA (scallops shrink with depth at FIXED recipe — DREM 2018; production recipes ramp
  so Park saw the opposite — don't gate on sign without the recipe); undercut monotone in etch:dep
  ratio + pressure (McVittie trend table).
- Scallop-depth measurement protocol: Park 2020 Fig 1c (tangent-line valley-to-crest, avg >=5 ripples).
- Independent geometry stats: Park 2020 (SPTS, 70µm/230µm holes): s=130 nm top / 230 nm bottom.

## Transport check (for the physical tier beyond emulation)
Ertl Eq 20 analytic bottom-center ion flux vs AR x: F = F_src[1 − (2x/√(1+4x²))^(κ+1)] (<1% vs MC).

## Sources (verified URLs)
- McVittie NNIN deck (Ayon slides, firm numbers): https://people.eecs.berkeley.edu/~pister/147fa14/Resources/BoschProc-STS.pdf
- Ayon et al. JES 146,339 (1999): https://iopscience.iop.org/article/10.1149/1.1391611
- Park et al. 2020 (open, protocol + geometry stats): https://doi.org/10.1186/s40486-020-00116-x
- Tillocher et al. 2021 (open, Config S): https://pmc.ncbi.nlm.nih.gov/articles/PMC8537062/
- Laermer & Schilp US 5,501,893 (process definition): https://patents.google.com/patent/US5501893A/en
- Ertl & Selberherr 2010 (open PDF above); ViennaPS example:
  https://github.com/ViennaTools/ViennaPS/tree/master/examples/boschProcess
