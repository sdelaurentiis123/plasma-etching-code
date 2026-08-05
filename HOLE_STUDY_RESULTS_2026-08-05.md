# 200:1 HAR hole study — phase-1 results (2026-08-05)

Series 1, 2 and 4 of `HOLE_STUDY_PLAN_2026-08-05.md`, plus the cascade
characterisation the plan lists as series 2b. Everything here is
frozen-geometry, deterministic and quadrature-exact — no Monte Carlo, no fitted
constant, no profile evolution (see *Not run* below). Reproduce with
`scripts/hole_study_phase1.py`; figures via `scripts/plot_hole_study.py`;
raw output `results/curated/hole_study/phase1.json`.

Geometry convention: hole diameter 1, radius 0.5, depth = AR, entrance at the
top, etch front at the bottom. Delivery is quoted as a fraction of the flux
entering the hole mouth. The etch-front area does not change with AR, so
delivered *fraction* is proportional to delivered flux per unit front area —
which makes its trend with AR the ARDE observable.

## 1. Geometry import and routing (receipts, not assumptions)

Each aspect ratio was written to a binary STL, re-read from disk, diagnosed and
measured. Every case: watertight, consistently oriented, out-of-roundness
`6.67e-5` against a faceting bound of the same order → **routed to the
axisymmetric operator** by the measured deviation, not by assertion.

## 2. Transport reference — the Clausing benchmark, reproduced in-study

Free-molecular (thermal radical) transmission, exact band algebra:

| AR | transmission | Santeler closed form | relative |
|---|---|---|---|
| 50 | 0.025287 | 0.025434 | 5.8e-3 |
| 100 | 0.012950 | 0.013013 | 4.8e-3 |
| 150 | 0.008710 | 0.008745 | 4.0e-3 |
| **200** | **0.006563** | 0.006585 | 3.3e-3 |

0.656 % at 200:1 — the study reproduces its own quoted benchmark, and agrees
with Santeler inside his stated ≤0.7 % error at every point.

**Documented limit (matters for phase 2).** The *general* body-of-revolution
operator — the one a tapered profile needs — does not certify its own
self-pair quadrature on a straight wall. Measured residuals against its 1e-4
tolerance: 9.29e-4 (band 0.20, az 32, gen 8), 4.39e-4 (band 0.10), 1.18e-4
(band 0.10, az 48, gen 16), then 4.28e-2 as finer bands move the failure to the
adjacent pair (0,1). It plateaus a factor ~1.2 above tolerance and cannot be
refined through. The refusal is **reported, never bypassed**; straight-hole
numbers come from the gated exact-algebra path. Tapered-profile work needs an
analytic self-pair kernel first — the same treatment the cylinder path already
has.

## 3. Line-of-sight acceptance vs the measured tail

Cone acceptance of the two-component beam at the declared reference energy
(1500 eV, inside the measurement's 1.4–2.0 keV band). Acceptance half-angle at
AR 200 is 0.2865°.

| tail fraction | cone @ AR 100 | cone @ AR 200 | sidewall share @ AR 200 |
|---|---|---|---|
| 0.00 (core only) | 0.9669 | 0.5736 | 0.4264 |
| 0.35 | 0.7095 | 0.3951 | 0.6049 |
| 0.50 | 0.5992 | 0.3186 | 0.6814 |
| 0.65 (declared reference) | 0.4888 | 0.2421 | **0.7579** |

The measured tail multiplies AR-200 sidewall interception by **1.78×** relative
to a core-only beam. (The plan's table quotes 0.140 → 0.607 for the same model
at 3465 eV; widths scale as `E^-1/2`, so the lower reference energy here gives a
wider beam and a larger core-only sidewall share. Same physics, declared
energy — the energy must always be quoted with an acceptance number.)

## 4. Where the etch front actually gets its energy

Deterministic specular cascade, production reaction rule verbatim
(`react = clip(0.9·kress(cos), 0, 1)`, continuing weight `1 − react`, Eq. 2.34
retention), production bounce cap 8. Flux closure — direct + cascaded + wall +
thermalised — holds to **≤2e-14** at every point.

AR 200:

| tail fraction | direct ions | cascaded hot particles | total | cascaded share |
|---|---|---|---|---|
| 0.00 | 0.1754 | 0.7707 | 0.9461 | 81.5 % |
| 0.35 | 0.1197 | 0.6886 | 0.8083 | 85.2 % |
| 0.50 | 0.0958 | 0.6534 | 0.7492 | 87.2 % |
| 0.65 | 0.0720 | 0.6182 | 0.6901 | **89.6 %** |

**Headline.** At 200:1, line-of-sight ions supply only 7–18 % of the energetic
flux reaching the etch front; **82–90 % arrives by wall reflection.** A model
without the reflection cascade is not wrong by a correction factor at this
aspect ratio — it is missing the dominant delivery channel.

**Second headline — the front is energy-rich and radical-poor.** Thermal
delivery at AR 200 is 0.656 % (section 2) against 69–95 % energetic delivery:
a ratio of **105–144×**. Any chemistry that needs a neutral flux at the bottom
of a 200:1 hole is starved by two orders of magnitude relative to the energetic
channel, which is the regime the surface model must be trusted in — and it is
outside the validated envelope (see appendix).

**Bounce-cap sensitivity, declared.** Raising the cap from the production 8 to
64 moves AR-200 cascaded delivery 0.6182 → 0.6283 (+1.6 %) at tail 0.65, and
the generation counter *saturates* at 64 — the cascade is not converged in
bounce count at extreme AR. The production number is the conservative one;
the +1.6 % is a floor on the truncation error, not a bound.

## 5. ARDE — and why the measured tail is not a correction but a sign change

Total energetic delivery vs aspect ratio (production cap):

| AR | tail 0.00 | tail 0.35 | tail 0.50 | tail 0.65 |
|---|---|---|---|---|
| 1 | 0.9997 | 0.9986 | 0.9981 | 0.9976 |
| 4 | 0.9989 | 0.9944 | 0.9924 | 0.9905 |
| 16 | 0.9956 | 0.9775 | 0.9697 | 0.9619 |
| 64 | 0.9826 | 0.9190 | 0.8918 | 0.8645 |
| 200 | 0.9461 | 0.8083 | 0.7492 | 0.6901 |
| **decline 1 → 200** | **5.4 %** | 19.1 % | 24.9 % | **30.8 %** |
| decline 1 → 16 | **0.41 %** | 2.1 % | 2.8 % | 3.6 % |

Every configuration is monotone declining, so none trips a formal anti-ARDE
flag — but that test is weak here by construction: delivered fraction starts at
≈1 and is bounded above by 1, so it cannot rise the way an *etch rate* can. The
informative quantity is the magnitude, and it is decisive:

- **Core-only beam: 0.41 % decline to AR 16, 5.4 % to AR 200.** Essentially no
  aspect-ratio dependence. This is the hole-geometry counterpart of the trench
  arc's anti-ARDE finding (`074068e`: delivery falling only 3 % to AR 16 while
  floor recession *rose* 27 % over AR 0–4). A beam this narrow lets the cascade
  return nearly everything to the front at any depth.
- **Measured two-component beam: 3.6 % to AR 16, 30.8 % to AR 200** — a factor
  **5.7× more ARDE** at 200:1 than the core-only beam.

The tail is what makes deep-hole delivery physical: it puts 60–76 % of the beam
outside the acceptance cone, so shadowing has something to shadow. This is the
study's strongest single statement about why the beam model matters.

**Gate on quoting profile numbers:** a configuration whose delivery curve is
flat (core-only, ≤5 % over the whole range) must not be used to quote depth or
profile predictions — flat delivery means the depth trend is set by chemistry
alone with no transport contribution, which is the failure mode the trench arc
diagnosed. Only the tail-bearing configurations pass this gate.

## Honesty appendix

**Validated.**
- Transport geometry at 200:1 against the analytic Clausing solution — 0.656 %,
  agreeing with Santeler inside his ≤0.7 % error at every aspect ratio here.
- STL import chain: exact-to-tessellation SDF, watertightness diagnostics, and
  a *measured* out-of-roundness driving the routing decision.
- Beam widths against a 0.1°-resolution measurement (Kim 2025), and the
  digitised Krüger Figure-4 marginal round trip.
- Cascade bookkeeping: flux closure to 2e-14; production reaction rule reused
  verbatim rather than re-derived.

**The ARDE number in section 5 is transport-only, and is a LOWER BOUND on
etch-rate ARDE.**  Added after the first-principles cascade audit
(`RESULTS_CASCADE_FIRST_PRINCIPLES_2026-08-05.md`, `fd61bb7`); it corrects a
reading this document previously invited.

- What the 30.8 % decline *is*: the fall in delivered **energetic** flux to the
  front, from shadowing plus cascade attenuation. Independently cross-checked —
  the 1-D analytic cascade law reproduces the 3-D gather to ~2 % at AR 200
  (0.679 vs 0.690 at tail 0.65), two different implementations and geometries.
- What it is **not**: a predicted etch-rate ARDE. Converting delivery to rate
  requires the surface model, and that model is currently *exactly ion-limited*
  (rate linear in ion flux, <1 % sensitivity across a 50× neutral range —
  `RESULTS_FLOOR_DELIVERY_2026-08-05.md`). An ion-limited rate law fed by a
  nearly AR-flat energetic supply cannot produce much ARDE by construction.
- The magnitude of what is missing: the one clean literature anchor in the local
  corpus is Huang's MCFPM, whose **etch rate falls 80 % by AR 40**
  (`huang_thesis.txt` L5578-5580) under his conditions. Our transport-only
  decline at AR 40 is ≈10 % (tail 0.65). The gap is not a transport error — the
  transport is graded against exact theory — it is everything that converts
  delivery into rate at depth: neutral conductance collapse (we *measure*
  thermal delivery at 0.656 % of source at AR 200, a ~150× starvation, which an
  ion-limited chemistry ignores), redeposition, tapering, and charging.
- Consequence for use: quote section 5 as **energetic delivery vs aspect ratio**,
  not as etch-rate ARDE. An etch-rate prediction at high AR is blocked on the
  limiting-regime question, which is the single defect that also owns the
  validation case's +29 % depth overshoot.
- Caveat on Huang's numbers as a reference at all: his published *ion* decay
  (2.0→0.3e15) sits below the fixed-beam straight-wall acceptance bound
  (13/53 = 0.245 vs his 0.150), so his curve embeds profile evolution, mask
  erosion, via geometry and re-arrival counting. It anchors an order of
  magnitude here, nothing finer.

**Extrapolated.**
- Surface chemistry beyond AR ≈ 17, where no matched profile data exists. The
  validation case additionally runs +29 % on depth, located in the ion-energy
  channel magnitude — a channel the source mechanism's own author flags as
  probably overestimated in his thesis.
- The energetic:thermal delivery ratio of 105–144× at AR 200 puts the surface
  model far outside the flux regime it was validated in.

**Not modelled.**
- In-feature charging. For a conventional plasma source this bounds every
  deep-AR claim here; the preregistered gate ladder is
  `RESEARCH_CHARGING_DEEP_AR_VALIDATION_2026-08-04.md` (hard ceiling: in-feature
  potential never exceeds the maximum ion energy and is nearly AR-independent
  above AR ≈ 17). A neutral-beam source is largely exempt from this term.
- Redeposition of etch products inside the hole.

**Declared open.**
- Tail fraction — swept as a band {0.00, 0.35, 0.50, 0.65}, never fitted. Its
  derivation from a sheath collision operator is the S2 rung of
  `RESEARCH_IADF_SUBDEGREE_AND_REACTOR_2026-07-29.md`.
- Cascade bounce cap: production 8, unconverged at AR 200 (+1.6 % at cap 64,
  generation counter saturating).
- The general body-of-revolution operator's self-pair quadrature (section 2) —
  blocks tapered-profile transport until an analytic kernel replaces it.

**Not run.**
- Profile evolution (the plan's series 3). The axisymmetric path is
  transport-only — there is no evolution driver on it — and the 3-D engine at
  200:1 with ≥10 cells across the hole needs ≥2000 cells of depth, which is not
  affordable at the resolution the mouth region demands. Phase 2 needs either an
  axisymmetric evolution driver or a declared reduced-resolution study; the
  transport and delivery results above are inputs to it, not substitutes.
