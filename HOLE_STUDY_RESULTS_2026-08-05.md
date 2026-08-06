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

**Gate on quoting profile numbers (re-worded 2026-08-06).** A flat delivery
curve is *not* unphysical and does not disqualify a configuration: Lam's own
published operating point is a sub-0.5° 1σ beam, and their Fig. 6(b) shows it
giving a flat ion-flux-vs-AR curve to AR ≈ 60, matching their experiment
(`RESEARCH_EXTREME_AR_FIELD_2026-08-06.md` §3). Our core-only branch (0.41 % to
AR 16, 5.4 % to AR 200) is the same qualitative statement, so the earlier
wording would have disqualified the best experimental reference in the field.

What *is* disqualified is the inference: **a flat-delivery configuration must
not be used to attribute a depth trend to transport.** When delivery is flat the
depth trend is set by chemistry, and quoting it as an aspect-ratio *transport*
effect is the failure mode the trench arc diagnosed. Tail-bearing
configurations carry a transport contribution that can legitimately be
attributed.

**Acceptance framing, sourced.** The PPPL analysis states the conversion
directly: "the tolerance parameter for the angular divergence of the beam is set
by the value of the aspect ratio of the etched features. For the value of 100,
it is roughly 0.5°… resulting in the angle of **0.25° in the laboratory frame**"
(Khrabrov & Kaganovich, arXiv:2604.04214v2, 8 Apr 2026). That is an independent
statement of the acceptance scale this section sweeps.

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
- The magnitude of what is missing (**anchor set corrected 2026-08-06**). Only
  three deep rate-vs-AR numbers exist, and they span a factor of two: Huang's
  MCFPM *simulation* falls **80 % by AR 40** (`huang_thesis.txt` L5578-5580)
  under his conditions; the deepest published *measurement* — Nguyen, Shkondin,
  Jensen, Hübner, Leussink & Jansen, JVST A **38**, 053002 (2020), the CORE
  silicon SF₆/O₂ series — falls **≈43 % by AR 54**; a de Boer/Blauw lineage
  series sits between. Nobody has reconciled them
  (`RESEARCH_EXTREME_AR_FIELD_2026-08-06.md` §2.1), and this document previously
  leaned on Huang alone, which overstated the size of the gap. Our
  transport-only decline at AR 40 is ≈10 % (tail 0.65) — below the whole
  measured band, as it must be, since transport is one term among several. The
  gap is not a transport error (the transport is graded against exact theory);
  it is everything that converts delivery into rate at depth: neutral
  conductance collapse (we *measure* thermal delivery at 0.656 % of source at
  AR 200, a ~150× starvation, which an ion-limited chemistry ignores),
  redeposition, tapering, and charging.
- Consequence for use: quote section 5 as **energetic delivery vs aspect ratio**,
  not as etch-rate ARDE. *(Superseded for the rate question:
  `HOLE_STUDY_RESULTS_PHASE2_2026-08-06.md` §1 supplies the chemistry-coupled
  rate-ARDE curve, which lands inside the 43–80 % measured/simulated band. The
  limiting-regime defect that blocked it was closed by the Gray-anchored ion
  laws, `2f1e218`.)*
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

## Phase-2 note (2026-08-05): E8 is available, the evolution driver is not

`thermalized_radical_return` is now complete and reachable from the feature step
(`RESULTS_E8_COUPLED_2026-08-05.md`, commits `543d8ae`/`f83ede8`): thermalized
cascade weight enters the same per-species neutral ledger that the diffuse
radiosity solve consumes as its direct term, so returned radicals re-emit and
diffuse at their own published sticking rather than depositing where they
thermalized. Measured effect at a mask-dominated trench: delivery to the front
roughly **doubles** versus local deposition (2.2x), reaching a few percent of
floor radicals.

**Why it matters more here than it did there.** Phase 1 measured that at AR 200
line-of-sight ions supply only 7-18 % of etch-front energy while wall-reflection
cascades supply 82-90 %, and that the front is radical-poor by 105-144x. A
feature whose energy arrives almost entirely by bouncing is exactly the regime
where the *thermalized tail of those same bounces* is the candidate radical
source — Huang's ">95 % above AR 10" statement is made for a feature of this
kind, not for a 0.85 µm mask over a 0.09 µm opening. The trench measurement does
not bound the hole result; it only shows that Krüger's stack is the wrong
geometry to test the claim.

**What phase 2 still needs, in order.**
1. **The axisymmetric evolution driver** — unchanged blocker. The axisymmetric
   path remains transport-only; delivery numbers (with or without E8) are inputs
   to a profile study, not substitutes for one.
2. **E8 on, fraction swept** over its physical band [0, 1] — the CFx+ share is
   unpublished for any of these reactors and must stay a declared sweep, never a
   fitted constant. Report the band, not a point.
3. **The two open transport items from phase 1**: the analytic self-pair kernel
   (blocks tapered profiles) and the bounce cap, unconverged at AR 200
   (+1.6 % at cap 64).
