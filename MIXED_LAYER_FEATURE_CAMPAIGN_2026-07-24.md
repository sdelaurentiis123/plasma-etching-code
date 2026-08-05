# Mixed-layer v1 at feature scale — first campaign record (2026-07-24)

Box 45735262 (RTX 3090), epoch through e4c5115, 10 nm grid, authority operator
(deterministic extruded exchange + common_refinement + continue_gas_cavity).
Audits: `results/curated/mixed_layer_feature_v1/`. No fitted knobs anywhere in
the chemistry; the single anchorable constant (k_v) was probed, never fitted.

## What the probes discovered (each starvation → a published-mechanism completion)

| probe result | diagnosis | completion (commit) |
|---|---|---|
| 0.85 nm/s at k_v=40 (40x capacity bought 2.8x) | no atomic-F in boundary; film-only F delivery caps at 4.4 nm/s | chemisorption channel, Krüger complex-formation p=0.2729/0.2 (a1cc77b) |
| 3.4 nm/s ceiling, k_v-insensitive | uniform sticking blankets etch front with film | substrate-dependent split: on-polymer 0.1 vs on-substrate 0.002 (7f774a5) |
| opening frozen at 90 | film growth not moving the boundary | growth velocity = -d(film)/dt via router growth field (dcfff45, 06b1136) |
| seam tear at step 4 | gas sheet nucleating at exactly-on-grid mask/oxide interface | interior-gas-nucleation invariant + receipt (b4fcef2) |
| refusal at overhang | near-grazing quadrature atom samples zero landings | drop immaterial (<=1e-3) unsampled measure with receipt (e4c5115) |

## Results (published constants only, k_v as labeled)

| run | t (s) | depth (nm) | opening (nm) |
|---|---|---|---|
| probe k_v=1 | 6 | 17.0 | **45.66** |
| probe k_v=4 | 6 | 19.0 | **44.07** |
| base k_v=5 | 33.4 | 33.4 | **0.000 (clogged)** |

Experimental base (Krüger fig. 7, 60 s): opening 45, depth 825.

## Findings

1. **Mouth narrowing 90 → ~45 nm emerges** from published constants (fresh-film
   sticking on the growing mask film). The R5 model needed the fitted
   crosslinked-growth-fraction knob for this number.
2. **But it does not equilibrate**: the mouth passes through 45 at ~6 s and
   seals at ~33 s. Deposition at the lip outruns film sputter ~50:1 at
   published fluxes, so nothing stops the closure.
3. **The missing mechanism is named and was already on the roadmap**:
   ion-induced crosslinking (Bruce/Graves ion-processed-skin dose
   competition). The lip receives heavy ion dose → crosslinked film → radical
   attachment drops ~3-5x (Krüger's own crosslinked-attachment law — the
   physics his resolved mask lattice carries and the retired 0.8765 knob
   lumped). Ion-shadowed deep sidewalls keep fresh-film attachment. That
   differential is what makes 45 nm an equilibrium rather than a waypoint.
4. Depth at these k_v values is mouth-transient-dominated; k_v anchoring is
   deferred until the mouth equilibrium physics is in.

## Next (in order)

1. Crosslink state in the film reservoir: crosslinked fraction x_c evolved by
   local ion dose (dose competition per RESEARCH_SYNERGY_CEILING /
   Bruce-Graves), attachment = blend(fresh 0.1, crosslinked ~0.02-0.03 — both
   published in Krüger's mechanism). No new fitted constants.
2. Rerun base 60 s: gate = mouth equilibrates in [40, 50] nm AND stays open to
   60 s. Then k_v Newton on depth, then the 8-condition scorecard.
3. If the mouth still misbehaves, the declared suspect is the ion-spectrum
   compression omission (grazing sputter enhancement at the lip).

## Campaign 2 (2026-07-25, box 45757468) — audit-corrected published-constant runs

Method change (user directive): stopped one-mechanism-at-a-time iteration; an Opus
completeness audit (RESEARCH_MECHANISM_COMPLETENESS_2026-07-25.md) tabled Krüger's
full mechanism against the implementation and caught a mis-lifted constant — film
sputter used 0.1384 (the COMPLEX sputter probability); the published polymer law is
p0=0.9, eth=20 eV, q=0.5, e0=500 eV, Kress B=9.3, at incident energy.

| run | mechanism state | 60 s endpoint |
|---|---|---|
| ml4 | + crosslinking (weak 0.1384 film sputter) | depth 66.7, opening 30.2, OPEN |
| ml5 | + Kress angular on the weak law | depth 59.2, opening 0.0 (reclogged) |
| ml6 | audit-corrected: published film law + two-state oxide, volatilization_yield=1.0 | **depth 811.1, opening 18.8, OPEN** |

**ml6 is the headline: depth within 1.7% of experiment (825) with ZERO anchored
constants** — every chemistry number as-published (the R5 campaign needed a fitted
0.587 yield scale + calibrated mask fraction for the same observable). Remaining
misses localize to the mouth-region energy budget: opening 18.8 vs 45, mask-top
erosion 133 nm vs ~0. The audit's P0.2 names the unwired mechanism: grazing-ion
specular reflection into hot neutrals (0.95 reflection, 0.90 energy retention) —
reflected flux bombards the opposite lip/upper sidewalls (opens the mouth, spares
the mask top by removing normal-incidence overkill... to be verified), and its
absence also currently overfeeds lip-face energy. Next single change: wire
GrazingSpecularIonReflection3D into the pilot ion path; then the scorecard.

## Campaign 3 (2026-07-26, box 45874720) — atom-resolved chemistry

| run | config | 60 s endpoint |
|---|---|---|
| ml9a | atoms, no reflection | depth 790.8, opening 22.5, OPEN — **checkpoint gate [780,870] PASSES** |
| ml9b | atoms + reflection | depth 170.9, opening 0.00 (sealed) ≈ ml8 |

Findings: (1) per-event chemistry preserves the zero-anchor depth result
(Jensen correction 811→791, −2.5%, physically sensible) and slightly improves
the opening (18.8→22.5). (2) The reflection regression is NOT a compression
artifact — two independent formulations agree. The redistribution itself
seals the mouth: P=0.95 at grazing strips ~all sputter from the mouth film
while secondaries fly deep. Pattern across ml4 (no kress, full grazing
sputter, mouth 30) / ml5+ml9b (grazing sputter suppressed, sealed) / ml9a
(kress only, 22.5): mouth survival is governed by near-grazing wall-film
sputter. Krüger holds 45 nm with BOTH kress and reflection active, so a
compensating mouth-region removal channel is still missing — candidates from
the audit queue: low-energy activation states (P2 #5: activated polymer
5-30 eV, p0=0.3), energetic FC-ion polymer sputter contribution, or
reflection's angular threshold applying differently on polymer. Next step is
receipt-level analysis of the archived checkpoints (ml9a vs ml9b per-face
energy budgets at the mouth), not further mechanism guessing. Reflection
stays OFF; ml9a is the new checkpoint of record.

## Scorecard 1 (2026-07-26, box 45918618) — ml9a config, 8 conditions, zero anchors

All runs: atoms, reflection off, volatilization_yield=1.0 (pure published).
Audits: results/curated/mixed_layer_scorecard_1/.

| condition | depth (nm) | opening (nm) | criterion | verdict |
|---|---|---|---|---|
| base (6 kW, O2 1.0) | 790.8 | 22.5 | depth 825±5% / opening 45±5 | depth **PASS** (−4.1%) / opening MISS |
| O2 0.5 | 243.1 | **0.000 sealed** | clog preserved | **HIT (exact)** |
| O2 1.5 | 475.7 | 15.2 | depth rank max at 1.5 | MISS (dip; mouth-starved) |
| O2 2.5 | 672.8 | 17.5 | necking absent (open) | **HIT** |
| 4 kW | 531.3 | 18.4 | r(4/6) ∈ [0.84, 0.94] | 0.672 — MISS |
| 8 kW | 858.0 | 22.4 | r(8/6) ∈ [0.97, 1.06] | 1.085 — near-miss (blind R5 was 1.21) |

Reading: the zero-anchor chemistry wins the absolutes and the topology
(base depth −4.1% with NO calibration; clog and necking-absent both emerge
from element ledgers with the O-saturation knob refused by construction) and
loses the graded transfers — and every graded miss carries the same
signature: the over-narrowed mouth (22.5 vs 45) modulates feature flux, so
depth trends inherit mouth error on top of true chemistry response (4 kW:
softer ions → even narrower mouth → double-counted depth loss; O2 1.5:
narrowest mouth of the sweep → dip). One residual, five symptoms.
Comparison: K24-DEKNOB-1 (reduced + derived law, 2 calibrated constants) won
the power ratios (0.928/1.019) but needed the calibrated base pair and a
declared O-saturation constant. The two models fail in complementary ways
that both point at mouth-region physics (activation states / FC-ion polymer
sputter / polymer reflection threshold — audit P1/P2 queue) as the next
mechanism, with receipt-level checkpoint analysis (ml9a vs ml9b) first.

## Campaign 4 close (2026-07-29) — verbatim-complete run and the honest verdict

| run | config | 60 s endpoint |
|---|---|---|
| ml13 | exact cascade, paper-set constants | **depth 852.1 (IN GATE), opening 24.8** |
| ml14 (partial, t=40) | + mask-AC armor | mouth ~19-20 trajectory; mask-armor != mouth lever |
| ml15 | COMPLETE verbatim converged set + activated-SiO2 + de-crosslink | depth 648.6 (below gate), opening 7.2 |

Verdict: with every Appendix-B row implemented and the thesis's own converged
constants, our engine UNDERPERFORMS our best partial configuration. The
chemistry transcription is now closed (nothing absent); the residual is
engine-fidelity class: per-particle MC with fully resolved IEAD azimuth and
voxel surfaces (his) vs deterministic compressed-azimuth level-set (ours),
plus declared omissions (activation energy windows, EX1 angular form).
CONFIGURATION OF RECORD: ml13 (depth in gate with reflection active, zero
added constants; mouth 24.8 vs 45 = documented known limitation). The mouth
gap moves from "missing chemistry" to "transport/representation fidelity" —
a categorically different problem, to be attacked with engine-level work
(azimuth resolution, voxel-vs-levelset comparison), not more constants.

## Campaign 5 (2026-08-02) — corrected beam (P1a sqrt-2 lift): the deconfounding

| run | config | 60 s endpoint |
|---|---|---|
| ml16a | fig-7 single-feature constants + corrected lift, dx 10nm | depth 590.5, opening 11.1, mask 850.25 (**mask exact**) |
| ml16b | ml13 paper-set constants + corrected lift, dx 10nm | depth 638.8, opening 12.25, mask 850.2 |
| ml16c | fig-7 constants at dx 5nm | OOM-killed on 62GB box (needs bigger RAM or domain redesign) |

Honest verdict: the geometrically-exact sqrt-2 beam correction REGRESSES both
headline metrics for every constant set — ml13's depth 852 was partly riding the
too-narrow beam (over-delivered floor flux), and the wider beam feeds the lip
film cascade more than it sputters it (mouth seals late). TRUE state at correct
transport: depth ~640 vs 825, mouth seals vs 39-45 neck. Mask physics exact.
Research verdicts (committed bd58ffb): mouth = polymer dep/removal balance at
lip (facet/charging/redep refuted); measured NER runs ABOVE cosine to 50-60 deg
(You 2023) where our lip apparently nets deposition from added flux; Izawa
sidewall CFx sticking 0.004 vs our ~0.1 class [VERIFY]; Krueger's 45 was a
FITTED optimizer target, his sim/SEM NECK minima 38.8/39.0 nearly agree; his
1nm voxels vs our 10nm remains open (ml16c redo needed). NEXT (in order): lip
removal-vs-angle audit against the cosine reference (local, free); sidewall
sticking definitions check; 5nm rerun on high-RAM box; adopt Top/Neck/z_neck
metrics.

## Campaign 6 (2026-08-04) — dx=5nm resolution verdict

| run | config | result |
|---|---|---|
| ml17a-dx5-short | fig-7 constants + corrected lift, **dx 5nm** | reached t=0.386s (500s/step); graded by matched-simulated-time vs ml16a |

Verdict (RESULTS_DX5_RESOLUTION_VERDICT_2026-08-04.md): **resolution is not the
driver of the top-band defect.** Floor etch rate grid-converged to 0.05%
(22.0 nm/s both grids); mouth closure only 8-11% slower at 5nm, Richardson
dx->0 limit 15-22% below the 10nm value, against a 10.8x (+980%) discrepancy —
refinement recovers under 3% of the gap. Throat pins at the mask top on both
grids from step 1. Krueger's 1 nm voxels are therefore acquitted as the
explanation. The ml16 depth undershoot is likewise not resolution: the floor
rate agrees to 0.05%, so depth is inherited from the throttled aperture.

Remaining unexplained: the late-time **throat reversal** — the 10nm reference
descends normally to 180 nm by t=30 s (opening 23.5) and then reverses to
130 nm while the aperture collapses to 11.1. At t=12 s petch reads 38.4 nm vs
the experimental 39.0 nm neck, i.e. the early trajectory is correct and the
defect is a late-time takeover by a second, higher constriction. That is a
profile-evolution question, not a per-face budget question.

Deployment lesson: three prior 5nm launches died on the extrusion guard *after*
it was fixed — stale module on the box from patching files onto a live tree.
Clean `git archive` deploy passed on the first attempt. Never patch a live tree.

Method: matched-simulated-time comparison against an archived trajectory cost
~2.5 h box time (~$0.50) for what a 12 s endpoint would have needed 17 h.

## Campaign 7 (2026-08-04) — deposition-driven crosslinking closes the early transient

`ml18-depxl-12s`, 12 s (not 60: 88 % of the defect completes by t = 8 s),
graded on `closure/etch` by window vs Krüger's run-average 0.0310.
Full record: `RESULTS_ML18_DEPOSITION_CROSSLINK_2026-08-04.md`.

| window | ml18 | × Krüger | ml16a baseline | × Krüger |
|---|---|---|---|---|
| 1–4 s | 0.0569 | **1.83×** | 0.1572 | 5.07× |
| 4–8 s | 0.0495 | **1.59×** | 0.1091 | 3.52× |
| 8–12 s | 0.0367 | **1.18×** | 0.0526 | 1.69× |

Aperture at t = 12: **64.89 nm vs baseline 38.25** (+26.6); closure budget spent
**49 % vs 101 %**; depth **+10.7 %** (246.0 vs 222.1) with no floor-chemistry
change — the wider mouth delivers more flux, as the neck regrade predicted;
mask 850.21 (exact).  All three windows land inside the declared 0.031–0.1 pass
band.  The fix is zero-constant: Krüger creates crosslinks *at deposition*
(thesis §2.2.3, Table 6.2 `P(s)+P(s)→PC(s)+PC(s)`, no `M(g)`) and breaks them
with ions; petch had creation riding the ion channel, which collapses ~200× on
a near-vertical lip — so our lip was the *least* crosslinked surface where his
is the most.  Residual 1.83× is bounded by the untabulated per-material bond
multiplicity (`[VERIFY]`; m = 3 → ≈1.4×, m = 8 → ≈1.0×).

## Campaign 8 (2026-08-05) — ml19 endpoint: the mouth equilibrates

Configuration of ml18 run to 60 s; stopped at t = 46.231 s on a topology
refusal (2-cell mask sliver; robustness, not physics). Full grade:
`RESULTS_ML19_ENDPOINT_2026-08-05.md`.

| observable | target | ml19 | prior best |
|---|---|---|---|
| mask constriction | 45 (his w_m) | **50.9 nm, EQUILIBRIUM** (drift −17 pm/s over 12 s) | sealed (11.1) |
| constriction depth | 200 (SEM) / 271 (sim) | 170.2 nm | 130 (ml16a) |
| mask remaining | 850 | **850.2 — exact** | 850.2 |
| etch depth | 825 @ 60 s | 835.9 @ 46.2 s; rate 1.22× → ~1066 @ 60 s | 590.5 (−28 %) |
| closure/etch (6 windows) | 0.0310 | 1.83× / 1.58× / 0.83× / 0.49× / 0.52× / 0.12× | 5.07× / 3.52× / 1.69× / 1.46× / 1.19× / 0.91× |

**First equilibrium mouth in the campaign.** Every previous run sealed
monotonically; the deposition-crosslinking correction (`3a931b1`, Krüger
§2.2.3 — crosslinks form during deposition, ions *break* them) makes lip growth
self-limiting. Residual +13 % on aperture is bounded by the untabulated bond
multiplicity (`[VERIFY]`).

**The compensating-error pair is now split.** ml16a undershot depth (590) because
its sealed mouth throttled flux; with the mouth open at ~51 nm the floor
over-etches at 1.22× Krüger's average. The mouth error is fixed; the floor rate
is the newly exposed, independent next target (normal-incidence oxide removal).

Follow-up: material-component sliver policy — a 2-cell mask fragment stalls the
adaptive controller (dt → 7.6e-05 s, time frozen). 60 s runs at dx = 10 nm are
not reliably completable until fragments below a declared size are reabsorbed.

## Campaign 9 (2026-08-05) — ml21, the full final mechanism: runs clean, endpoint not bought

First run carrying every landed mechanism: deposition-driven crosslinking with
the published per-material bond multiplicity (`914bb8d`, `924cc04`), the
Appendix-B angular classes on the oxide/mask ion rows (`4c66df1`), sub-resolution
sliver dissolution (`a467fe3`), and the sqrt-2 axisymmetric lift (`6e97ef3`).
Box 46904015, clean `git archive` deploy with all fixes verified post-extraction.

**Robustness: passes.** Zero errors, **zero topology events**, extrusion
projection max deviation **6.9e-18 mesh units**. The stall modes that ended ml19
(t = 46.2 s sliver) and blocked every 5 nm attempt (guard) are both gone.

**Matched-simulated-time vs ml19** (the only comparison a partial run supports):

| t (s) | opening ml21 | opening ml19 | Δ | depth ml21 | depth ml19 |
|---|---|---|---|---|---|
| 0.5 | 87.81 | 87.15 | **+0.66** | 12.06 | 13.27 |
| 1.0 | 86.83 | 85.91 | **+0.92** | 23.05 | 22.74 |

Less mouth closure at equal time, depth unchanged — the direction the bond
multiplicity predicts (harder lip), though the magnitude at t ≤ 1.4 s is small
and not yet a gate.

**Endpoint NOT obtained, deliberately.** ml21 integrates at **dt ≈ 0.057 s
against ml19's 0.15 s** — a 2.4× smaller stable timestep that persists past the
startup transient (still 0.057 at step 20). A 60 s endpoint therefore needs
~1000 steps ≈ 14–17 h, against ml19's 305. Two reasons not to spend it:

1. The timestep collapse is most likely *caused by* the angular classes, whose
   class-1 form (Kress B = 9.3) is **3–4× more peaked than the only in-chemistry
   measurements** (peak/normal 4.17 vs Cho 2000 / Schaepkens 1998 ≈ 1.3) and now
   drives the depth-setting SiO₂ row —
   `RESULTS_CASCADE_FIRST_PRINCIPLES_2026-08-05.md` §7.
2. Independently, ARDE is structurally impossible in the current configuration
   (ibid. §5): an exactly ion-limited chemistry fed by an AR-flat energetic
   supply cannot produce it, so a 60 s endpoint would re-measure the known +29 %
   depth miss at 14 h of box time.

**Next, in order:** (a) free frozen-geometry discrimination of B = 9.3 vs a
measurement-bounded B ≈ 1.7 on the oxide rows; (b) the limiting-regime question
(ion-limited vs neutral-limited), which owns both the missing ARDE and the depth
overshoot; (c) only then the 60 s endpoint and the 8-condition scorecard.

Box destroyed; `vastai show instances` returns NONE. Spend this pass ≈ $0.14.

## Campaign close (2026-08-05) — the stop rule fires

Definitive record: `VALIDATION_DOSSIER_KRUEGER_2026-08-05.md`.

| # | date | what it tested | endpoint |
|---|---|---|---|
| 1 | 07-24 | published constants, k_v probed | opening emerges 90→45.7 @ 6 s, then clogs @ 33 s |
| 2 | 07-25 | audit-corrected film law + two-state oxide | **depth 811.1 (−1.7 %), zero anchored constants**; mouth 18.8 |
| 3 | 07-26 | atom-resolved (per-event) ion chemistry | depth 790.8 (gate PASS), mouth 22.5; reflection regression proven real |
| S1 | 07-26 | 8-condition scorecard on ml9a | absolutes + topology PASS (clog exact); graded transfers miss on one mouth residual |
| 4 | 07-29 | exact cascade, verbatim-complete appendix set | ml13 depth 852 in gate; ml15 verbatim set 648.6/7.2 — appendix ≠ fig-7 set |
| 5 | 08-02 | √2 axisymmetric lift (`6e97ef3`) | deconfounding: ml13's 852 partly rode the narrow-beam bug; true state depth ~640, mouth seals |
| 6 | 08-04 | dx = 5 nm resolution | **resolution acquitted** — floor rate converged to 0.05 %, closure recovers < 3 % of the gap |
| 7 | 08-04 | deposition-driven crosslinking (`3a931b1`) | early transient closes: 5.07× → **1.83×** Krüger, all windows in band |
| 8 | 08-05 | ml19 endpoint | **mouth EQUILIBRATES 50.9 nm** (first in campaign), mask 850.2 exact, depth +29 % exposed |
| 9 | 08-05 | ml21, full final mechanism | robustness passes (zero topology events); endpoint not bought — flagged configuration |

**Stop rule.** The joint ion-channel solve (`953f01c`, `f322cb4`) enumerated the
entire discrete sourced space — 4 axes, 64 combinations — against nine measured
constraints and returned an **empty survivor set**. Two constraints select
assignments uniquely (Gray's beam dynamic range → oxide peak-normalisation;
the measured FC-film curve → the Barklund yield reading); four cannot be
satisfied by *any* sourced combination, and the required depth factor
(0.735–0.812) falls in the **gap between adjacent published options**
(0.628 | 0.868). The depth term is triple-attributed — cascade audit ~1.4× high,
joint solve 1.24–1.35× short, and the source author's own statement that "the
effect of ion energy … might be overestimated in the mechanism" (thesis
L4884-4888) — to a channel petch inherits from Appendix B row-for-row — and the capstone
measurement (`74eb2fa`) decomposes it per channel against Gray's *absolute*
yields: physical sputter **1.22x too strong**, chemically-enhanced **2.8x too
weak**. Wrong in opposite directions, which is why every single-parameter
attempt over-corrected, and why neutral-limited behaviour was structurally
unreachable. Receipted, not applied — both are 350 eV measurements against a
3406 eV front, pending the ZBL-vs-linear scaling question.

Closing it requires *departing* from the published mechanism with new physics
and its own gates (Kwon/Sawin site-limited adsorption,
`RESEARCH_NEUTRAL_LIMITED_REGIME_2026-08-05.md`), not further transcription.
Continuing to permute published options against a defect the source declares
would be fitting by search. **The campaign therefore closes with the bound
declared: depth in band at reference-energy conditions, +29 % at feature keV,
attributed and receipted.**

Delivered: mask exact, mouth equilibrating, neck-band agreement 41.8 vs 41.1,
clog and necking-absent topology exact, four permanent engine defects found by
measurement, and ~$1 of total compute spent.
