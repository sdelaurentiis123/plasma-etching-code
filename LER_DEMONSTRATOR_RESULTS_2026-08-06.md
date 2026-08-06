# LER demonstrator — results (2026-08-06)

Executes `LER_DEMONSTRATOR_PLAN_2026-08-05.md`. Two rungs: the free structural
gate against the published ion-shadowing transfer rule, and the end-to-end
measurement chain run on physics-generated edge ensembles.

Scripts: `scripts/ler_gate2_shadowing.py`, `scripts/ler_demonstrator.py`,
`scripts/plot_ler_demonstrator.py`. Data and figures:
`results/curated/ler_demonstrator/`. Box spend **$0** — everything local.

---

## Rung A (Gate 2) — the Constantoudis rule graded against exact shadowing

**The rule** (`RESEARCH_LER_EXPERIMENTAL_GATES_2026-07-29` §2.2, from
Constantoudis, Kokkoris & Gogolides, JM3 12(4), 041310 (2013) and the authors'
SPIE Newsroom summary): substrate σ_S vs resist σ_R is linear with **slope ≈
0.5** above a threshold **σ_R\* ≈ ξ_R/(c·tan θ_R), c ≈ 2.0–2.5**; equivalently
reduction requires (σ_R/ξ_R)·tan θ_R > 1/c. Stated mechanism: pure ion
shadowing by the rough, tapered resist sidewall, the substrate inheriting the
envelope of the resist edge. Their Fig. 4 conditions — ξ = 30 nm, α = 0.6,
resist 150 nm, sidewall 86.2°, selectivity 3, etch depth 150 nm — are used
verbatim here.

**What we computed.** A rough mask edge u(y) (Palasantzas synthesis, our Rung-0
module) on a tapered wall, with ions carrying transverse tangents (t_x, t_y) ~
N(0, s²) — the same tangent-plane Gaussian each component of the two-component
IADF uses. A ray from substrate point (x, y) clears the mask iff
`x + z·t_x ≥ u(y + z·t_y) + (h − z)·cot θ_R` for every height z, so the blocking
value is `B(y; t_x, t_y) = max_z [u(y + z·t_y) + (h − z)·cot θ_R − z·t_x]` and the
received flux fraction at (x, y) is the weighted CDF of B over the beam. The
etched edge is the level-quantile of B — exact, closed-form, no bisection and no
fitting. Ensemble: 4 seeds × 1024 points × 1 nm, σ_R swept 0.15 → 9 nm.

### Result: the static transfer is a rigid copy

| beam σ_θ | overall slope dσ_out/dσ_in | min σ_out/σ_in |
|---|---|---|
| 2° | 1.0000 | 1.0000 |
| 10° | 0.9925 | 0.9916 |
| 25° | 0.9885 | 0.9873 |
| 45° | 0.9897 | 0.9886 |

Against the published slope of **0.5**, static shadowing gives **≈1.0**. The
level sweep is equally flat: at σ_R = 2.7 nm, ξ = 30 nm the ratio is 1.000 at
every flux contour tested (0.02, 0.05, 0.15, 0.3, 0.5) and 0.997–0.998 in the
envelope tail (0.7 → 1.0).

**Why, mechanically.** At z = 0 the blocking value is `u(y) + h·cot θ_R` for
*every* direction, because the tilt terms vanish there. A positively tapered
resist line (foot protruding into the space, which is what an 86.2° sidewall
means) therefore presents its foot as the binding aperture, and the foot is
strictly y-local: the shadow boundary copies u(y) with a constant offset. The
y-coupling that could smooth the edge lives only in the upper tail of the
blocking distribution, where a single extreme direction dominates and merely
*translates* the roughness rather than enveloping it.

### The structural scalings do carry the published sign

| sweep | value | min σ_out/σ_in |
|---|---|---|
| correlation length ξ (beam 25°) | 15 nm | 0.9521 |
| | 30 nm | 0.9873 |
| | 60 nm | 0.9981 |
| sidewall θ_R (beam 25°) | 84.0° | 0.9899 |
| | 86.2° | 0.9873 |
| | 88.0° | 0.9846 |

All three dependencies point the way the rule says: reduction grows with the
roughness slope σ_R/ξ_R, grows with wall steepness tan θ_R, and **vanishes
identically without angular spread** (2° beam → ratio exactly 1.0000). The
magnitude is ~1/10 of the published effect (4.8 % maximum reduction vs 50 %).

### Verdict

**The published slope-0.5 rule is not reproducible from static shadowing, and
the shortfall is mechanistic rather than numerical.** Our operator reproduces
the rule's three scaling directions and none of its magnitude. The missing
ingredient is visible in their own stated conditions: **selectivity 3 over a
150 nm etch, i.e. 50 nm of resist erosion** — a moving, eroding mask whose edge
is itself smoothed and whose foot sweeps laterally by 50·cot θ_R ≈ 3.3 nm during
the etch. A frozen geometry cannot express that. Recorded as a bounded negative:
the rule encodes resist erosion coupled to the taper, not shadowing alone.

Consequence for the modality: `analytic_occlusion` transfers mask LER **1:1**,
so any LER reduction petch predicts must come from an explicitly modelled
mechanism (mask erosion, film smoothing, chemistry) — never from geometry
implicitly. That is a useful constraint on the next rung, not a defect.

---

## Rung B (demonstrator) — the measurement chain on engine-generated edges

Ensembles of **16 realizations** × 2048 points × 1 nm pushed through the same
exact shadowing operator, then handed to the H1 cross-spectral estimator
(`ler_transfer.estimate_transfer`) with no knowledge of the geometry.

| case | σ_in → σ_out | ⟨\|T\|²⟩ low band | ⟨\|T\|²⟩ high band | min γ² | intrinsic share |
|---|---|---|---|---|---|
| ξ=15 nm, beam 25° | 2.975 → 2.970 nm (0.9983) | 0.9979 | **0.9781** | 0.9198 | 5.8e-4 |
| ξ=30 nm, beam 25° | 2.972 → 2.972 nm (0.9999) | 0.9999 | 0.9989 | 0.9984 | 5.2e-6 |
| ξ=15 nm, beam 2° | 2.975 → 2.975 nm (1.0000) | **1.0000** | **1.0000** | **1.0000** | **8.0e-17** |

Band split at the correlation frequency 1/(2πξ).

**Three things this establishes.**

1. **Exact null control.** The narrow-beam case returns \|T\|² = 1.0000,
   coherence 1.0000 and an intrinsic floor of 8e-17 — the chain invents no
   transfer and no added roughness when the physics adds none. This is the
   engine-data analogue of the estimator's synthetic unity-limit gate.
2. **The shadowing transfer is a low-pass, resolved spectrally.** In the ξ=15 nm
   case the σ reduction is only 0.17 %, but frequency-resolved it is **10×
   stronger above the correlation frequency than below** (2.2 % vs 0.21 %
   attenuation). A σ-only metric would have missed the shape entirely — which is
   precisely the argument for the PSD-transfer architecture over 3σ ratios.
3. **Coherence separates correctly.** γ² falls to 0.92 exactly where the
   transfer acts (high band, short ξ) and stays at 1.0 where it does not.

Caveats recorded, not smoothed: \|T\|² at the Nyquist bin reads 1.02–1.04, a
finite-ensemble upward bias at the last bin; the ξ/α recovered by
`fit_edge_statistics` is unchanged between input and output at this attenuation
level, so this demonstrator cannot yet adjudicate the ξ/α literature
contradiction (§2.1 of the gates document) — that needs the chemistry-coupled
rung where the attenuation is large.

---

## Honesty appendix

**Validated.** The metrology and transfer modules (15 gates, `63fcdac` and the
Rung-0 set). The shadowing geometry is exact ray-blocking, closed-form in the
quantile, with no fitted content. The null control is exact to 1e-16 on
engine-generated data.

**Demonstrator scale, and what that excludes.** The transfer measured here is
**geometric only**: static mask, no chemistry coupling, no deposited film, no
feature evolution, no mask erosion. It is *not* an etch transfer function. The
plan's Rung B — chemistry-coupled wide-y trench runs with the extrusion guard
off, y ≥ 2.2 µm at Δy = 2 nm (N_y ≈ 1100), ~150 solves — remains unexecuted and
is the next rung. Its shakedown risk is unchanged: every trench run to date used
the quasi-2-D slab where y-structure is suppressed by construction, so that path
has never been exercised at scale.

**Not modelled.** Intrinsic (etch-added) roughness. The design anchor stays
zero-knob — Kushner's counting numbers are 1/√14 and 1/√10, i.e. Poisson, no
free scalar — and the sign constraint from the same literature (plasma alone
smooths; roughening requires micromasking) is unimplemented and ungated. Also
absent: mask erosion, which this pass identifies as the load-bearing mechanism
behind the published rule.

**Declared open.** Whether the ξ/α literature contradiction (Demokritos: ξ, α
almost unaffected; LTM: both rise whenever LWR drops) is chemistry-dependent —
unresolvable at this attenuation, flagged in the gates document as the highest
-value single result the modality could produce.
