# LER demonstrator — assembly plan (2026-08-05)

Predict how a rough mask edge transfers into an etched feature edge, with the
transfer function **measured** rather than assumed. Assembled from pieces
already built and gated; no new physics.

## What already exists (and its gate)

| piece | module | gate |
|---|---|---|
| roughness metrology | `src/petch/ler_metrology.py` | Palasantzas synthesis ↔ periodogram Parseval-consistent; HHCF; Mack noise floor; σ/ξ/α recovery (6 gates) |
| spectral transfer | `src/petch/ler_transfer.py` | H1 cross-spectral estimator; \|T\|² recovered to 6 % max / 2 % mean where coherence > 0.9; intrinsic floor to 2 % band-mean; unity/zero limits exact; Parseval closure 1e-3 (9 gates, `63fcdac`) |
| ion shadowing | `analytic_occlusion` exchange | exact, receipted 1e-9 |
| chemistry | deck + mixed layer | as in the validation record |

The estimator's substantive choice: the naive PSD ratio attributes the entire
etch-added floor to transfer and overstates it by ~270× in the gate's own
measurement; the cross-spectral estimator is asymptotically unbiased for
input-uncorrelated noise.

## Two rungs, cheapest first

### Rung A — structural gate, no digitisation, no box

Constantoudis/Gogolides publish closed-form transfer rules rather than curves
(`RESEARCH_LER_EXPERIMENTAL_GATES_2026-08-05.md` § 2): substrate vs resist 3σ is
linear with **slope ≈ 0.5** above a threshold **σ\*_R ≈ ξ_R/(c·tan θ_R)** with
**c ≈ 2.0–2.5**. That is pure ion-shadowing geometry — exactly what
`analytic_occlusion` computes exactly — so it is gradeable with a static
shadowing calculation over a synthesised rough edge ensemble. No fitting, no
digitisation, no feature evolution.

Deliverable: measured slope and threshold vs the published rule.

There is also a live contradiction in that literature worth settling: one camp
reports ξ and α *almost unaffected* by anisotropic transfer, the other reports
both rising whenever LWR drops. Our estimator returns all three, so the
demonstrator can adjudicate — a publishable result inside the deliverable.

### Rung B — measured |T(k)| gate (one box campaign)

The one place in the open literature that plots a before/after PSD **ratio** for
a named plasma etch — i.e. a measured \|T(k)\|² — with a free null control (a
treated resist whose before/after PSDs superpose, so T ≡ 1)
(`RESEARCH_LER_EXPERIMENTAL_GATES_2026-08-05.md` § 1). Three signed predictions
are preregisterable, and the digitisation trap is documented: the published
curves were noise-subtracted then re-aligned to a synthetic 3σ floor, so a naive
ratio biases toward 1 at high k — subtract the floor or rebuild from the fitted
(σ, ξ, α).

**Campaign sizing** (from that document § 5): y ≥ 2.2 µm at Δy = 2 nm →
N_y ≈ 1100; ~150 solves = 40 probe + 3 × 50-seed ensembles.

**Wide-y is the point.** Every trench run to date uses the extruded quasi-2-D
slab, where y-structure is *suppressed by construction* (the extrusion guard
certifies it). LER runs invert that: the roughness along y is the measured
quantity, so these are full 3-D runs and the guard must be off. This is the
first campaign that exercises that path at scale — budget a shakedown run.

## Intrinsic roughness (what the etch *adds*)

Not yet modelled, and the design anchor is zero-knob: Kushner's counting
numbers (~14 sticking radicals/site/cycle → 25 % thickness variation; ~10
ions/site → 30 %) are 1/√14 and 1/√10 — Poisson statistics, no free scalar.
Constraint from the same literature: plasma alone *smooths* (hills receive more
radical flux than valleys — the same view-factor asymmetry we already compute);
roughening requires micromasking. Any intrinsic model must reproduce that sign.

## Reporting

σ, ξ, α for input and output, \|T(k)\|² with coherence, and the intrinsic PSD —
with the α caveat stated: a white intrinsic floor destroys the self-affine slope,
so α must be read after floor subtraction, never off the raw fit.

## Order

Rung A first (free, and it exercises the whole metrology chain end to end),
then the Rung-B shakedown, then the ensemble campaign. Rung B benefits directly
from the GPU stages in `GPU_PORT_PLAN_2026-07-29.md` — 150 solves is where the
4.8× already landed (`067ed45`) pays for itself.
