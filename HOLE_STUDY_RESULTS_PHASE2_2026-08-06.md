# 200:1 HAR hole study — phase-2 results (2026-08-06)

The blocker phase 1 named in its own *Not run* section is cleared: the
axisymmetric path now has an evolution driver
(`src/petch/axisymmetric_evolution.py`), so the study can convert delivered
flux into an etch rate instead of stopping at delivery. Reproduce with
`scripts/hole_study_phase2.py`; figures via `scripts/plot_hole_study_phase2.py`;
raw output `results/curated/hole_study/phase2.json`. Everything is
deterministic, ran on a laptop in minutes, and spent no box time.

Geometry: 90 nm opening (radius 45 nm), aspect ratio = depth / 90 nm, Deck 1
(fluorocarbon/SiO₂, Gray-anchored ion laws), Krüger Table-6.1 mouth fluxes,
Kim-2025 two-component beam at 1500 eV. Rates are read after a 0.2 s
relaxation, so they are the mixed layer's steady response to the delivered
fluxes, not the bare-surface transient.

## The driver, and what it is allowed to claim

**Representation.** The front is a single-valued generator `r(depth)` on a fixed
axial band grid plus a scalar floor depth — chosen because it *is* the input the
Clausing-gated operator already consumes, so nothing sits between the moving
front and the benchmarked transport. Bowing is representable; true undercut is
not, and is declared, not guarded.

**Gates** (`tests/test_axisymmetric_evolution.py`, 15 tests):

| gate | result |
|---|---|
| enclosure in the Clausing limit (reflecting walls, absorbing floor) vs the gated benchmark | **exact**: rel 0 / 2.8e-14 / 8.9e-14 at AR 10 / 50 / 100 |
| enclosure row closure (exchange + escape = 1) | 2.2e-16 |
| band-resolved cascade vs phase-1 bottom delivery | **bitwise**, rel 0.0 at tail {0, 0.65} × AR {10, 50}; AR 200 tail 0.65 = 0.690134 reproduces the phase-1 table |
| cascade weight budget (floor + reacted + thermalised = entering) | < 1e-12 |
| diffuse solve conservation with an E8 birth term | 1e-12 |
| per-step neutral balance during evolution | ≤ 3.2e-15 |
| degenerate limit: 200×-wide short hole vs the 0-D blanket | 2.8264 vs 2.8263 nm/s, **rel 5.6e-5** (gate 1 %) |
| AR-10 smoke evolution | floor advances monotonically, walls passivate, ledgers close every step |

The two transport channels are exact for a straight cylinder, so the driver
measures its own straightness each step and **stops with a receipt** when the
profile leaves the declared envelope rather than extrapolating an operator
outside its gated regime. That stop is a reported observable — see §3.

## 1. Coupled rate-ARDE — the number phase 1 was blocked on

Floor etch rate (nm/s) vs aspect ratio, swept over the declared tail band:

| AR | tail 0.00 | tail 0.35 | tail 0.50 | tail 0.65 |
|---|---|---|---|---|
| 1 | 6.631 | 6.625 | 6.622 | 6.620 |
| 4 | 4.613 | 4.596 | 4.589 | 4.582 |
| 16 | 3.001 | 2.947 | 2.924 | 2.901 |
| 50 | 2.919 | 2.765 | 2.699 | 2.634 |
| 100 | 2.876 | 2.614 | 2.501 | 2.388 |
| 200 | 2.797 | 2.389 | 2.215 | **2.040** |
| **rate(200)/rate(1)** | 0.4218 | 0.3607 | 0.3344 | **0.3082** |

**The rate falls 69 % by AR 200 against a 31 % fall in energetic delivery.**
Phase 1 §5 measured delivery declining 30.8 % and forbade reading it as an
etch-rate prediction because the surface model was then exactly ion-limited.
With the Gray-anchored ion laws (`2f1e218`) the rate responds to radical supply,
and the conversion runs the other way: **chemistry more than doubles the ARDE
that transport alone produces.** The extra decline is neutral conductance —
thermal delivery falls as the Clausing transmission (0.656 % at AR 200), two
orders below the energetic channel.

Two secondary readings:

- **The knee is early.** Most of the fall is done by AR 16 (0.44 of the AR-1
  rate); AR 16 → 200 costs another 0.13. Neutral conductance collapses fastest
  at low aspect ratio.
- **The core-only beam's transport contribution saturates.** Its rate flattens
  past AR 16 (3.001 → 2.797, 7 % over a 12× depth range) while the measured
  two-component beam keeps declining (2.901 → 2.040, 30 %). A flat delivery
  curve is not unphysical — a sub-0.5° 1σ beam genuinely produces one, and that
  is Lam's own published operating point (`RESEARCH_EXTREME_AR_FIELD_2026-08-06.md`
  §3). What it means is that under a narrow beam the *residual* depth trend
  belongs to chemistry, and must not be attributed to transport. The tail is
  what keeps a transport-attributable aspect-ratio dependence alive at extreme
  depth.

**Scale check against the measured band.** Only three deep rate-vs-AR numbers
exist and they span a factor of two: Huang's MCFPM *simulation* falls 80 % by
AR 40 (`huang_thesis.txt` L5578-5580); the deepest published *measurement*,
Nguyen, Shkondin, Jensen, Hübner, Leussink & Jansen, JVST A **38**, 053002
(2020) (the CORE silicon SF₆/O₂ series), falls **≈43 % by AR 54**; a
de Boer/Blauw lineage series sits between. Nobody has reconciled them. Ours
falls **60 % by AR 50** at tail 0.65 — inside that 43–80 % band, and a different
chemistry from either, so this is a scale check, not a gate. The point that
survives the spread: the transport-only ≈10 % decline at AR 40 that phase 1
flagged as the open gap sat *below* every published number, and the coupled
curve does not.

**Acceptance scale, independently sourced.** Khrabrov & Kaganovich (PPPL,
arXiv:2604.04214v2, 8 Apr 2026) give the conversion directly: the angular
tolerance is set by the aspect ratio, "for the value of 100, it is roughly
0.5°… resulting in the angle of **0.25° in the laboratory frame**" — the same
scale as the 0.286° acceptance half-angle this study sweeps at AR 200.

## 2. Thermalised return (E8) is immaterial in the hole too

Sweeping the declared fluorocarbon share of the thermalised cascade over its
entire physical band changes the floor rate in the fourth decimal:

| AR | fraction 0.00 | 0.30 | 0.65 | 1.00 |
|---|---|---|---|---|
| 50 | 2.6336 | 2.6336 | 2.6336 | 2.6336 |
| 100 | 2.3884 | 2.3885 | 2.3885 | 2.3886 |
| 200 | 2.0401 | 2.0403 | 2.0404 | **2.0405** (+0.02 %) |

This is the geometry where Huang's ">95 % of floor radicals above AR 10 are
thermalised CFx" statement was made for, and the deep-hole result agrees with
the trench measurement (`RESULTS_E8_COUPLED_2026-08-05.md`): the source is real
but its *delivery to the floor* is not. A radical reborn on the wall re-emits at
its own published sticking and is re-absorbed by the wall long before it reaches
the front — the same conductance collapse that starves the plasma neutrals
starves the reborn ones. **E8 is complete, conserved, gated, and not the
deep-hole radical source at these conditions.** It stays default-off.

## 3. The straight-wall envelope — and why the tapered operator is not optional

Evolving from a straight hole until the measured straightness deviation reaches
the declared 2 % tolerance:

| initial AR | stop | process time | depth advance | mouth radius lost |
|---|---|---|---|---|
| 10 | envelope | 1.25 s | **+3.99 nm** | 1.24 nm |
| 50 | envelope | 1.50 s | **+3.95 nm** | 0.95 nm |

**The profile stops being straight after ~4 nm of etching.** Wall passivation
near the mouth outruns floor advance by roughly two orders of magnitude in
relative terms: the mouth loses ~1 nm of radius out of 45 while the floor gains
4 nm out of a 0.9–4.5 µm depth. A straight-wall transport model therefore cannot
carry a HAR hole etch to any interesting depth — not by a factor of two, by a
factor of thousands.

That converts phase 1's third declared-open item from a nice-to-have into the
critical path: the general body-of-revolution operator's **analytic self-pair
kernel** (phase 1 §2: the general operator plateaus ~1.2× above its own 1e-4
tolerance and cannot be refined through) is the blocker for every profile
result at depth. The driver is built to consume that operator the moment it
certifies — the geometry, chemistry, ledgers and stopping discipline above are
independent of which transport path supplies the fluxes.

## Honesty appendix

**Validated.**
- Enclosure assembly against the Clausing benchmark in the limit where the two
  problems coincide — exact to machine precision at AR ≤ 100.
- Band-resolved cascade against phase 1, bitwise.
- Degenerate limit against the 0-D blanket, 5.6e-5.
- Conservation: enclosure closure 2.2e-16, per-step neutral balance ≤ 3.2e-15,
  cascade weight budget < 1e-12.

**Extrapolated.**
- Surface chemistry beyond AR ≈ 17, where no matched profile data exists. The
  validation case (`VALIDATION_DOSSIER_KRUEGER_2026-08-05.md`) carries a bounded
  depth-channel residual, decomposed per channel and attributed with receipts.
- The rate-ARDE curve in §1 is a *model* prediction with no matched experiment
  at these aspect ratios; the Huang comparison is an order-of-magnitude scale
  check, not a gate.

**Not modelled.**
- In-feature charging (bounds every deep-AR claim for a conventional plasma
  source; a neutral-beam source is largely exempt).
- Redeposition of etch products.
- A mask layer: v1 is single-material, so §3's necking is the model's own wall
  passivation, not mask-shouldered necking.

**Declared open.**
- Tail fraction — swept, never fitted. There is now a *sourced route to closing
  it*: Kim, Kawamura, Fujitani, Naito, Iino, Fukumizu, Kurihara, Suzuki &
  Toyoda, Jpn. J. Appl. Phys. **64**, 096002 (2025) measure the main-to-total
  component ratio in a dual-frequency Ar CCP and find it "decreased
  exponentially with pressure", so the fraction can be quoted as a *measured
  function of pressure* with a stated reactor caveat rather than a free sweep.
  What is still missing is `f_tail` for a fluorocarbon HARC recipe at the wafer,
  which nobody has published (`RESEARCH_EXTREME_AR_FIELD_2026-08-06.md` gap 7).
- Thermalised-return fraction — swept over its whole physical band (§2), shown
  immaterial at any value.
- The general body-of-revolution self-pair kernel (§3) and the cascade bounce
  cap (phase 1: unconverged at AR 200, +1.6 % at cap 64).
- Undercut geometry: outside the representation, declared not guarded.

**Not run.**
- Evolution to full 200:1 depth. Not a cost limit — a physics-envelope limit
  (§3). It needs the tapered-profile operator first.
