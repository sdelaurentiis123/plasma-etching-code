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
