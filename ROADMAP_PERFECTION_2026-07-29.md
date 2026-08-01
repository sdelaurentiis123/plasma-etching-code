# Roadmap to "perfect" — synthesis of the 2026-07-29 research wave

Inputs: RESEARCH_IADF_SUBDEGREE_AND_REACTOR, RESEARCH_CHARGING_DEEP_AR_VALIDATION,
RESEARCH_BEAM_CONSTANTS_ATLAS, RESEARCH_LER_EXPERIMENTAL_GATES, GPU_PORT_PLAN
(all 2026-07-29, committed aee7c14). Baseline: ml13 config of record (depth in
gate, mouth 24.8 vs 45), chemistry transcription closed, deck framework landed,
LER transfer + STL importer landed, suite 1090+.

## The unifying discovery

The three "transport" gaps are one problem: **angular representation of the ion
population**. (a) Our IADF is collisionless mono-Gaussian; measured IADFs are
bi-Gaussian with a sheath-collision tail carrying 65% of ions at Krüger
conditions (his sheath is 30–65% collisional). We produce σ=0.15° where his
tool measures 0.83°. (b) Our transport compresses azimuth (quantified 1.5×
wall-flux deficit). (c) At AR 200 the acceptance half-angle is 0.29°, so both
errors are fatal there (30–70× sidewall under-delivery) — and the missing
wide tail plausibly explains much of the remaining mouth gap at AR ~9.
Charging then rides on top as a *variance* problem: ~0.5% azimuthal asymmetry
in the potential deflects AR-200 ions out of the cone, so deep-AR charging
validation must be statistical (twisting ensembles), not ensemble-mean.

## Ordered build plan

| # | Item | Cost | Gate | Source doc |
|---|---|---|---|---|
| P0 | **Angular-convergence harness (S0)**: error vs polar/azimuth/GH order at AR 30/100/200; test whether ml13 mouth gap = lost wide-angle tail | 0.5–1 d, no physics | self-report | IADF §C |
| P1 | **Two-component IADF (S1)**: (E,θ,φ) object, core+tail, analytic erf acceptance | 2–3 d | Kim 2025 *measured* IADF (0.1° res) | IADF §C |
| P2 | **LER Gate 2 (analytic ion-shadowing)**: petch vs Constantoudis slope/threshold rule, zero fitting, zero digitization | days | tan ω ≈ 0.5, σ* threshold | LER §2 |
| P3 | **GPU Stage 1** (hoist invariant sums; Stage 0 = bincount landed) then Stage 2 CSR+fused kernel after box spike | 1–2 wk staged | bitwise (S1) / 1e-14+run-repro (S2) | GPU plan |
| P4 | **Sheath collision tail derived (S2)**: anisotropic Born–Mayer cross-section; emits fast-neutral NEAD | 1–2 wk, frontier | tail fraction predicted vs Kim/Krüger | IADF §C |
| P5 | **Azimuth fix + angular AMR (S3)**: polar bin ≤ acceptance/5 inside cone only | 1 wk, concurrent w/ P4 | harness convergence; mouth re-run | IADF §C |
| P6 | **Deck 2 (SF₆/O₂, relabeled L3) + schema v2** (provenance_level required; energy-dep sticking; per-ion rows; redep + Arrhenius blocks) | days | bitwise vs belen arm | Atlas |
| P7 | **Charging gates D1→D2** (potential ceiling vs Huang&Kushner 2026; electron-shading vs Kamata measured) on existing solver | 1 wk | preregistered bands in doc §5 | Charging §5 |
| P8 | **LER Gate 1** (measured |T|² w/ turnover): ~150 wide-y solves, needs P3 speed | 2–3 wk + box | PSD ratio + null control | LER §1 |
| P9 | **Deck 3 Cl₂/poly-Si** (86% beam-buildable, 1 fitted mask constant) + blind profile gates (Mahorowala/Levinson) | 2–3 wk | blind profile | Atlas |
| P10 | **Reactor Tier-1 module B** (sheath closure vs Miller-Riley/Edelberg), then module A + gas-temperature balance (T_gas sets core width: 300↔1000 K = 16× AR-200 sidewall flux) | staged | design doc + corrections | IADF §B |
| P11 | **Charging D3–D5** (notch ladder, twisting *ensemble* 49%→12% ablation, deep-AR energy budget + Matsui negative control) — needs GPU MC tracer (first-in-class; K-SPEED refuted as prior art) | months | doc §5 bands | Charging |
| P12 | **AR-200 demonstration**: tail predicted not fitted, angular AMR, axisym operator, charging variance | after P4/P5/P11 | composed gates | all |

## Corrections applied to our own record (from the wave)
- Belen SF₆/O₂ constants are L3 (profile-fitted per the paper's abstract) → de
  Boer comparison is a transfer test, not blind. Docs to reword (P6).
- "Korean CUDA charging MC" prior-art claim refuted (K-SPEED has no charging).
- Krüger's reactor inputs are unvalidated by his own thesis ("treated as the
  ground truth") — Tier-1 gates must anchor on measurements (Kim 2025, Kamata).
- Charging ≠ deep-AR etch stop (redeposition is, per Ohiwa); Matsui AR>7 stop
  is a negative control we must NOT reproduce.
- Citation fixes listed per-doc (Chang/Sawin, Raja&Linne, Kawamura, Krüger PoP,
  Wang&Kushner, Zhai DOIs).

## Standing doctrine unchanged
Zero fitted knobs (provenance ladder makes this precise); receipts and
preregistered gates; one engine, decks as data; symmetry reductions per-case
with certifying guards; primary sources before model invention.
