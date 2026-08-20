# mack-2006-pattern-collapse

**Mack, "Pattern Collapse" — verbatim restatement of Tanaka's line criterion**

- **Citation:** C. A. Mack, "Pattern Collapse," *The Lithography Expert* column,
  November 2006 (Tutor55, version 8/10/06), Austin, Texas.
- **Primary record / retrieval:**
  `https://www.lithoguru.com/scientist/litho_tutor/Tutor55%20(Nov%2006).pdf`
  (direct curl, no auth)
- **Status:** PRIMARY FULL TEXT ARCHIVED —
  `research_sources/thesis_extracts/mack_2006_pattern_collapse_tutor55.txt`
  (6 pp.). SECONDARY with respect to Tanaka 1993, which it cites as ref. 1.
- **Topic:** the closed-form line pattern-collapse criterion, symbol definitions,
  and the numerical regime for resists

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | Eq. (1)–(2), verbatim: `R = ws/(2 cosθ)` and `F = σH/R = 2σH cosθ/ws` — "For a resist thickness H, the resulting capillary force per unit length of line". | Gives the Laplace/meniscus load. `ws` is the **space width** (gap), not the pitch. Force is per unit line length. |
| Q2 | Eq. (3), verbatim: `δ = (3/2)(F/E)(H/wl)^3`. | Cantilever sway of a line of width `wl` under total lateral load `F` per unit length distributed over the height — algebraically identical to δ = qH⁴/(8EI) with I = wl³/12 per unit length. |
| Q3 | Eq. (4), verbatim: `E/σ ≤ (4√3 Al³/ws)[3 As cosθ + sinθ + √(9 As² cos²θ + 6 As cosθ sinθ)]`, where `As = H/ws` (space aspect ratio) and `Al = H/wl` (line aspect ratio), attributed to "Tanaka calculated this critical point to occur when…". | The collapse criterion actually consumed. Sense: collapse when E/σ is **below** the right-hand side. |
| Q4 | Eq. (5), verbatim small-angle form: `E/σ ≤ (8√3 Al³/ws)[3 As cosθ + sinθ]`, valid "If the contact angle is less than about 80º and the aspect ratio of the space is high". | Note the prefactor doubles (4√3 → 8√3) because the square root is replaced by its leading term. |
| Q5 | Verbatim: "the surface tension of water, which at room temperature is about 0.072 N/m (72 dyne/cm)". | Consistent with IAPWS; IAPWS is the citable source. |
| Q6 | Verbatim: "Resists have been measured to have Young's modulus values in the range of 2 – 6 GPa"; "contact angles between water and resist are often in the 50 − 70º range". | Resist-regime context only; not transferable to TiO2. |
| Q7 | Verbatim: "The worst case (maximum value of 3Ascosθ + sinθ) occurs at an angle of tan-1(1/3As) which is generally less than 10º. Thus, hydrophilic resists produce greater pattern collapse." | Establishes θ→0 as (near-)worst case — the basis for treating θ = 0 as the conservative bound. |
| Q8 | Verbatim: "The two-dimensional model used here is really a worst case since it assumes very long lines and spaces." | Applicability caveat for line→pillar transfer. |

## Use in petch

Carries the Tanaka line criterion into
`research_sources/RESEARCH_PATTERN_COLLAPSE_CRITERION_2026-08-20.md`. Because
Tanaka 1993 itself is paywalled, this is the verbatim source of Eqs. (1)–(5)
there; the algebraic identity of its bracket with Chandra thesis Eq. 5.11 is an
independent cross-check that the transcription is faithful.
