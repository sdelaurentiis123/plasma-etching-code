# dupuy-2015

**Dupuy et al. (2015) SADP LWR**

- **DOI/URL:** Dupuy et al. (2015) SADP
- **Retrieval route:** open
- **Status:** FETCHED (demoted: PSDs shifted to 1)
- **Topic:** ler — Line-edge roughness: metrology, transfer, experimental gates

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | CD-SEM protocol (verified): Hitachi CG4000, rectangular scan, 1024×1024 px, magnification 300 000 (x) / 49 000 (y); pixel 0.44 nm (x) / 2.69 nm (y); smoothing 25 px in x for resist, 7 px otherwise; IA = 800 px, N = 400 measurement points, S = 2 → **L ≈ 2152–2200 nm, Δy ≈ 5.4 nm**; ≥200 images averaged (Dupuy uses N* > 200). ξ values encountered are "generally below 50 nm." | `RESEARCH_LER_EXPERIMENTAL_GATES_2026-07-29.md`:80 |
| Q2 | - **Protocol facts (verified, Lorusso et al. JM3 17(4) 041009 (2018)):** iRP 2018 mandates pixel size **0.8 nm in both x and y**, **no filtering**, unbiasing mandatory, ≥50 uncompressed images; the box-length requirement is tied to ξ — convergence to σ_inf goes as L/ξ, and because ξ "dropped from 35 nm down to 7 nm in the last 10 years" the classic **2 µm box can be relaxed to ~400 nm**; SEM-vs-AFM accuracy offset for *biased* settings varies 1–5 nm, for *unbiased* only 0.6–1.4 nm. - **LWR/LER algebra:** σ_LWR² = σ_L² + σ_R² + 2c σ_L σ_R with c the left/right edge correlation; c ≈ 0 gives the familiar LWR = √2 LER. Jiang, Li, Wang, Jiang, Huang, ICSICT 2012, DOI 10.1109/icsict.2012.6466733 (abstract verified) derive the **LWR ACF analytically from the two LER ACFs plus their cross-correlation**, and find **ξ_LWR decreases as the edge-edge correlation coefficient increases**. - **Load-bea | `RESEARCH_LER_EXPERIMENTAL_GATES_2026-07-29.md`:321 |
| U3 | [unquoted — verify on next use] \| Dupuy et al., SPIE 9428 (2015), hal-01869175 \| HAL full text (CC-BY) \| **read in full** \| | `RESEARCH_LER_EXPERIMENTAL_GATES_2026-07-29.md`:32 |
| U4 | [unquoted — verify on next use] ### 1.5 Companion open dataset: Dupuy et al. SPIE 9428, 94280B (2015) — full SADP ladder | `RESEARCH_LER_EXPERIMENTAL_GATES_2026-07-29.md`:156 |

_Harvested 2 quoted + 2 unquoted mentions across the repo's docs._
