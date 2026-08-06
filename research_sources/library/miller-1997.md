# miller-1997

**Miller & Riley (1997) sheath model**

- **DOI/URL:** Miller & Riley (1997)
- **Retrieval route:** publisher
- **Status:** not-fetched
- **Topic:** reactor-sheath — Reactor-scale and sheath closure models

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| U1 | [unquoted — verify on next use] ### S4 — Tier-1 module B (sheath-voltage waveform) — *unchanged from the design doc, revised inputs* `reactor_tier1/sheath_voltage.py` as designed (CCP voltage division → V_dc + harmonic `PeriodicSheathVoltage`), but seeded with the **complete** Krüger recipe from §B.3 (P_lf = 8.0 kW, P_hf = 2.5 kW, −500 V on the top electrode, r = 15 cm, gap 4 cm) and with a **sheath-edge** density rather than n_e,max in the Child-law closure. Benchmark the sheath dynamics against Miller & Riley (1997) and Edelberg & Aydil (1999) rather than only against Krüger. **Gate S4 = the revised Gate 1a:** score against the digitized Fig-4 **for energy only** (mean energy 15%, horn split 25%) and against **Kim 2025 for angle** (§S1 gate). Explicitly declare that Fig-4's 0.25° grid cannot grade angle. | `RESEARCH_IADF_SUBDEGREE_AND_REACTOR_2026-07-29.md`:518 |

_Harvested 0 quoted + 1 unquoted mentions across the repo's docs._
