# raja-linne

**Raja & Linne**

- **DOI/URL:** 10.1063/1.1524020
- **Retrieval route:** publisher
- **Status:** not-fetched (DOI corrected from 1.1519941)
- **Topic:** reactor-sheath — Reactor-scale and sheath closure models

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | \| "Analytical model for ion angular distribution functions…" DOI 10.1063/1.1519941 [VERIFY] \| **Wrong DOI.** Correct: **10.1063/1.1524020**, Raja & Linne, *J. Appl. Phys.* **92**, 7032–7040 (2002). \| | `RESEARCH_IADF_SUBDEGREE_AND_REACTOR_2026-07-29.md`:440 |
| U2 | [unquoted — verify on next use] ### S1 — Replace the angular law of the boundary (2–3 days) New `IonAngularEnergyDistribution` object replacing the transverse Gauss–Hermite tensor in `collisionless_sheath_boundary_state` (`src/petch/boundary_state.py:853`). Requirements: - explicit (E, θ, φ) representation with **arbitrary polar resolution** and support out to ≥ 5°; - **two-component** structure (thermal core + non-thermal tail) so each can be integrated against the geometric acceptance **analytically** — an erf per component, exactly as in the §A.7 table — giving exact bottom/wall splits at any AR with no quadrature noise; - keeps the existing phase↔energy pairing (§A.4) so the energy–angle correlation survives; - `T_perp` **derived** from `T_gas` (G1), never a keyword default; - azimuthal dependence allowed (Raja & Linne drift term) so S3 has something to consume. **Gate S1 (measured, not HPEM):** reproduce Kim 2025  | `RESEARCH_IADF_SUBDEGREE_AND_REACTOR_2026-07-29.md`:480 |

_Harvested 1 quoted + 1 unquoted mentions across the repo's docs._
