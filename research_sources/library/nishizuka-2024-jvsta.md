# nishizuka-2024-jvsta

**Nishizuka et al. (TEL), JVST A 42, 043003 (2024)**

- **DOI/URL:** 10.1116/6.0003515
- **Retrieval route:** publisher (paywalled)
- **Status:** ABSTRACT-ONLY
- **Topic:** harc-field — HARC / extreme-AR field practice and ARDE measurements

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | - **CD-SAXS**: *CDSAXS study of 3D NAND channel hole etch pattern edge effects and etched hole pattern variance*, Proc. SPIE **12955**, 1295539 (2024), DOI 10.1117/12.3010927. Measures "CD, CD profile, tilt and distortion, specifically for holes at pattern edge (outer holes) and holes inside pattern (inner holes)" and "hole etch behavior, especially versus etch depth, and its impact on the final hole pattern variance" [Q-relay]. Companion: *Inline metrology of high aspect ratio hole tilt using small-angle x-ray scattering*, Proc. SPIE **12053**, 1205312 (2022) — "sub-nanometer precision" on tilt without a structural model [Q-relay]. - **FIB-SEM tomography**: Zhang, Lee, Klochkov, Korb, Sorkhabi, Lan, Pichumani, Tekleyohannes, Wang, Sallis, Ningen, Kim, Teoh, Polubotko, Pirkle, Foca — *3D reconstruction of 3D NAND memory etch profiles using FIB-SEM: identifying variances in etched hole pa | `RESEARCH_EXTREME_AR_FIELD_2026-08-06.md`:420 |
| Q2 | ### B1. Nishizuka, Igosawa, Yokoyama, Sako, Moki, Honda (TEL) 2024 — 3D HARC topography - **JVST A 42(4), 043003 (2024)**, DOI 10.1116/6.0003515. All authors Tokyo Electron. - Measured: **experimental 3D HARC hole profiles** (vertical *and* horizontal cross-sections, i.e. 3D-imaged), reproducing "distortion" and "twisting" (non-axisymmetric features along the hole axis). A cell-based Particle Monte Carlo simulator is fitted to the SEM/3D data via a model-optimization algorithm. Data lives in the profile figures (digitize vertical + azimuthal cross-sections) — this is genuine *experimental profile* data from a tool vendor. - Chemistry: fluorocarbon HARC (dielectric) — **match**. - Blind suitability: **high value but boundary-limited.** Real experimental profiles (E-grade), but the TEL reactor boundary is not fully published → needs a Tier-1 CCP surrogate. Best used as a held-out *profile- | `RESEARCH_INDUSTRY_DATASETS_2026-07-23.md`:94 |
| Q3 | 1. **Nishizuka, Igosawa, Yokoyama, Sako, Moki, Honda (Tokyo Electron), *J. Vac. Sci. Technol. A* 42(4), 6.0003515 (2024), DOI 10.1116/6.0003515.** Full text **not obtained** (AIP paywall). *[abstract, verbatim]*: "we created the models for HARC etch with a **cell-based Particle Monte Carlo topography simulator** by **fitting both vertical and horizontal cross-sectional profiles carefully to the experimental results**. Moreover, we attempted to apply a **model optimization algorithm**. By collaboration of human and the algorithm, modeling engineers can minimize a try-and-error approach… **the distortion and twisting profiles were reproduced very well**." → **The closest published analogue of our optimization-against-a-profile workflow.** `[VERIFY]` whether they report a mouth/top-CD residual. **Worth a paid/ILL retrieval.** 2. **Shen/Lill (Lam) 2023 (§2.5)** — "top CD, bow CD, and taper"  | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:657 |
| U4 | [unquoted — verify on next use] \| 12 \| Nishizuka, Igosawa, Yokoyama, Sako, Moki, Honda (Tokyo Electron Miyagi) — *Precise and practical 3D topography simulation of high aspect ratio contact hole etch by using model optimization algorithm* — **J. Vac. Sci. Technol. A 42, 043003 (2024)**, DOI 10.1116/6.0003515 \| rare experimental *horizontal* cross-sections of HARC holes; twisting/distortion target \| paywalled \| | `RESEARCH_EXTREME_AR_FIELD_2026-08-06.md`:757 |
| U5 | [unquoted — verify on next use] \| **Nishizuka *et al.* (TEL), *JVST A* 42, 6.0003515 (2024)**, DOI 10.1116/6.0003515 \| AIP paywall \| Closest published analogue of our fit-to-profile workflow; may report a mouth residual. \| | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:810 |

_Harvested 3 quoted + 2 unquoted mentions across the repo's docs._
