# ohiwa-1998

**Ohiwa et al., JJAP 37, 5060 (1998)**

- **DOI/URL:** JJAP 37, 5060 (1998)
- **Retrieval route:** publisher
- **Status:** ABSTRACT-ONLY
- **Topic:** harc-field — HARC / extreme-AR field practice and ARDE measurements

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | 2. **Charging alone does not produce the de Boer-style deep floor collapse — literature agrees with our own STEP-2 refutation.** In the best-documented deep dielectric case (Huang *et al.*, JVST A **37**, 031304, 2019 [FULL], AR 40, tri-frequency CCP, Ar/C₄F₈/O₂), charging reduces mean ion energy at the etch front from 1940 → 1050 eV (−46 %) and lengthens the etch by **only ~30 %** (36 min → 48 min, *"which are in reasonable agreement with experiments"* quoting 40–50 min HVM times for AR 40–50). Etch *stop* appears only in the low-power arm (2.5 kW, before AR = 40). Separately, Ohiwa *et al.* (JJAP **37**, 5060, 1998) [ABS] attribute experimental HARC etch stop to **redeposition of sputtered fluorocarbon**, not charging: *"the redeposition of sputtered species from the fluorocarbon polymer on the hole sidewall induces the etch stop at the bottom of the high-aspect hole"*, and note ions * | `RESEARCH_CHARGING_DEEP_AR_VALIDATION_2026-07-29.md`:31 |
| Q2 | - **Ohiwa, Kojima, Sekine, Sakai, Yonemoto, Watanabe, *Jpn. J. Appl. Phys.* 37, 5060 (1998)** (Toshiba), DOI 10.1143/jjap.37.5060. *[abstract]*: "**the redeposition of sputtered species from the fluorocarbon polymer on the hole sidewall induces the etch stop** at the bottom of the high-aspect hole… the etch stop in a high-aspect-ratio hole is determined by **the balance between the effects of high-energy-species bombardment and etch inhibition of carbon species**." → the canonical statement that **redeposition of FC polymer is a real transport channel** inside the feature, and that the clog is a *balance*, not a threshold. - **Izawa, Negishi, Yokogawa, Momonoi, *Jpn. J. Appl. Phys.* 46, 7870 (2007)** (Hitachi), DOI 10.1143/jjap.46.7870. *[abstract]*, and these are **directly liftable constants**: "Sticking coefficients of radicals on the sidewall have been estimated by comparing the obse | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:487 |
| Q3 | **Supporting/independent:** You *et al.* 2023 (§2.1, measured NDR/NER), Shen/Lill 2023 (§2.5, "polymer deposition on the sidewalls of the mask, so-called necking"), Ohiwa 1998 (§2.6, balance). | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:551 |
| Q4 | - **Lee 2010 (§2.2), Faraday-cage isolated:** "The redeposition of particles sputtered from the mask slope on the contact-hole sidewall **resulted in sidewall necking**." - **Ohiwa 1998 (§2.6):** FC-polymer redeposition induces the etch stop. → If we add redeposition without simultaneously adding grazing removal, **our mouth gets worse.** | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:648 |
| U5 | [unquoted — verify on next use] **Pass criteria:** - **D5a:** `ΔE/E₀` at AR 40 = **46 % ± 15 pp**. - **D5b (channel separation):** hot-neutral flux to the etch front changes by **< 10 %** between charging-on and charging-off. If petch's charging suppresses the neutral/hot-neutral channel, the neutralisation bookkeeping is wrong. - **D5c (etch-time budget):** charging lengthens the etch by **20–50 %**, not by 2×+. - **D5d (NEGATIVE control — must FAIL to reproduce):** petch must **not** reproduce Matsui/Nakano/ Petrović/Makabe's AR > 7, 300 eV etch stop (APL **78**, 883, 2001). At 300 eV incident and AR ≥ 7 in an oxide trench, petch must still deliver non-zero power to the etch front. Reproducing that etch stop is a **failure**, because the manufacturing record (AR > 100 holes) refutes it. - **D5e (scope discipline):** petch must **not** claim charging as the mechanism for the de Boer SF₆/O₂ cryo-Si AR > 20 floor collap | `RESEARCH_CHARGING_DEEP_AR_VALIDATION_2026-07-29.md`:586 |
| U6 | [unquoted — verify on next use] \| **5** \| **Mask-slope redeposition is a real necking channel** (Lee 2010) and **FC-polymer redeposition drives etch stop** (Ohiwa 1998). \| **Measured** (Faraday cage) / measured. \| Do **not** add redeposition alone. If added, it must be paired with #1 or the mouth closes faster. \| | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:763 |

_Harvested 4 quoted + 2 unquoted mentions across the repo's docs._
