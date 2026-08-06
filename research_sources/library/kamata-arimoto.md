# kamata-arimoto

**Kamata & Arimoto, JAP 80, 2637 / JVST B 14, 3688**

- **DOI/URL:** JAP 80, 2637 / JVST B 14, 3688
- **Retrieval route:** publisher
- **Status:** ABSTRACT-ONLY
- **Topic:** charging — Feature charging, notching, electron shading

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | 4. **Direct experimental in-feature potential measurements exist, but only up to AR ≈ 6–7 in published detail.** The two instrument families are (a) Samukawa/Ohtake/Jinnai **on-wafer monitoring chips** with top and bottom electrodes inside real SiO₂ contact structures — Shimmura *et al.* JVST A **22**, 433 (2004) [ABS, AR = 5.7], Ohtake *et al.* JVST A **24**, 2172 (2006) [ABS, pulsed vs CW], Ohtake *et al.* JVST B **25**, 400 (2007) [ABS, *IEDF and electron energy at the hole bottom*], Jinnai *et al.* JVST B **25**, 1808 (2007) [ABS, *"the charge accumulation potential between the top and bottom of the contact-hole structures increased with the aspect ratio of the contact holes"*]; and (b) Kamata & Arimoto's **current-through-the-dielectric** measurements, JAP **80**, 2637 (1996) and JVST B **14**, 3688 (1996) [ABS], which give *"The dc self-bias potential difference reached about 100 V | `RESEARCH_CHARGING_DEEP_AR_VALIDATION_2026-07-29.md`:55 |
| Q2 | \| **Kamata & Arimoto**, *Suppression of electron shading effect by a counter radio frequency bias in plasma etching*, **JVST B 14, 3688–3691 (1996)**, DOI 10.1116/1.588648 [ABS] \| Hole patterns, AR ≥ 2, vs substrate rf bias to 400 V \| dc self-bias potential *difference* between HAR hole pattern and open area. Verbatim: *"The dc self-bias potential difference reached about 100 V with a substrate rf bias voltage of 400 V for an aspect ratio of 2."* Also: difference *"increased, independent of the hole pattern aspect ratio, at the lower substrate rf bias voltage and they tended to saturate"* at higher bias \| **Yes** — V(bias) saturation curves. The *saturation* is itself the ceiling physics of §0.1 measured experimentally. \| | `RESEARCH_CHARGING_DEEP_AR_VALIDATION_2026-07-29.md`:91 |
| Q3 | ## Corrections applied to our own record (from the wave) - Belen SF₆/O₂ constants are L3 (profile-fitted per the paper's abstract) → de Boer comparison is a transfer test, not blind. Docs to reword (P6). - "Korean CUDA charging MC" prior-art claim refuted (K-SPEED has no charging). - Krüger's reactor inputs are unvalidated by his own thesis ("treated as the ground truth") — Tier-1 gates must anchor on measurements (Kim 2025, Kamata). - Charging ≠ deep-AR etch stop (redeposition is, per Ohiwa); Matsui AR>7 stop is a negative control we must NOT reproduce. - Citation fixes listed per-doc (Chang/Sawin, Raja&Linne, Kawamura, Krüger PoP, Wang&Kushner, Zhai DOIs). | `ROADMAP_PERFECTION_2026-07-29.md`:41 |
| U4 | [unquoted — verify on next use] \| **Kamata & Arimoto**, *Charge build-up in Si-processing plasma caused by electron shading effect*, **J. Appl. Phys. 80, 2637–2642 (1996)**, DOI 10.1063/1.363179 [ABS] \| Si substrate + 500 nm SiO₂ line-and-space, varied pattern size; Ar ICP 2–40 mTorr \| **Electron and ion currents through the dielectric structure**; floating potential vs pattern size; floating-potential difference vs T_e (2 → 4 eV) \| **Yes** — figures are current/potential vs pattern size and vs pressure. This is the *only* published measurement that separates the electron and ion current channels. \| | `RESEARCH_CHARGING_DEEP_AR_VALIDATION_2026-07-29.md`:90 |
| U5 | [unquoted — verify on next use] **References:** Kamata & Arimoto, JAP **80**, 2637 (1996), DOI 10.1063/1.363179 (electron/ion currents through the dielectric structure vs pattern size; floating-potential difference vs T_e 2 → 4 eV); and JVST B **14**, 3688 (1996), DOI 10.1116/1.588648 (dc self-bias potential difference vs substrate rf bias to 400 V, AR = 2; saturation behaviour; counter-rf-bias suppression). | `RESEARCH_CHARGING_DEEP_AR_VALIDATION_2026-07-29.md`:482 |
| U6 | [unquoted — verify on next use] \| P7 \| **Charging gates D1→D2** (potential ceiling vs Huang&Kushner 2026; electron-shading vs Kamata measured) on existing solver \| 1 wk \| preregistered bands in doc §5 \| Charging §5 \| | `ROADMAP_PERFECTION_2026-07-29.md`:34 |

_Harvested 3 quoted + 3 unquoted mentions across the repo's docs._
