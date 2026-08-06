# hashimoto-shading

**Hashimoto, electron shading theory**

- **DOI/URL:** Hashimoto, electron shading
- **Retrieval route:** publisher
- **Status:** not-fetched
- **Topic:** charging — Feature charging, notching, electron shading

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | \| **Jinnai, Orita, Konishi, Hashimoto, Ichihashi, Nishitani, Kadomura, Ohtake**, *On-wafer monitoring of charge accumulation and sidewall conductivity in high-aspect-ratio contact holes during SiO₂ etching process*, **JVST B 25, 1808–1813 (2007)**, DOI 10.1116/1.2794050 [ABS] \| 8-in wafer, HARC structures *"comparable with the practical interconnect structures of recent DRAM devices"* \| Verbatim: *"the charge accumulation potential between the top and bottom of the contact-hole structures increased with the aspect ratio of the contact holes"* — i.e. **the experimental V(AR) curve** \| **Yes, and it is the target curve.** Need the figure; AR range not stated in the abstract. **[VERIFY]** the exact AR span. \| | `RESEARCH_CHARGING_DEEP_AR_VALIDATION_2026-07-29.md`:95 |
| Q2 | - **Hashimoto**, *Charge Damage Caused by Electron Shading Effect*, **JJAP 33, 6013 (1994)**, DOI 10.1143/jjap.33.6013 [ABS] (and JJAP **32**, 6109, 1993). Antenna-MOS capacitors under photoresist with HAR openings. Verbatim: *"This damage increased with the pattern's aspect ratio, and occurred even when the test wafer was cut into chips about 5 mm square and mounted on a wafer with insulation."* And: *"The damaging current from this mechanism increased by a factor of more than ten with a decrease in the gate oxide thickness only from 8 nm to 6 nm, implying that the degree of shading depends on the gate charging voltage."* Digitizable: breakdown/damage fraction vs AR; F–N current vs oxide thickness. **This is the founding electron-shading dataset.** - **Arita, Akamatsu, Asano**, *Reduction of Charge Build-Up during RIE by Using SOI Structures*, **JJAP 36, 1505 (1997)**, DOI 10.1143/jjap. | `RESEARCH_CHARGING_DEEP_AR_VALIDATION_2026-07-29.md`:101 |
| Q3 | Hashimoto's model (JJAP **32**, 6109, 1993; **33**, 6013, 1994) [ABS]: photoresist with HAR openings geometrically shades the underlying conductor from *obliquely incident* electrons while transmitting the *normally incident* ion flux; local flux imbalance charges the bottom positive **without any wafer-scale potential difference** (proved by the cut-chip-on-insulated-wafer control). The non-linearity — damage current up >10× for an 8 → 6 nm oxide — implies the shading fraction itself depends on the accumulated gate voltage, i.e. **the shading is self-consistent, not geometric**. Kamata & Arimoto (1996) [ABS] then measured the two currents separately and showed the floating potential difference *grows with T_e* (2 → 4 eV) — the same T_e lever Hwang & Giapis modelled (PRL **79**, 845, 1997; JAP **81**, 3433, 1997 [ABS]: *"Larger values of T_e cause the potential of the upper photoresist s | `RESEARCH_CHARGING_DEEP_AR_VALIDATION_2026-07-29.md`:428 |
| U4 | [unquoted — verify on next use] ### 4.4 Electron shading theory (Hashimoto) in one line | `RESEARCH_CHARGING_DEEP_AR_VALIDATION_2026-07-29.md`:426 |
| U5 | [unquoted — verify on next use] **[ABS] — abstract fetched verbatim (Crossref JATS or OSTI)** 12. K. Hashimoto, *Charge Damage Caused by Electron Shading Effect*, Jpn. J. Appl. Phys. **33**, 6013 (1994). DOI 10.1143/jjap.33.6013. (Also JJAP **32**, 6109, 1993.) 13. T. Kamata and H. Arimoto, *Charge build-up in Si-processing plasma caused by electron shading effect*, J. Appl. Phys. **80**, 2637–2642 (1996). DOI 10.1063/1.363179. 14. T. Kamata and H. Arimoto, *Suppression of electron shading effect by a counter radio frequency bias in plasma etching*, J. Vac. Sci. Technol. B **14**, 3688–3691 (1996). DOI 10.1116/1.588648. 15. N. Fujiwara, T. Maruyama, M. Yoneda, *Profile Control of poly-Si Etching in Electron Cyclotron Resonance Plasma*, Jpn. J. Appl. Phys. **34**, 2095 (1995). DOI 10.1143/jjap.34.2095. 16. N. Fujiwara, T. Maruyama, M. Yoneda, *Pulsed Plasma Processing for Reduction of Profile Distortion Induced by Charg | `RESEARCH_CHARGING_DEEP_AR_VALIDATION_2026-07-29.md`:670 |

_Harvested 3 quoted + 2 unquoted mentions across the repo's docs._
