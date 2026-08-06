# kwon-2006-jvsta

**Kwon et al., JVST A 24, 1906/1914/1920 (2006)**

- **DOI/URL:** JVST A 24, 1906/1914/1920 (2006)
- **Retrieval route:** MIT DSpace thesis (bad OCR on 2004 copy)
- **Status:** PARTIAL (Fig 3.4 replot of Gray 1993)
- **Topic:** beam-yields — Beam-measured yields, thresholds and sticking (the L0/L1 provenance floor)

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | \| Kwon & Sawin 2006 "atomic F flux identified as the dominant etch driver" \| `RESEARCH_BEAM_CONSTANTS_ATLAS` entry + abstract, `[Q-relay]` \| body not read for this note; qualitative use only \| | `RESEARCH_F_SUPPLY_BAND_2026-08-06.md`:784 |
| U2 | [unquoted — verify on next use] The depth MISS changed sign when the inherited two-row energy-law anomaly was replaced with Gray's measured laws. It is decomposed, not merely attributed: the physical channel read 1.22x too strong and the chemical channel 2.8x too weak against Gray's absolute yields, and the residual now sits on fluorine delivery to the front. Two declared-open supply items remain (the Kwon/Sawin adsorption element beyond the scalar; the unpublished CFx+ fraction, swept and shown immaterial). E8 was built, gated, and measured immaterial at this geometry (+0.02% at AR 200 over the whole physical band). | `BENCHMARK_CERTIFICATION_2026-08-06.md`:114 |
| U3 | [unquoted — verify on next use] Closing it requires *departing* from the published mechanism with new physics and its own gates (Kwon/Sawin site-limited adsorption, `RESEARCH_NEUTRAL_LIMITED_REGIME_2026-08-05.md`), not further transcription. Continuing to permute published options against a defect the source declares would be fitting by search. **The campaign therefore closes with the bound declared: depth in band at reference-energy conditions, +29 % at feature keV, attributed and receipted.** | `MIXED_LAYER_FEATURE_CAMPAIGN_2026-07-24.md`:339 |
| U4 | [unquoted — verify on next use] **F4. There is an L1 alternative to the Krüger fluorocarbon deck.** Kwon's translating-mixed-layer series (JVST A 24, 1906/1914/1920, 2006) + Yin's beam yields (JVST A 26, 161, 2008) + Guo's 20-reaction SiO₂ table (thesis Table 4.1, 2009) form a fluorocarbon/SiO₂ deck whose coefficients were fitted to **blanket QCM yields**, never to a feature. Retiring Krüger's five profile-fitted constants therefore does not require new measurements — it requires **swapping to the beam-regressed deck and re-running the Krüger feature as a blind gate**. This is a far cheaper retirement path than the one currently written in `chemistry_deck.py`'s docstring. | `RESEARCH_BEAM_CONSTANTS_ATLAS_2026-07-29.md`:61 |
| U5 | [unquoted — verify on next use] Kwon thesis §3.3 (pp. 72-75), the model that reproduces Gray's beam data in Fig. 3.4: | `RESEARCH_NEUTRAL_LIMITED_REGIME_2026-08-05.md`:731 |
| U6 | [unquoted — verify on next use] \| **E4** \| **Stoichiometric enforcement of F (or C) per removal** — 4 F per SiF4, or Huang's 1 C + 1-3 F per SiO2 formula unit \| closes the ledger; makes the F budget, not a rate constant, the limiter at depth \| Kwon p. 76 `SiO2 + 4F -> SiF4 + O2`; Huang App. E rows (`SiO2CF -> SiF + CO2` etc., L10229-10236) \| | `RESEARCH_NEUTRAL_LIMITED_REGIME_2026-08-05.md`:796 |
| U7 | [unquoted — verify on next use] Source: Gray, Tepermeister & Sawin, *JVST B* **11**, 1243 (1993), 350 eV Ar+ on SiO2, replotted as Kwon (ScD, MIT DMSE 2004) Fig. 3.4 p. 76 — floor **0.28** (F/Ar+ -> 0, pure physical sputter), plateau **1.10** (F-saturated). | `RESULTS_ABSOLUTE_YIELD_2026-08-05.md`:7 |
| U8 | [unquoted — verify on next use] \| C1a \| Gray beam dynamic range `Y(F/Ar+ -> 0)/Y(sat)` \| 0.20-0.30 \| Gray, Tepermeister & Sawin, *JVST B* **11**, 1243 (1993), 350 eV Ar+ on SiO2; replotted Kwon (ScD, MIT 2004) Fig. 3.4 p.76 — floor 0.28, plateau 1.10 \| | `RESULTS_ION_CHANNEL_SOLVE_2026-08-05.md`:35 |
| U9 | [unquoted — verify on next use] 1. **The half-rise position is set by `s` and nothing else in this space.** Gray's measured knee implies `s ~ 0.06` on bare SiO2 — a factor **~17 below** petch's 1.0, and not reachable by either enumerated option (1.0 or 0.0, the latter removing the rise altogether).  So the survivor set is empty *by construction*, and the constraint that empties it names its own fix: the site-limited adsorption coefficient of Kwon/Sawin element **E1**, which is a constant petch does not carry.  It is **not** invented here — `s ~ 0.06` is what inverting Gray's published curve yields, and it is recorded `[VERIFY]` against the Gray/Kwon body until the number is read directly rather than inverted. 2. **`s` does not fix the dynamic range** (0.873 -> 0.878 across a 33x sweep). C1a and C1b are controlled by different terms: the knee position by the adsorption coefficient, the floor-to-plateau ratio by the bare | `RESULTS_ION_CHANNEL_SOLVE_2026-08-05.md`:68 |
| U10 | [unquoted — verify on next use] \| C1b Gray half-rise \| **NO** — needs `s ~ 0.06`, space holds {1.0, 0.0} \| the site-limited adsorption coefficient (Kwon/Sawin E1), a constant petch does not carry \| | `RESULTS_ION_CHANNEL_SOLVE_2026-08-05.md`:134 |

_Harvested 1 quoted + 9 unquoted mentions across the repo's docs._
