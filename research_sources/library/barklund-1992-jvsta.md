# barklund-1992-jvsta

**Barklund & Blom, JVST A 10, 1212 (1992)**

- **DOI/URL:** 10.1116/1.578229
- **Retrieval route:** reprinted as Chang thesis Fig. 4.16 (digitized there)
- **Status:** VIA CHANG THESIS FIG 4.16
- **Topic:** angular-laws — Angular yield laws (class-1 physical, class-2 chemical)

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | **(a) Magnitude (soft — convention-dependent).** Barklund's ordinate is labelled "Normalized Etch Rate" and his abstract says "angular dependence of the etch rate". Whether the ion flux was tilt-corrected is **not stated in the abstract and I could not obtain the body** (AIP 403; no OA copy). Two readings: - already flux-corrected (Chang calls it "etching yield" in his prose) → measured FC-polymer enhancement is **1.45× at 60-65°**; petch is **2.9× too high**. - raw rate, `R = Y·cosθ` → implied `Y(60°)/Y(0°) = 1.348/0.5 = 2.70`, `Y(65°) = 3.43`; petch is ~1.2-1.5× too high. In both readings petch is **above** the only FC-polymer measurement, which is the same direction `RESEARCH_LIP_CERTAINTY` reached from the SiO₂ proxies (Cho 2000 ≈ 1.3, Schaepkens 1998 ≈ 1.33). `[VERIFY]` the Barklund flux convention against JVST A 10, 1212 body before quoting a single number. | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:354 |
| Q2 | 1. **The class label compresses two different curves.** Krueger's legend gives one "class 1" for oxide, mask and polymer alike, but the measured angular yields differ by material: SiO2 in fluorocarbon peaks 1.30-1.33 (Cho 2000, Schaepkens 1998) while an FC *film* peaks **1.448 at 65 deg** under Ar+ (Barklund & Blom, JVST A 10, 1212 (1992), via RESEARCH_VERIFY_HUNT).  If each row carries its own measured curve, the polymer row is nearly flat and the oxide row is strongly peaked — which is exactly the split the data above wants.  This is the most likely resolution and it is a *measurement* change, not a convention change. 2. **The depth channel has a second defect of its own**, so no single normalisation can satisfy all three columns.  The depth wants a factor ~0.78 where the oxide-only convention delivers 0.54; the residual would sit in the complex-row energy form (the ZBL-vs-published-li | `RESULTS_ANGULAR_CONVENTION_2026-08-05.md`:105 |
| U3 | [unquoted — verify on next use] > **A. M. Barklund and H.-O. Blom**, *Influence of polymer formation on the angular > dependence of reactive ion beam etching*, **J. Vac. Sci. Technol. A 10(4), 1212–1216 > (1992)**. DOI **10.1116/1.578229** (resolved via Crossref this session). | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:311 |
| U4 | [unquoted — verify on next use] Barklund & Blom publisher abstract, verbatim (OpenAlex), final sentence: | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:315 |
| U5 | [unquoted — verify on next use] **Ar⁺ on CHF₃-deposited fluorocarbon polymer** (normalised etch rate, Barklund): | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:333 |
| U6 | [unquoted — verify on next use] Reading 1 is the next change to try, and it is cheap: digitize Barklund & Blom for the polymer row, keep the oxide row peak-normalised, and re-forecast the coupled rate before any run. | `RESULTS_ANGULAR_CONVENTION_2026-08-05.md`:126 |
| U7 | [unquoted — verify on next use] \| C3 \| FC-film angular peak/normal \| 1.30-3.50 \| Barklund & Blom, *JVST A* **10**, 1212 (1992): 1.448 at 65 deg (yield reading) to 2.70 (raw-rate reading); the band spans both because the flux convention is `[VERIFY]` \| | `RESULTS_ION_CHANNEL_SOLVE_2026-08-05.md`:38 |
| U8 | [unquoted — verify on next use] \| **Barklund yield reading** \| **1.448** \| **PASS** \| | `RESULTS_ION_CHANNEL_SOLVE_2026-08-05.md`:95 |
| U9 | [unquoted — verify on next use] \| Barklund raw-rate reading \| 5.131 \| FAIL \| | `RESULTS_ION_CHANNEL_SOLVE_2026-08-05.md`:96 |
| U10 | [unquoted — verify on next use] **48 of 64 combinations are eliminated here, at zero cost, on a measurement.** The measured FC-film curve selects the Barklund yield reading uniquely, which also settles the `[VERIFY]` flux-convention question in the only way consistent with the data: the raw-rate reading is *more* peaked than the Kress form it was meant to bound. | `RESULTS_ION_CHANNEL_SOLVE_2026-08-05.md`:98 |
| U11 | [unquoted — verify on next use] Best row: **`barklund_yield \| kress9.3_peaknorm \| zbl \| unity`**, 3 failures (C2, C5, C6).  Every other combination is worse.  Thermal-F sticking `zero` additionally breaks C1a (no F uptake, so the yield curve is flat), which is the sweep's own check on that axis. | `RESULTS_ION_CHANNEL_SOLVE_2026-08-05.md`:123 |
| U12 | [unquoted — verify on next use] \| C3 FC-film shape \| **YES** — uniquely by the Barklund yield reading \| settles the flux-convention `[VERIFY]` \| | `RESULTS_ION_CHANNEL_SOLVE_2026-08-05.md`:133 |
| U13 | [unquoted — verify on next use] Flagged follow-up with a receipt already in hand: the polymer row has a measured in-chemistry counterpart too — Barklund & Blom, JVST A 10, 1212 (1992), Ar+ on a fluorocarbon film, peak 1.448 at 65 deg — 2.9x below the Kress form now carried there. Changing it needs its own graded run against the lip results. | `RESULTS_LIMITING_REGIME_2026-08-05.md`:111 |

_Harvested 2 quoted + 11 unquoted mentions across the repo's docs._
