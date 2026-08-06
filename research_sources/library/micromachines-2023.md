# micromachines-2023

**TU Wien/ViennaPS, Micromachines (2023)**

- **DOI/URL:** Micromachines (2023), TU Wien/ViennaPS
- **Retrieval route:** open
- **Status:** FULL TEXT: research_sources/thesis_extracts/mask_geometry_micromachines_2023.txt
- **Topic:** sf6-si — SF6/O2 on silicon (the partner-relevant arm)

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | ## Gate Config S — short-cycle "smooth" (Tillocher 2021, Micromachines 12,1143, open access, firm) Ultrafast switching: 500 ms SF6 etch / 50 ms passivation, 1000 cycles, 10 µm trench: 1. D = 60.8 µm ± 10%                   (firm text, Fig 4) 2. p = 60.8 nm ± 10%                   (D/N) 3. s ≤ 30 nm                           (upper bound, "residual roughness") 4. ARDE: depth(4µm)/depth(10µm) = 49.8/60.8 = 0.82 ± 0.05 (both firm) | `BOSCH_BENCHMARK_SPEC.md`:51 |
| Q2 | \| **7** \| **Facet angle is an OUTPUT of the local deposition/angular-etch balance, not an input** (Mahorowala & Sawin 2002); **primary facet has little effect on the oxide profile** (Kim/Hudson 2007); a **large** facet angle is equivalent to **no** facet (Micromachines 2023). \| Simulation + abstract-level + open-access sweep. \| Deprioritize "add mask faceting" as the mouth fix. If added, expect it to move the **bow**, and expect weak faceting to put max width *just under the mask* and strong faceting to move it down. \| | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:765 |
| U3 | [unquoted — verify on next use] \| `mask_geometry_micromachines_2023.pdf` (+`.txt`) \| Micromachines 14, 665 — systematic mask taper/facet-angle sweep (ViennaPS). \| 3,340,692 \| | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:797 |
| U4 | [unquoted — verify on next use] \| `mask_geometry_micromachines_2023.txt` \| 18 \| — \| | `RESEARCH_SF6_RELEVANCE_2026-08-06.md`:124 |
| U5 | [unquoted — verify on next use] ### 3.4 The Micromachines 2023 citation, re-verified `[Q]` | `RESEARCH_SF6_RELEVANCE_2026-08-06.md`:192 |
| U6 | [unquoted — verify on next use] At the arm's `Emean = 100 eV`: **R = 150 sits 2.5x ABOVE half-rise**, implying `theta_F = 0.715` — the saturated, neutral-rich branch. Micromachines' independent flux table (§3.4) gives R = 300-550, further above. | `RESEARCH_SF6_RELEVANCE_2026-08-06.md`:262 |
| U7 | [unquoted — verify on next use] **Bosch** - Ayon et al., J. Electrochem. Soc. 146, 339 (1999) via the McVittie NNIN deck — Config R gates (28.2 µm / 434 nm / 140 nm / 250 nm). Deck: https://people.eecs.berkeley.edu/~pister/147fa14/Resources/BoschProc-STS.pdf - Tillocher et al., Micromachines 12, 1143 (2021) — Config S (ultrafast, 60.8 nm/cycle, smooth). Open: https://pmc.ncbi.nlm.nih.gov/articles/PMC8537062/ - Park et al., Micro Nano Syst. Lett. 8, 14 (2020) — scallop measurement protocol. - Ertl & Selberherr, Microelectron. Eng. 87, 20 (2010) — the academic 3D reference (validation-free; our bar). Open: https://www.iue.tuwien.ac.at/pdf/ib_2009/hashed_links/ep4PPErIJjnr4Y_us.pdf - Laermer & Schilp, US 5,501,893 — the process definition. - VLSet-AE, Microsyst. Nanoeng. (2026) — inverse SEM-measurement model (cited + distinguished: we are forward-predictive); its 16-run dataset = future Config T. | `STATE_OF_THE_PHYSICS.md`:98 |

_Harvested 2 quoted + 5 unquoted mentions across the repo's docs._
