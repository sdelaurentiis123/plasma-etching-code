# chae-2003-jvsta

**Chae, Vitale & Sawin (2003)**

- **DOI/URL:** https://doi.org/10.1116/1.1539085
- **Primary full text:** https://www.researchgate.net/publication/249508706_Silicon_dioxide_etching_yield_measurements_with_inductively_coupled_fluorocarbon_plasmas
- **Retrieval route:** author-uploaded full text
- **Status:** PRIMARY FULL TEXT ONLINE + VERIFIED EXCERPT
- **Extract:** `research_sources/thesis_extracts/chae_2003_jvsta_primary_excerpt.txt`
- **Topic:** beam-yields — Beam-measured yields, thresholds and sticking (the L0/L1 provenance floor)

## Verified claims added 2026-08-06

| claim | implication |
|---|---|
| QCM directly measures SiO2 etch/deposition yield versus energy, ion-neutral ratio, and angle | independent surface-validation topology |
| at 500 V, Ar sputter is about 0.25 SiO2/ion and C2HF5 plasma etching about 2/ion | reactive CFx+ and mixed-plasma synergy cannot be replaced by Ar+ alone |
| CF+ is the dominant positive ion for most reported C2HF5 conditions | species identity matters |
| neutral-only deposition is about 30% of total; ion-associated deposition is roughly 2-3 times ion flux | surface film is a dynamic ion-neutral balance |
| PTFE density 2.1 g cm^-3 is assumed for film-mass conversion | not an atomically exact density |
| fitted active-site creation yield is sensitive to monomer mass and omitted energy/flux/film-removal physics | topology clue only; fitted number is not importable |

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | \| **S-E** \| S. A. Vitale, H. Chae, H. H. Sawin, *Silicon etching yields in F₂, Cl₂, Br₂, and HBr high density plasmas*, **J. Vac. Sci. Technol. A 19, 2197–2206 (2001)**, DOI [10.1116/1.1378077](https://doi.org/10.1116/1.1378077) \| Si in F₂, Cl₂, Br₂, HBr \| `Y` vs **ion energy, ion angle, plasma composition**. Stated explicitly as "a database of experimental values needed for feature profile evolution modeling". Key results: `Y ∝ √E` for all four; **Cl₂ and HBr yields nearly identical** (HBr's lower Si etch rate is a *flux* effect, not a yield effect); Cl₂ `Y(θ)` falls rapidly above **60°**, HBr `Y(θ)` starts falling at small off-normal angles \| | `RESEARCH_BEAM_CONSTANTS_ATLAS_2026-07-29.md`:123 |
| Q2 | **Full citation:** S. You, H. S. Yang, D. Jeon, H. Chae, C.-K. Kim, "Controlling Bowing and Narrowing in SiO2 Contact-Hole Etch Profiles Using Heptafluoropropyl Methyl Ether as an Etchant with Low Global Warming Potential," *Coatings* **13**(8), 1452 (2023). DOI 10.3390/coatings13081452. Ajou University + SKKU (SAINT). **Open access; fetched:** `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/tmp/pdfs/coatings2023_bowing_narrowing.pdf` (+ `.txt`). ICP, 13.56 MHz source + 13.56 MHz bias, ACL mask, HFE-347mcc3/O2/Ar. | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:253 |
| U3 | [unquoted — verify on next use] \| Logue, Michael — *Control of Electron and Ion Energy Distributions in ICPs…*, PhD, U. Michigan \| `.../tmp/pdfs/logue_michael_phd_thesis.pdf` \| NOT RELEVANT. Fetched and ruled out. \| | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:88 |
| U4 | [unquoted — verify on next use] ### 2.1 The single best quantitative mouth/neck dataset found: You, Yang, Jeon, Chae, Kim, *Coatings* 13, 1452 (2023) | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:251 |
| U5 | [unquoted — verify on next use] \| `logue_michael_phd_thesis.pdf` (+`.txt`) \| Logue (Michigan) — ruled out. \| 9,601,217 \| | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:792 |
| U6 | **Historical claim retracted:** the old rate-gap memo said Karahashi superseded this source and independently settled the magnitude question. It did not: Karahashi resolves pure CFx+ ions, while Chae/Vitale/Sawin concerns fluorocarbon-plasma surface synergy. Keep this source unimportable until full text and provenance are archived. | `RESULTS_DEPTH_IDENTIFIABILITY_2026-08-06.md` |
| U7 | [unquoted — verify on next use] 1. the common engine completes ordinary topology, state transfer, and recovery unattended; 2. a clean authoritative operator calibrates only declared physical closures; 3. that frozen operator predicts a complete untouched feature-profile transfer set inside honest uncertainty, or produces a decomposed miss with no hidden retuning; 4. at least one charged-profile and one HAR/chemistry campaign subsequently repeat that discipline; 5. petch is no worse than Vienna on shared truth-based cases at matched error; 6. the exact path has a measured bounded runtime and the fast path is certified against it; 7. a partner can supply a documented boundary/material/geometry deck and receive replayable results without repository archaeology. | `UNIFIED_ENGINE_VALIDATION_EXECUTION_PROGRAM_2026-07-17.md`:697 |

_Harvested 2 quoted + 5 unquoted mentions across the repo's docs._
