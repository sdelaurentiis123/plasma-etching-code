# pargon-2013-jvstb

**Pargon et al. (LTM), JVST B 31, 012205 (2013)**

- **DOI/URL:** JVST B 31, 012205 (2013)
- **Retrieval route:** publisher; peer twin of Azarnouche
- **Status:** ABSTRACT/relay
- **Topic:** ler — Line-edge roughness: metrology, transfer, experimental gates

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | ### 1.4 Sandia / imec / Leti resist-LER-transfer studies - **imec ADI→AEI EUV PSD pairs** — single-step transfer functions on real stacks; the industry design rule *"litho owns low-frequency roughness, etch owns high-frequency (grows ξ)"* is the qualitative shape any `T(f)` must reproduce. Fractilia/Intel SPIE 2024: **material ranking flips between ADI and AEI** → you cannot predict AEI roughness without the etch model. (`RESEARCH_LER_EXPERIMENTAL_SOURCES` §2/§4). - **CEA-Leti / Pargon gate-patterning series** — unbiased PSDs at *every* stack step (resist→BARC→ hard-mask→poly-Si) across ≥5 pre-treatments. The richest per-step measured transfer functions (`RESEARCH_LER_EXPERIMENTAL_SOURCES` §2, ranked target #1). | `RESEARCH_LER_MODALITY_DESIGN_2026-07-24.md`:94 |
| U2 | [unquoted — verify on next use] ### Stage D — LER/PSD stochastic modality (design NOW, build after B) The continuum engine is smooth by construction; LER is fluctuation physics. Design doc in flight (Opus agent): seed measured mask-edge PSDs, add discrete fluctuation sources (polymer clusters, ion shot noise), validate PSD *transfer functions* against the Leti/Pargon series. Nobody has a PSD-validated roughness simulation — open lane. Gate: reproduce one published PSD-transfer measurement blind. | `PROGRAM_ROADMAP_2026-07-24.md`:46 |
| U3 | [unquoted — verify on next use] ## 1. Pargon / LTM-CNRS (Grenoble): the only published *measured* PSD ratio through a named etch | `RESEARCH_LER_EXPERIMENTAL_GATES_2026-07-29.md`:60 |
| U4 | [unquoted — verify on next use] \| **1** \| **Incoming mask-edge roughness (PSD transfer)** \| Full band; dominates low-f (< ~1/50nm) \| `PSD_in(f)` = *measured* resist PSD (σ,ξ,α); transferred by engine `T(f)`. Zero fitted content. \| Naulleau/Gallatin AO 42 (2003); imec ADI→AEI; Pargon per-step PSDs \| | `RESEARCH_LER_MODALITY_DESIGN_2026-07-24.md`:154 |
| U5 | [unquoted — verify on next use] **Rung 2 — Blind protocol on CEA-Leti/Pargon per-step series (sources target #1).** Unbiased PSDs at every gate-stack step across ≥5 pre-treatments. Freeze the engine, seal a held-out subset of treatments/steps, calibrate any *declared* boundary inputs (fluxes/IADF for that recipe) on the non-held-out subset, then **blind-predict the held-out per-step output PSDs**. Boundary data: per-step input PSD + recipe. **Gate (Krüger-style):** predicted held-out `PSD_out(f)` inside a preregistered band across the frequency window, revealed once. Weakness to declare: Cl2/HBr Si chemistry (matches our Si capability, not SiO2/FC) — so run Rung 2 on Si, keep FC for Rung 3. | `RESEARCH_LER_MODALITY_DESIGN_2026-07-24.md`:264 |
| U6 | [unquoted — verify on next use] - **Translational-invariance assumption in (a):** valid for straight lines; an imported STL geometry may have corners/curvature → `T̂(k)` diagonal breaks, need block operator. Fine for v1 (straight line), flag for Rung 4. - **Level-set + stochastic forcing (b) numerics:** band-limited noise vs reinit cadence is unsolved; prototype on a 1-D edge before 3-D. KPZ scaling gives the *target*, not the discretization. - **Metrology self-consistency:** our simulated edge has no SEM noise; the datasets' *unbiased* PSDs already subtracted it. Decide once: compare in the unbiased domain (cleaner) and only forward-model SEM noise when comparing to *raw* (pre-2017) data. Prefer unbiased-domain comparison throughout. - **[VERIFY] all digitized PSD values** from Demokritos/Leti JM3 8, 043004 (2009) and Pargon figures before any gate — the sources doc flags these as figure-digitized, not tabulated. - ** | `RESEARCH_LER_MODALITY_DESIGN_2026-07-24.md`:332 |

_Harvested 1 quoted + 5 unquoted mentions across the repo's docs._
