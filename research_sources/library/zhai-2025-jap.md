# zhai-2025-jap

**Zhai, Filipović, Chen, JAP 137, 063302 (2025)**

- **DOI/URL:** JAP 137, 063302 (2025)
- **Retrieval route:** publisher
- **Status:** not-fetched
- **Topic:** modeling-sota — Modeling state of the art and competitor codes

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | - **Yook, You, Park, Chang, Kwon, Yoon, Yoon, Shin, Yu, Im**, *Fast and realistic 3D feature profile simulation platform for plasma etching process*, **J. Phys. D: Appl. Phys. 55, 255202 (2022)**, DOI 10.1088/1361-6463/ac58cf [ABS via publisher]. Jeonbuk National University + KRISS + Samsung-adjacent. GPU-parallelised neutral and ion transport, hash-map 3-D level set, adaptive surface meshing, two-layer surface reaction model. Verbatim: *"the speedup ratio, as compared to a single central processing unit (CPU), is approximately 200"*. Applied to HAR oxide etching and MOSFET SiN spacer. Predecessor code is referred to as **K-SPEED**; earlier GPU work presented at GEC 2013 (*GPU based 3D feature profile simulation of high-aspect ratio contact hole etch process under fluorocarbon plasmas*, Bull. APS 2013, GEC MR1.020) and ICOPS 2016 (DOI 10.1109/plasma.2016.7534064). **I could NOT verify th | `RESEARCH_CHARGING_DEEP_AR_VALIDATION_2026-07-29.md`:284 |
| Q2 | ### ViennaPS (TU Wien) — general geometry + GPU flux, but NO charging Sparse-field **level-set** surface (ViennaLS/HRLE) advected by fluxes from **top-down Monte-Carlo ballistic ray tracing** (ViennaRay, Embree on CPU / **NVIDIA OptiX on GPU**). General geometry, no hardcoded features. - **No charging module exists**: its 27 process-model headers contain zero charging/potential/ Laplace/electron code; the SF6O2 paper states the model is kinetic chemistry "without electrostatic charging effects." - **Why (structural):** a level-set is a scalar signed-distance field that "does not hold directional information" and carries no volumetric state (charge, potential). A charging field solve needs a *separate* volumetric representation — hence ViennaCS (a voxel-like cell set) exists alongside it. Level-sets are great at moving surfaces (topology change, pinch-off), bad at holding volumetric charg | `RESEARCH_SOTA.md`:46 |
| U3 | [unquoted — verify on next use] \| **6** \| **Charging.** \| Real for mask faceting/bowing/microtrench (Zhai 2025) and for twisting; **no source found claiming charging sets the mouth width.** \| Low priority for the mouth. \| | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:55 |
| U4 | [unquoted — verify on next use] \| **Zhai *et al.*, *J. Appl. Phys.* 137, 063302 (2025)**, DOI 10.1063/5.0243470 \| AIP 403 (despite OA flag) \| Charging + reflection → mask faceting/bowing/microtrench, validated vs SEM. \| | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:811 |

_Harvested 2 quoted + 2 unquoted mentions across the repo's docs._
