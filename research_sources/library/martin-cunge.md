# martin-cunge

**Martin & Cunge, plasma smoothing**

- **DOI/URL:** Martin & Cunge — plasma smoothing
- **Retrieval route:** publisher
- **Status:** ABSTRACT/relay
- **Topic:** ler — Line-edge roughness: metrology, transfer, experimental gates

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | **Martin & Cunge, JVST B 26, 1281 (2008)**, DOI 10.1116/1.2932091 (abstract verified): high-density plasmas "do not generate roughness during silicon etching; on the contrary they tend to smooth existing roughness"; 20 nm-high/20 nm-wide Si pillars are rapidly smoothed by Cl₂ and SF₆; the smoothing mechanism is **radical-starved transport — hills receive more etchant flux than valleys**; reported roughening in F-based plasmas is due to **AlF_x micromasking** from chamber walls. | `RESEARCH_LER_EXPERIMENTAL_GATES_2026-07-29.md`:369 |
| Q2 | ### 1.5 kMC / etch-front-roughening theory (the intrinsic-noise anchor) - **Drotar, Zhao, Lu, Wang, "Mechanisms for plasma and reactive ion etch-front roughening," Phys. Rev. B 61, 3012 (2000)** ([APS](https://link.aps.org/doi/10.1103/PhysRevB.61.3012) · [RG PDF](https://www.researchgate.net/publication/235508642_)). (2+1)-D flux-redistribution MC: etch-front roughening with re-emission gives **universal KPZ exponents α ≈ β ≈ z ≈ 1**, matching experiment. This fixes the *scaling* of the stochastic-forcing term — the intrinsic-roughness PSD slope and its temporal growth law are not free; they are the KPZ/re-emission class. Enumerated roughening mechanisms: stochastic noise, shadowing, etchant re-emission, micromasking, ion scattering. - **Kuboi et al., "Insights into different etching properties of CW and ALE for SiO2 and Si3N4 using voxel-slab model," JVST A 37(5), 051004 (2019)** ([AIP] | `RESEARCH_LER_MODALITY_DESIGN_2026-07-24.md`:103 |
| Q3 | \| **5** \| **Redeposition granularity** \| Discrete sputter-product landing; contributes micromasking bumps, ~1–3 nm \| Redep flux already computed (`surface_product_redeposition_3d.py`); granularity variance = Poisson on landed-product count per cell. Gated OFF when redep flux → 0 (respects Martin-Cunge "no unconditional roughening") \| Martin & Cunge JVST B 26 (2008); Nakazaki/Ono two-mode roughening \| | `RESEARCH_LER_MODALITY_DESIGN_2026-07-24.md`:158 |

_Harvested 3 quoted + 0 unquoted mentions across the repo's docs._
