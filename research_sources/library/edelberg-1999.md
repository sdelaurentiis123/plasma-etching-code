# edelberg-1999

**Edelberg & Aydil (1999)**

- **DOI/URL:** Edelberg & Aydil (1999)
- **Retrieval route:** publisher
- **Status:** not-fetched
- **Topic:** reactor-sheath — Reactor-scale and sheath closure models

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | - **Lam Research (Kim, Hudson, Cooperberg, Edelberg, Srinivasan, *Thin Solid Films* 515, 4874 (2007)):** necking ← net polymer deposition rate on the sidewall; bowing ← ion scattering off the **secondary** facet; and — the load-bearing negative result for us — **the *primary* photoresist facet angle "showed only a small influence on the SiO2 etch profile."** - **Lam Research 2023 (Shen, Lill *et al.*, JJAP 62, SI0801):** "Polymer deposition on the sidewalls of the mask, so-called **necking**, is considered the main root cause of sidewall roughness." | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:32 |
| Q2 | - **Benoit-Cattin & Bernard, J. Appl. Phys. 39, 5723 (1968)** — the collisionless bimodal IED: constant sheath width, uniform field, sinusoidal sheath voltage, ions entering at Bohm speed give the classic **two-horn (bimodal) distribution** whose splitting ΔE ∝ (1/ωτ_ion), narrowing with frequency and ion mass. Its high-frequency limit is the **arcsine distribution** (density ∝ 1/√(1−((E−Ē)/ΔE)²), horns at ±ΔE/2) — **this is exactly the form petch already implements** in `chemistry.py::_ied_yield` (`'bimodal'` branch, cited to Kawamura 1999). - **Kawamura, Vahedi, Lieberman, Birdsall — "Ion energy distributions in rf sheaths: review, analysis and simulation" — Plasma Sources Sci. Technol. 8, 313 (1999)**, DOI [10.1088/0963-0252/8/3/202](https://iopscience.iop.org/article/10.1088/0963-0252/8/3/202) — the definitive review spanning collisionless→collisional and low→high ωτ_ion; the analyti | `RESEARCH_REACTOR_TIER1_DESIGN_2026-07-24.md`:145 |
| U3 | [unquoted — verify on next use] **Primary claim: Kim, Hudson, Cooperberg, Edelberg, Srinivasan (Lam Research), *Thin Solid Films* 515(11), 4874–4878 (2007), DOI 10.1016/j.tsf.2006.10.023.** Full text **not obtained** (Elsevier 403). The claim is preserved in two independent, *fetched* secondary sources: | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:525 |
| U4 | [unquoted — verify on next use] \| **Kim, Hudson, Cooperberg, Edelberg, Srinivasan, *Thin Solid Films* 515, 4874 (2007)**, DOI 10.1016/j.tsf.2006.10.023 \| Elsevier 403 \| **The** mouth-mechanism paper. Its depositor-flux sweep is the exact experiment we want to reproduce. **Highest-value retrieval.** \| | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:807 |
| U5 | [unquoted — verify on next use] ### Gate 2 — independent global-model measured densities - Reproduce an **independent** paper's *measured* densities/T_e, not HPEM output. Best targets: Gudmundsson CF₄ (2019) or O₂/Ar (2007) density-vs-pressure curves (the paper validates its own set), and **Efremov** Ar/CF₄/O₂ / CHF₃ Langmuir-probe densities (2021) — these give F, CF₂, CF₃, ion-density vs mixing ratio. SF₆/O₂ arm: the ICP SF₆/O₂/Ar 2014 model's density-vs-%O₂. Score density trends (E2/E3) within factor-of-2 absolute, correct monotonic trends. - Also fold in a **measured-IEDF** check against Sobolewski/Edelberg-Aydil CF₄ IEDs (10 mTorr) to grade the sheath module against experiment, not just HPEM. | `RESEARCH_REACTOR_TIER1_DESIGN_2026-07-24.md`:361 |
| U6 | [unquoted — verify on next use] \| P10 \| **Reactor Tier-1 module B** (sheath closure vs Miller-Riley/Edelberg), then module A + gas-temperature balance (T_gas sets core width: 300↔1000 K = 16× AR-200 sidewall flux) \| staged \| design doc + corrections \| IADF §B \| | `ROADMAP_PERFECTION_2026-07-29.md`:37 |

_Harvested 2 quoted + 4 unquoted mentions across the repo's docs._
