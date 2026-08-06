# zhang-thesis

**Zhang (Yiting), PhD thesis, Univ. Michigan (2015)**

- **DOI/URL:** PhD thesis, Univ. Michigan (2015)
- **Retrieval route:** cpseg theses index
- **Status:** FULL TEXT: research_sources/thesis_extracts/zhang_yiting_phd_thesis.txt
- **Topic:** modeling-sota — Modeling state of the art and competitor codes

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | - **Nishita–Nakamae (SIGGRAPH '85; and TOG "Shading models for point and linear sources")**: first exact penumbra/umbra decomposition for polygonal area sources — shadow volumes per (source, blocker) pair, penumbra = convex-hull volumes, umbra = intersection volumes; illumination at a point computed by clipping the source polygon against blocker projections and applying the exact contour-integral (Lambert) formula. This is the per-point EXACT visible-polygon primitive we need, published 1985. [dl.acm.org/doi/10.1145/282918.282938](https://dl.acm.org/doi/10.1145/282918.282938), [history.siggraph.org](https://history.siggraph.org/learning/the-visibility-skeleton-a-powerful-and-efficient-multi-purpose-global-visibility-tool-by-durand-drettakis-and-puech/) - **Heckbert 1992 (3rd Eurographics Workshop on Rendering)**: *incomplete* discontinuity meshing — builds only the EV wedges (source edge | `RESEARCH_EXACT_3D_OCCLUSION_2026-07-21.md`:43 |
| Q2 | **Has our exact failure been diagnosed before? Partially, twice** — see §4. The closest match is **Zhang (Michigan, 2015)**: MCFPM-3d reproduced necking + bowing *phenomenology* but **"does not precisely reproduce the positions of the necking and bowing effect"**, and he attributes it to the **ion angular distribution**, not to chemistry: *"a slight change in angular distribution may contribute significantly to different shape evolutions."* No source I found reports a simulated mouth that seals when the experiment stays open, so our specific magnitude failure appears to be **undiagnosed in the open literature** — which is consistent with it being an implementation-side angular-delivery defect rather than a missing chemical channel. | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:57 |
| Q3 | Verbatim (`zhang_yiting_phd_thesis.txt`, §7.3 Model Validation, He/Cl2 Si trench, **Lam Research ICP**, line/pitch 50/100 nm, 60 nm oxide + 60 nm nitride mask, ∆x=∆y=∆z=1.25 nm): > "The masks show erosion with increasing etch time, which can be seen in the measured SEMs … > With the thickness of mask continuing to decrease, ions with large horizontal velocities will > bombard the sidewalls of the feature, and thus causes sidewall etching. After ions strike on the > surface, there will be high energy neutrals reflecting back to the plasma and bombarding the > surface again. **This high energy particle reflection brings about the necking and bowing effect** > as observed in the third column of Fig. 7.2. **There is a difference of necking and bowing > positions between the experimental measurements and the predicted simulation results.** This is > mainly due to the absence of reactor scale  | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:150 |
| U4 | [unquoted — verify on next use] \| **Zhang, Yiting — *Low Temperature Plasma Etching Control through Ion Energy Angular Distribution and 3-Dimensional Profile Simulation*, PhD, U. Michigan, 2015** \| `.../tmp/pdfs/zhang_yiting_phd_thesis.pdf` (+`.txt`) \| **HIGH.** MCFPM-3d vs Lam-ICP SEM time series; explicit admission of necking/bowing **position** mismatch. \| | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:81 |
| U5 | [unquoted — verify on next use] ### 1.3 Zhang Yiting thesis (Michigan 2015) — the only found statement of a necking/bowing model-vs-experiment MISS | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:144 |
| U6 | [unquoted — verify on next use] \| `zhang_yiting_phd_thesis.pdf` (+`.txt`) \| **Zhang (Michigan 2015) — MCFPM-3d vs SEM; published necking/bowing position miss.** \| 15,915,085 \| | `RESEARCH_MOUTH_LITERATURE_BROAD_2026-08-02.md`:785 |

_Harvested 3 quoted + 3 unquoted mentions across the repo's docs._
