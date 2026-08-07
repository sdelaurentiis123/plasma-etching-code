# mahorowala-1998-thesis

**Mahorowala, fixed-time Cl2/HBr poly-Si feature corpus and simulator**

- **Citation:** A. P. Mahorowala, *Feature profile evolution during the high
  density plasma etching of polysilicon*, PhD thesis, MIT Department of
  Chemical Engineering (1998).
- **Handle:** `http://hdl.handle.net/1721.1/50514`
- **MIT DSpace item:** `9af24424-813c-401d-9ca1-33cb863d86c3`
- **Bitstream:** `42415621-MIT.pdf`,
  `c04d2da3-cbc2-4898-97ed-c83c58293e05`
- **PDF SHA-256:**
  `e561dfe9780c0e27439b2b7288788c18333b8228152b6cc40baf72ba2edf4b6a`
- **Status:** PRIMARY FULL THESIS READ + TABLE 2.2/FIGURE 2.4 PIL-AUDITED
- **Rights:** MIT thesis viewed locally; source pixels are not redistributed.
- **Topic:** Cl2/HBr poly-Si absolute fixed-time rates, designed feature
  profiles, mask effects, deposition, microtrenching, and charging

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | The experiment used a Lam TCP 9400SE ICP with a separately RF-biased wafer electrode at 60 °C, 10 mTorr, pure Cl2 or pure HBr, and optical endpoint monitoring. | Establishes a documented commercial-reactor feature corpus. Recipe knobs do not by themselves identify the wafer species/IEAD boundary. |
| Q2 | The structures used 500 nm poly-Si on gate oxide, either about 475 nm photoresist after ARC etch or a 200 nm oxide mask, 250 nm lines, and spacings from 250 nm to 2.5 µm. Figure 2.4 is the oxide-mask 310 nm-spacing montage. | Supplies physical geometry for direct profile replay and multiple same-condition transport widths. |
| Q3 | Table 2.2 holds pressure and etch time fixed at 10 mTorr and 75 s while varying inductive power (250–550 W), RF bias power (20–140 W), and Cl2 flow (25–175 sccm). Eleven runs report rates; runs 8 and 12 are marked overetched. | A designed 11-point absolute fixed-time rate/depth board with corresponding SEM profiles; no target-time adjustment. |
| Q4 | The reported poly-Si rates span 900–3675 Å/min, equivalent to 112.5–459.375 nm removed in 75 s. Oxide-mask rates and selectivities are reported separately. | Direct dimensional rate/time evidence. The conversion is algebraic and contains no depth fit. |
| Q5 | Photoresist-masked Cl2 samples were deliberately timed to approximately equal depth, whereas oxide-masked Table 2.2 samples were not; the thesis states that their etch depth therefore varied. | Only Table 2.2 is promoted as a fixed-time absolute-depth board. The equal-depth photoresist montage cannot be used as a blind depth test. |
| Q6 | At the center condition the thesis estimates about 2 mA/cm2 ion current, 100 eV sheath gain (about 120 eV including plasma potential), a neutral/ion ratio near 100, and a 10° ion angular FWHM; the profile simulator actually used 35 eV and neutral/ion ratio 10 after declaring profile insensitivity in its chosen regime. | These are model estimates/reductions, not measured per-run boundary data. They cannot be relabeled as species-resolved IEAD/IAD measurements or used to grant first-principles depth authority. |
| Q7 | The experimental profiles show mask-dependent deposition, bowing, faceting, and microtrenching; oxide masking minimizes product deposition/redeposition and preserves its initial mask profile better than photoresist. | Causal morphology targets for reflection, product return, mask erosion, and sidewall-state ablations. |

## Vision audit and use

Table 2.2 was inspected on a checksum-bound 600 dpi render of PDF/print page
39. The full table grid and every retained value were reconciled against the
source. Figure 2.4 on PDF/print page 42 was inspected at 600 dpi and contains
all thirteen run-numbered SEM panels. The numerical table and reproducible
PIL/checksum audit live in
`data/experimental/mahorowala_1998_cl2/`; source pixels do not.

This source repairs the old literature statement that the chlorine profile
lineage lacked fixed exposure times: Levinson Figure 11 still does, but
Mahorowala Table 2.2 does not. It does **not** solve the reactor boundary.
The clean validation architecture is now beam-derived Chang/Sawin surface
kinetics plus an independently validated Cl2 reactor provider, followed by
the untouched 75 s rate/profile matrix. A facility-conditioned mode may use
independent same-reactor rate information to set a boundary before scoring
profile shapes, but the scored Table-2.2 rate itself may not select that
boundary and then count as a prediction.
