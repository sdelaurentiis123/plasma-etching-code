# cagomoc-2023-thesis

**Charisse Marie Donato Cagomoc, Osaka University doctoral dissertation
(2023), “Molecular Dynamics Simulation of SiO2 and SiN Etching for 3D NAND
Memory Device Applications”**

- **DOI:** 10.18910/91922
- **Handle:** https://hdl.handle.net/11094/91922
- **Retrieval route:** official Osaka University OUKA dissertation PDF
- **Status:** FULL TEXT + PIL-AUDITED FIGURES; PDF not redistributed
- **PDF SHA-256:** `f1eb74b1b42bc12cc89fe60426b327c873fd75f68c2012c2f04e3c54adcade9f`
- **Extracted text:** `research_sources/thesis_extracts/cagomoc_2023_dissertation.txt`
- **Extracted-text SHA-256:** `904fc9293090f6cf666e2d9ce5723a74473aaf6a7c730ff842738052327fcc01`
- **Topic:** atomistic surface physics; fluorocarbon/SiO2 radical-ion synergy;
  nanohole product escape and redeposition

## Claims table

| # | source-backed claim | consumed by |
|---|---|---|
| Q1 | Figure 5.10 reports steady-state Si removal from flat SiO2 for normal-incidence 2000 eV CF3+ with CF3 radical/ion ratios 0, 25, 50, 100, 200, and 300. The digitized yields are approximately 1.81, 3.94, 3.82, 2.27, 0, and 0 Si atoms per ion. | `data/surface_interactions/cagomoc_2023/` |
| Q2 | The source states that 25:1 nearly doubles the ion-only yield, 50:1 is similar, 100:1 decreases, and excessive fluorocarbon accumulation at 200:1 and 300:1 causes etch stop and film deposition. | surface-law topology constraint |
| Q3 | The SiO2/SiN simulations show a C/F/substrate mixing layer typically a few nanometers thick; its thickness increases with ion energy while retained C and F can decrease. | finite mixed-layer architecture |
| Q4 | In nanoholes, the source attributes lower nominal Si yield than flat-wafer yield to redeposition on sidewalls and ion interception by tapered surfaces; it also states that slowly traveling products may be removed prematurely by its finite MD cycle. | product-escape/redeposition claim boundary |
| Q5 | For the radical arm, 0.5 eV normal CF3 is a computational surrogate for roughly 0.026 eV thermal radicals. Appendix D obtains similar approximately 2 nm, F/C=3 films after equal dose, but this does not validate angular delivery into a feature. | radical-transfer limitation |
| Q6 | The source does not simulate simultaneous radical-ion nanohole etching and explicitly defers it; the Figure 5.10 curve therefore cannot identify Krueger's feature boundary or absolute depth. | no-fit/no-transfer boundary |

## Visual audit

Figure 5.10 (PDF page 121) was rendered at 400 dpi. Original-resolution vision
inspection and independent Pillow RGB/dimension checks were followed by
dark-pixel axis localization and saturated-blue marker-core localization. The
render SHA-256 is
`e8888593f82d107a15b329e91acf20af1bb60e59ac3fb39aebc4d11f52adacf0`.
Figures 5.7--5.9 were also visually inspected and their hashes are frozen in
`data/surface_interactions/cagomoc_2023/source_manifest.json`.

## Evidence ceiling

This is classical molecular dynamics with an empirical interatomic potential.
It is valuable evidence for response topology and atom/product bookkeeping,
not a direct measurement, reactor boundary, DFT-trained event kernel, or
species transfer from CF3+/CF3 at 2000 eV to Krueger's unresolved
C4F6/Ar/O2 ion population.
