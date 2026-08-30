# miao-2016-cryo-grating

**H. Miao, L. Chen, M. Mirzaeimoghri, R. Kasica, and H. Wen,
"Cryogenic Etching of High Aspect Ratio 400 nm Pitch Silicon Gratings," 2016**

- **DOI:** 10.1109/JMEMS.2016.2593339
- **Primary full text:** https://pmc.ncbi.nlm.nih.gov/articles/PMC5084849/
- **NIST record:** https://www.nist.gov/publications/cryogenic-etching-high-aspect-ratio-400-nm-pitch-silicon-gratings
- **Status:** PRIMARY OPEN FULL TEXT READ ONLINE; NUMERICAL RECIPE, FIGURE
  CAPTIONS, AND REPORTED OUTCOMES VERIFIED
- **Topic:** cryogenic SF6/O2 silicon etching, mask selectivity, undercut,
  high-aspect-ratio pattern transfer

## Claims table

| # | verified claim | use boundary |
|---|---|---|
| P1 | The optimized two-minute condition was -110 +/- 2 C, 8 mTorr, 1000 W ICP, 10 W RF, 52 sccm SF6, and 8 sccm O2 on an Oxford Plasmalab 100 ICP tool. | Direct recipe evidence for this tool and wafer family; powers do not imply the same wafer flux or ion energy on another chamber. |
| P2 | Polymer, Cr, SiO2, and Cr-on-polymer masks produced trench depths of 3.9, 3.8, 3.7, and 3.5 um, with half-depth widths near 185, 200, 190, and 205 nm. | Strong process-window anchor; not an absolute-rate transfer to another mask polymer, pattern loading, or tool. |
| P3 | Reported Si-to-mask selectivity was about 15 for polymer, greater than 500 for Cr, and about 60 for SiO2. Cr-on-polymer retained Cr-like selectivity until Cr loss. | Supports mask-budget feasibility and the Cr-on-polymer option. It does not identify an unknown printed polymer's selectivity. |
| P4 | At fixed 52 sccm SF6, oxygen below 8 sccm sharply increased undercut and reduced the neck; oxygen above 8 sccm reduced undercut but also reduced depth and increased passivation. | Fixes the sign of the first O2 derivative board around the selected center. |
| P5 | At 8 mTorr, the 1000 W ICP / 10 W RF / 52:8 recipe reached about 10.6 um after ten minutes with Cr-on-polymer; higher ICP or RF variants etched faster but became undercut-limited. | Selects the low-RF branch when fidelity outranks throughput. |
| P6 | The authors state that feature size materially changes the usable cryogenic recipe. | Requires local calibration/uncertainty propagation for a different layout even when it is geometrically easier. |

## Consumed decision

The recipe is the primary evidence anchor for an initial fidelity-first silicon
transfer board. It is not used as a universal surface coefficient or as proof
that a different etcher will reproduce the reported absolute depth.
