# yoshie-2023-apsusc

**Yoshie et al., Applied Surface Science 638, 157981 (2023)**

- **DOI/URL:** 10.1016/j.apsusc.2023.157981
- **Retrieval route:** open-access publisher article and author-posted full text;
  official Elsevier figure rasters are checksum-pinned in
  `data/experimental/yoshie_2023/digitization_manifest.json`
- **Status:** FULL TEXT ONLINE + FIGURES PIL-AUDITED; 7 blanket and 49 feature
  observations machine-digitized after preregistration
- **Topic:** cyclic C4F8/SF6 silicon depth, transient reactor boundary, ARDE

## Verified claims

| # | claim | consumed by |
|---|---|---|
| Q1 | The experiment used continuous 80 sccm Ar, alternating 15 sccm C4F8 and 15 sccm SF6 one-second pulses, 4 Pa, 400 W ICP power, and a 50 W bias applied for 0.25 s per cycle. | `data/experimental/yoshie_2023/README.md` |
| Q2 | The blanket sample was a 370 nm poly-Si film; the patterned sample used an approximately 500 nm SiO2 mask on a Si substrate with 80--140 nm openings. | `scripts/audit_yoshie_2023_blanket_transfer.py` |
| Q3 | Blanket exposures used 75 cycles. Pattern exposures used 675 cycles for Figure 5 and 450 cycles for Figure 6. The reported rate is depth divided by cumulative bias-on time. | `scripts/digitize_yoshie_2023_figures4_6.py`; `scripts/audit_yoshie_2023_blanket_transfer.py` |
| Q4 | The paper describes FC-film deposition from CFx, Ar+-assisted bottom-film removal, and Si removal by four F atoms with or without Ar+ assistance, yielding volatile SiF4. | future cyclic-Si surface closure |
| Q5 | The plasma discussion identifies Ar+ and SFx+ (x=1--5), F radicals, and CF/CF2 deposition precursors, but reports no species-resolved wafer flux or measured ion energy-angle distribution. | `data/experimental/depth_cross_chemistry_v1/README.md` |
| Q6 | OES shows CF/CF2 peaking during C4F8 introduction and F peaking during SF6 introduction, with F persisting for roughly 500 ms. The article also warns that emission ratios depend on electron temperature. | future reactor-boundary closure |
| Q7 | XPS Table 1 reports different post-cycle C/F/S surface compositions for timings II and III, directly requiring a persistent surface inventory rather than a memoryless rate law. | future cyclic-Si surface closure |
| Q8 | Full-resolution raster replay gives the published Figure-4 maxima (466 and approximately 591 nm per bias-minute) and all Figure-5/6 marker/error-bar locations. | `scripts/digitize_yoshie_2023_figures4_6.py` |

## Transfer warning

The Figure-4 blanket rates are a same-reactor boundary observable, not an
automatic multiplicative feature-depth calibration. Material state and exposure
history both change between blanket and patterned samples. The read-only audit in
`results/yoshie_2023_blanket_transfer/audit.json` quantifies the resulting
non-transfer: the 8 s timing-I feature/blanket rate ratio remains above 2.5 even
after digitization allowances, and the timing rank changes.
