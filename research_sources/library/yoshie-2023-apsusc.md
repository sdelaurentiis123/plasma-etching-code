# yoshie-2023-apsusc

**Yoshie et al., Applied Surface Science 638, 157981 (2023)**

- **DOI/URL:** 10.1016/j.apsusc.2023.157981
- **Retrieval route:** open-access publisher article and author-posted full text;
  official Elsevier figure rasters are checksum-pinned in the depth and
  reactor-state manifests under `data/experimental/yoshie_2023/`
- **Status:** FULL TEXT ONLINE + FIGURES PIL-AUDITED; 7 blanket, 49 feature,
  7 electron-density-window, 60 OES observations digitized, and 19 published
  XPS Table-1 rows transcribed
- **Topic:** cyclic C4F8/SF6 silicon depth, transient reactor boundary, ARDE

## Verified claims

| # | claim | consumed by |
|---|---|---|
| Q1 | The experiment used continuous 80 sccm Ar, alternating 15 sccm C4F8 and 15 sccm SF6 one-second pulses, 4 Pa, 400 W ICP power, and a 50 W bias applied for 0.25 s per cycle. | `data/experimental/yoshie_2023/README.md` |
| Q2 | The blanket sample was a 370 nm poly-Si film; the patterned sample used an approximately 500 nm SiO2 mask on a Si substrate with 80--140 nm openings. | `scripts/audit_yoshie_2023_blanket_transfer.py` |
| Q3 | Blanket exposures used 75 cycles. Pattern exposures used 675 cycles for Figure 5 and 450 cycles for Figure 6. The reported rate is depth divided by cumulative bias-on time. | `scripts/digitize_yoshie_2023_figures4_6.py`; `scripts/audit_yoshie_2023_blanket_transfer.py` |
| Q4 | The paper describes FC-film deposition from CFx, Ar+-assisted bottom-film removal, and Si removal by four F atoms with or without Ar+ assistance, yielding volatile SiF4. | future cyclic-Si surface closure |
| Q5 | The plasma discussion identifies Ar+ and SFx+ (x=1--5), F radicals, and CF/CF2 deposition precursors, but reports no species-resolved wafer flux or measured ion energy-angle distribution. | `data/experimental/depth_cross_chemistry_v1/README.md` |
| Q6 | OES shows CF/CF2 peaking during C4F8 introduction and F peaking during SF6 introduction, with F persisting after the pulse. The article also warns that emission ratios depend on electron temperature. | `scripts/digitize_yoshie_2023_reactor_state.py`; future reactor-boundary closure |
| Q7 | XPS Table 1 reports timing-II elemental S/C/F = 32.2/39.5/28.3% and timing-III = 18.5/45.8/35.7%, with separately normalized S 2p, C 1s, F 1s, and Si 2p bonding-state partitions. These are direct surface-memory constraints, not layer thicknesses or absolute site inventories. | `data/experimental/yoshie_2023/table1_xps_surface_composition.csv`; future cyclic-Si surface closure |
| Q8 | Full-resolution raster replay gives the published Figure-4 maxima (466 and approximately 591 nm per bias-minute) and all Figure-5/6 marker/error-bar locations. | `scripts/digitize_yoshie_2023_figures4_6.py` |
| Q9 | Figure 12 measures a strongly phase-dependent bulk electron-density trace. Full-resolution replay retains three samples over each bias window; the 8 s timing-II window crosses the fast post-SF6 density collapse. Electron density alone is not a measured positive-ion flux in this changing electronegative plasma. | `data/experimental/yoshie_2023/figure12_bias_window_electron_density.csv` |
| Q10 | Figure 14 reports CF, CF2, and F optical-emission intensities normalized to Ar, not absolute ground-state densities. The digitization therefore preserves the ratios and background-subtracted negative points without converting them to flux. | `data/experimental/yoshie_2023/figure14_phase_resolved_oes.csv` |

## Transfer warning

The Figure-4 blanket rates are a same-reactor boundary observable, not an
automatic multiplicative feature-depth calibration. Material state and exposure
history both change between blanket and patterned samples. The read-only audit in
`results/yoshie_2023_blanket_transfer/audit.json` quantifies the resulting
non-transfer: the 8 s timing-I feature/blanket rate ratio remains above 2.5 even
after digitization allowances, and the timing rank changes.
