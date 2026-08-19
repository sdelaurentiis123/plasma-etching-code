# depla-2024-oxide-sputter

**Depla and Van Bever, semi-empirical oxide sputter-yield calculation**

- **Citation:** D. Depla and J. Van Bever, “Calculation of oxide sputter
  yields,” *Vacuum* **222**, 112994 (2024).
- **DOI:** `10.1016/j.vacuum.2024.112994`
- **Primary record:**
  `https://biblio.ugent.be/publication/01HN0EVQ4K9X7ZKEMNFYB2R2RP`
- **Status:** PRIMARY AUTHOR-REPOSITORY FULL TEXT; FIGURES 1/2
  CHECKSUM-PINNED, VISUALLY AUDITED, AND TiO2 FIT CURVE DIGITIZED
- **Receipt:**
  `data/experimental/depla_2024_tio2_sputter/digitization_manifest.json`

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | The paper compiles oxide sputter yields as total target atoms emitted per incident Ar ion and converts molecular yields using stoichiometry. | The TiO2 formula-unit scale is therefore the plotted atom yield divided by three; silently treating the plotted axis as TiO2 units would overstate removal by 3x. |
| Q2 | Its TiO2 curve is a Seah–Nunney semi-empirical fit with oxygen surface binding energy as the fitted free parameter; Figure 2 gives `3.52 ± 1.14 eV` at 95% confidence. | This supplies a bare-oxide physical-sputter reference and prior structure, not the reactive fluorocarbon law. |
| Q3 | The low-energy Crist TiO2 input spans `0.2–3 keV`, but the source reports no ion current density; Depla and Van Bever infer current density from metal reference yields. | The resulting low-energy yield and fitted curve are not a direct current-normalized TiO2 beam measurement. |
| Q4 | The digitized fitted curve gives approximately `0.128`, `0.192`, and `0.279` stoichiometric TiO2 formula units per Ar ion at `200`, `276`, and `400 eV`. | This is far below the Oxford effective-removal requirement of roughly `1.34–1.78` formula units per delivered ion, showing that reactive ion-assisted chemistry must carry most removal if the reactor dose bracket is right. It is not a ceiling or a transferable coefficient. |

## Use in petch

The curve is retained as a model-form and unit-conversion constraint. It may
seed a broad sensitivity prior for the bare physical-sputter channel, but it
cannot certify ALD TiO2 under CHF3/SF6/O2, mixed molecular ions, a fluorinated
surface, or the Oxford NPG80. The source PDF and rendered pixels are not
redistributed; their SHA256 receipts and the replayable curve are committed.
