# sandia-2001-fluorocarbon-mechanisms

**C2F6 and CHF3 plasma-reactor and surface mechanisms**

- **Citation:** P. Ho, J. E. Johannes, R. J. Buss, and E. Meeks,
  *Modeling the Plasma Chemistry of C2F6 and CHF3 Etching of Silicon Dioxide*,
  Sandia National Laboratories report SAND2001-1292 (2001).
- **DOI:** `10.2172/782704`
- **Official full text:** `https://www.osti.gov/servlets/purl/782704`
- **PDF SHA256:**
  `3e401d0d5c5ffb0308767bffb2d4d952b1c6110f79f95c03b781e927342a4d8f`
- **Local focused extraction:**
  `research_sources/thesis_extracts/sandia_2001_chf3_table9.txt`
- **Status:** PRIMARY OSTI FULL REPORT + TABLE 9 TRANSCRIBED AND
  CONSERVATION-CHECKED

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| S1 | The CHF3 gas mechanism adds the 38 reactions in Table 9 to 122 non-C2F6 reactions from Table 2. | Table 9 is a hydrogen-bearing extension, not the full standalone fluorocarbon network. |
| S2 | The source rate law is `k = A T^B exp(-C/T)`, with temperature in kelvin, `C` in kelvin, and electron energy loss in eV. | The executable deck converts the source law exactly to the solver's electron-temperature-in-eV convention and tests equality against the source expression. |
| S3 | Rows 1--19 are mostly fits to the Kushner/Zhang CHF3 cross-section set; rows 20--25 copy CHF3 rate shapes to CHF2; the remaining ion-ion and dissociative-recombination rows are estimates. | Petch retains `24 regressed` and `14 estimated` labels. Conservation passing does not upgrade their empirical status. |
| S4 | The silicon and oxide surface additions include adsorption/abstraction coefficients chosen as reasonable values or adjusted for good experimental fits, and mass-weighted ion yields extrapolated from the C2F6 mechanism. | Those surface parameters are not imported into the TiO2 target model; the report is gas-topology evidence here. |
| S5 | The report compares its mechanisms with diagnostics in specific Sandia/GEC reactor configurations and discusses machine-specific power-coupling choices. | Neither its coupling efficiency nor its surface fits transfer to the Oxford NPG80 without an independent equipment closure. |

## Executable decision

All 38 Table-9 rows are installed in
`petch.reactor_global.sandia_chf3_mechanism`.  Atoms, charge, electron-energy
ledger coverage, source-law conversion, and evidence counts are unit-tested.
The deck is not yet a complete target-reactor solve: Table-2 CFx/O/F chemistry,
SF6 daughters, transport, wall loss, exhaust, and power balance must be coupled
before it can produce wafer fluxes.
