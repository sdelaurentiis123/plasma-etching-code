# christophorou-olthoff-2000-sf6

**NIST-evaluated SF6 electron interactions and swarm constraints**

- **Citation:** L. G. Christophorou and J. K. Olthoff, “Electron Interactions
  with SF6,” *Journal of Physical and Chemical Reference Data* **29**,
  267--330 (2000).
- **DOI:** `10.1063/1.1288407`
- **Official NIST PDF:**
  `https://srd.nist.gov/jpcrdreprint/1.1288407.pdf`
- **PDF SHA256:**
  `b13fac820570646e6b72a59e573fb85ba99814f72cbb0962193b969036a54ad0`
- **Focused extraction:**
  `research_sources/thesis_extracts/christophorou_olthoff_2000_sf6_verified.txt`
- **Executable tables:**
  `data/experimental/christophorou_olthoff_2000_sf6/`
- **Status:** PRIMARY NIST FULL TEXT + TABLES 9, 14--17, 20, 24--25,
  27--28, AND 35--37 VISUALLY AUDITED AT 400 DPI

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| S1 | Tables 9, 14, 15, 17, 20, and 28 separately provide evaluated total scattering, momentum transfer, vibration, total ionization, neutral dissociation, and total attachment curves with different recommendation classes. | The aggregate deck deconvolves low-energy total scattering before re-adding vibration and attachment, preventing collision double counting. Table 20 remains a deduced approximation. |
| S2 | Tables 35--37 provide pure-SF6 effective ionization, drift velocity, and an assessed attachment-rate product versus reduced field. | These are independent swarm gates for the deterministic EEDF solve; they do not specify a reactor state or wafer flux. |
| S3 | Table 16 gives nine partial positive-ion values at 100 eV. Their preferred-analysis sum is `6.51e-20 m2`, while the independently recommended Table-17 total is `6.53e-20 m2`. | The 0.31% mismatch is source rounding/analysis spread. The executable branch split preserves the Table-17 total and the Table-16 fractions. Table 16 is one-energy evidence, not nine measured curves. |
| S4 | Tables 24, 25, and 27 resolve SF6-, SF5-, SF4-, SF3-, SF2-, F2-, and F- attachment products; Table 28 is the sum of those recommended/suggested partials. | The product-resolved deck preserves Table 28 exactly, uses the reported partial shapes, and reconciles only the source-table rounding. Blank Table-27 cells remain unreported rather than measured zero. |
| S5 | The review states that SF5- dominates dissociative attachment below about 1.5 eV and F- above it; beyond about 0.3 eV total attachment is dissociative attachment. | This gives a direct physics check on the resolved branches and materially changes electronegativity/product supply relative to an aggregate `SFx-` sink. |
| S6 | Table 3 reports appearance/formation energies for the fragment positive ions, but the review plots rather than tabulates their energy-dependent partial curves. | Positive-ion energy dependence currently uses a declared threshold/100-eV-anchor closure renormalized to the measured total. Figure 15 digitization or original partial data is the upgrade path; the closure is never described as direct curve measurement. |

## Executable decision

`derive_nist_product_resolved_sf6_replay()` retains the evaluated aggregate
momentum, vibration, and neutral-dissociation inputs, replaces aggregate
attachment by seven source-resolved products, and replaces aggregate
ionization by nine explicitly graded product branches. Every energy-grid knot
reconstructs the evaluated attachment and ionization totals to floating-point
tolerance. Double-ion branches create two electrons. This closes primary SF6
charged-product bookkeeping, but heavy-particle reactions, walls, open-flow
balance, absorbed-power coupling, sheath transport, and TiO2 surface yields
remain required before a Zhu depth can be predicted.
