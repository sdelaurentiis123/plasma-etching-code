# woo-2024-c4f6-thesis

**Woo 2024 CF4/C4F6/He ICP reactor, patterned-rate, and profile board**

- **Citation:** Byungjun Woo, “A Study on the Etching Characteristics of High
  Aspect Ratio Oxide Etching using C4F6 Plasma in Inductively Coupled Plasma
  with Low Frequency Bias Power,” M.S. thesis, Korea University (2024).
- **DOI:** `10.23186/korea.000000288984.11009.0001569`
- **Repository:** `https://dcollection.korea.ac.kr/srch/srchDetail/000000288984`
- **Archived PDF:** `research_sources/woo_2024_c4f6_thesis.pdf`, SHA-256
  `16bdf4843d0218fe6801ed418fcdd342ef2f825d29a7074d0a1e334b77229523`
- **Extracted full text:**
  `research_sources/thesis_extracts/woo_2024_c4f6_thesis.txt`
- **Status:** PRIMARY FULL TEXT + FIGURE 4.1 PIL-DIGITIZED + PROFILE AND
  DIAGNOSTIC PAGES VISUALLY AUDITED
- **Topic:** C4F6-specific reactor diagnostics, absolute patterned SiO2/ACL
  rates, and feature-profile identifiability.

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | The planar ICP used 10 mTorr, 140 sccm total CF4/C4F6/He flow, 100 sccm He, a 17 C lower electrode, 13.56 MHz source coupling, and 2 MHz bias. The gas-ratio sweep fixed source/bias power at 150/400 W. | Defines a same-reactor knobs/diagnostics/rate board. It is not Krüger's dual-frequency C4F6/Ar/O2 CCP. |
| Q2 | The sample was a 2 x 2 cm center coupon patterned at nominal 100 nm line / 500 nm pitch, with 2400 nm SiO2, 1400 nm ACL, and 50 nm SiON. Cross sections were measured by FE-SEM. | Figure 4.1 is a patterned absolute-rate board, not a blanket-wafer rate. |
| Q3 | Figure 4.1 reports five SiO2 and five ACL rates versus C4F6 fraction. The body gives endpoint SiO2 rates 48.155 and 31.922 nm/min and ACL rates 25.12 and 12.902 nm/min. | The 600 dpi PIL replay produces all ten points and reconciles every printed endpoint within 0.126 nm/min. These rates can grade an independently conditioned reactor/surface pipeline; they may not set its missing fluxes. |
| Q4 | The same sweep measured electron temperature, aggregate ion-current density, self-bias, and relative CFx/F optical emission. Increasing C4F6 lowered Te and patterned rate, increased the plotted ion-current density, increased CFx emission, and lowered F emission. | A reactor model must match several observables together. OES intensity is not an absolute radical flux and aggregate probe current is not an ion-species vector or IEAD. |
| Q5 | Figure 4.6 compares features at 2100-2200 nm only after exposure time was adjusted between gas fractions to obtain similar depth. | The SEM profiles can grade width/bowing shape conditionally. They are not value-blind absolute-depth targets, and target depth may not select petch time. |
| Q6 | The text says ion current rises from 0.0168 to 0.05255 mA/cm2 but labels that a 21% increase; the arithmetic is 212.8%. It also labels the power-sweep composition 56.25% while printing 15/25 sccm CF4/C4F6, which is 62.5%. | Both conflicts are quarantined in the machine-readable audit. Neither the percentage nor the power-sweep mixture may become a model constant until resolved from raw logs or an author correction. |
| Q7 | Source-power sweeps change both ion current and relative reactive-species emission; bias-power sweeps leave Te/current nearly fixed while self-bias changes strongly. | A future reduced reactor provider has a useful qualitative separation test: source power controls plasma production while low-frequency bias mainly changes ion energy. The source does not publish the IEAD needed for a quantitative feature boundary. |

## Numerical and claim gate

The authoritative data receipt is
`data/experimental/woo_2024_c4f6/digitization_manifest.json`; the
classification receipt is
`results/curated/woo_2024_c4f6_board/audit.json`. The Figure 4.1 rates are
allowed as experimental outputs. No Woo value is currently imported as a
surface-law coefficient, neutral flux, species mixture, or Krüger boundary.
