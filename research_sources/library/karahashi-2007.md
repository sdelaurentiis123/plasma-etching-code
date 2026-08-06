# karahashi-2007

**Karahashi, ion-beam SiO2/CFx**

- **DOI/URL:** Hyomen Kagaku 28, 60 (2007)
- **Retrieval route:** fetched
- **Status:** FULL TEXT: research_sources/thesis_extracts/karahashi_2007_sio2_cfx_ionbeam.txt; PDF in research_sources/
- **Topic:** beam-yields — Beam-measured yields, thresholds and sticking (the L0/L1 provenance floor)

## Claims table

| # | source-backed claim | consumed by / boundary |
|---|---|---|
| Q1 | The experiment uses energy-controlled, mass-analyzed single-species ions under ultrahigh vacuum without gas-phase reaction or incident neutral radicals. | `data/experimental/karahashi_2007/`; isolates reactive-ion identity and energy, not a plasma mixture. |
| Q2 | At 1000 eV the text reports about `0.3 SiO2/ion` for F+ and `1.5 SiO2/ion` for CF3+. | Text cross-check for the visually audited Figure-4 digitization. |
| Q3 | The abstract states that yield increases with ion energy and with the number of fluorine atoms in CFx+, gradually saturates above 1000 eV, and drops into fluorocarbon-film growth below 500 eV. | Species and energy ordering; “gradually saturated” is not a hard 1.5 ceiling. |
| Q4 | Figure 4 digitization at 1000 eV gives F+ `0.3232`, CF+ `0.6751`, CF2+ `1.1957`, and CF3+ `1.4703 SiO2/ion`. | Required direct-beam ladder; forbids species-agnostic validation. |
| Q5 | The same CF3+ series reaches `1.8736` at 1500 eV and `1.7549` at 2000 eV. | Retracts the former 1.5 universal-ceiling claim. |
| Q6 | Per-fluorine yield approximately follows the square root of energy allocated to each F atom over the measured series. | Supports the measured-domain energy trend; does not authorize extrapolation above 2000 eV or to larger ions. |
| Q7 | Angular yields rise, peak near 60 degrees, and fall toward grazing incidence; the reported 60/0 ratio depends on species. | Angular-class constraint at the source condition; no energy-wide angular table was digitized. |

## Corrected consumption

The former library rows that described a “4.7% independent validation” and a
`1.5` universal ceiling were invalid. The end-to-end default mechanism
discarded energetic ion identity and returned the same yield for every CFx+
name. The source now feeds only the opt-in
`Karahashi2007ReactiveIonYieldTable`, which reproduces/interpolates the
digitized data inside measured normal-incidence support and refuses
extrapolation. Matching its own table is reproduction evidence, not validation.
