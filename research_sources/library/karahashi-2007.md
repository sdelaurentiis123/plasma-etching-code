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
| Q7 | Angular yields rise, peak near 60 degrees, and fall toward grazing incidence; the reported 60/0 ratio depends on species. | `figure6_angular_yields_unknown_energy.csv`; condition-unknown angular-class constraint, not an energy-resolved production law. |
| Q8 | Figure 10 resolves CF3+-driven desorbed SiF, SiF2, and SiF4 fractions at 500, 1000, and 2000 eV. SiF2 dominates at or below 1000 eV; SiF increases and SiF4 decreases with energy. The Figure-10 incidence angle is not reported, while the same source says product composition changes with angle. | `figure10_cf3_product_fractions.csv`; product-branching validation only. It cannot be condition-matched to the explicitly normal-incidence Figure-4 yield board. |
| Q9 | The text assigns SiF primarily to prompt collision-cascade ejection, whereas SiF2 and SiF4 include thermally activated desorption of collision-generated precursors; the SiF4 delay is about 0.5 ms. | Requires the event architecture to distinguish prompt removal from delayed product transport. No diffusion coefficient or escape length is reported. |

## Corrected consumption

The former library rows that described a “4.7% independent validation” and a
`1.5` universal ceiling were invalid. The end-to-end default mechanism
discarded energetic ion identity and returned the same yield for every CFx+
name. The source now feeds only the opt-in
`Karahashi2007ReactiveIonYieldTable`, which reproduces/interpolates the
digitized data inside measured normal-incidence support and refuses
extrapolation. Matching its own table is reproduction evidence, not validation.
Figure 6 is now independently digitized and visually audited, but its energy
remains unreported; Arts et al. (2021) explicitly label it unknown. Its normal
points strongly identify 1000 eV by cross-figure consistency, recorded only as
an inference.

Figure 10 now adds a checksum-bound product-identity board. Its ordinate is
`SiFx / sum(SiFx)`, so it cannot be multiplied into a total removal rate as if
it were an escape probability. It instead grades whether a condition-matched
event kernel routes removed silicon into the measured product identities while
a separate conserved transport state handles delayed escape or reincorporation.
No such condition-matched total-yield join is currently available: Figure 4 is
normal incidence, but the Figure-10 incidence angle is unreported.
