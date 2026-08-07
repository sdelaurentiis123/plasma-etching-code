# an-2026-apsusc

**An et al., Applied Surface Science 716, 164601 (2026)**

- **DOI/URL:** 10.1016/j.apsusc.2025.164601; arXiv:2508.00327v1
- **Retrieval route:** open arXiv preprint + author data/code repository
- **Status:** FULL TEXT READ; VERIFIED EXCERPTS + RELEASED FIGURE DATA PINNED
- **Topic:** atomistic-surface — DFT-trained reactive-ion event physics

## Claims table

| # | source-backed claim | consumed by / boundary |
|---|---|---|
| Q1 | Separate SiO2 and Si3N4 neural potentials are trained on DFT energy/force/stress labels across bulk, slab, gas, mixed-layer, and fluorocarbon-film environments; ZBL supplies the short-range collision repulsion. | Atomistic evidence candidate; stronger than a yield-regressed surface coefficient, but still an approximate PBE-trained interatomic potential. |
| Q2 | Figure 3 reports normal-incidence Si-removal yields for CF+, CF2+, CF3+, and CHF2+ on SiO2 up to 1000 eV, plus CF2+/CHF2+ on Si3N4. | `data/surface_interactions/an_2026_nnp/`; no-fit transfer audit against petch's separately digitized Karahashi board. |
| Q3 | Experimental etch yields are comparison observables, not NNP loss targets; the released work was nevertheless developed and published with those comparisons, so petch does not label the result blind. | Calibration firewall in `AN-NNP-KARAHASHI-TRANSFER-R1`. |
| Q4 | The MD protocol deletes designated stable products after each impact and allows them to escape regardless of formation depth. | Explicit nonfundamental transport/escape operation; prevents treating atomistic resolution as complete physical accuracy. |
| Q5 | For low-energy CF+ on SiO2, the model delays deposition and etches more than 10 nm before film growth versus about 2 nm experimentally; an arbitrary 2 nm product-removal-depth rule improves agreement. | Direct evidence that finite product escape/diffusion and evolving film state are required closures. |
| Q6 | The present calculations omit incident neutrals and radicals and use normal-incidence ions. | Does not identify the Krüger reactor boundary or an angular response law. |
| Q7 | Code, potential, and figure data are released in the author repository, but no license file was found at pinned commit `4bcd0350`. | Numerical data may be audited with provenance; do not import implementation or weights without permission/license. |
| Q8 | At pinned commit `4bcd0350`, `SM_codes/PlasmaEtchSimulator/calc/byproduct.py` deletes every disconnected cluster matching the byproduct list without testing formation depth; its optional depth-gated diagnostic hard-codes `20.0` Å and labels it an “arbitrary value.” | Code-level confirmation that the reported low-energy improvement is a sensitivity study, not a measured product-transport law. Source SHA-256 `a0471f888bd885a84b7188a4111aa18a4c7c6a68d8061e0b7d5e6193204babd0`. |
| Q9 | The released SiO2 transient archive and the CF 300 eV default/2 nm histories are checksum-pinned in `AN-NNP-KARAHASHI-TRANSFER-R1`. | Reproducibility receipt only; the arbitrary 2 nm gate is not consumed as a production parameter. |

## Corrected consumption

This source does not become a table lookup in the production surface model.
Its high-energy SiO2 transfer is a strong candidate for the prompt
species/energy event kernel. Its low-energy failure identifies the necessary
coupling to a finite mixed layer, fluorocarbon-film evolution, and
depth-/time-dependent product escape. The result is therefore a physics
architecture constraint, not a new depth-fitting knob.

The code audit also rules out a superficially attractive shortcut: copying
the 2 nm diagnostic threshold into petch would merely reproduce an arbitrary
author sensitivity. A production escape closure must instead conserve the
formed product until an independently constrained transport/reaction process
either releases it or reincorporates it.
