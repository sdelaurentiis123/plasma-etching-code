# casey-2021-pt-foundations

**Physical definitions for pulsed-Townsend transport properties**

- **Citation:** M. J. E. Casey et al., "Foundations and interpretations of
  the pulsed-Townsend experiment," *Plasma Sources Science and Technology*
  **30**, 035017 (2021).
- **DOI:** [10.1088/1361-6595/abe729](https://doi.org/10.1088/1361-6595/abe729)
- **Retrieval:** [Ulster institutional repository publisher PDF](https://pure.ulster.ac.uk/ws/files/91167940/Casey_et_al_2021_PSST_Foundations_and_interpretations_of_the_pulsed_Townsend_experiment_30_035017.pdf)
- **Local source PDF SHA-256:**
  `21acda40d9bba7fad36a97b1ba614010327846ea225914b8a55eac46249e74c1`
- **Extract:**
  `research_sources/thesis_extracts/casey_2021_pt_verified_excerpt.txt`
- **Status:** PRIMARY CC-BY PUBLISHER FULL TEXT READ; EQUATIONS 24 AND 30-32
  VERIFIED
- **Topic:** flux/bulk/legacy pulsed-Townsend transport definitions

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | Flux and bulk transport coefficients are distinct when electron-number-changing collisions are active. | Any swarm board must name the measurement-equivalent coefficient; generic "drift" is insufficient in attaching/ionizing gases. |
| Q2 | The legacy Brambring PT representation is not the physical continuity equation under nonconservative collisions. Its fitted `W_tilde`, `alpha_tilde_T`, and `D_tilde_L` are not universal kinetic coefficients. | Legacy PT values cannot directly grade flux or modern bulk solvers. |
| Q3 | Equations 30-32 give `R_net = alpha_tilde_T W_tilde`, `W_B = W_tilde + alpha_tilde_T D_tilde_L`, and `D_B,L = D_tilde_L`. | Transform a legacy PT result only with all same-study co-observables, or re-fit the original current transient using the physical equation. |
| Q4 | The current amplitude can provide a flux-drift measurement if the initial electron number is measured accurately. | This special calibrated-current route does not reclassify a transit-time drift value as flux drift. |

## Executable decision

Lan-Jeon reports `Wv` without the same-study Townsend and longitudinal-
diffusion values needed for the Casey transformation. Therefore neither petch
flux drift nor BOLSIG density-gradient bulk drift is a measurement-equivalent
target for that curve. Retuning the collision set to force agreement is
rejected. This decision has no reactor-state, wafer-flux, or depth authority.
