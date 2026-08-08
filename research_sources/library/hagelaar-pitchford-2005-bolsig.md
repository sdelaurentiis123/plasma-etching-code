# Hagelaar & Pitchford (2005) — deterministic electron Boltzmann solver

- **Citation:** G. J. M. Hagelaar and L. C. Pitchford, “Solving the
  Boltzmann equation to obtain electron transport coefficients and rate
  coefficients for fluid models,” *Plasma Sources Science and Technology*
  **14**, 722–733 (2005), DOI
  [10.1088/0963-0252/14/4/011](https://doi.org/10.1088/0963-0252/14/4/011).
- **Accessed:** 2026-08-08 from a University of Antwerp plasma-school course
  packet containing the complete article; article PDF pages 14–19 were
  rendered and visually checked against the extraction.
- **Local text:**
  `research_sources/thesis_extracts/hagelaar_pitchford_2005_bolsig.txt`
- **Topic:** deterministic two-term electron kinetics, conservative
  energy-space discretization, transport/rate moments, and solver evidence
  boundaries.

## Verified source statements

| ID | Source statement | Executable consequence |
|---|---|---|
| H1 | Equations 3 and 9 use the two-term expansion `f = f0 + f1 cos(theta)` and normalize the isotropic energy distribution with `integral sqrt(epsilon) F0 d epsilon = 1`. | The reactor kinetic state stores an EEPF, not an arbitrary probability density, and normalization is a bordered equation of the solve. |
| H2 | Equations 5–7 distinguish the isotropic and anisotropic equations and eliminate `F1` to obtain a convection–diffusion continuity equation in energy with a nonlocal collision source. | A Maxwellian rate table is not relabeled as a Boltzmann solution; physical coefficient and collision-source assembly remain explicit later gates. |
| H3 | Equations 36–38 define the reaction-rate moment, mean energy, mobility, and diffusion from the converged distribution. | Distribution moments are implemented as linear operators with analytic JVP/VJP contracts. Transport definitions remain distinct from the later pulsed-Townsend and steady-state-Townsend validation observables. |
| H4 | Equations 44–45 discretize the energy-space convection–diffusion flux with an exponential Scharfetter–Gummel scheme. The article states that it is exact when local coefficients are constant and preserves nonnegative solutions in the regimes considered. | The reference backend uses conservative face fluxes and must reproduce the constant-coefficient exponential equilibrium before collision physics can be trusted. |
| H5 | The low- and high-energy boundaries impose zero energy-space electron flux, while normalization closes the otherwise singular steady system. | Both boundary fluxes are represented explicitly; a compatibility multiplier and residual ledger expose inconsistent manufactured sources rather than hiding them. |
| H6 | Section 3 reports fast convergence for many fluid-model uses but also documents failure of the two-term approximation when the distribution is strongly anisotropic. | This solver layer can supply reactor EEDF/rate infrastructure only. It cannot pass the direct chlorine swarm board or claim multi-term angular accuracy. |
| H7 | Equations 22--25 derive a high-frequency two-term heating operator with a momentum-transfer factor `sigma_tilde/(2(sigma_tilde^2+q^2))` for the printed peak field amplitude, where `q = (omega/N)/(gamma sqrt(epsilon))`. The current BOLSIG+ manual instead defines AC `E/N` as RMS so that the AC equation tends exactly to the DC result as frequency tends to zero. | The implementation uses the equivalent RMS-field form, so the explicit factor `1/2` is absorbed by converting peak amplitude to RMS amplitude. A manufactured constant-cross-section gate checks the exact frequency-suppression factor and the zero-frequency DC limit. |
| H8 | The article warns that the quasi-stationary high-frequency ansatz is appropriate when energy loss during a cycle is small, gives GHz-range microwave/optical fields as the practical regime, and says intermediate frequencies require a more complete time-dependent Boltzmann solution. | Applying the operator to Malyshev's 13.56 MHz ICP is evidence-bounded as a sensitivity. It cannot by itself certify reactor state or replace the future time-periodic/nonlocal electron-heating rung. |

## Executable decision

The article is the primary authority for the fixed-grid conservative kinetic
foundation. The implementation is an original MIT-repository implementation;
neither BOLSIG+ source/binaries nor third-party LXCat data are copied into the
repository. BOLSIG+ and LGPL BOLOS may be used locally as independent numeric
oracles, never as redistributed dependencies or evidence substitutes.

This entry does not upgrade reactor, wafer-flux, or feature-depth status. The
next physics gates are collision-source assembly, particle/power closure, and
transport-definition-correct multi-term density-gradient/SST validation.
