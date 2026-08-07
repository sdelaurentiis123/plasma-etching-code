# economou-1995-2002-cl-transport

**Chlorine neutral diffusion and wall boundary in two-dimensional ICP models**

- **1995 citation:** D. P. Lymberopoulos and D. J. Economou,
  “Two-Dimensional Simulation of Polysilicon Etching with Chlorine in a High
  Density Plasma Reactor,” *IEEE Transactions on Plasma Science* **23**,
  573–580 (1995), DOI `10.1109/27.467977`.
- **1995 author-hosted full text:**
  `https://www.chee.uh.edu/sites/chbe/files/faculty/economou/ieee_lymberopoulos_icp_95.pdf`
- **1995 PDF SHA-256:**
  `142bdb74f969ddc494a9e659e911f56a66fa12994966a18ea1ee2bc4588934b6`
- **2002 citation:** B. Ramamurthi and D. J. Economou, “Two-Dimensional
  Pulsed-Plasma Simulation of a Chlorine Discharge,” *Journal of Vacuum
  Science & Technology A* **20**, 467–478 (2002),
  DOI `10.1116/1.1450581`.
- **2002 author-hosted full text:**
  `https://www.chee.uh.edu/sites/chbe/files/faculty/economou/jvst_02_pulsed_cl2_2d.pdf`
- **2002 PDF SHA-256:**
  `b2f1450df1bac12f7b405254bc4b03900bf926f0ee327a0f043c6b37441060e8`
- **Status:** TWO PRIMARY FULL TEXTS READ; TRANSPORT EQUATIONS/TABLES AND
  TEMPERATURES VISUALLY AUDITED AT 400–500 DPI
- **Local extraction:**
  `research_sources/thesis_extracts/economou_chlorine_neutral_transport.txt`

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| E1 | The 1995 model solves neutral diffusion separately from charged transport and uses a thermal wall flux for Cl recombination plus ion-neutralization return. | Neutral wall loss and ion return must remain separate ledger terms. |
| E2 | The 1995 paper says neutral diffusivities came from Chapman–Enskog theory using Lennard–Jones parameters from a handbook; it reports `N D_Cl = 6.21e18 cm^-1 s^-1`. | Exact SI conversion is `6.21e20 m^-1 s^-1`; it is a published model coefficient, not a measurement. |
| E3 | The 1995 base condition explicitly uses a gas temperature of 500 K. | The printed coefficient is scoped to 500 K because the collision parameters needed to reconstruct its temperature law are not printed. |
| E4 | The 2002 model uses the Robin boundary `-D grad(n) = [gamma/(2(2-gamma))] n vbar - J_Cl+`, includes a 3.8 ms convective residence time, and again prints `N D_Cl = 6.21e18 cm^-1 s^-1`. | Confirms the partial-reflection boundary coefficient and that pumping is a separate loss. |
| E5 | The 2002 base condition explicitly uses a gas temperature of 300 K and cites the 1995 paper for the transport table. | The identical 300 K/500 K coefficient is a temperature-provenance conflict; both reproductions are quarantined from prediction. |
| E6 | The 1995 simulation uses 10 mTorr, 3560 W plasma power, 200 sccm Cl2, and a nonuniform two-dimensional reactor; the 2002 base case uses 20 mTorr, 320 W peak pulsed power, and a GEC ICP. | Neither coefficient is a measured Lam-tool boundary or a feature-depth calibration. |

## Executable decision

The two printed reduced diffusivities are represented as separate,
single-temperature `published_model` providers. Both carry the conflict in
their provenance and return `supports_prediction = False`. Cross-temperature
use is rejected rather than supplied with an assumed power law.

The exact cylindrical transport solver may use either provider for
source-reproduction checks. A predictive chlorine reactor still needs a
measured or evaluated Cl–Cl2 transport coefficient with uncertainty and
temperature dependence.
