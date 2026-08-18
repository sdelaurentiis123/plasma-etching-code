# voloshin-2007-chf3-kinetics

**Validated gas-phase kinetics in CHF3:H2:O2 mixtures**

- **Citation:** D. G. Voloshin, K. S. Klopovskiy, Y. A. Mankelevich,
  N. A. Popov, T. V. Rakhimova, and A. T. Rakhimov, "Simulation of
  Gas-Phase Kinetics in CHF3:H2:O2 Mixtures," *IEEE Transactions on Plasma
  Science* **35** (2007), 1691--1703.
- **DOI:** `10.1109/TPS.2007.906780`
- **Author-posted full text read online:**
  `https://www.researchgate.net/publication/3168274_Simulation_of_Gas-Phase_Kinetics_in_hboxCHF_3hboxH_2_hboxO_2_Mixtures`
- **Status:** PRIMARY AUTHOR-POSTED FULL TEXT READ ONLINE; TABLE III PIXEL
  TRANSCRIPTION OPEN

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| V1 | The modeled neutral chain includes `CHF3 + F -> CF3 + HF`, electron-impact HF dissociation that regenerates F and H, and H-driven conversion of CFx into CFx-1 plus HF. | This topology is necessary to avoid treating feed fractions as steady plasma composition. |
| V2 | The article text gives `1.82e-12 cm3/s` for `CHF3 + F -> CF3 + HF` at the modeled `Tg=350 K` state. | Landed as the default 350 K branch in `lim_2014_chf3_oxygen_chemistry.py`; it is not averaged with Lim's conflicting 700 K compilation. |
| V3 | Increasing an assumed radical wall-loss probability from `0.05` to `0.1` changes calculated radical densities by roughly 15%. | Wall recombination is an explicit uncertainty direction, not a hidden fit knob. |
| V4 | Validation used measured electron density and electron temperature as model inputs because generator power did not uniquely close the electron state. | The paper validates chemistry at a diagnosed plasma state; it does not justify converting the Zhu `150 W` forward demand directly into electron density or field. |

## Executable decision

The source supports topology, the text-resolved R13 coefficient, and a
declared wall-loss sensitivity.  Its full Table III remains pixel-only and is
not transcribed.  Sandia Table 9 and Lim Table I provide the executable
daughter extension; the Voloshin R13/Lim R26 disagreement is an explicit
source branch.  The validation used measured electron-state inputs and does
not by itself close an equipment power balance.
