# wang-olthoff-1999-ion-flux

**Absolute, mass-resolved chlorine ion fluxes in an inductively coupled GEC
reactor**

- **Citation:** Y. Wang and J. K. Olthoff, “Ion Energy Distributions in
  Inductively Coupled Radio-Frequency Discharges in Argon, Nitrogen, Oxygen,
  Chlorine, and Their Mixtures,” *Journal of Applied Physics* **85**,
  6358–6365 (1999).
- **DOI:** `10.1063/1.370138`
- **Official NIST PDF:**
  `https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=1506`
- **PDF SHA-256:**
  `17702f0ffb904cca42760867c693c382b174c68f02ef7b5064ec95933adea460`
- **Local text:**
  `research_sources/thesis_extracts/wang_olthoff_1999_ion_energy.txt`
- **Figure-9 data:**
  `data/experimental/wang_olthoff_1999/figure9_chlorine_ion_flux.csv`
- **Digitization manifest:**
  `data/experimental/wang_olthoff_1999/digitization_manifest.json`
- **Status:** PRIMARY NIST FULL TEXT + FIGURE 9 PIL-AUDITED AT 600 DPI

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| W1 | The experiment is a five-turn, 13.56 MHz planar ICP in a GEC reference cell, with a 4.13 cm quartz-window-to-grounded-electrode gap and 10 sccm flow. | Defines a reproducible ICP geometry; it is neither a Lam proprietary chamber nor a dielectric-etch CCP. |
| W2 | Reported RF power is net power to the coil matching network; the authors state that actual plasma-dissipated power is approximately 80% of the reported value. | Keeps the 300 W instrument node separate from the approximate 240 W plasma boundary. |
| W3 | Total absolute ion flux is derived from sampled current divided by the 10 µm orifice area; mass-resolved IED intensities are scaled to that independently measured total. | Separates the absolute-current observable from the species-fraction observable. |
| W4 | Figure 9 gives total, `Cl+`, and `Cl2+` flux versus pressure at 300 W for pure `Cl2`, and total, `Cl+`, `Cl2+`, and `Ar+` for 20% `Cl2` / 80% Ar. | Supplies held-out reactor/sheath targets for magnitude, mixture transfer, and species composition. |
| W5 | In pure `Cl2`, the total flux increases from about 7.46 to 17.04 mA/cm² over 2.7–6.7 Pa, while the `Cl+` fraction rises from about 0.83 to 0.96 and `Cl2+` falls. | Stronger gate than electron density alone; tests both total production/loss and the species network. |
| W6 | The pure-`Cl2` `Cl+` dominance differs from Woodworth et al.’s nearly equal `Cl+`/`Cl2+` intensities; the authors identify surface conditions as the likely cause. | Wall state/recombination must remain an explicit measured or bounded boundary. The disagreement cannot be used to fit an electron-impact branching fraction. |
| W7 | The paper reports ±1 eV energy-scale uncertainty for IEDs and describes total-current repeatability for argon, but Figure 9 has no chlorine error bars. | The committed flux allowance covers digitization only; no chlorine measurement uncertainty is invented. |

## Executable decision

The 24 Figure-9 markers are frozen before a mixed Ar/Cl2 reactor solver exists.
The current evaluated molecular-ionization table remains aggregate: reactor
species fluxes may validate a complete chemistry/wall/sheath model, but they do
not identify the elementary `Cl2+` versus `Cl+` ionization branch by themselves.
