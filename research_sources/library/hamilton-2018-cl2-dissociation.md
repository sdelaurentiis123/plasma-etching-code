# hamilton-2018-cl2-dissociation

**State-resolved electron-impact dissociation of molecular chlorine**

- **Citation:** J. R. Hamilton, J. Tennyson, J.-P. Booth, T. Gans, and
  A. R. Gibson, “Calculated electron impact dissociation cross sections for
  molecular chlorine (Cl2),” *Plasma Sources Science and Technology* **27**,
  095008 (2018).
- **Article DOI:** `10.1088/1361-6595/aada32`
- **Data DOI:** `10.15124/b11c65cf-2913-4c63-a522-2f57006cfb8a`
- **Article PDF SHA-256:**
  `7073912aeeeacedb82c75bcb3018e721f1b6ea7d596c721e2fd81d1808d57add`
- **Official OPJ SHA-256:**
  `1e53fd091c1685f38326a9da5a1e78ba97387bc9c375783e8bb6c2a1d7fa0272`
- **Status:** PRIMARY FULL TEXT + OFFICIAL CC-BY DATASET EXTRACTED; FIGURES
  4--5 VISUALLY AUDITED AT 600 DPI; 235-POINT RATE REPRODUCTION PASS
- **Local text:**
  `research_sources/thesis_extracts/hamilton_2018_cl2_dissociation.txt`
- **Extraction manifest:**
  `research_sources/digitized/hamilton_2018_cl2_dataset_manifest.md`

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| H1 | Industrial low-temperature-plasma electron temperatures are stated as 0.3--5 eV, making near-threshold energy resolution material. | Freezes the runtime rate-table domain; no extrapolation. |
| H2 | Table 2 gives calculated vertical excitation energies of 3.252, 4.348, 6.498, 7.537, 7.790, 7.257, 8.228, and 9.219 eV for the eight retained dissociative states. | Event-specific electron-energy losses replace Lee's nonphysical use of one rate-fit exponent. |
| H3 | All dissociative states retained in the results lead to two ground-state Cl atoms. | Makes each executable reaction exactly atom- and charge-conserving: `e + Cl2 -> e + 2Cl`. |
| H4 | The R-matrix calculation is fixed-nuclei `Cl2(v=0)` and uses transition-specific scaling above the ionization potential. | Evidence class is `semi_empirical`, not measured or first-principles exact. |
| H5 | Cosby's experiment contains a distribution of vibrational states and is not directly comparable to the `v=0` calculations. | Forbids relabeling broad graphical consistency as an independent experimental pass. |
| H6 | Figure 5 defines `x=1` as Maxwellian and supplies a total dissociation-rate curve versus effective electron temperature. | Independent numerical target for the analytic integration, not a coefficient-selection target. |
| H7 | The article data-management section points to the official York dataset DOI; the archive contains the exact arrays behind the figures. | The 50,000-point cross sections supersede plot digitization. |

## Executable decision

Eight state-resolved Maxwellian providers and reactions land from the official
dataset. Their summed rate reproduces the authors' separately supplied
Maxwellian total to `0.4443%` maximum and `0.1855%` mean over 235 temperatures.
The runtime reduction evaluates all eight channels in roughly 22 microseconds
on the development CPU and fails outside 0.3--5 eV.

This closes the neutral-dissociation particle source and its excitation-energy
ledger. The evaluated particle-deck builder replaces exactly the one legacy
neutral-dissociation row with these eight channels while leaving the original
Lee deck available for source reproduction. It does not close molecular-ion
branching, attachment/detachment power exchange, wall recombination, absorbed
RF power, or a Lam reactor.
