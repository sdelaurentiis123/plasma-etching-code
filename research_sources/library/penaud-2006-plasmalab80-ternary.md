# penaud-2006-plasmalab80-ternary

**Oxford PlasmaLab 80 SF6/O2/CHF3 self-bias and anisotropy table**

- **Citation:** Julien Penaud, *Contributions à la conception et à la
  réalisation de transistors MOS à grille multiple*, PhD thesis, Université
  des Sciences et Technologies de Lille (2006).
- **Official full text:**
  `https://pepite-depot.univ-lille.fr/LIBRE/Th_Num/2006/50376-2006-Penaud.pdf`
- **PDF SHA256:**
  `fceed7b5bd6d2ed43a7fe5fe3b1f1b8eeb8ae338ad3d1144326cbc41d2ba3e46`
- **Status:** PRIMARY FULL THESIS; TABLES 3.2--3.5 AND FIGURES 3.13--3.17
  VISUALLY AUDITED

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| P1 | The patterns are transferred with an Oxford PlasmaLab 80 RIE. | Establishes the equipment family for the tabulated process measurements. |
| P2 | Table 3.5 gives `SF6/O2/CHF3 = 15/10/10 sccm`, `10 mTorr`, `50 W`, `276 V` DC bias, and anisotropy angle `9.2 +/- 2 degrees`. | Closest primary bias anchor: same active-gas set and exactly the target recipe's `P/p = 5 W/mTorr`. Gas fractions, target material, loading, and exact machine differ. |
| P3 | Tables 3.2--3.4 give `256 V` for CF4/O2, `267 V` for SF6/O2, `333 V` for SF6/O2/Ar, and `270 V` for SF6/CHF3 at the same `50 W`, `10 mTorr`. | Direct evidence that chemistry changes self-bias even at fixed reported power and pressure. |
| P4 | The stated mechanism assigns chemical removal to F radicals, sidewall passivation to oxygen-derived SiOxFy, and horizontal passivation removal to positive CFx+/SFx+ ions. | Species/topology evidence for the deterministic chemistry deck; not quantitative TiO2 surface kinetics. |

## Executable decision

The `276 V` point is selected mechanically as the same-active-gas,
same-reduced-drive family anchor. It is not relabeled as a measurement on the
Zhu NPG80 condition and cannot certify absolute depth.
