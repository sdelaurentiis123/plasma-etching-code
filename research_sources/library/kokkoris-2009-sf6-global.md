# kokkoris-2009-sf6-global

**SF6 daughter and heavy-particle global chemistry**

- **Citation:** G. Kokkoris, A. Panagiotopoulos, A. Goodyear, J. Cooke, and
  E. Gogolides, “A global model for SF6 plasmas coupling reaction kinetics in
  the gas phase and on the surface of the reactor walls,” *Journal of Physics
  D: Applied Physics* **42**, 055209 (2009).
- **DOI:** `10.1088/0022-3727/42/5/055209`
- **Focused extraction:**
  `research_sources/thesis_extracts/kokkoris_2009_sf6_table1_verified.txt`
- **Status:** PRIMARY FULL TEXT + TABLE 1 VISUALLY AUDITED IN THE READ-ONLY
  REACTORLAB SOURCE MIRROR; EXECUTABLE SUBSET CONSERVATION-CHECKED

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| K1 | Equation 2 fits each electron rate as `ln(k)=A+B ln(Te)+C/Te+D/Te^2+E/Te^3`, with separate Druyvesteyn and Maxwellian coefficient rows. | Both shapes are executable sensitivity branches. Neither is called an exact rate for Petch's solved non-Maxwellian EEPF. |
| K2 | Table 1 rows G4--G7 and G11--G16 supply daughter-SFx/F2 dissociation and ionization; G19 supplies F2 dissociative attachment. | These rows extend the measured parent-feed deck without repeating parent SF6 rows G1--G3, G8--G10, or G17--G18. |
| K3 | G20--G21 give neutral-collision detachment; G38--G41 neutral rearrangement; G42--G49 ion-ion neutralization; and G50 ion-molecule conversion. | These close important daughter/charged paths while preserving explicit atoms, charge, and electrons. |
| K4 | Table-1 footnote f says the printed G35--G37 `F + SFx` fall-off coefficients are valid at 2 Pa. | The Zhu NPG80 recipe is 3e-2 Torr = 3.9997 Pa. Those rows are quarantined until the original pressure law is recovered; the 2-Pa values are not silently copied. |
| K5 | The paper's wall coefficients were fitted/selected for its Oxford ICP380 reactor and surface state. | Wall topology can inform the model, but numerical wall probabilities are not transferred to the different single-source NPG80 CCP. |

## Executable decision

`build_zhu_supplemental_chemistry()` installs 36 Kokkoris daughter/heavy
reactions plus ten explicitly selected Sandia CHF2/recombination rows. It
exposes the two Kokkoris EEDF assumptions, excludes all measured-parent
duplicates, and excludes the three pressure-specific fall-off rows. Oxygen
heavy chemistry, target-machine wall probabilities, and true daughter
cross-section integrations remain open physics layers.
