# pateau-2014-sf6-o2

**O2 daughter and SFx/O titration chemistry**

- **Citation:** A. Pateau, A. Rhallabi, M.-C. Fernandez, M. Boufnichel, and
  F. Roqueta, “Modeling of inductively coupled plasma SF6/O2/Ar plasma
  discharge: Effect of O2 on the plasma kinetic properties,” *Journal of
  Vacuum Science & Technology A* **32**, 021303 (2014).
- **DOI:** `10.1116/1.4853675`
- **Focused extraction:**
  `research_sources/thesis_extracts/pateau_2014_sf6_o2_verified.txt`
- **Status:** PRIMARY FULL TEXT + TABLES I, III, V, VII VISUALLY AUDITED IN
  THE READ-ONLY REACTORLAB SOURCE MIRROR; EXECUTABLE SUBSET
  CONSERVATION-CHECKED

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| P1 | Table V gives O/O(1d), O2+, O+, and O- daughter reactions, including O2+ dissociative recombination and associative detachment. | Parent O2 dissociation/ionization/attachment rows are replaced by the newer measured Song deck; only daughter/closure rows R75--R85 are installed. |
| P2 | Table VII gives SFx/O reactions that form SOxFy species; several release F. The paper attributes the O2-driven F increase to these titration channels and increased dissociation. | The full conserved R119--R137 topology is installed except R138. This is the key mechanistic bridge from a small O2 feed to fluorine delivery. |
| P3 | Table III R38 applies detachment to seven SFx/F negative ions over nine neutral partners; R50 applies a common ion-ion neutralization coefficient over eight positive and seven negative ions. | Only pairs absent from the Kokkoris block are added, preserving every parent negative ion's volume loss path without double counting. |
| P4 | R138 is printed as unimolecular but carries the table's bimolecular units. | R138 is quarantined; no arbitrary first-order reinterpretation is imported. |
| P5 | The model is an SF6/O2/Ar ICP at 1500 W, 20--40 mTorr, 200 sccm, R=18 cm and L=17.5 cm. | Reaction topology and sourced coefficients are useful; its reactor geometry, coupling efficiency, and wall state are not transferred to the 150 W NPG80 CCP. |
| P6 | The paper assumes a Maxwellian EEDF and uses rate fits; its calculated F and electronegativity trends are compared with independent Pessoa measurements. | It is a validation target and closure source, not evidence that a Maxwellian represents the NPG80 EEPF. |

## Executable decision

The Zhu supplemental network adds 11 oxygen daughter reactions, 32 SFx/O
titration reactions, and the non-duplicated R38/R50 charge-closure pairs.
Together with the Kokkoris and Sandia blocks it contains 188 conserved volume
reactions over a 56-species reactor basis (including neutral CH for conserved
CH+ wall return). Cross-family O-/fluorocarbon-ion recombination and
CHF3/O/S neutral chemistry remain incomplete and are exposed as limitations.
