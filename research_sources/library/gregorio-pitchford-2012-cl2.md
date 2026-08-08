# gregorio-pitchford-2012-cl2

**Updated compilation of electron--Cl2 scattering cross sections**

- **Citation:** J. Gregorio and L. C. Pitchford, “Updated compilation of
  electron--Cl2 scattering cross sections,” *Plasma Sources Science and
  Technology* **21**, 032002 (2012).
- **DOI:** `10.1088/0963-0252/21/3/032002`
- **Author-posted full text:**
  `https://www.researchgate.net/publication/254499582_Updated_compilation_of_electron-Cl2_scattering_cross_sections`
- **Status:** PRIMARY FULL TEXT READ ONLINE; NOT LOCALLY ARCHIVED; FIGURE 4
  NUMBERS NOT IMPORTED
- **Topic:** molecular-chlorine electron-collision compilation, swarm
  inversion, and momentum-transfer provenance

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | The recommended set is an updated compilation, not a single direct measurement. The authors compare a two-term Boltzmann calculation with the then-available swarm data and make small modifications to improve that comparison. | Every channel must retain its individual measured/calculated/regressed provenance; “recommended” does not mean measured. |
| Q2 | The molecular momentum-transfer cross section is Rescigno's 1994 close-coupling calculation with a low-energy extrapolation suggested by Christophorou and Olthoff. Agreement between Rescigno's calculated total-elastic cross section and Gote--Ehrhardt measurements is used as indirect support for the calculated momentum-transfer result. | The COMSOL/SIGLO Cl2 momentum-transfer lineage is theory plus indirect consistency evidence, not a momentum-transfer measurement. |
| Q3 | Set 2 leaves the momentum-transfer row unchanged but multiplies the electronic-excitation and dissociation cross sections by 1.3 to improve agreement with older swarm ionization data. | Set 2 is explicitly swarm-adjusted. A reactor deck may replay it, but cannot label all rows ab initio or independently measured. |
| Q4 | The authors obtain generally poor agreement with the old drift-velocity and characteristic-energy measurements. They decline to tune momentum transfer or vibration to those data because the measurements are not reliable enough, and call for new measurements and a new analysis. | The 2012 drift/energy comparison is an open validation gate, not a successful transport board. |
| Q5 | The complete cross-section set was distributed through SIGLO/LXCat. | The database is the implementation route, but its current terms and dataset metadata are a separate evidence/licensing layer. Raw database bytes are not redistributed here. |

## Executable decision

No numeric cross section is imported from this paper. The official COMSOL
table remains eligible only for hash-pinned implementation replay until its
nodes are mapped to a primary/evaluated source and checked against the newer
swarm evidence. The exact elastic-energy operator can consume a future
approved momentum-transfer table without changing this provenance decision.
