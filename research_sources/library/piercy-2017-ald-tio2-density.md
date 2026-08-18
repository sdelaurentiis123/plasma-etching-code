# piercy-2017-ald-tio2-density

**Piercy, Leng & Losego, ALD TiO2 density versus growth temperature**

- **Citation:** B. D. Piercy, C. Z. Leng, and M. D. Losego, “Variation in
  the density, optical polarizabilities, and crystallinity of TiO2 thin films
  deposited via atomic layer deposition from 38 to 150 C using the titanium
  tetrachloride-water reaction,” *J. Vac. Sci. Technol. A* **35**, 03E107
  (2017).
- **DOI:** `10.1116/1.4979047`
- **Primary record:** `https://doi.org/10.1116/1.4979047`
- **Status:** PRIMARY ARTICLE RECORD + NUMERICAL ABSTRACT READ; FULL PDF NOT
  LOCALLY ARCHIVED
- **Topic:** ALD TiO2 mass density and phase

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | X-ray reflectometry and optical ellipsometry give mass densities increasing from `3.25` to `3.68 g cm-3` between `38` and `125 C` for TiCl4/H2O ALD TiO2. | Low endpoint of a Zhu film-density sensitivity only. The deposition precursor, temperature, and phase of Zhu's 700 nm film are still unknown, so this is not a selected target coefficient. |
| Q2 | Raman and XRD show amorphous films below `150 C` and anatase at or above `150 C`. | Supports keeping film phase explicit rather than assigning rutile bulk density from the word “TiO2.” |

## Use in petch

The `3.25 g cm-3` endpoint participates only in the target-free atom/dose
clearance gate in `audit_zhu_npg80_tio2_depth_gate.py`.  It is not installed in
the reactor or TiO2 surface model.  A condition-specific density measurement
supersedes this sensitivity immediately.
