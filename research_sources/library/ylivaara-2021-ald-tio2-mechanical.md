# ylivaara-2021-ald-tio2-mechanical

**Ylivaara et al. — measured elastic modulus and hardness of ALD TiO2 (TiCl4/H2O)**

- **Citation:** O. M. E. Ylivaara, A. Langner, X. Liu, D. Schneider, J. Julin,
  K. Arstila, S. Sintonen, S. Ali, H. Lipsanen, T. Sajavaara, S.-P. Hannula and
  R. L. Puurunen, "Mechanical and optical properties of as-grown and thermally
  annealed titanium dioxide from titanium tetrachloride and water by atomic
  layer deposition," *Thin Solid Films* **732**, 138758 (2021).
- **DOI:** `10.1016/j.tsf.2021.138758`
- **License:** CC BY (publisher's PDF / version of record).
- **Retrieval route that worked:** Aaltodoc DSpace API,
  `https://aaltodoc.aalto.fi/server/api/core/bitstreams/bbf794b3-ad46-4621-9bbe-1588075509e3/content`
  (ScienceDirect returns 403; JYX record is metadata-only).
- **Status:** PRIMARY FULL TEXT ARCHIVED —
  `research_sources/thesis_extracts/ylivaara_2021_ald_tio2_mechanical.txt`
  (Table 3 numeric rows extracted).
- **Topic:** Young's modulus, hardness, residual stress, density and phase of
  ALD TiO2 vs growth temperature

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | Table 3 (verbatim rows, ALD temperature / thickness / nanoindentation E ± st.dev / hardness / LSAW E ± st.dev): 80 °C, 97 nm → 156 ± 19 GPa, 6.7 ± 0.3 GPa, 107.4 ± 0.4 GPa. 110 °C, 98 nm → 152 ± 5 GPa, 6.9 ± 0.1 GPa, 126.7 ± 0.5 GPa. 150 °C, 98 nm → 149 ± 4 GPa, 7.3 ± 0.1 GPa, 129.0 ± 0.5 GPa. | **The amorphous-film values.** Two independent methods disagree by ~20–45 % (nanoindentation high, LSAW low); carry the whole band, not one number. Films are TiCl4/H2O ALD on **silicon**, ~100 nm thick. |
| Q2 | Table 3, crystalline rows: 250 °C, 102 nm → 159 ± 7 / 10.5 ± 0.8 / 153.4 GPa; 300 °C, 102 nm → 165 ± 16 / 9.7 ± 1.0 / 166.5 GPa; 350 °C, 98 nm → 165 ± 15 / 9.6 ± 1.4 / 167.5 GPa. | Anatase branch. Phase matters more for LSAW than for nanoindentation. |
| Q3 | Sec. 3.1, verbatim: "The as-grown TiO2 was amorphous at low growth temperatures, 80, 110 and 150 ◦C, and at lower film thicknesses… A polycrystalline film with anatase was observed for about 50 nm TiO2 film, grown at 200 ◦C. For films grown at 300 ◦C, anatase phase was already detected for about 10 nm film. No other phases were detected." | Phase assignment for Q1/Q2. |
| Q4 | Sec. 2, verbatim: "The elastic modulus values were calculated using Poisson's ratio of 0.27 [46] for ALD TiO2." Ref [46] = L. Borgese, M. Gelfi, E. Bontempi, P. Goudeau, G. Geandier, D. Thiaudière, L. E. Depero, "Young modulus and Poisson ratio measurements of TiO2 thin films deposited with atomic layer deposition", *Surf. Coat. Technol.* **206**, 2459–2463 (2012), DOI `10.1016/j.surfcoat.2011.10.050`. | ν = 0.27 for ALD TiO2, sourced. Borgese itself **not** read — its reported E (~151 GPa, per third-party summary) is `[unquoted — verify on next use]`. |
| Q5 | Sec. 2, verbatim method: "Nanoindentation was performed with a Hysitron Triboindenter for the elastic modulus and hardness using a cube-corner indenter with a 90◦ total induced angle and 40 nm tip radius. The indent depth was kept less than 10 % of the film thickness." | Substrate-effect control; supports using these as film-intrinsic values. |
| Q6 | Sec. 3.2 / Table 3, verbatim: "ALD TiO2 films were under tensile stress in the scale of hundreds of MPa" (values 230–967 MPa across the table). | Residual tensile stress is present but is ~3 orders below E; not a first-order term in the collapse criterion. Relevant if the fracture branch is ever activated. |

## Use in petch

Sole sourced Young's modulus for amorphous ALD TiO2 in
`research_sources/RESEARCH_PATTERN_COLLAPSE_CRITERION_2026-08-20.md`.
Boundary: measured on planar films on Si, TiCl4/H2O chemistry, ~100 nm thick.
A different precursor (e.g. TDMAT), a different temperature, or the 700 nm
device stack supersede these immediately.
