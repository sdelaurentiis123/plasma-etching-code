# stafford-2010-cl-wall

**Measured, state-dependent O and Cl wall recombination**

- **Citation:** L. Stafford, J. Guha, R. Khare, S. Mattei,
  O. Boudreault, B. Clain, and V. M. Donnelly, “Experimental and
  modeling study of O and Cl atoms surface recombination reactions in O2
  and Cl2 plasmas,” *Pure and Applied Chemistry* **82**, 1301–1315
  (2010).
- **DOI:** `10.1351/PAC-CON-09-11-02`
- **Official full text:**
  `https://iupac.org/publications/pac/pdf/2010/pdf/8206x1301.pdf`
- **PDF SHA-256:**
  `2f74388576d435d9b3ae3843d5fa14f6a941ef61124b406f3ee7a7496e464b08`
- **Status:** PRIMARY FULL TEXT READ; PAGES 1301/1309–1312 VISUALLY
  CHECKED; FIGURE 8 PIL-AUDITED AT 600 DPI; 39 MARKERS FROZEN
- **Local extraction:**
  `research_sources/thesis_extracts/stafford_2010_cl_wall_recombination.txt`
- **Digitization:**
  `data/experimental/stafford_2010/figure8_chlorine_wall_recombination.csv`
  and adjacent `digitization_manifest.json`

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| S1 | Reactor Cl atoms are lost primarily through heterogeneous wall/substrate recombination, and the wall is dynamically conditioned by plasma exposure. | Wall state is a reactor state variable, not a universal material constant. |
| S2 | The spinning-wall experiment separates delayed Langmuir–Hinshelwood recombination and corrects the plasma-on Cl2 signal using the measured plasma-off adsorption/desorption response. | Establishes the measurement basis behind Figure 8; the points are not generic model outputs. |
| S3 | Conditioned stainless steel and anodized aluminum acquire silica-rich oxychloride layers due to discharge-tube erosion. | Bare SS, bare Al, smooth quartz, and seasoned reactor walls cannot share one gamma. |
| S4 | Figure 8 spans 1.25–20 mTorr and 100–600 W; marker shape encodes pressure but individual point powers are not published. | The digitized table leaves per-point power unidentified rather than inventing it. |
| S5 | Conditioned-SS gamma_Cl is reported as 0.004–0.03 and rises with n_Cl/n_Cl2; conditioned anodized-Al values follow the same trend at roughly twice the SS values. | Supplies a measured envelope and material ordering within the published domain. |
| S6 | Cl2 competes with Cl for adsorption sites; absolute molecular impingement density also matters, so n_Cl/n_Cl2 alone is not sufficient. | Forbids a one-coordinate empirical law from being labeled fundamental closure. |
| S7 | In Figure 9 no constant gamma_Cl curve fits the full pressure-dependent dissociation dataset; the required gamma rises with pressure. | Negative validation gate for any constant-gamma reactor solver. |
| S8 | In a conditioned/roughened quartz tube the inferred gamma_Cl is 0.02–0.04, versus less than 0.005 for smooth annealed quartz under an atomic beam. | Quantifies the magnitude of reactor seasoning and roughness effects. |
| S9 | Figure 8 is reproduced from the authors’ refs. 31 and 34. | The 2010 figure is suitable for an evidence board; retrieve the two original articles before promoting an interpolation to a production provider. |

## Executable decision

The 39 direct Figure-8 markers are frozen as a no-fit wall-model validation
board. They replace the idea of selecting one chlorine recombination
coefficient. They do **not** yet define an executable wall law: pressure and
power are partly confounded, individual powers are absent, absolute impingement
density is not tabulated, and the original Figure-8 articles have not yet been
independently archived.

The next native reactor implementation must therefore expose wall-condition
state and fail outside a declared evidence domain. It may be graded against
this board, but feature depth must not choose its wall parameters.
