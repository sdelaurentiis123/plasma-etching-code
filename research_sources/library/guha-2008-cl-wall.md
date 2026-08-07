# guha-2008-cl-wall

**Atomic and molecular chlorine interaction with conditioned anodized Al**

- **Citation:** J. Guha, V. M. Donnelly, and Y.-K. Pu, “Mass and Auger
  electron spectroscopy studies of the interactions of atomic and molecular
  chlorine on a plasma reactor wall,” *Journal of Applied Physics* **103**,
  013306 (2008).
- **DOI:** `10.1063/1.2828154`
- **PDF SHA-256:**
  `ae2bc2d1016c01920a1a17cb377786aa96991185f999fadd3d7f6318e659b668`
- **Local extraction:**
  `research_sources/thesis_extracts/guha_2008_cl_wall_interactions.txt`
- **Visual audit:** PDF pages 13--14 were rendered at 450 dpi and visually
  inspected. Render SHA-256 values are
  `f1ceb69d302c191b1ebec4d7507421a853a7ebe4cec25e75af98af97553b4607`
  and `1f3c722295b1737cba90ec5990619bcc4b25276f0ee92e3827571f0f9a712867`.
- **Status:** PRIMARY FULL TEXT + FIGURES 13--14 VISUALLY AUDITED

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| G1 | In a 600 W, 5 mTorr Cl2 plasma, the anodized-Al wall acquired an approximate `Al2Si2O10Cl3` surface layer; silicon came from quartz-tube erosion. | The relevant boundary is conditioned oxychloride, not bare anodized Al. |
| G2 | Direct time-dependent desorption measurements give `gamma_Cl` from about 0.01 to 0.1 and an increase with `n_Cl/n_Cl2`, consistent with competition between Cl and Cl2 for sites. | Independent support for ratio-dependent wall kinetics over the measured ratio range. |
| G3 | Measurements over `n_Cl/n_Cl2 = 0.1--0.8` collapse onto a common trend, and the authors explicitly warn against extrapolating it beyond the investigated densities. | Strict evidence domain for any interpolation; no low- or high-ratio extrapolation. |
| G4 | At 1.25 mTorr the estimated Cl mean free path is about 4 cm, comparable with the wall distance. Franck--Condon Cl born near 0.5 eV can reach the wall before thermalizing; evaluating flux at 300 K then underestimates incident flux and overestimates gamma. At 5 mTorr and above, the paper considers thermalization reasonable because the mean free path is below 1 cm. | A Maxwellian 300 K wall-flux law must fail closed unless thermalization is supported. This is a transport-validity condition, not extra gamma uncertainty. |
| G5 | Ion-bombardment effects on chlorine recombination were not investigated. | Prevents importing these data into an ion-conditioned wall model without a separate study. |

## Use decision

Guha is the direct anodized-Al complement to the Stafford stainless-steel
dataset. Its strongest new contribution is the non-thermal low-pressure
warning: direct surface probabilities remain conditional on how incident atom
flux was inferred. The reactor code therefore requires a provenance-bearing
`ChlorineIncidentVelocityState`; a Maxwellian state must carry its reference
temperature, and assumed thermalization cannot support a predictive wall-loss
claim even when the wall-response measurement itself is direct.
