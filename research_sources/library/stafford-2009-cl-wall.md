# stafford-2009-cl-wall

**Direct chlorine recombination on plasma-conditioned stainless steel**

- **Citation:** L. Stafford, R. Khare, J. Guha, V. M. Donnelly,
  J.-S. Poirier, and J. Margot, “Recombination of chlorine atoms on
  plasma-conditioned stainless steel surfaces in the presence of adsorbed
  Cl2,” *Journal of Physics D: Applied Physics* **42**, 055206 (2009).
- **DOI:** `10.1088/0022-3727/42/5/055206`
- **PDF SHA-256:**
  `39a33b7c677f46750851ed243bfc152edf334c067e894a8e111889ccc39e4b70`
- **Local extraction:**
  `research_sources/thesis_extracts/stafford_2009_cl_wall_recombination.txt`
- **Visual audit:** PDF pages 7--8 were rendered at 450 dpi and visually
  inspected. Render SHA-256 values are
  `5261ae653dff646ee808a7369a0e92468f2890efc8b43a699a6254c449e32d83`
  and `b44335589d2ecf09ec1e0f6f57e3d76e094f92f69f3f4eaee9789b4649294490`.
- **Status:** PRIMARY FULL TEXT + FIGURES 8--10 VISUALLY AUDITED

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| S1 | After hours in the Cl2 plasma, electropolished 304L stainless steel was coated by a Si-oxychloride-rich layer with approximate Fe:Si:O:Cl composition 1:13:13:3. | Wall material labels must describe the plasma-conditioned surface, not the unloaded alloy. |
| S2 | Spinning-substrate measurements give `gamma_Cl = 0.004--0.03` on conditioned stainless steel and an increase with `n_Cl/n_Cl2`. | Direct local wall evidence; forbids a universal scalar recombination probability. |
| S3 | Under matched plasma conditions, conditioned-stainless values are about two to three times lower than conditioned anodized-Al values. | Surface state materially changes reactor composition; area fractions and conditioning require explicit evidence. |
| S4 | Re-analysis of chlorine dissociation versus pressure shows the best effective gamma falls from about 0.03 at 0.2 mTorr to 0.01--0.015 at 10 mTorr. Treating gamma as a function of Cl/Cl2 improves agreement over a constant. | Mechanistic cross-check only: those effective values came from reactor-model comparison, while Figure 8 is the direct measurement. |
| S5 | The measured recombination probability is defined from delayed Cl2 desorption divided by the incident Cl flux; the incident thermal speed was evaluated at an assumed 300 K gas/wall temperature. | The probability is conditional on the incident-flux model. Non-thermal incident Cl requires a separate boundary. |
| S6 | The surface mechanism is interpreted as Langmuir--Hinshelwood recombination with adsorbed Cl2 blocking Cl adsorption/recombination sites. | Supports state dependence on Cl/Cl2 ratio; does not by itself provide a first-principles site-kinetic parameter set. |

## Use decision

This primary paper independently confirms the state-dependent wall law already
frozen from Stafford et al. 2010 and explains why the `gamma_Cl = 0.035` Lam
fit is not portable. It also exposes the thermal-incident-flux assumption that
must be explicit before a local gamma is called predictive.
