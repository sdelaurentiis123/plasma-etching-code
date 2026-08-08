# malyshev-1998-lam-cl2

**Measured chlorine dissociation in a commercial Lam Alliance TCP reactor**

- **Citation:** M. V. Malyshev, V. M. Donnelly, A. Kornblit, and
  N. A. Ciampa, “Percent dissociation of Cl2 in inductively coupled,
  chlorine-containing plasmas,” *Journal of Applied Physics* **84**, 137--146
  (1998).
- **DOI:** `10.1063/1.368010`
- **PDF SHA-256:**
  `569ab180bb8cab71dda0860e0350f03d24b391e084e19116fefb74f1c719789c`
- **Local extraction:**
  `research_sources/thesis_extracts/malyshev_1998_lam_chlorine_dissociation.txt`
- **Digitized data:** `data/experimental/malyshev_1998_lam/`
- **Status:** PRIMARY FULL TEXT + FIGURES 3, 7--8, 11 NATIVE-PIXEL PIL AUDIT
  + FOOTNOTE 16 EQUATION VISUALLY AUDITED
- **Footnote-16 visual audit:** publisher PDF page 10 was rendered at 500 dpi;
  render SHA-256
  `e5eb5218302f592c2ce01dcd8ad002cb6563077e7fabf06dd5a3a5380c94754d`.
  The equation, units, masses, mixing rules, Lennard--Jones constants, `1.25`
  multiplier, and rounded `0.15 cm2/s` anchor were read from the pixels.

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| M1 | The apparatus is a commercial Lam Research Alliance metal etcher with a high-flow chamber, a 12-inch TCP coil, a 6-inch wafer, anodized-Al walls at 333 K, and selectable 11 or 6.5 cm TCP-window-to-wafer gaps. | Supplies a commercially relevant chlorine reactor board. It is not the geometry of Krüger's fluorocarbon reactor. |
| M2 | TCP source power spans 10--900 W and was measured into the matching network with inline meters. | This is forward-power evidence, not absorbed-plasma-power evidence; a coupling boundary remains mandatory. |
| M3 | The Cl2/rare-gas flows were 114/6 sccm at 2 and 10 mTorr, 57/3 sccm at 1 mTorr, and 38/2 sccm at 0.5 mTorr. Figures 7--8 use a SiO2-coated wafer with substrate RF bias off. | Defines the open-system inputs for a source-condition reproduction. Flow and pressure must not be collapsed into an invented residence time. |
| M4 | Absolute Cl2 density is derived independently with Ar and Xe actinometry; the authors estimate the averaged value is accurate to about +/-25%. Figure-bar endpoints are the two reductions. | The endpoints are not statistical sigma. The dataset preserves the authors' uncertainty statement separately from pixel uncertainty. |
| M5 | Dissociation reaches about 70% at 900 W and 1--2 mTorr and decreases as pressure increases at fixed power. | Primary reactor-validation target; no surface or feature-depth fit is permitted. |
| M6 | Reducing the gap from 11 to 6.5 cm lowers dissociation at 2 mTorr and below but slightly raises it at 10 mTorr. | Strong geometry/transport discriminator that a well-mixed rate fit alone should not receive credit for. |
| M7 | At 6.5 cm and 10 mTorr, the apparent Cl2 density above 100% near 90 W is attributed to enhanced emission during a discharge-mode transition. | Preregistered exclusion from grading, not a negative dissociation datum. |
| M8 | In the paper's model, wall recombination probability was the sole adjusted parameter and `gamma_Cl = 0.035` gave the reported pressure dependence. | A source-model fit, not an independent constant. It is excluded from the reactor validation inputs; direct wall measurements supersede it. |
| M9 | Figure 3 reports OES electron-temperature measurements versus TCP power for both 11 and 6.5 cm gaps; the article states Te rises roughly 20--30% between 20 and 900 W. The 11 cm values were reported in source ref. 3 and the 6.5 cm values were previously unpublished. | The 62 visible markers supply measured-Te conditioning across 0.5--20 mTorr. The article supplies no Te uncertainty here, TCP power is not absorbed power, and the assumed Maxwellian EEDF is not independently validated by this figure. |
| M10 | The paper uses a 21.5 cm chamber radius, an 11 cm active-plasma gap giving 16.0 L, a 43,000 cm3 chamber volume, and effective `V/A` lengths of 3.6 cm (11 cm gap) and 2.5 cm (6.5 cm gap). | Requires distinct active-plasma and neutral-control volumes. The exact cylinders reproduce both quoted `V/A` values; a one-volume residence/source ledger is not a faithful Lam reproduction. |
| M11 | Figure 11 reports volume-average electron densities derived from Langmuir-probe analysis. Measurements along a line 1.35 cm above the wafer/chuck were converted assuming radial symmetry and an axial `sin(pi*h/gap)` distribution. Reducing the gap decreases average electron density by about a factor of two. | The 27 resolved markers condition/grade the volume-average electron state. They are not local sheath-edge density or wafer flux; the article reports no density uncertainty and points to the probe analysis elsewhere. |
| M12 | Footnote 14 uses `ne = 1e11 cm^-3`, `kd = 7e-9 cm3/s`, and therefore `kd*ne = 700 s^-1`. But the article's printed `kdis = 4.52e-8 exp(-7.40/Te)` law, evaluated at the highest measured 11 cm/10 mTorr Figure-3 Te and augmented by the article's maximum stated 1/7 attachment contribution, gives at most `1.13e-9 cm3/s` or `113 s^-1`. | The footnote arithmetic is internally consistent but its nominal `kd` is 6.18x above the maximum supported by the article's own printed rate law and measured Te board. The footnote rate and `200x` residence-time statement are quarantined from calibration; the measured boards and explicit equations remain usable with their stated uncertainty limits. |
| M13 | Footnote 16 prints the complete Chapman--Enskog Cl-in-Cl2 diffusion equation, masses 35/70 g/mol, `sigma_Cl=3.548 A`, `sigma_Cl2=4.115 A`, `epsilon_Cl/k=75 K`, `epsilon_Cl2/k=357 K`, and a final factor of 1.25 tied to a reported 0.15 cm2/s room-temperature, one-atmosphere measurement. | Reconstructs temperature-dependent diffusion without importing the temperature-conflicted Economou constant. The factor 1.25 is source-declared and is never retuned to reactor or depth data. Neufeld's later evaluated collision integral explicitly upgrades the Hirschfelder table cited by Malyshev, so the provider is source-parameterized rather than a source-identical replay. Missing physical uncertainty keeps it non-predictive. |
| M14 | The gas is initially at 333 K, while the article explicitly says gas temperature increases with power and publishes no powered gas-temperature board. | A 333 K transport replay is an initial-condition sensitivity, not a condition-resolved prediction or temperature extrapolation of the Stafford wall data. |
| M15 | Eq. 11 uses `n0_Cl2 = n_Cl2 + n_Cl/2` and says the in-reactor pressure/number density rises even though wall recombination in the gauge line leaves the pressure reading unchanged. | For each measured relative-Cl2 row, the transport diagnostic combines Eq. 11 with the reported Cl2/rare-gas flows: `p_particles/p_gauge = 1 + x_Cl2,0(1-relative_Cl2)`, with `x_Cl2,0=0.95`. It does not dissociate the rare-gas inventory or equate gauge pressure with powered particle density. |

## Use decision

The 62 clean Figure-3 markers condition measured electron temperature, the 27
clean Figure-11 markers condition volume-averaged electron density, and the 38
clean Figure 7--8 markers are the first Lam-equipment dissociation board for
the native chlorine reactor. Together they separate observed electron state
from the dissociation test. They do not measure absorbed power, species-
resolved wafer flux, ion energy, or etched depth.

The reactor model may either calibrate on a preregistered subset and grade the
held-out conditions, or use independently measured boundaries and grade all
markers. Feature depth may never select the reactor parameters used here.

The runtime Figure-3 provider is hash-locked to the audited CSV. It permits
only an exact marker or an explicitly labeled linear interpolation between two
unambiguous markers in one fixed-gap, fixed-pressure series. It refuses
pressure/gap extrapolation and overlapping marker clusters.

The measured-state Eq.-7 inversion is frozen under
`results/curated/reactor_global_chlorine/`. It replaces the paper's legacy
neutral-dissociation fit with Hamilton's eight-state rate, retains the
published attachment channel, and solves only for the wall-return frequency
required by each supported marker.

The separate neutral-transport diagnostic reconstructs a source-parameterized
temperature-dependent Cl/Cl2 diffusivity with the Neufeld collision integral,
applies Eq. 11's particle-density multiplier, and inverts the exact cylindrical
Robin mode. It reports only a model-conditioned effective `gamma` at an
explicitly declared gas temperature. Because the powered gas temperature and
distributed wall state are unmeasured, it does not supply a predictive local
wall law, wafer flux, or feature depth.
