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
