# tian-2018-vuv-thesis

**Spectrum- and trapping-resolved ion/VUV equipment-model benchmark**

- **Citation:** P. Tian, *Controlling Photon and Ion Fluxes in Low Pressure
  Low Temperature Plasmas*, PhD dissertation, University of Michigan (2018).
- **Official full text:** `https://cpseg.eecs.umich.edu/pub/theses/tian_peng_phd_thesis.pdf`
- **Local extraction:** `research_sources/thesis_extracts/tian_peng_phd_thesis.txt`
- **Status:** PRIMARY FULL THESIS READ; FIGURES 5.3, 5.6, 5.10(a,c), AND
  5.12(c) QUANTITATIVELY/VISUALLY PIL-AUDITED
- **Digitized data:**
  `data/experimental/tian_2017_vuv_figures/digitized_figures_5_10_5_12.csv`
- **Reproduction script:** `scripts/digitize_tian_2017_vuv_figures.py`
- **Spatial-moment script:**
  `scripts/digitize_tian_2017_base_spatial_moments.py`

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| T1 | The HPEM/radiation calculation resolves Ar 104.8/106.7 nm and Cl 139 nm emission, line trapping, wall state, pulse dynamics, and wafer angular flux. | Spectrum/trapping benchmark for a reduced deterministic model, not a production Monte Carlo prescription. |
| T2 | For the 22.5-cm diameter, 12-cm high, 150 W, 20 mTorr, 200-sccm case, Ar/Cl2=80/20 gives about 6.2e14 VUV and 4.5e15 ion cm^-2 s^-1; the 95%-Cl2 endpoint is about 2.3e13 VUV and 1.3e15 ion cm^-2 s^-1. | Independent magnitude and mixture-trend comparator; it is simulated, not a Lam measurement. |
| T3 | Pressure, mixture, wall recombination, and pulsing move the VUV/ion ratio by orders of magnitude. | Proves that one fixed beta cannot bridge knobs to depth. |
| T4 | Atomic-Cl 139-nm emission uses explicit excitation/radiative kinetics and optical trapping. | Kemaneci reaction-18 emissivity without escape/trapping is only a transparent upper sensitivity. |
| T5 | Figure 5.10(a,c) resolves the substrate 104.8, 106.7, and 139 nm fluxes, their three-line total, total photon/ion ratio, and 139-nm spectral fraction over the full 5--95% Cl2 sweep. | The 300-dpi Poppler render is SHA-256 pinned and 64 visible markers are replayed from PIL/NumPy-audited pixel coordinates. These are source model outputs, not measurements. |
| T6 | Figure 5.12(c) predicts strong, composition-dependent resonance trapping: the digitized 106.7-nm factor falls from about 534 to 12, the 104.8-nm factor from about 230 to 51, while the actual 139-nm factor (the printed curve is multiplied by ten) lies around 4--11. | Quantitative regression gate for deterministic escape/trapping. It cannot validate atomic branching, surface photon yield, or feature depth. |
| T7 | In the base argon case the text gives substrate fluxes near 3.7e14 cm^-2 s^-1 at 106.7 nm and 1.9e14 cm^-2 s^-1 at 104.8 nm, with trapping factors about 368 and 216. | The deterministic zonal partial-frequency-redistribution audit reproduces the two trapping factors with 3.46% combined MAPE and correct ordering. This is a source-model regression, not an experimental validation. |
| T8 | Figure 5.3 resolves the Ar(1s4) field associated with 106.7 nm; the corresponding Ar(1s2) spatial field for 104.8 nm is not printed. | Three source-derived spatial zones are digitized for 106.7 nm. Reusing that field for 104.8 nm remains an explicit sensitivity and prevents a formal two-line spatial validation claim. |

## Digitization boundary

The committed 97-marker table uses only visible square markers; plotted spline
segments are not recast as data.  Four pixels is carried as placement
uncertainty.  The source PDF, rendered pages, axes, and marker neighborhoods are
checksum/vision checked by the reproduction script.  Because the thesis curves
come from HPEM, their correct role is an independent implementation regression
target.  Treating them as reactor measurement or as a depth calibration is
explicitly forbidden by the manifest.

The separate base-case spatial manifest contains source-derived zone moments,
not raw plasma measurements or pixel-perfect field values.  Its page hashes
and QA overlays are retained, and the missing 104.8-nm upper-state field is a
named source limitation rather than silently copied as data.
