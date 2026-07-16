# Jeong et al. 2023 fixed-duration SiO2 energy/flux controls

Primary source: W. Jeong et al., “Contribution of Ion Energy and Flux on High-Aspect-Ratio
SiO2 Etching Characteristics in a Dual-Frequency Capacitively Coupled Ar/C4F8 Plasma:
Individual Ion Energy and Flux Controlled,” *Materials* **16**, 3820 (2023), DOI
[10.3390/ma16103820](https://doi.org/10.3390/ma16103820), PMCID
[PMC10222222](https://pmc.ncbi.nlm.nih.gov/articles/PMC10222222/).

Unlike the 2022 Jeon/Jeong pulse study, this article explicitly states that every scored coupon was
etched for 20 minutes.  Figure 7 is therefore suitable for absolute depth comparison.  It separates
an ion-energy sweep at approximately fixed electron density from an electron-density/ion-flux sweep
at approximately fixed self-bias, and reports 200, 100, and 60 nm trenches.  The 60 nm series is a
held-out charging-sensitive test: increasing energy does not release the approximately 500 nm etch
stop, while increasing ion flux and the accompanying heavy-radical population relax it slightly.

`digitized_figure7_depths.csv` records the colored marker centres and both linear axis maps.  A
conservative 35 nm digitization interval covers marker thickness, JPEG compression, and several
source pixels.  The publication does not state the statistical meaning of experimental uncertainty,
so that interval is never presented as total measurement uncertainty.  Exactly one point—the 200 nm,
890 V energy-control reference—is designated as a magnitude calibration anchor.  The remaining 17
points are frozen held-out transfer observations.

`digitized_figure6_radicals.csv` is deliberately a different evidence class.  Figure 6 contains
radical densities calculated by the authors' in-house volume-averaged plasma model, not direct
measurements.  These values may provide boundary-state development inputs, but may not be scored as
experimental validation.  The heavy (`C4F7`, `C3F6`, `C2F4`) and light (`CF3`, `CF2`, `CF`) labels
follow the article's mechanism discussion.

## Reproducible sources

- Europe PMC XML: `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10222222/fullTextXML`
  (SHA-256 `249045f4e77a47fd4e01fe77e7beb05b413d3468051f780600d4a8dbef86507c`).
- Figure 6: `https://cdn.ncbi.nlm.nih.gov/pmc/blobs/cdcd/10222222/25483ca8cce4/materials-16-03820-g006.jpg`
  (SHA-256 `3e4ea56418343dbc13bf3109e778a852181f2a473d169d35ad5ccfef0baf6d53`).
- Figure 7: `https://cdn.ncbi.nlm.nih.gov/pmc/blobs/cdcd/10222222/884a7130302a/materials-16-03820-g007.jpg`
  (SHA-256 `9c8acd0e9a7219ea5f99e097f7977fcb3ed490635fca3a0f9f69cb0a15a6508c`).

The control labels use the values stated in the Figure 4/5 captions: energy-control self-biases of
450/890/1270 V at electron density about 2.0e15 m^-3, and flux-control electron densities of
1.1/1.9/3.1e15 m^-3 at self-bias about 740 V.  Self-bias is an energy-scale proxy, not a measured
IEDF; electron density supports an ion-flux closure, not a direct ion-flux measurement.
