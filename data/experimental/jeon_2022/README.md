# Jeon et al. 2022 SiO2 trench-depth transfer data

Primary source: W.-n. Jeon et al., “Investigation into SiO2 Etching Characteristics Using
Fluorocarbon Capacitively Coupled Plasmas: Etching with Radical/Ion Flux-Controlled,”
*Nanomaterials* **12**, 4457 (2022), DOI
[10.3390/nano12244457](https://doi.org/10.3390/nano12244457), PMCID
[PMC9781520](https://pmc.ncbi.nlm.nih.gov/articles/PMC9781520/).

The article and figures are CC BY 4.0. The CSV is a transcription of the plotted experimental marker
centres in Figures 4b, 7b, and 9b; it is not an author-supplied raw table. Published error bars are
visible, but their statistical meaning is not stated in the article, so they are not silently converted
to standard deviations. `digitization_uncertainty_nm=35` covers roughly six vertical pixels (marker
thickness, JPEG compression, and manual centre selection) and is separate from the unspecified
measurement error bars.

## Reproducible source artifacts

The Europe PMC full-text XML was downloaded from
`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9781520/fullTextXML`; SHA-256:
`4402d3ccf3ad876bba2a28a08ec6ef1dc0744162defdb8ac7f5886d1a88a1e6d`.

The plotted source images can be downloaded from the following CC-BY CDN URLs:

| Figure | URL | SHA-256 |
|---|---|---|
| 4 | `https://cdn.ncbi.nlm.nih.gov/pmc/blobs/70a5/9781520/8664a589e5c1/nanomaterials-12-04457-g004.jpg` | `d211901e102023a9a320d38fc7ff02dc89798d56e93e6f7a94501dd225f5ff4c` |
| 7 | `https://cdn.ncbi.nlm.nih.gov/pmc/blobs/70a5/9781520/e8e33a0ab345/nanomaterials-12-04457-g007.jpg` | `4cd8712f248eca6a91a57e660d5bd7b0fd3ea3e023d03f8b01d9cd06c62e7855` |
| 9 | `https://cdn.ncbi.nlm.nih.gov/pmc/blobs/70a5/9781520/ed8b0c340ed0/nanomaterials-12-04457-g009.jpg` | `df1733a5e03d50201d1fee80309342189e3b57afc9108a0f7f4340ca8651215f` |

The physical-control ratios in `digitized_plasma_controls.csv` come from Figures 3c, 6c, and 8c:

| Figure | URL | SHA-256 |
|---|---|---|
| 3 | `https://cdn.ncbi.nlm.nih.gov/pmc/blobs/70a5/9781520/6d5f58d93bda/nanomaterials-12-04457-g003.jpg` | `3781dd55d34ac737c46ebd61725c229813cfbc07615faf4523117c4c1021850a` |
| 6 | `https://cdn.ncbi.nlm.nih.gov/pmc/blobs/70a5/9781520/ac8fd03386dd/nanomaterials-12-04457-g006.jpg` | `7992bcc9f8fa7ea26b34640861e71de290f385677fc18474cc0dae6fdc7602ae` |
| 8 | `https://cdn.ncbi.nlm.nih.gov/pmc/blobs/70a5/9781520/2e501630da98/nanomaterials-12-04457-g008.jpg` | `bd69da4595a144bb936d58b2b71c8daed9c6a92606056e53a5c60286400e1f1d` |

Those neutral-to-ion ratios are **diagnostic-derived**, not directly measured fluxes. Jeon et al.
integrated measured radical densities (excluding stable Ar and C2F4), assumed 300 K radical thermal
flux, estimated Bohm ion flux from measured electron density using 3 eV electrons and Ar ion mass, and
used a simplified on/off duty model for pulsed operation. Published error bars are visible but their
statistical semantics are not stated. The CSV therefore retains `not_specified` rather than inventing a
standard deviation.

`digitized_electron_bias_controls.csv` independently replays the black electron-density markers and
blue self-bias markers in Figures 3b, 6b, and 8b. It stores each marker's y pixel, the linear or log-axis
map, and conservative digitization bounds (4e14 m^-3 and 20 V). Electron density supports replay of the
paper's Bohm-flux estimate using its assumed 3 eV electron temperature and Ar ion mass. Self-bias is an
energy-scale diagnostic only: it is not mislabeled as a measured IEDF, and the ion composition remains
unmeasured.

The vertical pixel-to-depth maps were least-squares fits to the four labeled y-axis ticks
(400, 800, 1200, 1600 nm) in each source image. Every CSV row retains its marker-centre pixel,
slope, and intercept so the reported depth can be replayed rather than trusted as an opaque number.

## Experimental conditions and scope

- 2.4 um PECVD SiO2 beneath an amorphous-carbon mask; trench widths 60–200 nm.
- Ar/C4F8 capacitively coupled plasma; 13.56 MHz; 300 W; 20 mTorr; total flow 100 sccm;
  electrode temperature 10 degC.
- Continuous-wave C4F8 fractions are varied in Figure 4. Figures 7 and 9 fix pulse-on time at 1 ms
  and vary pulse-off time at 20% and 80% C4F8, respectively.
- The paper does not report etch duration, a measured IEDF, or species-resolved surface fluxes.
  Consequently these observations constrain relative depth/profile transfer, not absolute etch rate
  from a complete physical boundary state.
- Depth-versus-width shape ratios remain usable without inventing the missing duration because every
  width within one plotted condition shares one coupon exposure.  Pulsed/CW depth ratios do **not**
  enjoy that cancellation: the paper does not report whether conditions shared wall-clock duration or
  cumulative RF-on duration.  At 1 ms on / 1 ms off those hypotheses differ by a factor of two in wall
  time and neutral fluence.  The related experiment from the same group (Cho et al., *Materials* 14,
  5036 (2021), DOI 10.3390/ma14175036) used 10 min for CW and 20 min for a 50%-duty pulse, demonstrating
  why “same etch conditions” cannot establish Jeon's exposure protocol.  It is context, not evidence
  that Jeon used the same schedule.
- The loader therefore omits cross-exposure pulse/CW targets by default.  A development analysis may
  request `common_wall_time` or `common_rf_on_time` explicitly; every resulting target records that the
  hypothesis was not reported by the source.  Neither hypothesis supports a validation claim without
  author data or another primary record.  All intervals propagate digitization budgets only and are not
  mislabeled as total experimental uncertainty.

## Preregistered split

Only Figure 4b at 20% C4F8 continuous wave is marked `calibration`. The 40% and 80% gas-fraction
curves remain held-out transfer observations.  The raw pulse rows retain their originally registered
labels, but their cross-condition use is blocked on the unresolved exposure protocol.  Within-condition
width shapes remain usable development evidence.  The observed 1 ms regime reversal—greater depth at
20% and smaller depth at 80%—must not be scored as a predictive pulse/CW ratio until exposure is known.
