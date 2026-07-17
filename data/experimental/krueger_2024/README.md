# Krüger et al. 2024 — SiO2 calibration and transfer evidence

Source: Krüger et al., *Autonomous hybrid optimization of a SiO2 plasma etching mechanism*,
J. Vac. Sci. Technol. A 42, 043008 (2024), DOI
[10.1116/6.0003554](https://doi.org/10.1116/6.0003554). The authors' manuscript is available from
the [University of Michigan Computational Plasma Science and Engineering Group](https://cpseg.eecs.umich.edu/pub/articles/JVSTA_42_043008_2024.pdf).
Source checked 2026-07-11.

The article is published under an exclusive license to AVS. This directory redistributes no figures,
SEM pixels, article text, or supplementary mechanism. It contains only a small, attributed transcription
of numerical facts and qualitative findings needed to define reproducible validation targets.

## Evidence contract

- `base_case_metrics.csv` contains the six experimental SEM-derived targets in Table IV. These are the
  **calibration** observations. The source used one base-case SEM after a 60 s C4F6/Ar/O2 etch of SiO2
  through an amorphous-carbon mask.
- `base_case_boundary_fluxes.csv` contains Table I wafer fluxes. They are **HPEM simulation outputs**, not
  measurements. The paper explicitly treats them as ground truth for its proof of concept while noting
  their model uncertainty. petch must preserve that distinction.
- `transfer_observations.csv` contains only claims stated in Sections IX.A–B. Experimental trends are
  **held out** from calibration. Three MCFPM etch depths are retained only as `reference_only` values and
  can never be scored as experimental validation.
- `digitized_figure4_iead.csv` is a probability-weighted numerical projection of the combined positive-ion
  energy/angular distribution in Figure 4(a). `digitized_figure4_iead_metadata.json` records the source
  checksum, pixel calibration, logarithmic color mapping, binning error, and digitization sensitivity.
  The copyrighted raster is not redistributed.

The profile images themselves are not numerically recoverable from the paper without digitization. These
tables therefore gate scalar/trend validation only. Pixel- or contour-level matching requires source SEMs,
an explicit digitization uncertainty, or new measurements; it must not silently use the printed figures.

The paper's Table V probabilities are **experiment-calibrated mechanism inputs**, not first-principles
constants. The four-feature optimization reports 0.0909 for bare-SiO2 sputtering, 0.1384 for
SiO2-fluorocarbon-complex sputtering, 0.2729 for complex formation, 0.0628 for O-based polymer etching,
and 0.0842 for polymer deposition on the amorphous-carbon mask. Petch now exposes these values through
`krueger_2024_reduced_projection()` factories for oxide and mask, but only as a **development replay**.
Krüger's source mechanism additionally resolves species-specific complexes, SiO2 redeposition, excited
sites, polymer identities, crosslinking/bond breaking, and a species-resolved ion/hot-neutral mixture.
Those omissions and the assumed amorphous-carbon density remain in machine-readable provenance; the
projection cannot earn a predictive claim merely by reproducing the calibration case.

The Figure 4(a) IEAD is HPEM output, not a measured wafer distribution. Its digitized mean energy is
approximately 3.465 keV, its central 90% interval is approximately 0.98--4.62 keV, and the signed-angle
standard deviation is approximately 0.84 degrees. Reasonable pixel/colorbar perturbations move the mean
by less than 15 eV peak-to-peak. The table retains the joint energy-angle correlation; it is not replaced
by a monoenergetic or factorized closure.

## Base process context

- CCP pressure: 10 mTorr.
- C4F6/Ar/O2 flows: 140/100/105 sccm.
- Low-/high-frequency powers: 8.0/2.5 kW at 1/40 MHz.
- Initial mask thickness/opening: 850/90 nm; feature grid: 1 nm; etch time: 60 s.
- Transfer cases vary O2/C4F6 at 6 kW low-frequency power or vary low-frequency power at otherwise
  related conditions. They are intentionally outside the single-SEM calibration case.

Checksums are verified by `load_krueger_2024_evidence`; update them only after independently checking a
source correction:

- `base_case_metrics.csv`: `5d51d124a93e1f942a9b999649b8adcf217662967d9ea2a40089f72940992351`
- `base_case_boundary_fluxes.csv`: `ad50b6099a52d2c2cc00eb4eade496b9d75c41d19881c5fec9e905f9dfd3808b`
- `transfer_observations.csv`: `85cef607f20ab5e56e606666aa7e0e6241abb546d0369277b21833542e04d425`
- source author-manuscript PDF: `65b7750b2b773c3725d8f09f778b5b728ce9974a4548a5d522d19256f6bf9a51`
- `digitized_figure4_iead.csv`: `783f7084b5ba6dc71eb89efeb70871cabdb7378f87db44f8fa9dcb1f3adb6ce4`
- `digitized_figure4_iead_metadata.json`: `7904a700afdcd116c6f57ef35aeb5555661ffed0d004304d815d20f79840ca55`
