# Hwang–Giapis 1997 source-backed boundary and notch-profile evidence

Primary modeling source: Hwang and Giapis, “On the origin of the notching effect during etching in
uniform high density plasmas,” *Journal of Vacuum Science & Technology B* **15**, 70–87 (1997),
DOI `10.1116/1.589258`.

This directory contains three distinct evidence products:

- `fig4a_ion_energy_distribution.csv`: the ion-energy distribution used as a plasma-to-feature
  boundary condition.
- `fig4b_electron_energy_distribution.csv`: the electron-energy distribution used by the
  source-faithful 2-D charging replay. Its companion electron angular law is the paper's explicit
  Figure 5(b) fit, `cos(theta)^0.6`.
- `fig13_notch_profile.csv`: the open-circle experimental profile shown in Figure 13, transcribed
  from the paper’s comparison with the Nozawa experiment.

The tables are not interchangeable. Figures 4(a), 4(b), and the analytic Figure 5(b) angular fit are
engine inputs; the Figure 13 contour is an experimental score target and is never used to modify the
incoming particles, the etch-yield law, the oxide-scattering law, the 50-collision removal threshold,
or the charging boundary condition.

The source calculation is explicitly two-dimensional. The Hwang--Giapis replay therefore samples the
published 2-D energy/angle laws directly. It does not sample a three-dimensional Maxwellian and fold
the unmodeled out-of-plane energy back into the trajectory plane. Three-dimensional consumers retain
a separately declared half-Maxwellian electron closure; the two representations are not silently
interchanged.

## Figure 4(b) EEDF extraction

`scripts/digitize_hwang_giapis_1997_eedf.py` traces the black Figure 4(b) curve in the 240-dpi source
render and writes 0.25 eV bins. The curve is below the rendering's resolved height above 12 eV, so the
script applies a declared 15-pixel resolution threshold and does not invent an analytic tail. The
result has mean energy `3.59999 eV` and probability mass `0.74533` below `5 eV`. Calibration,
checksums, uncertainty, and the tail policy are recorded in `fig4b_digitization_manifest.json`.

## Figure 13 pixel mapping

The source profile crop is
`tmp/pdfs/hwang_giapis_1997/figure_crops/fig13_profile_source_3x.png` (2160 × 1350 pixels,
SHA-256 `abe6d851259c4b6f27bb8632ddfa86dd7127e586373462509ffa182fcaf82167`).
The crop comes from PDF page 12 / printed page 81. Its physical anchors are:

- original poly-Si sidewall: `x = 756.0 px` → notch depth `0.0 um`
- far poly-Si edge: `x = 1744.0 px` → line width `0.5 um`
- poly-Si top: `y = 604.0 px` → height above oxide `0.3 um`
- oxide interface: `y = 1204.0 px` → height above oxide `0.0 um`

Therefore:

```text
notch_depth_um = 0.5 * (pixel_x - 756.0) / (1744.0 - 756.0)
height_above_oxide_um = 0.3 * (1204.0 - pixel_y) / (1204.0 - 604.0)
```

The 22 rows are the open-circle centers. They were identified from connected dark rings and visually
reconciled against the original crop. Point 7 overlaps the printed simulation curve and was manually
centered at `(875.0, 956.0)`; the uppermost ring overlaps the top interface and was centered on
`y = 604.0`. Raw pixels and both linear transforms remain in every row so the physical coordinates
can be replayed without trusting the transcription.

## Uncertainty and claim boundary

A conservative `0.005 um` digitization bound is used in both axes. This is wider than the typical
marker-center ambiguity and equals one cell of the source’s 5 nm local profile grid.

The paper reports no statistical measurement uncertainty for this contour. Agreement may therefore
be described quantitatively as a source-reproduction or experiment-comparison result, but the contour
alone cannot satisfy a strict uncertainty-based experimental-validation claim. It is also development
evidence rather than held-out evidence because the profile was inspected while implementing the
source-faithful replay. No target value may be used to tune the model.

The maximum digitized undercut (`0.2157 um`) independently agrees with the approximately `0.2144 um`
large-open-area notch depth digitized from the Nozawa primary experiment near the same width ratio.
That cross-check is diagnostic only; it does not create an unreported measurement uncertainty.
