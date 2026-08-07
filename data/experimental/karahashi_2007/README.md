# Karahashi 2007 mass-selected reactive-ion beam yields

This directory contains an auditable digitization of Figure 4 in K. Karahashi,
*Surface Science* (Japan) **28**, 60–65 (2007). The review reports the
mass-selected SiO2 beam experiments from Karahashi *et al.*, *Journal of Vacuum
Science & Technology A* **22**, 1166 (2004), DOI `10.1116/1.1761119`.

The experiment isolates surface physics unusually cleanly: energy-controlled
single-species ions impinge on SiO2 under ultrahigh vacuum, with no gas-phase
reactions and no incident neutral-radical flux. Figure 4 therefore tests
reactive-ion identity and energy directly. It does **not** represent a
fluorocarbon-plasma ion mixture or ion–neutral co-incidence.

## Reproduction and visual audit

The source PDF is archived at
`research_sources/karahashi_2007_hyomen_kagaku_28_60.pdf`. Render PDF page 4
with Poppler at 600 dpi:

```bash
mkdir -p tmp/pdfs/karahashi_figures
pdftoppm -f 4 -l 4 -r 600 -png -singlefile \
  research_sources/karahashi_2007_hyomen_kagaku_28_60.pdf \
  tmp/pdfs/karahashi_figures/page4_600dpi
```

Then verify the source, page-render checksum, axis pixels, committed table, and
manifest, and make a QA overlay:

```bash
python3 scripts/digitize_karahashi_2007_fig4.py \
  --overlay tmp/pdfs/karahashi_figures/fig4_digitization_overlay.png
```

`figure4_reactive_ion_yields.csv` retains every numerical yield together with
the raw marker-center and error-cap pixels. `digitization_manifest.json`
records the checksums and linear axis transform. Marker centers and caps were
located at full 600-dpi resolution, cross-checked with dark-pixel contours,
and reconciled on the overlay. A conservative ±0.011 yield digitization
allowance covers 3.5 vertical pixels, which is wider than the line thickness.

The source plots error bars but the accessible text does not define their
statistical meaning. They are consequently named `plotted_lower_yield` and
`plotted_upper_yield`, not confidence bounds.

## Scientific result and limits

At 1000 eV, the digitized CF3+ yield is `1.4703 SiO2/ion`, consistent with the
text's rounded value of `1.5`. It is not a universal ceiling: the same series
reaches `1.8736` at 1500 eV and is `1.7549` at 2000 eV. The full figure also
shows large species dependence at fixed energy. At 1000 eV the digitized
ladder is:

```text
F+    0.3232 SiO2/ion
CF+   0.6751 SiO2/ion
CF2+  1.1957 SiO2/ion
CF3+  1.4703 SiO2/ion
```

That ladder is a required direct-beam validation target for any
species-resolved surface law. Agreement with only the CF3+ point cannot
validate a species-agnostic ion mechanism.

Figure 4 also places filled CF+ and CF2+ markers on the zero-yield baseline at
low energy. The accompanying text says that below 500 eV yield drops abruptly
and fluorocarbon film grows. Those baseline points are not entered as exact
numeric zeros because the plot cannot distinguish zero etch from negative net
removal clipped at the axis. They remain a qualitative deposition/etch-stop
constraint.

This evidence supports interpolation only within each species' positive-yield
measurement range. It does not support extrapolation beyond 2000 eV, does not
bound molecule-assisted plasma etching, and cannot identify the unpublished
ion-species mixture in the Krüger reactor.

## Figure 10 product branching

Figure 10 is independently replayed from PDF page 6 at 600 dpi:

```bash
python3 scripts/digitize_karahashi_2007_fig10.py \
  --overlay tmp/pdfs/karahashi_figures/fig10_product_overlay.png
```

The nine-point table resolves the fractions of detected SiF, SiF2, and SiF4
for CF3+ at 500, 1000, and 2000 eV. SiF2 dominates at or below 1000 eV,
whereas the SiF fraction rises and the SiF4 fraction falls with energy. The
source text identifies prompt collision-cascade ejection for SiF and a delayed
precursor route for SiF2/SiF4 (about 0.5 ms delay for SiF4).

The ordinate is normalized only over detected SiFx products. It supplies a
product-identity/branching gate, not an absolute product yield, a
prompt-versus-delayed numerical partition, or a diffusion/escape length.
The source does not report the ion incidence angle for Figure 10. Because the
same source shows that the product fractions change with incidence angle,
Figure 10 cannot be condition-matched to the explicitly normal-incidence
Figure-4 yield board without an additional assumption.

## Figure 9 product timing

Figure 9 is independently replayed from PDF page 5 at 600 dpi:

```bash
python3 scripts/digitize_karahashi_2007_fig9.py \
  --overlay tmp/pdfs/karahashi_figures/fig9_product_peak_overlay.png
```

The source condition is a mass-selected 1000 eV CF3+ beam at 30 degrees from
normal, pulsed with 100 microsecond FWHM. The highest measured marker in each
arbitrary-unit panel occurs at `0.2391 ms` for SiF, `0.2993 ms` for SiF2, and
`0.5403 ms` for SiF4, each with a `0.006 ms` digitization allowance. The
original-resolution PIL overlay was inspected with all three crosshairs
centered on the intended source markers.

This is a timing/topology constraint, not a transport parameter. The measured
clock convolves the beam pulse, product flight time and velocity distribution,
detector geometry, and the surface reaction. The source does not report the
instrument-response deconvolution, a diffusion coefficient, a formation-depth
distribution, an escape probability, or a prompt/delayed branching fraction.
The panel ordinates are unlabeled arbitrary units and are deliberately not
digitized or compared across products. A predictive product-transport model
must carry a conserved delayed inventory, but Figure 9 alone cannot assign its
kinetic coefficients.
