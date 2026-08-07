# Yoshie et al. 2023 cyclic SF6/C4F8 silicon depth board

Primary source: T. Yoshie et al., “Bias-supply timing tailored to the
aspect ratio dependence of silicon trench etching in Ar plasma with
alternately injected C4F8 and SF6,” *Applied Surface Science* **638**,
157981 (2023), DOI
[10.1016/j.apsusc.2023.157981](https://doi.org/10.1016/j.apsusc.2023.157981).
The article is CC BY 4.0.

This dataset was machine-digitized only after the target IDs, complete
timing/width grid, source URLs/checksums, and calibration/held-out split were
frozen by the `depth-cross-chemistry-v1` preregistration (commit
`dec49c464af246b04f0168588dc436a76442e64c0294248b324ceb96cf9d0c62`).
The figures had been visually inspected to choose panels and establish
scientific relevance; this is value-blind with respect to machine
digitization, not a claim that no human had ever seen the plots.

## Experimental meaning

The reactor used continuous Ar and alternating 1 s C4F8 / 1 s SF6 injections.
The substrate bias was on for 0.25 s per cycle.  Figure 5 reports 675 cycles
in 45 min; Figure 6 reports 450 cycles in 60 min.  The paper defines the
plotted rate as depth divided by cumulative bias-on time, so absolute depths
are:

```text
Figure 5 depth = plotted rate * 168.75 s / 60 s/min
Figure 6 depth = plotted rate * 112.5 s / 60 s/min
```

`figure4_blanket_poly_si_rates.csv` contains the same-reactor blanket poly-Si
rates at the selected bias timings.  Those seven values may condition a
facility-specific boundary before feature values are used.
`figures5_6_feature_depths.csv` contains every preregistered feature marker:
three timings by seven widths for the 4 s cycle and four timings by seven
widths for the 8 s cycle.  All 49 are held-out transfer targets.

This creates a possible Tier-B test—independent blanket conditioning followed
by feature depth prediction—but only for a mechanism that resolves the
experimental transfer.  The blanket is a 370 nm poly-Si film exposed for 75
cycles; the patterned target is bulk Si exposed for 450 or 675 cycles.  A
read-only audit (`scripts/audit_yoshie_2023_blanket_transfer.py`) shows that a
single multiplicative blanket scale is invalid: for the 8 s timing-I condition,
every feature rate is more than 2.5 times the blanket rate even after the
digitization allowances, and the timing rank changes.  The blanket observable
may therefore condition a reactor boundary only inside a material- and
cycle-history-resolved mechanism.

The source does not report species-resolved wafer fluxes or a measured IEAD.
The plotted error bars are retained, but their statistical semantics are not
stated, so they are not silently called standard deviations or confidence
intervals.

## Pixel replay and visual audit

Download the official publisher rasters without renaming them:

```bash
mkdir -p tmp/sources/yoshie_2023
curl -L -o tmp/sources/yoshie_2023/yoshie_fig4.jpg \
  https://ars.els-cdn.com/content/image/1-s2.0-S0169433223016604-gr4_lrg.jpg
curl -L -o tmp/sources/yoshie_2023/yoshie_fig5.jpg \
  https://ars.els-cdn.com/content/image/1-s2.0-S0169433223016604-gr5_lrg.jpg
curl -L -o tmp/sources/yoshie_2023/yoshie_fig6.jpg \
  https://ars.els-cdn.com/content/image/1-s2.0-S0169433223016604-gr6_lrg.jpg
```

Then verify checksums, dimensions, axis pixels, the two committed CSVs, and
the manifest, and optionally draw QA overlays:

```bash
python scripts/digitize_yoshie_2023_figures4_6.py \
  --source-dir tmp/sources/yoshie_2023 \
  --overlay-dir tmp/yoshie_2023_overlays
```

The Figure-4 y transforms reproduce the paper’s printed maxima: 464.959
versus 466 nm/bias-min for the 4 s cycle, and 592.683 versus approximately
591 nm/bias-min for the 8 s cycle.  Full-resolution color-component centers,
error-cap pixels, and a conservative digitization allowance are retained in
the committed tables and manifest.

The initial-pattern SEMs in Figures 5(a) and 6(a) were also inspected at full
resolution with a pixel grid against their 500 nm scale bars.  The apparent
vertical openings are the trenches through the SiO2 mask; their bottoms
coincide with the Si/mask interface to within the raster resolution.  There is
no hidden initial Si depth large enough to explain the super-blanket 8 s
timing-I rates.
