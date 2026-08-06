# Takada 2005 molecule–ion co-incidence yields

This directory contains an auditable digitization of the SiO2 series in
Figure 3 of N. Takada, H. Toyoda, and H. Sugai, “Evidence of Radical-free
Etching of SiO2 by Fluorocarbon Molecule under Ion Bombardment,” *Transactions
of the Materials Research Society of Japan* **30**[1], 319–322 (2005). The
closely related journal article is N. Takada *et al.*, *Journal of Applied
Physics* **97**, 013534 (2005), DOI `10.1063/1.1829400`.

The source is archived and text-extracted at:

- `research_sources/takada_2005_radical_free_etching.pdf`
- `research_sources/thesis_extracts/takada_2005_radical_free_etching.txt`

## Reproduction and visual audit

Render source PDF page 3 at 600 dpi:

```bash
mkdir -p tmp/pdfs/takada_2005
pdftoppm -f 3 -l 3 -r 600 -png -singlefile \
  research_sources/takada_2005_radical_free_etching.pdf \
  tmp/pdfs/takada_2005/page3_600dpi
```

Then replay the PIL/NumPy axis audit, verify both source checksums and the
committed files, and produce a QA overlay:

```bash
python3 scripts/digitize_takada_2005_fig3.py \
  --overlay tmp/pdfs/takada_2005/figure3_overlay.png
```

The Figure-3 x axis is logarithmic. The CSV therefore retains both each
experimental flux-ratio setpoint and the marker-center pixel used as an
independent placement check. The red C5F8 and blue CF2 overlay centers were
visually reconciled against the 600-dpi source after dark-pixel axis
localization. A conservative ±0.015 yield digitization allowance corresponds
to eight vertical pixels and is wider than the scanned marker edge.

## What the experiment establishes

At 400 eV Ar+ energy, the digitized C5F8/Ar+ co-incidence series is
non-monotone:

```text
C5F8/Ar+ ratio     SiO2/Ar+ yield
0.25               0.6697
0.50               1.0162
1.00               1.1969
2.50               1.0591
10.0               0.7945
```

The source text independently rounds the ratio-0.25 and ratio-1 yields to
`0.67` and `~1.2`. It also reports `2.5` at 900 eV and ratio 1, outside
Figure 3. The related *Journal of Applied Physics* article's publisher
abstract reports `2.4` for that nominal condition, so downstream audits retain
the source discrepancy as `2.4–2.5`. Thus a stable fluorocarbon parent
molecule can participate directly in ion-assisted oxide removal, while
excessive molecule/ion ratio turns the balance back toward fluorocarbon-layer
formation.

This is not a C4F6 law. The measured molecule is C5F8, whose ring and double
bond are identified by the authors as a potentially unique chemisorption
route. The table is therefore an analog mechanism constraint and a proof that
`1.5 SiO2/ion` is not a universal plasma-surface ceiling—not a Krüger boundary
condition, not a feature prediction, and not permission to transplant the
numeric curve across molecular identity or energy.
