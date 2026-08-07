# Krüger thesis Figure 6.17 boundary visual audit

This board answers one narrow completeness question: did a careful visual
inspection of the thesis reveal the two reactor-boundary quantities missing
from the absolute-depth calculation?

No. At full rendered resolution Figure 6.17 labels nine curves: `CF2`,
`C3F4`, `O`, `C2F3`, `CF`, `Ions`, `CF3`, `CO`, and `C`. It contains no
stable-`C4F6` parent curve, and the positive-ion population remains one
aggregate `Ions` curve rather than a species-resolved flux/IEAD.

The receipt is checksum-bound to the 400 dpi render that was inspected.
`scripts/audit_krueger_thesis_boundary_figure.py` uses PIL to verify the exact
image and its dimensions. The label transcription is a full-resolution human
vision audit. No numerical curve values were digitized because curves for the
missing variables do not exist in the figure.

The figure is also a transfer condition (`P_lf = 6.0 kW`,
`P_hf = 2.5 kW`, varied `O2/C4F6`), not the 8 kW base-case boundary. The
result therefore cannot be used to transplant any plotted flux into the base
case. It confirms underidentification; it does not imply zero parent flux.

Replay against the locally archived source and page render:

```bash
python scripts/audit_krueger_thesis_boundary_figure.py \
  --source-pdf tmp/pdfs/krueger_thesis_2024.pdf \
  --rendered-page /private/tmp/krueger_thesis_pdf228_fig6_17.png
```
