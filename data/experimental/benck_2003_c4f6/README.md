# Benck 2003 C4F6/Ar mass-resolved ion-current board

This directory digitizes the five unambiguous reactor-boundary series in
Figures 9 and 10 of Benck, Goyette, and Wang, *Journal of Applied Physics*
**94**, 1382–1389 (2003), DOI `10.1063/1.1586978`: total positive-ion current
and the Ar+, CF+, CF2+, and CF3+ components. Figure 9 sweeps 25–100% C4F6 at
1.33 Pa; Figure 10 sweeps 0.67–2.66 Pa at 50% C4F6.

The source pixels are not redistributed because the archived scan carries an
AIP redistribution notice. The numerical table retains every full-page
600-dpi marker center, the three logarithmic-axis calibration ticks, the PDF
and render hashes, a conservative 10.1% digitization bound, and the source's
separate 20% corrected-transmission uncertainty.

With an authorized local source copy, reproduce the render and visual audit:

```bash
pdftoppm -f 6 -l 6 -singlefile -png -r 600 \
  /path/to/benck_2003_c4f6.pdf /tmp/benck_page6
python scripts/digitize_benck_2003_figure9.py \
  --source-pdf /path/to/benck_2003_c4f6.pdf \
  --render /tmp/benck_page6.png \
  --overlay /tmp/benck_figure9_overlay.png
python scripts/digitize_benck_2003_figure10.py \
  --source-pdf /path/to/benck_2003_c4f6.pdf \
  --render /tmp/benck_page6.png \
  --overlay /tmp/benck_figure10_overlay.png
```

Both overlays were inspected at original resolution. Overlapping C/F/Cx and
SiFx/COFx glyphs were excluded instead of guessed. The retained curves agree
with the source text: total and Ar+ current fall as C4F6 replaces Ar; CF+ and
CF2+ rise nearly proportionally through 75% and then level; CF3+ continues to
rise toward pure C4F6. At 50/50 C4F6/Ar, total, CF+, and CF2+ current decrease
with increasing pressure while CF3+ is nonmonotonic.

The independently drawn Figure-9 50%-C4F6 column and Figure-10 1.33-Pa column
represent the same condition. Their five retained currents reconcile within
3.53%, below either panel's conservative 10.2% digitization bound. This checks
the pixel extraction, not a second physical experiment, so it is not counted
as an independent validation case.

Together these are the strongest quantitative ion-side targets in the library
for a C4F6 reactor provider. It is not transplanted into Krüger: Benck used a
200 W ICP with a grounded surface, no oxygen, and a different sampling
location, while Krüger used a high-power dual-frequency biased CCP. Neither
figure measures stable neutral C4F6. A model can earn held-out composition and
pressure-response grades here without gaining authority to normalize
Krüger's 825 nm endpoint.
