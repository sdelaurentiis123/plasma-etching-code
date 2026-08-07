# Woo 2024 C4F6 ICP board

This directory contains an original-pixel digitization of Figure 4.1 from
Byungjun Woo's 2024 Korea University M.S. thesis,
`woo-2024-c4f6-thesis`.

- `figure4_1_page70_600dpi.png` is the checksum-bound full PDF-page render.
- `figure4_1_patterned_etch_rates.csv` contains five SiO2 and five ACL
  patterned-rate points.
- `digitization_manifest.json` records the source hash, axis transform,
  line-segment fits, body-text reconciliation, reactor conditions, and claim
  boundary.
- `results/curated/woo_2024_c4f6_board/figure4_1_overlay.png` is the
  full-resolution visual receipt.

The plotted rates are useful absolute patterned-rate observables for a future
CF4/C4F6/He reactor-to-feature model. They are not a species-resolved kinetic
boundary. The thesis supplies aggregate ion-current density, electron
temperature, self-bias, and relative OES, but no species-resolved ion flux,
IEAD, or absolute neutral flux.

The Figure 4.6 SEM sweep is not a blind depth board: the body text explicitly
states that exposure time was adjusted between conditions to obtain similar
2100-2200 nm depths. The source also contains two unresolved internal
inconsistencies retained in the manifest and audit:

1. `0.0168 -> 0.05255 mA/cm2` is a `212.8%` increase, not the printed `21%`.
2. the power sweep is called `56.25%` C4F6 while the printed
   `15/25 sccm CF4/C4F6` pair is `62.5%` C4F6.

Replay:

```bash
python scripts/digitize_woo_2024_c4f6_board.py
```
