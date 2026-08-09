# Mahorowala 1998 axisymmetric ion-transport audit

This is a conservative spatial lift of a solved 0-D inventory, not an independent feature-depth prediction. No etch depth selected any reactor or transport parameter.

- conditions: `9`
- maximum particle-ledger residual: `4.888e-16`
- 14x10 to 20x14 center-flux change: `1.63%` (3% gate: `True`)
- formal feature-depth passes: `0`

| source moment | raw wafer / Lee-global range | after center-current renormalization |
|---|---:|---:|
| uniform | 1.4824--1.8426 | 0.9090--1.1299 |
| top_center_broad | 1.4126--1.7574 | 0.9125--1.1353 |
| top_center_compact | 1.5769--1.9810 | 0.9122--1.1460 |
| top_annular_wise_class | 0.5843--0.5962 | 0.9977--1.0180 |
| inductive_em_annular | 0.6023--0.6369 | 0.9874--1.0442 |

The raw 1.4--2.0x correction is not an absolute-depth gain in the existing diagnostic-conditioned board: its independent 400 W/100 sccm center-current anchor renormalizes that common offset away. What remains is a roughly +/-15% trend correction. The spread is a source-field uncertainty, not permission to choose the mode that best matches depth. The discriminating measurement is a radially resolved, species-resolved ion-current profile at the wafer plane.
