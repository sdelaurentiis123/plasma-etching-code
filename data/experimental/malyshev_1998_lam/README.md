# Malyshev 1998 Lam Alliance chlorine reactor data

This directory freezes the measured markers in Figures 7 and 8 of Malyshev
et al., *J. Appl. Phys.* **84**, 137--146 (1998), DOI
`10.1063/1.368010`.

The measurements are Cl2 number density relative to its plasma-off value in a
commercial Lam Research Alliance metal etcher. The two figures provide 11 cm
and 6.5 cm TCP-window-to-wafer gaps at 0.5, 1, 2, and 10 mTorr. The committed
CSV also reports dissociation as `100 - relative Cl2 percent`.

Important boundaries:

- TCP source power was measured into the matching network; it is not measured
  absorbed plasma power.
- The plotted bars span independent Ar- and Xe-based actinometry reductions;
  they are not statistical sigma.
- The paper estimates absolute Cl2 density accuracy at about +/-25%.
- The documented 10 mTorr, approximately 100 W emission anomaly is excluded.
- Smooth curves in the figures are model outputs and are not digitized.
- These data validate reactor dissociation, not wafer flux or feature depth.

Replay the native-pixel extraction, checksum verification, CSV, manifest, and
QA overlays with:

```bash
python scripts/digitize_malyshev_1998_lam_dissociation.py \
  --source-pdf /path/to/137_1_online.pdf \
  --overlay-directory tmp/pdfs/malyshev_1998_qa \
  --write
```
