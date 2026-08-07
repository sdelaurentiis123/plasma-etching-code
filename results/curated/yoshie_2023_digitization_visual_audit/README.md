# Yoshie 2023 full-resolution digitization replay

Status: **all five publisher rasters and all digitized marker boards visually
reconciled at original resolution on 2026-08-06**.

This receipt is an independent human-vision replay of the existing
checksum-bound PIL digitizers. It changes no extracted value and is not a
model calibration.

## Replayed evidence

| figure | committed observations | visual assertion |
|---|---:|---|
| 4 | 7 blanket-rate markers | every selected marker center lies on the plotted symbol for its declared timing |
| 5 | 21 held-out feature markers plus error caps | all centers and upper/lower caps lie on the plotted board |
| 6 | 28 held-out feature markers plus error caps | all centers and upper/lower caps lie on the plotted board |
| 12 | 7 bias windows, each sampled at start/mid/end | window endpoints and midpoints align; the short timing-II collapse in the 8 s cycle is preserved |
| 14 | 60 OES markers | every CF, CF2, and F center lies on its color-coded plotted marker |

The replay used the official publisher JPEG dimensions and hashes already
pinned in the two committed manifests. Overlays were generated locally and
inspected with original-resolution image rendering; publisher rasters and
derived overlays are not redistributed.

## Commands

```bash
python scripts/digitize_yoshie_2023_figures4_6.py \
  --source-dir tmp/sources/yoshie_2023 \
  --overlay-dir tmp/yoshie_2023_overlays
python scripts/digitize_yoshie_2023_reactor_state.py \
  --source-dir tmp/sources/yoshie_2023 \
  --overlay-dir tmp/yoshie_2023_reactor_state_overlays
```

The machine-readable receipt records the exact source and overlay SHA-256
values. The science boundary remains unchanged: these measurements provide an
excellent held-out profile board and cycle-state diagnostics, but Yoshie does
not publish species-resolved wafer fluxes or an IEAD.
