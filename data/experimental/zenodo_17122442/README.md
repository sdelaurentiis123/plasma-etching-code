# Bosch wafer measurements — Zenodo 17122442

Source: Sayyed et al., *A Multi-Model Dataset for BOSCH Plasma-Etching: Optical Emission Spectra,
Process Parameters, and Wafer Measurements for Data-Driven Plasma Modeling*, version 1 (2025),
DOI [10.5281/zenodo.17122442](https://doi.org/10.5281/zenodo.17122442).

License: Creative Commons Attribution 4.0. The authors and source must be credited on reuse.

Included here:

- `Si_Oxide_etch_9_points.csv` — 684 measurements from 76 processed 200-mm silicon wafers, at nine
  nominal positions per wafer (some source measurements are unavailable).
- Nine rows forming one complete nine-position wafer have blank experiment, lot, and wafer identifiers
  in the source CSV. The physical measurements are retained and the loader exposes those identifiers as
  missing; identified-wafer analyses must exclude or separately handle this record.
- Source MD5: `78515caf25e29e558e1859b92f8a4827`, verified at acquisition on 2026-07-11.
- All measurement values and coordinates are in micrometres, per the dataset README.

Experiment: SPTS Omega i2L DSi Rapier Bosch process, SF6 etch/C4F8 passivation, 100 cycles with 4.5 s
etch and 1.5 s passivation. Wafers have more than 99.5% exposed silicon and a nominal 1 µm SiO2 mask.
The study varies chamber conditioning and records wafer sequence, synchronized OES, and machine data.
Only the small measurement table is vendored here; the larger OES and process NetCDF files remain at
Zenodo and can be acquired by DOI when reactor-scale work begins.

Important scope: these are reactor/wafer-scale depth, selectivity, uniformity, and drift observations.
They do not validate feature-profile charging, ARDE, sidewall shape, or scallops.
