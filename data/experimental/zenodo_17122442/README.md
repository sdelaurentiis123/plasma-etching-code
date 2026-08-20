# Bosch wafer measurements — Zenodo 17122442

Source: Sayyed et al., *A Multi-Model Dataset for BOSCH Plasma-Etching: Optical Emission Spectra,
Process Parameters, and Wafer Measurements for Data-Driven Plasma Modeling*, version 1 (2025),
DOI [10.5281/zenodo.17122442](https://doi.org/10.5281/zenodo.17122442).

License: Creative Commons Attribution 4.0. The authors and source must be credited on reuse.

Included here:

- `Si_Oxide_etch_9_points.csv` — 684 measurements from 76 processed 200-mm silicon wafers, at nine
  nominal positions per wafer (some source measurements are unavailable).
- `Si_Oxide_etch_89_points.csv` — 7,832 measurements from 88 identified wafers. Its distinct source
  schema preserves 157 originally unavailable post-oxide measurements as `N/A` in
  `postox_thickness_nan`, alongside the processed `postox_thickness` used in the derived columns.
- Nine rows forming one complete nine-position wafer have blank experiment, lot, and wafer identifiers
  in the source CSV. The physical measurements are retained and the loader exposes those identifiers as
  missing; identified-wafer analyses must exclude or separately handle this record.
- Source MD5: `78515caf25e29e558e1859b92f8a4827`, verified at acquisition on 2026-07-11.
- 89-point source MD5: `446e75b040eea37b634eeb8f763a62fc`, verified at acquisition on 2026-07-12.
- All measurement values and coordinates are in micrometres, per the dataset README.
- `Process_data.nc` - all 96 available wafer process records, 31 common machine channels at 5 Hz
  (the first day exposes 13 additional channels); source MD5
  `4567d24ec2125102a2e5129203ba31fa`.
- `Dictionary_process.nc` - the shared lossless process decoder; source MD5
  `0dde5a3a913eb1fa8512ef2f8748fb34`.
- `Lot_status.xlsx` - experiment date, lot, conditioning class, and source missingness; source MD5
  `339beca13f321dc3af244bf2d2ce284c`.
- `Readme.pdf` - the four-page source documentation; source MD5
  `f9a9bd323bd9e227a486249a631e8468`.
- `process_wafer_summary.csv` and its manifest - a deterministic label-free extraction from only
  the two process NetCDF files. Rebuild/check it with
  `python scripts/extract_zenodo_bosch_process_features.py [--check]`.

Experiment: SPTS Omega i2L DSi Rapier Bosch process, SF6 etch/C4F8 passivation, 100 cycles with 4.5 s
etch and 1.5 s passivation. Wafers have more than 99.5% exposed silicon and a nominal 1 µm SiO2 mask.
The study varies chamber conditioning and records wafer sequence, synchronized OES, and machine data.
The compact process NetCDF and its dictionary are vendored here. The ten daily OES files total about
7.9 GB and remain at Zenodo; they are not required to reproduce the machine-channel extraction.

The source does not explicitly map the anonymized gas-channel numbers to chemical names. The process
extractor declares a waveform-role inference: `Gas5Flow` is the 600-unit, approximately 4.5-second
etch waveform and is mapped to SF6; `Gas4Flow` is the 300-unit, approximately 1.5-second passivation
waveform and is mapped to C4F8. This matches the source recipe and 100-cycle timing, but remains an
inference until the authors or tool export provides the channel dictionary.

Important scope: these are reactor/wafer-scale depth, selectivity, uniformity, and drift observations.
They do not validate feature-profile charging, ARDE, sidewall shape, or scallops.
