# Mahoney et al. 1994 Table I argon, 100 W transcription manifest

- **Primary source:** L. J. Mahoney, A. E. Wendt, E. Barrios, C. J.
  Richards, and J. L. Shohet, “Electron-density and energy distributions in a
  planar inductively coupled discharge,” *Journal of Applied Physics* **76**,
  2041–2047 (1994), DOI `10.1063/1.357672`.
- **Primary host:** University of Wisconsin MINDS author/institutional
  repository, record `1793/9564`.
- **Source PDF SHA-256:**
  `acd59d5def6373a81f1ec73248608b7c67a8769b010f2287fa05b04fd8cc61b7`.
- **Visual source:** PDF page 4 rendered at 200 dpi.
- **Render SHA-256:**
  `cc2694e820db4b38660666cbbcf5aac721b0d92cf4165f4fd5304e129943af74`.
- **Method:** direct transcription of the five `100 W` rows in printed
  Table I. This is not a curve digitization and no OCR value was accepted
  without native-resolution visual inspection.
- **Independent cross-check:** the `10` and `20 mTorr` peak-density values
  agree visually with the corresponding 100 W symbols in Figure 11. Figure
  11 was inspected from PDF page 6 at 200 dpi; render SHA-256
  `906940002ddc809e8c60dbc4aaed0c4c7aa225c487dd53b148857e64e75de3b9`.
- **Unit conversion for grading:** the printed density unit is
  `10^10 cm^-3`; multiply by `1e16` to obtain `m^-3`.
- **Duplicate 20 mTorr rows:** retained deliberately. They were acquired
  with different pumping paths and provide a direct reproducibility/system
  sensitivity check rather than a value to average away.
- **Copyright boundary:** only numerical facts and checksums are committed;
  source PDF pixels are not redistributed.

## Condition facts read from the same primary paper

- Grounded metal liner: `22.8 cm` inner diameter and `13.7 cm` length.
- Peak-density probe location: `r = 0 cm`, `z = 5.0 cm`.
- Argon pressure range in Figure 11: `2–20 mTorr`; Table I extends to
  `100 mTorr`.
- The paper reports **net rf power** as incident minus reflected power, with
  reflected power below 5% after matching. It explicitly warns that some
  further power can be dissipated in the induction coil and matching network.
  Therefore `100 W` is not silently relabeled as measured absorbed plasma
  power.
- The paper does not report neutral-gas temperature. Validation must bracket
  that missing boundary rather than tune it to the measured density.
- The authors state that their electron-density diagnostic can read two to
  five times below ion-density determinations even in electropositive
  discharges. That stated diagnostic interval is frozen as a comparison
  boundary, not inferred from the model result.
