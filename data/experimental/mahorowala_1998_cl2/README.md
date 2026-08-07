# Mahorowala 1998 fixed-time Cl2/poly-Si board

This directory transcribes Table 2.2 of A. P. Mahorowala’s MIT PhD thesis,
*Feature profile evolution during the high density plasma etching of
polysilicon* (1998), handle `1721.1/50514`. It is a designed 13-run matrix in
a Lam TCP 9400SE: inductive power, wafer RF-bias power, and Cl2 flow vary
while pressure (10 mTorr) and time (75 s) remain fixed. Eleven runs report
poly-Si and oxide rates; runs 8 and 12 are explicitly marked overetched.

The `derived_*_removed_nm` columns are unit conversions, not fitted values:

```text
depth_nm = rate_A_per_min * 75 / 60 * 0.1
```

They span 112.5–459.375 nm of poly-Si removal. Figure 2.4 supplies the
corresponding complete SEM montage for 250 nm lines and 310 nm spaces. The
source PDF and pixels are not redistributed under the thesis rights notice.
With a local copy, reproduce the checksum and original-pixel audit:

```bash
pdftoppm -f 39 -l 39 -singlefile -png -r 600 \
  /path/to/42415621-MIT.pdf /tmp/mahorowala_page39
pdftoppm -f 42 -l 42 -singlefile -png -r 600 \
  /path/to/42415621-MIT.pdf /tmp/mahorowala_page42
python scripts/audit_mahorowala_1998_cl2_fixed_time.py \
  --source-pdf /path/to/42415621-MIT.pdf \
  --table-render /tmp/mahorowala_page39.png \
  --profile-render /tmp/mahorowala_page42.png \
  --table-overlay /tmp/mahorowala_table2_2_overlay.png
```

This is a substantially stronger chlorine feature board than a knob-only
recipe: absolute time and rates are published, every condition maps to an SEM
panel, and the beam-derived Chang/Sawin surface law is independent. It still
does not publish species-resolved wafer fluxes or measured IEAD/IAD for each
run. The thesis’s 2 mA/cm2 and 100–120 eV center-condition values are estimates.
Accordingly the board can grade a reactor provider or a facility-conditioned
profile transfer, but it cannot yet grant a first-principles knobs-to-depth
pass.
