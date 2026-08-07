# Tinacba et al. 2021 SF5+ atomistic/beam depth board

Primary source: E. J. C. Tinacba, T. Ito, K. Karahashi, M. Isobe,
and S. Hamaguchi, “Molecular dynamics simulation for reactive ion etching of
Si and SiO2 by SF5+ ions,” *Journal of Vacuum Science & Technology B* **39**,
043203 (2021), DOI `10.1116/6.0001230`.

Figure 8 compares a DFT-informed modified-Stillinger–Weber MD calculation
against an identified beam experiment. The ion is mass selected, its energy is
checked by an energy-mass analyzer, dose is measured with a Faraday cup at the
sample position, and etched depth is measured by contact profilometry. No
radical beam was used. The authors convert depth `d`, number density `N`, and
dose `D` to yield with `Y=dN/D`; predicting yield at fixed `D,N` is therefore
exactly an absolute depth-per-dose prediction.

The MD potential was not fit to the beam depth or yield. Its S–F part was
conditioned on DFT and a published S–F bond energy, while the older Si/O/F
potentials were transferred. This makes the filled-circle experimental points
independent validation of the open-diamond atomistic calculation. It is
retrospective, not prospectively blinded.

The model has a material limitation that matters: sulfur carries five F atoms
and has the correct mass/radius, but S–S, S–Si, and S–O chemistry is
intentionally suppressed. The source argues the agreement indicates F-assisted
sputtering dominates on this board; the result must not be transferred to low
energies or chemistries where sulfur reactions matter.

## Reproducible vision audit

Download the author-hosted PDF:

```bash
curl -L -o tmp/sources/tinacba_2021/paper.pdf \
  https://researchmap.jp/satoshi_hamaguchi/published_papers/37820817/attachment_file.pdf
```

Its required SHA-256 is
`c0be3b475aa17b396c1f788baee14ba37b9026b264bb870dc2553055f27b31ad`.
Render PDF page 8 at 300 dpi:

```bash
pdftoppm -f 8 -l 8 -r 300 -png -singlefile \
  tmp/sources/tinacba_2021/paper.pdf \
  tmp/sources/tinacba_2021/page8_300dpi
python scripts/digitize_tinacba_2021_sf5.py \
  --overlay tmp/sources/tinacba_2021/figure8_overlay.png
```

The script verifies the 2550 × 3375 render checksum, dark axis pixels, the
committed CSV and manifest, and every full-resolution marker center. A
conservative ±0.04 Si/ion digitization bound covers 3.9 pixels; it is not
presented as experimental statistical uncertainty.

The board identifies a mass-selected SF5+ surface provider. It does not identify
an SF6 reactor, an SFx+ mixture, an angular law, neutral-F synergy, or the
Yoshie/Krüger wafer boundary.
