# Nozawa 1995 notching benchmark — primary evidence

Primary experiment: T. Nozawa et al., “The Electron Charging Effects of Plasma on Notch Profile
Defects,” *Japanese Journal of Applied Physics* **34**, 2107–2113 (1995), DOI
[`10.1143/JJAP.34.2107`](https://doi.org/10.1143/JJAP.34.2107).

Hwang and Giapis subsequently modeled these measurements in “On the origin of the notching effect
during etching in uniform high density plasmas,” *JVST B* **15**, 70–87 (1997), DOI
`10.1116/1.589258`. The CSV in this directory comes from the **Nozawa primary experiment**, not from
the later simulation.

## Source and extraction provenance

- source PDF SHA-256:
  `87500f53f0286aae0597b14168b0991791db73138ee8d67a5e5e3bc56b329c67`
- source: author/user-supplied lawful copy; the copyrighted PDF is not vendored
- source format: eight-page scan; quantitative pages contain embedded 4960 x 7008 px, one-bit images
- extraction: `pdfimages -png Nozawa_1995_Jpn._J._Appl._Phys._34_2107.pdf page`
- extracted PDF page 6 / printed page 2111 SHA-256:
  `efc83d088516a5b74869d2ed744d57062e57091153f1729c2d50514e21f71f61`
- extracted PDF page 7 / printed page 2112 SHA-256:
  `53f3b150e5d6d89fda0ef1f71a8b35097b464f60454d2d3052d72efda1b4ac7a`
- evidence table: `digitized_notch_curves.csv`
- evidence-table SHA-256:
  `2e472385e002aebf94f2f0ec299f877180786b6095791ec3079a80fc48a22ec2`

OCR/Markdown is used only as a searchable reading aid. Numerical values come from the 600-dpi source
pixels and are replayed from the recorded marker coordinates and linear axis maps. Pixel coordinates
are in the full extracted page-image coordinate system, not in an undocumented screenshot crop.

## Quantitative experiment selected before engine scoring

The table records 17 markers from three curves:

1. Figure 10b open-area width, open-circle series (`L=0.6 um`, `S=0.6 um`): seven calibration
   markers. This curve contains the sharp rise and large-open-area saturation used in the later
   Hwang–Giapis validation study.
2. Figure 8b space width with all connected lines sharing one pad: five held-out markers. Notch depth
   falls as the space widens.
3. Figure 9b space width with each connected line group tied to its own pad: five held-out markers.
   Notch depth rises as the space widens.

The opposite Figure-8/Figure-9 trends are intentional. They test whether the engine transports charge
and electron supply through the declared electrical topology; a generic width-dependent fit cannot
pass both by construction. Figure 7 (perimeter ratio) and the Figure-1/10 SEM topology observations
remain useful expansion evidence but are not silently mixed into this first frozen numeric split.

## Uncertainty and claim boundary

Marker centers were inspected against the original page pixels. A conservative `0.01 um` notch-depth
digitization bound is wider than the center-localization ambiguity and is comparable to half the
printed marker height. Control-axis bounds are `0.12 um` for Figures 8/9 and `0.30 um` for Figure 10.

The paper reports no statistical measurement uncertainty or error bars for these notch-depth markers.
Every row therefore says `measurement_uncertainty_semantics=not_reported`. The code may report error,
trend reproduction, and digitization-only coverage, but `score_notching_benchmark_3d` refuses the
headline “validated notch prediction” until an experimental uncertainty is independently justified.
Digitization uncertainty is never relabeled as measurement uncertainty.

The four values embedded in the legacy `scripts/notching_depth_gate.py` have no image checksum, raw
pixel coordinates, or replayable axis map. They remain a historical shape diagnostic and are not C4
evidence.

## Process conditions reported by the paper

- resist mask: 10000 A; poly-Si: 3000 A; SiO2: 1000 A on silicon
- poly-Si sheet resistance: 30 ohm
- 200% overetch for the SEM notch-depth measurements
- microwave power: 700 W
- 400 kHz RF bias: 10 W, 60 V peak-to-peak
- pressure: 3.0 mTorr
- Cl2 flow: 100 sccm
- magnetic field at wafer: 500 G
- wafer temperature: 0 degrees C
- reported poly-Si etch rate: 2325 A/min
- selectivity: 100 versus SiO2; 20 versus resist

These reactor settings do not uniquely determine the feature-boundary IEDF, IADF, electron-energy
distribution, or absolute charged-particle flux. Any engine replay must list those as inferred boundary
inputs with provenance and sensitivity bounds rather than pretending the paper measured them.
