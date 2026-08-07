# An et al. DFT-trained NNP/MD outputs

`figure3_reported_etch_yields.csv` is a lossless transcription of the
simulation entries in the authors' released Figure 3 `dat.yaml`, pinned in
`source_manifest.json`. It contains both active plot values and the zero-yield
deposition values that the authors retained as comments. The latter are marked
`not_plotted`; no uncommenting or regime inference is hidden.

The values are **model outputs, not measurements**. They are stored so the
DFT-trained atomistic model can be tested, without fitting, against petch's
independently digitized Karahashi mass-selected beam board. The table is not an
executable yield closure and does not authorize interpolation, energy
extrapolation, or transfer to a reactor ion mixture.

The authors' code and neural-potential weights are not redistributed here.
The inspected Git commit did not contain a license file. See the manifest for
the exact source commit, checksums, and scientific/implementation boundary.
