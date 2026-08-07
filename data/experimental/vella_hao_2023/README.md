# Vella–Hao et al. Si/Cl₂/Ar⁺ ALE absolute etch per cycle

Primary source: J. R. Vella, B. Hao, V. M. Donnelly, and D. B. Graves, “Comparison of
Molecular Dynamics Simulations and Experiments of Si–Cl₂–Ar⁺ Atomic Layer Etching,”
*Plasma Processes and Polymers* 20, 2200198 (2023), DOI
[10.1002/ppap.202200198](https://doi.org/10.1002/ppap.202200198), OSTI
[2248044](https://www.osti.gov/biblio/2248044).

The pinned DOE manuscript is `https://www.osti.gov/servlets/purl/2248044`, with SHA-256
`789bf50302fc2ed9175403c47f895e1fb8db5481be5fb980739cad97f33c2218`.

`figure8_epc.csv` digitizes the six black experimental square markers in Figure 8 from a Poppler
240-dpi render of manuscript page 17. The plot frame used for the affine axis map is:

- zero-energy x pixel: 657
- 250 eV x pixel: 1450
- zero-depth y pixel: 827
- 8 nm/cycle y pixel: 268

Thus `E = (x - 657) / (1450 - 657) * 250 eV` and
`EPC = (827 - y) / (827 - 268) * 8 nm/cycle`. Marker centers were recovered with a PIL black-pixel
mask and connected components rather than manually clicked. A 1.25-rendered-pixel extraction audit
corresponds to ±0.394073 eV and ±0.0178891 nm/cycle. Source error bars are not assigned numerical
values here because the black bars are smaller than or obscured by the markers at this resolution.

The experiment measured the positive-ion flux as `3.7e16 cm^-2 s^-1` and used a 3 s ion-bombardment
step. It did **not** measure a species-resolved IEAD. The plotted mean ion energy was inferred from
electrical measurements and assumptions about electron temperature/plasma potential. This board is
therefore a no-depth-fit cross-source validation with a facility-conditioned energy boundary, not a
blind Tier-A prediction.

The source uses “monolayer” both as a fluence shorthand (`1e15 impacts/cm²`) and, in the associated
atomistic work, as a material count (72 Si atoms in the simulation cell). The petch depth audit never
equates those labels. It converts the atomistic trajectory through its physical cell area and material
atomic density, then applies the measured experimental fluence in `ions/cm²`.
