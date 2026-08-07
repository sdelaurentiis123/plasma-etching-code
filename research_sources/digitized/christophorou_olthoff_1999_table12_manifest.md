# Christophorou--Olthoff 1999 Table 12 transcription manifest

## Source identity

- Primary source: L. G. Christophorou and J. K. Olthoff, “Electron
  Interactions With Cl2,” *Journal of Physical and Chemical Reference Data*
  **28**, 131--169 (1999).
- DOI: `10.1063/1.556036`
- Official NIST PDF:
  `https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=8765`
- Source PDF SHA-256:
  `a70bbed40bae014a551dd91fc96e322a258d959674ea53fe50ea2c6f4e020a6b`
- Source location: journal page 149, PDF page 19.

## Pixel audit

- Renderer: Poppler `pdftoppm`, 300 dpi PNG.
- Full audited page filename: `page-19.png`.
- Full page SHA-256:
  `46a3b5a8b9aa41ec69279cbf45112161139be9ec49a7bd0e964b048a4099d5e3`.
- Inspection: original-resolution visual review of all 29
  energy/cross-section pairs, the units `10^-20 m2`, the “suggested total
  ionization” label, and the text defining the evaluated curve.

The PDF and render are not redistributed. These hashes pin the exact pixels
used to validate
`christophorou_olthoff_1999_table12_cl2_total_ionization.csv`.

## Physical and evidence boundary

The executable provider uses the independently printed `11.481 eV`
ground-state ionization energy and linearly interpolates the evaluated table.
It does not assign a scalar uncertainty: the review explicitly says Table 12
averages the Kurepa--Belic and Stevie--Vasile measurements even though their
magnitudes differ by more than their combined quoted uncertainties.

The table is total ionization only. The same page states that no partial
electron-impact ionization data exist and that the relative production of
`Cl2+` and `Cl+` is unknown. The table can close total positive-ion production
and its associated minimum power sink, but cannot close the species-resolved
ion flux delivered to a wafer.

The table ends at 100 eV. The generic Maxwellian support check therefore
permits the industrial `0.3--5 eV` temperature band identified by Hamilton et
al. and fails closed when the unmeasured high-energy kernel becomes material.
