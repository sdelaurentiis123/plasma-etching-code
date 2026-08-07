# Christophorou--Olthoff 1999 Table 25 transcription manifest

## Source identity

- Primary source: L. G. Christophorou and J. K. Olthoff, “Electron
  Interactions With Cl2,” *Journal of Physical and Chemical Reference Data*
  **28**, 131--169 (1999).
- DOI: `10.1063/1.556036`
- Official NIST PDF:
  `https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=8765`
- Source PDF SHA-256:
  `a70bbed40bae014a551dd91fc96e322a258d959674ea53fe50ea2c6f4e020a6b`
- Source location: journal page 160, PDF page 28.

## Pixel audit

- Renderer: Poppler `pdftoppm`, 300 dpi PNG.
- Full audited page filename: `page-28.png`.
- Full page SHA-256:
  `3755914af7e970b369a1e7bdf3578bea65ff961858cf8e214759dde95263c940`
- PIL/ImageMagick-equivalent crop:
  `1350x1550+1200+250`, normalized output `1348x1550`.
- Audited crop SHA-256:
  `6a01d03172e2d49619998e0593d14e9b547ad01803f45a6677068768cf599c25`
- Inspection: original-resolution visual review of all 48 energy/cross-section
  pairs, the units `10^-20 m2`, and the Hayes attribution.

The PDF and render are not redistributed. These hashes pin the exact pixels
used to validate
`christophorou_olthoff_1999_table25_atomic_cl_ionization.csv`.

## Physical boundary

Table 25 contains rounded selected measurements, including tiny entries below
the current NIST ASD ground-state threshold. The executable provider preserves
the raw rows but forces the cross section to zero below the independently
sourced `12.967633 eV` threshold. No point was optimized.

The table ends at 200 eV. Maxwellian averaging therefore fails closed whenever
the fraction `(Emax/Te + 1) exp(-Emax/Te)` of a constant-cross-section kernel
above the measured support exceeds `1e-6`. This permits the 2--10 eV reactor
domain and rejects materially unsupported high-temperature use.
