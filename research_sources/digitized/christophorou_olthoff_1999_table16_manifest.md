# Christophorou--Olthoff 1999 Table 16 transcription manifest

## Source identity

- Primary source: L. G. Christophorou and J. K. Olthoff, “Electron
  Interactions With Cl2,” *Journal of Physical and Chemical Reference Data*
  **28**, 131--169 (1999).
- DOI: `10.1063/1.556036`
- Official NIST PDF:
  `https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=8765`
- Source PDF SHA-256:
  `a70bbed40bae014a551dd91fc96e322a258d959674ea53fe50ea2c6f4e020a6b`
- Source location: journal page 153, PDF page 23.

## Pixel audit

- Renderer: Poppler `pdftoppm`, 500 dpi PNG.
- Full audited page dimensions: `4247x5570` pixels.
- Full page SHA-256:
  `4dfa2863a0622a06a0a58ccdc4b35b1d3fd3907af5cb12e872742389260064ad`.
- Table crop: ImageMagick grayscale `1950x2300+300+240`.
- Audited crop SHA-256:
  `97f73d5fcb067bd86a1415a2ff8c4aa097b51da279b32f4a4e2d19cbb3274164`.
- Inspection: original-resolution visual review of all 42
  energy/cross-section pairs, the `10^-20 m2` units, and the “suggested total
  dissociative electron attachment cross section” label.
- Independent text-layer comparison: all 42 pairs agree with Poppler's
  layout-preserving extraction. The pixel image, rather than OCR, is the
  transcription authority.

The PDF and render are not redistributed. These hashes pin the exact pixels
used to validate
`christophorou_olthoff_1999_table16_cl2_attachment.csv`.

## Physical and evidence boundary

The table spans only `0.05--11.8 eV`. It neither samples the near-zero end of
the exothermic attachment channel nor closes the high-energy EEDF tail. The
executable object therefore integrates only the printed support and reports
the missing lower and upper Maxwellian kernel fractions. It does not set the
cross section to zero outside the table and is not a complete particle-rate
provider.

Dissociative attachment removes an incident electron from the electron
population. Its directly supported kinetic-energy sink is the collision
moment `<sigma v E>`, not the activation parameter in a fitted particle-rate
law. The support-only object evaluates that moment separately from
`<sigma v>`. Neither moment may close a reactor electron-power solve until
the unmeasured support is shown immaterial or supplied by primary data.

The review states that it adjusted the Kurepa--Belic cross section upward by
30 percent to align with electron-swarm measurements. It does not assign the
resulting suggested table a scalar uncertainty, so the implementation carries
none.
