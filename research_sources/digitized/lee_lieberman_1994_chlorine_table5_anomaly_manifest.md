# Lee--Lieberman 1994 chlorine Table 5 anomaly manifest

## Source identity

- Primary report: C. Lee and M. A. Lieberman, UCB/ERL M94/49, revised
  28 November 1994.
- Primary record:
  `https://digicoll.lib.berkeley.edu/record/134386`
- Source PDF SHA-256:
  `f2049e7041984d658d23688e8e8112a8d8e8a524172a8d2e335be8fde7fc2e23`
- Source location: report page 27, PDF page 31.

## Visual audit and finding

- Renderer: Poppler `pdftoppm`, 300 dpi PNG.
- Full page SHA-256:
  `77bb0df77d646a11855dc2ecffd3aeb1f1aae60eab60aa4f32399df03cb732bb`
- Audited 1500x900 grayscale crop SHA-256:
  `346a44af64dd2ad663e46e1a3bec20d2936b05d3e1980c386fa65e77b07acd40`
- The first row visibly prints
  `e + Cl(2P) -> Cl(3D) + 2e`.
- The remaining atomic rows are `4D`, `4P`, `4S`, `5D`, and `5P`.

The principal/orbital labels denote atomic excited states; the NIST critical
review independently identifies calculated Cl excitation channels including
3d, 4d, 4p, 4s, 5d, and 5p. An excitation event must retain one electron.
The printed first row creates net charge and cannot be imported literally.

## Executable decision

Table 5 is quarantined as a source defect until the original cited channel is
recovered. The predictive electron-energy ledger remains fail-closed. The row
is not silently repaired and is not interpreted as an ionization channel.
