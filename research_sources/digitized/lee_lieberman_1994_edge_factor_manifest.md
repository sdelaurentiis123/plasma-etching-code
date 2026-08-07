# Lee--Lieberman 1994 sheath-edge equations: visual audit manifest

## Source identity

- Primary report: C. Lee and M. A. Lieberman, UCB/ERL M94/49, revised
  28 November 1994.
- Record: `https://digicoll.lib.berkeley.edu/record/134386`
- PDF: `https://digicoll.lib.berkeley.edu/record/134386/files/ERL-94-49.pdf`
- DOI of journal version: `10.1116/1.579366`
- PDF SHA-256:
  `f2049e7041984d658d23688e8e8112a8d8e8a524172a8d2e335be8fde7fc2e23`

The source PDF is not redistributed in this repository.

## Render audit

Poppler `pdftoppm` rendered source PDF pages 11 and 22 at 500 dpi. The
temporary review images had these SHA-256 values:

- PDF page 11, report page 7, Eqs. 13--14:
  `ecbb5b681d7a81fef3d77195952f2f1b30d9df60fc778dc93ab1dbf1e3bb4e30`
- PDF page 22, report page 18, Eqs. A.10--A.11:
  `0f1b0436ecd3188f364e84f0026d856ac80ed66e9ba25cdf38824e5dd273d7ae`

Both images were visually inspected at native render resolution. This check
was necessary because text extraction displaced fraction bars and omitted
the temperature-ratio symbol in Eq. 13.

## Exact implementation boundary

For volume-average electronegativity
`alpha_avg = n_minus / n_e` and electron-to-ion temperature ratio
`gamma = T_e / T_i`, Eqs. 13--14 give the same multiplicative correction for
the axial and radial edge ratios:

`C_en = (1 + 3 alpha_avg / gamma) / (1 + alpha_avg)`.

The remaining denominator is the electropositive expression when it is
written in terms of the ambipolar diffusion coefficient `D_a`. The code does
not insert `gamma` into that denominator a second time.

The derivation assumes:

- uniform electron density away from the sheath;
- a parabolic negative-ion density that vanishes at the sheath edge;
- `alpha_0 = (3/2) alpha_avg` for the profile conversion; and
- a sheath-edge to bulk-density ratio independent of positive-ion species.

Consequently, the implementation is a faithful Lee--Lieberman closure. It is
not evidence that these spatial assumptions hold in an arbitrary reactor.
