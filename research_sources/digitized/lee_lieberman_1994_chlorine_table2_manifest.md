# Lee--Lieberman 1994 chlorine Table 2 transcription manifest

## Source identity

- Primary report: C. Lee and M. A. Lieberman, UCB/ERL M94/49, revised
  28 November 1994.
- Primary record:
  `https://digicoll.lib.berkeley.edu/record/134386`
- Source PDF SHA-256:
  `f2049e7041984d658d23688e8e8112a8d8e8a524172a8d2e335be8fde7fc2e23`
- Source location: report page 24, PDF page 28.

## Visual audit

- Renderer: Poppler `pdftoppm`, 300 dpi PNG.
- Audited render filename: `page-28.png`.
- Audited render SHA-256:
  `5e699b16498e718bb90fac5b7cdddbd64659157e5118828960c39984dffce67a`
- Inspection: original-resolution visual review of every reaction,
  coefficient, sign, exponent, superscript, and unit in Table 2.

The source PDF and render are not redistributed in the repository. The
checksum makes a future audit deterministic against the same primary pixels.

## Transcription boundary

This is a formula transcription, not curve digitization. The executable rows
are in `src/petch/reactor_global/chlorine.py`; the human-readable transcription
is in
`research_sources/thesis_extracts/lee_lieberman_1994_global_model.txt`.

The transcription preserves the printed rate fits in `cm^3 s^-1` with
electron temperature in eV, then performs the unit conversion in the generic
rate-law layer. No coefficient was optimized.

Two source-model boundaries remain explicit:

1. the printed `Cl*` mutual-neutralization product is lumped into tracked
   ground-state `Cl` because the report omits chlorine metastable balances;
2. physical electron energy losses are not inferred from rate-fit exponents.
   Electron-power evaluation remains disabled until Tables 4--5 and the cited
   primary channel sources close that ledger.

The printed atomic-chlorine ionization expression uses `log(Te/12.96)`.
The implementation currently uses natural logarithm, records that
interpretation, and requires direct confirmation from the cited Lennon source
before predictive promotion.
