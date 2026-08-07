# Gray 1993 Ar+/F silicon beam board

Primary source: D. C. Gray, *Beam simulation studies of plasma-surface
interactions in fluorocarbon etching of Si and SiO2*, MIT PhD thesis (1993),
handle [1721.1/13187](https://hdl.handle.net/1721.1/13187).

`table5_9_ar_f_si_model_parameters.csv` is a direct transcription of Table
5-9.  The values outside parentheses are independent `s0`/`beta2` regressions;
the parenthesized values are the thesis's constant-`s0 = 0.2` refit.  They are
not raw etch-yield measurements and are not assigned unstated uncertainties.

The source PDF SHA-256 used for the page-level visual audit was
`be6bce26b699b3172cf67bb68e4d12e039fd3ea775f73873ee1aaf25164c065b`.
At 300 dpi, PDF pages 246, 247, and 252 were inspected at original resolution
to verify Eqs. 5-30/5-31, Table 5-9, and Eq. 5-34.  The relevant OCR is
archived in
`research_sources/thesis_extracts/gray_thesis_1993_ocr_sections.txt`.

The direct closure in `src/petch/gray_argon_fluorine_si.py` uses no feature
depth.  It is valid only for atomic F, Ar+, normal incidence, and its declared
beam energy/flux-ratio board.  In particular, it does not license replacing
the unpublished SFx+ mixture in the Yoshie reactor with Ar+.
