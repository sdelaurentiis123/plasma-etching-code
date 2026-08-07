# Gray 1993 Ar+/F silicon and silicon-dioxide beam board

Primary source: D. C. Gray, *Beam simulation studies of plasma-surface
interactions in fluorocarbon etching of Si and SiO2*, MIT PhD thesis (1993),
handle [1721.1/13187](https://hdl.handle.net/1721.1/13187).

`table5_9_ar_f_si_model_parameters.csv` is a direct transcription of Table
5-9.  The values outside parentheses are independent `s0`/`beta2` regressions;
the parenthesized values are the thesis's constant-`s0 = 0.2` refit.  They are
not raw etch-yield measurements and are not assigned unstated uncertainties.

`table5_10_ar_f_sio2_model_parameters.csv` is the corresponding direct
transcription for SiO2.  Its parenthesized values are the constant-`s0 = 0.02`
refit.  The separate Si and SiO2 rows must not be pooled: their sticking and
ion-assisted removal laws differ by more than a scale.

The source PDF SHA-256 used for the page-level visual audit was
`be6bce26b699b3172cf67bb68e4d12e039fd3ea775f73873ee1aaf25164c065b`.
At 300 dpi, PDF pages 243--247 and 252 were inspected at original resolution
to verify the stated mechanism limitations, branching law, Tables 5-7/5-8,
Eqs. 5-30/5-31, Tables 5-9/5-10, and Eqs. 5-34/5-35.  Raster SHA-256 values
for the decisive page renders
were `3625197771d2dd2c7f90325b5c56cde2f5afb3c6d204097957a183be333`
(p. 243), `b7006792b10c6266ac3b6692556de3c68c2b12015e8c296ba66391400e4f0fac`
(p. 244), `6ef669220860f4852fb2f779f433f9388d0b3d83526bb39faced7bedf74b6b52`
(p. 245), `7e4bb5dd68102cfb1baf3def89aee438bc4df9649de71dd3c0430e8a67a5e5ee`
(p. 246), `04c4229dc7b62aac301c00941fa45aa9915bdac759b10890ccd81263875463d5`
(p. 247), and
`28ac695f8ff5343636b0b519c291198bf7877ca3ab731585c38b433b264e9b67`
(p. 252).  The relevant OCR is
archived in
`research_sources/thesis_extracts/gray_thesis_1993_ocr_sections.txt`.

The direct closures in `src/petch/gray_argon_fluorine_si.py` and
`src/petch/gray_argon_fluorine_sio2.py` use no feature depth.  They are valid
only for atomic F, Ar+, normal incidence, and their declared beam
energy/flux-ratio boards.  In particular, they do not license replacing the
unpublished SFx+ mixture in the Yoshie reactor with Ar+.

These are source-faithful *reduced* beam closures, not atomistic reaction
potentials.  Gray explicitly calls the common SiF2 surface representation a
gross simplification, calls the two-new-site SiO2 cascade count arbitrary,
and says the full mechanistic description could not be verified for lack of
in-situ SiFx measurements (PDF pp. 243--245).  The implementations conserve
the atoms implied by the printed product topology; that accounting does not
promote the topology itself above its published evidence grade.
