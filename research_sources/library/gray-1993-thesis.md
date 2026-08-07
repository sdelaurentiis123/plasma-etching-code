# gray-1993-thesis

**Gray, PhD thesis, MIT (1993)**

- **DOI/URL:** PhD thesis, MIT (1993)
- **Retrieval route:** MIT DSpace bitstream
- **Status:** FULL TEXT (OCR sections): research_sources/thesis_extracts/gray_thesis_1993_ocr_sections.txt
- **Topic:** beam-yields — Beam-measured yields, thresholds and sticking (the L0/L1 provenance floor)

## 2026-08-06 silicon and silicon-dioxide beam closure audit

The MIT DSpace PDF used for the audit has SHA-256
`be6bce26b699b3172cf67bb68e4d12e039fd3ea775f73873ee1aaf25164c065b`.
PDF pages 243--247 and 252 were rendered at 300 dpi and inspected at original
resolution.  The visual audit confirms:

- Tables 5-7/5-8 explicitly route a SiF2-like near-surface species through
  separate SiF2 and SiF4 volatile-product branches;
- Eqs. 5-27/5-28 give `b_Si = 0.009 sqrt(E_i)` and
  `b_SiO2 = 0.007 sqrt(E_i)`;
- Eq. 5-30:
  `theta_F = s0 R / [s0 R + 2 beta2 (1 + b)]`;
- Eq. 5-31 separates thermal, bare physical, and fluorinated
  `beta2(1+b)theta_F` removal;
- Tables 5-9/5-10 contain the Ar+/F/Si and Ar+/F/SiO2 rows now transcribed under
  `data/surface_interactions/gray_1993/`; and
- Eqs. 5-34/5-35 are
  `beta2,Si = 0.687 [sqrt(E_i) - sqrt(4)]` and
  `beta2,SiO2 = 0.053 [sqrt(E_i) - sqrt(4)]`.

The resulting implementations,
`src/petch/gray_argon_fluorine_si.py` and
`src/petch/gray_argon_fluorine_sio2.py`, are
species/energy/site/product beam closures with explicit elemental
bookkeeping.  They are not called first-principles, contain no feature-depth
fit, and refuse off-normal, reactive-ion, and out-of-board use in strict
mode.

### Evidence ceiling

The same page-level audit also fixes the ceiling on these claims.  Gray calls
the common SiF2 near-surface representation a "gross simplification" (p. 244),
the two-new-site count for a SiO2 cascade event an "arbitrary assumption"
(p. 243), and says that the full mechanism could not be verified because
in-situ SiFx concentrations versus energy and flux were unavailable
(p. 243).  Therefore:

- the printed yield and branching laws are direct controlled-beam
  benchmarks within their declared board;
- product-resolved atom accounting is stricter than the legacy SiF4-only
  approximation; but
- neither implementation is evidence for an elementary event network or for
  transfer to a fluorocarbon reactor boundary.

The decisive 300 dpi page-render SHA-256 values are
`3625197771d2dd2c7f90325b5c56cde2f5afb3c6d204097957a183be333`
(p. 243),
`b7006792b10c6266ac3b6692556de3c68c2b12015e8c296ba66391400e4f0fac`
(p. 244),
`6ef669220860f4852fb2f779f433f9388d0b3d83526bb39faced7bedf74b6b52`
(p. 245),
`7e4bb5dd68102cfb1baf3def89aee438bc4df9649de71dd3c0430e8a67a5e5ee`
(p. 246),
`04c4229dc7b62aac301c00941fa45aa9915bdac759b10890ccd81263875463d5`
(p. 247), and
`28ac695f8ff5343636b0b519c291198bf7877ca3ab731585c38b433b264e9b67`
(p. 252).

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | Gray 1993 **Table 5-10, "Ar+/F-SiO2 Model Parameters"** (p. 247; archive L1492-1500). Verbatim: | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:103 |
| Q2 | Gray 1993 **Table 5-1, "Si and SiO2 Low Energy Sputtering Model Parameters"** (p. 159; archive L132-135), verbatim: | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:141 |
| Q3 | Gray 1993 **Table 6-1, "Ar+/CF2 SiO2 Etching Model Parameters"** (p. 305; archive L1873-1878): | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:190 |
| Q4 | **Measured coverage.** - Chemically-enhanced SiO2: Gray's own data runs **to 2000 eV** and the sqrt law holds there (p. 252, quoted above). 3406 eV is only **1.70x** beyond his top measured point. - Physical sputter of SiO2: Gray's Figure 5-2 spans **0-1225 eV** (Gray 1992, Chapman 1980, Oostra 1986, Steinbruchel 1989). 3406 eV is 2.8x beyond. - Gray's beam itself ran 20-500 eV (p. 155, archive L9; and his Table 5-3 header, archive L1081: "E_i = Ion energy (20-500 eV)"); the 1000 eV and 2000 eV points come from Gerlach-Meyer 1981 and Tu et al. 1981 respectively (Figs. 5-5, 5-13). | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:446 |
| Q5 | - Gray p. 157 (archive L77): *"Sigmund [1969] has developed analytical solutions ... in the medium **(0.1-1 KeV)** and high **(>1 KeV)** energy ranges"* — 1 keV is the regime boundary. - Chang thesis L2249-2251: *"In the low ion energy regime **(< 1 kV)**, the binary particle interactions can be characterized by a Born-Mayer-type cross section, and the sputtering yield is linear to the square root of ion incident energy"*. | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:473 |
| Q6 | `BENCHMARK_CERTIFICATION_2026-08-06.md` already carries Gray's measured constants: SiO₂ ion-enhanced yield `β_e = 0.053(√E − √4)` "the number of SiF₄ molecules removed from **fluorine saturated** surface regions per incoming ion" and F adsorption `s0 = 0.02` (Gray MIT thesis 1993, Table 5-10 + p.246, co-regressed). The certification also records: "Gray's own printed parameters predict [a half-rise at F/Ar⁺ of] **96–127** across the plausible energy assignment". Those are the numbers §4 grades the band against. | `RESEARCH_F_SUPPLY_BAND_2026-08-06.md`:361 |
| Q7 | Against Gray's own printed half-rise of **96–127** in F/Ar⁺ [`BENCHMARK_CERTIFICATION` §(e) half-rise note], a floor ratio of 0.55–3.6 is **27–230× below half-rise**. **Atomic F at any point in the sourced band cannot fluorine-saturate this floor.** The floor's saturation, if it happens, must come from chemisorbed CFₓ (bare 0.278 / activated 0.85–0.9), which is 3–9× more sticky per strike and carries 1–3 bound F per arrival. Consistent with `RESULTS_ML23_GRAY_LAWS` §5: "the chemical channel alone carries the observed etch **if the floor is fluorine-saturated**". | `RESEARCH_F_SUPPLY_BAND_2026-08-06.md`:648 |
| Q8 | **Primary, fetched:** 5. O. Kwon, "Surface kinetics modeling of silicon oxide etching in fluorocarbon plasmas," ScD thesis, MIT DMSE, 2004, advisor H. H. Sawin. https://hdl.handle.net/1721.1/28353 (local copy: scratchpad `kwon_thesis.pdf`, PDF pages = printed page + 1). Key: **p. 76 Fig. 3.4 (SiO2 yield vs F/Ar+, Gray 1993 data, 350 eV)**, p. 36 Fig. 2.6 (Butterbaugh three-beam), p. 35 (three-regime text), pp. 64-65 (TML assumptions), **pp. 72-75 (SiO2/F TML rate laws)**, p. 40 Fig. 2.10 and p. 52 Fig. 2.20 (reactor N/I 70-230). | `RESEARCH_NEUTRAL_LIMITED_REGIME_2026-08-05.md`:875 |
| Q9 | 1. **Thermal re-emission of thermalized radicals** (the E8 completion) — the only path by which the source reaches the etch front, and the mechanism behind Huang's >95%. 2. **The (s0, B0) adsorption pair** (Kwon/Sawin E1). Gray's printed SiO2 sticking is 0.02 (thesis p.246: "setting s,=0.2 and 0.02 for the cases of silicon and SiO, etching respectively"), but it is half of a co-regressed pair: landing the scalar alone moves the measured Gray half-rise from 1.94 to 104.9 against his measured 27 +/- 8, and breaks four validated chemistry gates carrying s = 1 in their own closed forms. Left at 1.0 with a gate pinning the choice. 3. The CFx+ fraction for this reactor — unpublished; swept, never fitted. | `RESULTS_E8_THERMALIZED_RETURN_2026-08-05.md`:115 |
| U10 | [unquoted — verify on next use] Half-rise note: the previously quoted target of 27 +/- 8 came from a digitized *replot* (Kwon Fig. 3.4) that the harness itself flags `[VERIFY]`. Gray's own printed parameters predict 96-127 across the plausible energy assignment, so the primary source refutes the replot. Graded against the primary source, the pair transplant moved petch from 0.038x to 1.76x Gray's own model — a 26x improvement — and the residual 1.76x is recorded as an open MISS, not smoothed. | `BENCHMARK_CERTIFICATION_2026-08-06.md`:77 |
| U11 | [unquoted — verify on next use] The dossier's named path was run. Krüger's two-row `n = 1` anomaly on the SiO₂ channels replaced by Gray's own measured √E laws (`2f1e218`), then graded at feature scale (`ml23-gray-12s`). | `MIXED_LAYER_FEATURE_CAMPAIGN_2026-07-24.md`:353 |
| U12 | [unquoted — verify on next use] **Answer: it is measured, it is `sqrt(E) - sqrt(Eth)`, and the measurement is in the same source that gave us the absolute magnitudes.** Gray's own thesis tabulates the F-saturated SiO2 yield at **six ion energies from 20 eV to 2000 eV** and publishes the fitted law. Nothing here is inferred from a review or a modelling paper — it is the primary beam data. | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:9 |
| U13 | [unquoted — verify on next use] `research_sources/thesis_extracts/gray_thesis_1993_ocr_sections.txt` (2063 lines, committed to the archive, NOT to git yet). | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:19 |
| U14 | [unquoted — verify on next use] **(b) The thesis body** (Gray 1993 p. 252; archive L1612-1618): | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:56 |
| U15 | [unquoted — verify on next use] - `s_0 = 0.02` on SiO2 is the **site-limited F adsorption coefficient** that `RESULTS_ION_CHANNEL_SOLVE` §3 inverted from Gray's half-rise as `s ~ 0.06` and recorded `[VERIFY]`. Gray's own value is **0.02**, and the regressed column shows it rising with ion energy (0.013 -> 0.045 over 20-500 eV). The inverted 0.06 is the right order; the source number is **0.02-0.045**. This closes the C1b `[VERIFY]` on the low side and it is *smaller* than inverted, i.e. the knee should sit even further right. - `b` is exactly `0.007 * sqrt(E)` at all six energies (0.031/sqrt20 = 0.00693; 0.313/sqrt2000 = 0.00700). The auxiliary parameter is itself sqrt-scaled. | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:129 |
| U16 | [unquoted — verify on next use] The supporting data set is Figure 5-2 (p. 162), plotted twice — panel (a) yield vs `E` over **0-1200 eV**, panel (b) yield vs `E^0.5` over `0-35` (i.e. to 1225 eV) — with four sources: Gray's own 1992 data, Chapman 1980, Oostra 1986, Steinbruchel 1989. Its caption: | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:157 |
| U17 | [unquoted — verify on next use] **Constant-`s_0` column (0.13, 0.55, 0.60, 0.85, 1.10, 2.25) — Gray's own preferred fit:** | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:232 |
| U18 | [unquoted — verify on next use] Honest reading: on the raw column a free-parameter linear fit is nearly as good as sqrt (0.944 vs 0.954), because the 2 keV raw point (2.845) sits above the sqrt curve. On Gray's own preferred column the sqrt form wins decisively (0.994 vs 0.744). **In both columns the two forms petch actually runs — the Appendix-B linear anchored at its own threshold, and the ZBL-deposited-in-layer shape — score negative or near-zero R^2**, because both predict roughly 5x at 2 keV where the measurement is 2.25-2.85. | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:253 |
| U19 | [unquoted — verify on next use] \| **Steinbruchel sqrt, Gray's own sputter threshold** \| **18** \| **3.74** \| | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:389 |
| U20 | [unquoted — verify on next use] \| **Steinbruchel sqrt, Gray's own chem threshold** \| **4** \| **3.37** \| | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:390 |
| U21 | [unquoted — verify on next use] 1. `RESULTS_ABSOLUTE_YIELD` reads the Gray/Kwon floor as **0.28** at 350 eV. Gray's own sputtering model gives **0.201** at 350 eV, so the bare row is **1.69x** too strong against the primary source, not 1.22x against the replot. The 0.28 reading corresponds to `Y_sputter` at ~600 eV under Gray's own constants. Recommend re-reading Kwon Fig. 3.4's stated energy — the Kwon thesis DSpace text layer is unusable OCR, so this stays `[VERIFY]`. 2. The 1.10 plateau anchor is **confirmed and located**: Gray's total saturated SiO2 yield at 350 eV is `beta_e(1+b) + Y_sputter = 0.885*1.131 + 0.201 = 1.20`, and the constant-`s_0` table entry at 500 eV is literally `1.10`. The plateau is real; the 350 vs 500 eV attribution should be tightened. | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:423 |
| U22 | [unquoted — verify on next use] **Under Gray's own measured energy laws, the bare/complex balance at the feature front inverts from 1.04:1 to 0.18:1 — the chemical channel becomes ~5.6x the physical one.** The fork's forecast (0.30:1 under sqrt with Krueger's thresholds) was directionally right and conservative; with Gray's own constants the inversion is stronger. | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:434 |
| U23 | [unquoted — verify on next use] \| Steinbruchel sqrt (Gray's own thresholds), extrapolated \| 3.74 \| 3.37 \| | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:485 |
| … | 25 further unquoted mentions elided — `grep -n` the patterns above across `*.md` | — |

_Harvested 9 quoted + 39 unquoted mentions across the repo's docs._
