# chang-thesis

**Chang, PhD thesis, MIT (1721.1/50356)**

- **DOI/URL:** MIT DSpace 1721.1/50356
- **Retrieval route:** dspace.mit.edu/bitstreams/5285146e-c906-4f6b-ab72-1b12764d8011/download
- **Status:** FULL TEXT: research_sources/thesis_extracts/chang_thesis.txt; PDF research_sources/chang_thesis.pdf
- **Topic:** beam-yields — Beam-measured yields, thresholds and sticking (the L0/L1 provenance floor)

## Claims table

Rows relocated verbatim from the repo's research/results docs. `consumed by` = the
doc that recorded it (that doc states which constant/decision it fed).

| # | recorded claim (as written in our docs) | consumed by |
|---|---|---|
| Q1 | The polynomial **coefficients are not printed** — not in §5.3, not in Ch. 7 (the DSMC simulator chapter says only "The reaction probability of ions was calculated based on the measured ion angular dependence", `chang_thesis.txt` ~L6950). The Ch. 4 companion prints the *form* only: `Y(φ) = Σ_{i=1..6} a_i cos^i(φ) = c(φ)·Y(φ=0°)` (Eq. 4.1, p. 93) with the `a_i` unlisted. **`a_i` / `c(φ)` closed-form: NOT FOUND** (thesis is the fullest available rendering; the JVST paper itself is AIP-403 and its abstract quotes only the two percentages). | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:177 |
| Q2 | Remaining gap: I have Chang's **verbatim third-party attribution** of the 35 eV to Joubert, not Joubert's own printed number. The Joubert abstract (OpenAlex, verbatim, retrieved) does not contain it — it is a body number. `[VERIFY]` the Joubert 1994 body (AIP paywall; no OA copy found). Recommend re-labelling the deck row **"35 eV — Joubert/Oehrlein/Surendra 1994 via Chang 1998 p. 90"** rather than "Chang–Sawin", which is demonstrably wrong. | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:443 |
| Q3 | \| `RESEARCH_LIP_CERTAINTY_2026-08-04.md` §2.4(b) \| "There is **no measured FC-film angular sputter yield** in the literature" \| **False.** Barklund & Blom, JVST A 10, 1212 (1992), DOI 10.1116/1.578229 — Ar⁺ and CF₄ RIBE on a 5000 Å CHF₃-plasma polymer film. Reprinted as Chang Fig. 4.16, digitised in §3. \| | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:480 |
| Q4 | \| `RESEARCH_BEAM_CONSTANTS_ATLAS_2026-07-29.md` L693, L726 \| "The only 35 eV in Chang's corpus is the Ar⁺-on-Si physical-sputter threshold" \| **Incomplete.** Chang p. 90 also cites *"a 35 eV threshold energy reported for etching of silicon dioxide in a CHF₃ plasma [Joubert, 1994]"* = Joubert, Oehrlein & Surendra, JVST A 12, 665 (1994), DOI 10.1116/1.578850. Re-label the deck row. \| | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:482 |
| U5 | [unquoted — verify on next use] \| **T-1** \| J. P. Chang, *Study of plasma-surface kinetics and simulation of feature profile evolution in chlorine etching of patterned polysilicon*, PhD MIT ChemE, 1998 (adv. Sawin) \| [1721.1/50356](http://hdl.handle.net/1721.1/50356) \| **Good OCR** — tables recoverable \| | `RESEARCH_BEAM_CONSTANTS_ATLAS_2026-07-29.md`:138 |
| U6 | [unquoted — verify on next use] \| B5 substrate sputter bare \| Chang Table 3.1 (A=0.04, E_th=35 eV) + Guo Table 2.1 (Yamamura closed form) \| **L1/L2** \| | `RESEARCH_BEAM_CONSTANTS_ATLAS_2026-07-29.md`:463 |
| U7 | [unquoted — verify on next use] ``` <scratch>/chang1998.pdf   chang1998.txt   # T-1, 1721.1/50356, good OCR, all tables <scratch>/yin2007.pdf     yin2007.txt     # T-4, 1721.1/38973, clean text <scratch>/guo2009.pdf     guo2009.txt     # T-6, 1721.1/46600, clean text, Tables 2.1/2.2/2.3/3.2/4.1/4.2 <scratch>/kanarik2017.pdf kanarik2017.txt # open OSTI copy, purl/1376399 <scratch>/jin2003.pdf     jin2003.txt     # T-2, 1721.1/28357 — OCR UNUSABLE <scratch>/kwon2004.pdf    kwon2004.txt    # T-3, 1721.1/28353 — OCR UNUSABLE ``` | `RESEARCH_BEAM_CONSTANTS_ATLAS_2026-07-29.md`:704 |
| U8 | [unquoted — verify on next use] Chang (MIT ScD 1999, `research_sources/chang_thesis.pdf`) writes the same form explicitly (his Eq. 3.2-3.4, thesis text extract L2266-2280): | `RESEARCH_ENERGY_SCALING_2026-08-05.md`:86 |
| U9 | [unquoted — verify on next use] > **J. P. Chang**, *Study of plasma-surface kinetics and simulation of feature profile > evolution in chlorine etching of patterned polysilicon*, PhD thesis, MIT Chemical > Engineering, 1998 (advisor H. H. Sawin). MIT DSpace handle **1721.1/50356**. > Direct PDF (open, no auth): `https://dspace.mit.edu/bitstreams/5285146e-c906-4f6b-ab72-1b12764d8011/download` > — 21.8 MB, 1.6 dpi-clean scan, `pdftotext` OCR is good. > Local: `scratchpad/hunt/chang_thesis.pdf` + `.txt`. | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:17 |
| U10 | [unquoted — verify on next use] Chang thesis §5.3, **p. 115** (rendered page image `scratchpad/hunt/chang_p115-115.png`; text `chang_thesis.txt` L5691-5706): | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:155 |
| U11 | [unquoted — verify on next use] Source of the CF_x⁺ curve: **T. M. Mayer, R. A. Barker, L. J. Whitman**, *Investigation of plasma etching mechanisms using beams of reactive gas ions*, **J. Vac. Sci. Technol. 18(2), 349 (1981)** (Chang bibliography, `chang_thesis.txt` L7796). | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:289 |
| U12 | [unquoted — verify on next use] **That is false.** Chang thesis **Fig. 4.16, p. 104** (image `scratchpad/hunt/fig416b-104.png`) reprints one, and names it. | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:306 |
| U13 | [unquoted — verify on next use] Source, from Chang's bibliography (`chang_thesis.txt` L7656): | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:309 |
| U14 | [unquoted — verify on next use] **Confirmed, Table 3.1, p. 73** (`chang_thesis.txt` L3109-3131), Ar⁺ on Si physical sputtering, universal (Steinbrüchel) energy dependence: | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:406 |
| U15 | [unquoted — verify on next use] **THE MISSED ONE — Chang thesis §4.1, p. 90** (`chang_thesis.txt` L4125-4130), verbatim: | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:424 |
| U16 | [unquoted — verify on next use] Reference resolved (`chang_thesis.txt` L7761 + Crossref): **O. Joubert, G. S. Oehrlein, M. Surendra**, *Fluorocarbon high density plasma. VI. Reactive ion etching lag model for contact hole silicon dioxide etching in an electron cyclotron resonance plasma*, **J. Vac. Sci. Technol. A 12(3), 665–670 (1994)**, DOI **10.1116/1.578850**. | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:430 |
| U17 | [unquoted — verify on next use] \| `RESEARCH_BEAM_CONSTANTS_ATLAS_2026-07-29.md` §T-1 \| Chang thesis listed as a target \| **Acquired.** `1721.1/50356`, direct bitstream URL in §0. \| | `RESEARCH_VERIFY_HUNT_2026-08-05.md`:483 |
| Q18 | Chapter 3 explicitly calls the surface law a “simplified model,” collapses physisorption, chemisorption, desorption, and ion-induced desorption into one overall adsorption reaction, and assigns the resulting coefficient the name “surface chlorination coefficient, s.” | `src/petch/chang_sawin_chlorine_si.py` |
| Q19 | Chapter 3 Table 3.4 reports, for 100 eV Ar+ with atomic Cl, `Y0=0.07`, `s=0.30`, and `beta=3.59`; for 100 eV Ar+ with Cl2 it reports `Y0=0.07`, `s=0.07`, and `beta=0.83`. The text says beta and s are derived by linear regression of the reciprocal-yield relation. | `src/petch/chang_sawin_chlorine_si.py` |
| Q20 | Chapter 3 Eqs. 3.9--3.11 use the steady site balance `theta=sR/(sR+4 beta)` and the yield `Y=c(phi)[Y0(1-theta)+beta theta]`. The factor four comes from the simplified reaction `Si + 4Cl -> SiCl4 + 4*`. | `src/petch/chang_sawin_chlorine_si.py` |
| Q21 | The thesis states that SiCl4 is assumed as the major low-energy product and that unsaturated SiClx is omitted from the Chapter-3 model; Chapter 5 separately warns that omission of unsaturated chlorides causes deviation at low flux ratio. | `src/petch/chang_sawin_chlorine_si.py` |

_Harvested 8 quoted + 13 unquoted mentions across the repo's docs._
