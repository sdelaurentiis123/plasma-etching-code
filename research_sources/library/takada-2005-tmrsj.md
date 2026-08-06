# takada-2005-tmrsj

**Takada, Toyoda & Sugai, stable-parent molecule / ion co-incidence on SiO2**

- **Citation:** N. Takada, H. Toyoda, and H. Sugai, “Evidence of Radical-free
  Etching of SiO2 by Fluorocarbon Molecule under Ion Bombardment,”
  *Transactions of the Materials Research Society of Japan* **30**[1],
  319–322 (2005).
- **Related article:** N. Takada, H. Toyoda, I. Murakami, and H. Sugai,
  *Journal of Applied Physics* **97**, 013534 (2005),
  DOI `10.1063/1.1829400`.
- **Retrieval route:** open publisher PDF,
  `https://www.mrs-j.org/pub/tmrsj/vol30_no1/vol30_no1_319.pdf`
- **Status:** FULL TEXT:
  `research_sources/thesis_extracts/takada_2005_radical_free_etching.txt`;
  archived PDF:
  `research_sources/takada_2005_radical_free_etching.pdf`
- **PDF SHA-256:**
  `9034f445f575b85c0b9c95d79b81699c42eafa47a347a7e73a6de73cbe222e25`
- **Topic:** beam-yields — stable-molecule / ion co-incidence, surface-film
  balance, and the limits of a pure-ion yield ceiling

## Claims table

| # | verbatim source claim | use and boundary |
|---|---|---|
| Q1 | “The sample installed at the beam irradiation part is exposed to Ar+ beam and C5F8 flood simultaneously, in a high vacuum condition.” | Establishes controlled molecule/ion co-incidence. It is not a fluorocarbon plasma. |
| Q2 | “The beam energy (E) of Ar+ is controlled from 100 to 900 eV” and “The Ar+ beam is directed to a sample in parallel to the surface normal.” | Defines the energy and incidence support. |
| Q3 | “The C5F8 beam is directed from an angle of 45° to the sample” and “The flux ratio (R) of C5F8 molecule to Ar+ is varied from 0.25 to 250.” | Defines molecule incidence and the explored ratio domain. |
| Q4 | “SiO2 etching yield of 2.5 is obtained in the condition of E=900 eV and the flux ratio R=1.” | Directly refutes treating the 1.5 CF3+ single-ion value as a universal surface-physics ceiling. It does not identify C4F6 behavior. |
| Q5 | “In the case of C5F8/Ar+ co-incidence, an etching yield of ~1.2 is observed at R~1.” | Text cross-check for the digitized Figure-3 maximum (`1.1969`) at 400 eV. |
| Q6 | “the present result in fig. 3 shows an etching yield of 0.67 at R=0.25 and E=400 eV.” | Text cross-check for the digitized ratio-0.25 point (`0.6697`). |
| Q7 | “The contribution of C5F8 molecule to the SiO2 etching was estimated to be roughly 30-40% under typical plasma condition if we consider the survival C5F8 of a few % in plasma.” | Shows that even a small surviving parent fraction can matter. This is the authors’ C5F8 estimate, not a C4F6 or Krüger flux measurement. |
| Q8 | “The C5F8 molecule has a ring structure with a double bond and four single bonds.” The authors then say the double bond may bind to bombardment-created SiO2 dangling bonds and that this “unique process of C5F8 may cause results different from other fluorocarbon molecules.” | Mandatory prohibition on transplanting the C5F8 curve to C4F6. Use only as an analog mechanism envelope until C4F6 beam data exist. |
| Q9 | The related *Journal of Applied Physics* publisher abstract reports that the C5F8 yield “reaches the value of 2.4 at 900 eV,” whereas this open TMRSJ paper reports `2.5` at the same nominal energy and ratio. | Preserve `2.4–2.5` as a source discrepancy. Do not select `2.5` merely because it nearly equals the Krüger target normalization. |

## Consumed data

- `data/experimental/takada_2005/figure3_sio2_coincidence_yields.csv`
  digitizes the C5F8 and CF2 SiO2 series at 400 eV with retained pixels,
  checksum manifest, and a PIL/NumPy replay script.
- The 900 eV, ratio-1 value is textual, not a Figure-3 digitization. The open
  TMRSJ text says `2.5`; the related JAP publisher abstract says `2.4`.
- No value from this source is a default petch chemistry constant. It bounds a
  missing mechanism and invalidates the former universal-ceiling argument.
