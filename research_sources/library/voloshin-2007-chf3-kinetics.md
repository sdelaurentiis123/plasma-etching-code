# voloshin-2007-chf3-kinetics

**Validated gas-phase kinetics in a CHF3 discharge**

- **Citation:** D. G. Voloshin, A. M. Efremov, and K.-H. Kwon, "On the
  Reaction Kinetics in a CHF3 Plasma," *IEEE Transactions on Plasma Science*
  **35** (2007).
- **DOI:** `10.1109/TPS.2007.906780`
- **Author-posted full text read online:**
  `https://www.researchgate.net/publication/3454720_On_the_Reaction_Kinetics_in_a_CHF3_Plasma`
- **Status:** PRIMARY AUTHOR-POSTED FULL TEXT READ ONLINE; TABLE III PIXEL
  TRANSCRIPTION OPEN

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| V1 | The modeled neutral chain includes `CHF3 + F -> CF3 + HF`, electron-impact HF dissociation that regenerates F and H, and H-driven conversion of CFx into CFx-1 plus HF. | This topology is necessary to avoid treating feed fractions as steady plasma composition. |
| V2 | The article text gives `1.82e-12 cm3/s` for `CHF3 + F -> CF3 + HF`. | This one text-resolved constant may be sensitivity-tested; the remainder of Table III is not imported until the table image is independently transcribed. |
| V3 | Increasing an assumed radical wall-loss probability from `0.05` to `0.1` changes calculated radical densities by roughly 15%. | Wall recombination is an explicit uncertainty direction, not a hidden fit knob. |
| V4 | Validation used measured electron density and electron temperature as model inputs because generator power did not uniquely close the electron state. | The paper validates chemistry at a diagnosed plasma state; it does not justify converting the Zhu `150 W` forward demand directly into electron density or field. |

## Executable decision

The source currently supports topology and a declared wall-loss sensitivity.
It does not yet supply a landed full numeric mechanism because its Table III
was available only as page imagery in the accessible copy.  Sandia Table 9 is
the checksum-pinned executable CHF3 extension; any later Voloshin expansion
must preserve the distinction between measured electron-state inputs and a
closed equipment power balance.
