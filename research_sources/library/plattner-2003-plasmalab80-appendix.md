# plattner-2003-plasmalab80-appendix

**Oxford PlasmaLab RIE 80 CHF3/Ar process appendix with measured self-bias**

- **Citation:** Luca Plattner, *A study in biomimetics: nanometer-scale,
  high-efficiency, dielectric diffractive structures on the wings of butterflies
  and in the silicon chip factory*, PhD thesis, University of Southampton
  (2003), Appendix B.
- **Official record:** `https://eprints.soton.ac.uk/260031/`
- **Official appendix:** `https://eprints.soton.ac.uk/260031/5/AppendixB.pdf`
- **Appendix PDF SHA256:**
  `b42b9a357d347cc203ff1357fa6e420e73595c413484464b5905312bddbb37fc`
- **Status:** PRIMARY THESIS APPENDIX; TABLE B.5 VISUALLY AUDITED

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| L1 | Table B.5 identifies an Oxford PlasmaLab RIE 80 at `13.56 MHz`. | Establishes the equipment family. |
| L2 | The tabulated condition is `6 sccm CHF3`, `7.5 sccm Ar`, `25 mTorr`, `200 W`, and `-400 V` DC bias. | A high-voltage family witness demonstrating that the finite values in the Oxford-80 corpus extend to `400 V`; chemistry and loading differ from the target. |

## Executable decision

The magnitude `400 V` is carried as a sensitivity history. Its sign is stored
by the module-level convention (magnitude of negative electrode self-bias), not
lost or confused with a positive electrode potential.
