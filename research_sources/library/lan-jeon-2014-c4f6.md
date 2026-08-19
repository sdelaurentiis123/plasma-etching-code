# lan-jeon-2014-c4f6

**Swarm-inverted electron collision set for C4F6**

- **Citation:** P.-T. Lan and B.-H. Jeon, "Determination of the Electron
  Collision Cross-section Set for the C4F6 Molecule by Using an Electron
  Swarm Study," *Journal of the Korean Physical Society* **64**, 1320-1326
  (2014).
- **DOI:** [10.3938/jkps.64.1320](https://doi.org/10.3938/jkps.64.1320)
- **Local source PDF SHA-256:**
  `82d672f5b60611894a4584aa503e4c26d66aec3c66792550de46fb660f77aeb6`
- **Extract:**
  `research_sources/thesis_extracts/lan_jeon_2014_c4f6_verified_excerpt.txt`
- **Status:** PRIMARY FULL TEXT; TABLES 1-2 TRANSCRIBED AND
  CONSERVATION-CHECKED; FIGURE 7 PIL-AUDITED AT 600 DPI
- **Topic:** C4F6 electron collisions and pulsed-Townsend swarm inversion

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | Tables 1-2 print 37 momentum-transfer points and nine inelastic processes in units of `1e-16 cm2`. | The committed SI deck checksum-locks every printed row and unit conversion. |
| Q2 | The set was assembled from measured, borrowed, and trial-and-error-adjusted curves so two-term Boltzmann drift agreed with pure-C4F6 and C4F6/Ar PT measurements. The stated acceptable deviation is `+/-5%` over `0.35-1200 Td`. | This is a swarm-inverted working set. Its fitted rows are not independent validation data. |
| Q3 | The source explicitly identifies the comparison observable as average electron drift velocity `Wv` obtained by the pulsed-Townsend method. | `Wv` is not silently relabeled as flux or universal bulk drift; [`casey-2021-pt-foundations`](casey-2021-pt-foundations.md) governs the definition boundary. |
| Q4 | Figure 7 prints the pure-C4F6 markers; the archived seven-page source contains Figures 1-7 although the text refers to Figure 8. | Only visually recoverable Figure-7 markers are executable. The missing figure remains a source anomaly. |

## Executable decision

`load_lan_jeon_2014_c4f6_replay()` consumes the exact Table-1/2 arrays. The
official BOLSIG+ comparator independently reproduces petch flux drift to
`0.234%` without retuning, which corroborates the local two-term numerics. The
deck does not identify a reactor state, ion product branching, wafer flux,
surface response, or depth.
