# kawaguchi-2020-cl2

**Updated Cl2 electron-collision set and transport analysis**

- **Citation:** S. Kawaguchi, K. Takahashi, and K. Satoh, “Electron collision
  cross section set of Cl2 gas and electron transport analysis in Cl2 gas and
  Cl2/N2 mixtures,” *Japanese Journal of Applied Physics* **59**, SHHA09
  (2020).
- **DOI:** `10.35848/1347-4065/ab72ce`
- **Local publisher PDF SHA-256:**
  `624863578c10a1f0a5d06f3a1bc7a77344c009b27d8c6183899ffcd1fefd6f30`
- **Extract:**
  `research_sources/thesis_extracts/kawaguchi_2020_cl2_cross_sections.txt`
- **Status:** PRIMARY FULL TEXT; TABLE I AND FIGURES 1--12 VISUALLY AUDITED;
  NATIVE NUMERIC CROSS-SECTION ARRAYS NOT SUPPLIED
- **Topic:** updated molecular-chlorine cross sections and swarm validation

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | The constructed set contains one elastic momentum-transfer, three vibrational, eight dissociative-electronic, two electronic, one ion-pair, one dissociative-attachment, and four ionization cross sections. | This is a useful modern transport-calibrated candidate deck, not a uniformly measured or ab-initio deck. |
| Q2 | The Gregorio--Pitchford momentum-transfer cross section is retained below `1 eV`; its shape above `1 eV` is modified to reproduce a momentum-transfer curve obtained by integrating Gote--Ehrhardt elastic differential measurements. | The elastic row has mixed calculated/experiment-constrained provenance. It is not a direct momentum-transfer measurement. |
| Q3 | Hamilton dissociative-excitation cross sections are multiplied by `2.5`, and Rescigno Rydberg-excitation cross sections by `4.0`, to reproduce measured ionization coefficients. | These large, explicit swarm-calibration factors must remain visible and cannot silently replace the evaluated physical deck. |
| Q4 | Below `0.2 eV`, the calculated dissociative-attachment cross section is multiplied by `1.05` to reproduce a measured attachment-rate coefficient; above `0.2 eV`, its shape follows measured attachment data. | Attachment is experiment-constrained but still assembled from different evidence regimes. |
| Q5 | Basner--Becker measured ionization cross sections are reduced by `15%`, which the authors state is within their experimental uncertainty, to reduce the calculated ionization coefficient. | The provider should retain the original measurement and the applied scale as separate provenance fields. |
| Q6 | The sum of the eight dissociative-excitation channels only roughly reproduces the measured total neutral-dissociation cross section. | This is a useful aggregate check, not species/state-resolved validation of all eight rows. |
| Q7 | The set reproduces the 2018 drift markers well and the longitudinal-diffusion markers through about `330 Td`, then trends slightly low. The authors decline further modification because pulsed and time-of-flight diffusion definitions can differ and request new TOF data. | Drift and diffusion form the best available board; diffusion above `330 Td` remains an explicit residual, not something to tune away. |
| Q8 | The article infers rate-coefficient reliability from transport agreement over approximately `100--760 Td`, while also stating that no measured Cl2 rate-coefficient data are available for direct comparison. | Transport validation supports the EEDF/collision set conditionally; it is not direct validation of every reaction-rate curve. |
| Q9 | The Monte Carlo treatment calculates mean-arrival-time drift `Wm` and `NDL` from spatiotemporal EVDF development, while `alpha/N`, `eta/N`, and `(alpha-eta)/N` are defined by the steady-state-Townsend EVDF. The rate coefficients use the pulsed-Townsend EVDF. | Flux mobility, scalar diffusion, or one generic EEDF cannot be substituted for every plotted observable in attaching Cl2. The native solver and residual board must carry the experiment-specific transport definition. |
| Q10 | The transport calculation assumes isotropic electron scattering and uniform sharing of residual ionization energy between the scattered and ejected electrons. | A deterministic replay may use these as explicit source-reproduction closures. The momentum-transfer row alone does not supply higher angular collision moments, and the isotropic closure is not differential-scattering validation. |

## Executable decision

No numeric cross section is imported from the plotted curves. The set is
classified as a candidate **swarm-calibrated provider**, kept side by side with
the primary/evaluated deck. Native author arrays are preferred; plot
digitization is a lossy fallback. Any implementation must reproduce the 52-row
2018 board without fitting those same markers and must expose the five stated
modifications above.
