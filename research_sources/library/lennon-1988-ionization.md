# lennon-1988-ionization

**Lennon et al., evaluated electron-impact ionization of atoms and ions**

- **Citation:** M. A. Lennon, K. L. Bell, H. B. Gilbody, J. G. Hughes,
  A. E. Kingston, M. J. Murray, and F. J. Smith, “Recommended Data on the
  Electron Impact Ionization of Atoms and Ions: Fluorine to Nickel,”
  *Journal of Physical and Chemical Reference Data* 17, 1285--1363 (1988).
- **DOI:** `10.1063/1.555809`
- **Official NIST full text:**
  `https://srd.nist.gov/jpcrdreprint/1.555809.pdf`
- **PDF SHA-256:**
  `ef4ecdc0522a3cbd0608d8af9227d60261b5ff69e53ff4e6735bb2e7cc8d6d76`
- **Status:** PRIMARY NIST FULL TEXT READ + EQ. 5--6/TABLES 1--2 PIXEL-AUDITED
- **Topic:** evaluated atomic-chlorine ionization cross section, Maxwellian
  rate coefficient, analytic fit, and fit domain

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | Equation 5 defines the rate coefficient as the ionization cross section times electron velocity, averaged over a Maxwellian velocity distribution. | Establishes what the fitted rate represents; it is not an arbitrary Arrhenius coefficient. |
| Q2 | For `I/10 < kT < 10 I`, Eq. 6 is `exp(-I/kT) (kT/I)^(1/2) sum(a_n [log10(kT/I)]^n)`. | Resolves the ambiguous `log` printed in Lee--Lieberman Table 2 as base 10 and supplies a hard validity domain. |
| Q3 | Table 2 supplies `a0...a5` in `cm3 s^-1`; the Cl I row begins `1.4193e-7`, while Lee--Lieberman rounds that coefficient to `1.419e-7`. | The reactor deck preserves Lee's printed rounded row for faithful reproduction but inherits the primary functional form and domain from Lennon. |
| Q4 | The neutral-chlorine ionization potential used in the fit is approximately 12.96 eV. | Physical ionization energy for the electron-power ledger; it is also the fit's reference energy. |
| Q5 | The paper evaluates and recommends underlying cross-section data, then fits the recommended rates for convenient use. | Evidence classification is evaluated/regressed primary data, not ab initio prediction. |
| Q6 | The Cl I Table-2 coefficients are `1.4193e-7`, `-1.8637e-8`, `-5.4395e-8`, `3.3056e-8`, `-3.5393e-9`, and `-2.9152e-8 cm3 s-1`. | Lee's values are rounded transcriptions; the rounding is far too small to explain the measured-table/rate-fit difference. |
| Q7 | Section 6.17 explicitly recommends the Hayes et al. experimental Cl I data; Table 1 assigns the Cl I cross-section fit a 15% estimated error. | The later NIST Table-25 ±14% measurement and the Lennon analytic law are related evidence, not independent laboratories. They must not be counted twice. |

## Use decision

This source closes one real transcription ambiguity before the chlorine solver
is built. The executable law is base-10, not natural-log, and rejects
temperatures outside the published Eq.-6 domain. Lee's rounded coefficients
remain source-faithful; a later sensitivity can quantify whether replacing
them with Lennon's higher-precision row is material.
