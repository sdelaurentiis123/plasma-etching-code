# lennon-1988-ionization

**Lennon et al., evaluated electron-impact ionization of atoms and ions**

- **Citation:** M. A. Lennon, K. L. Bell, H. B. Gilbody, J. G. Hughes,
  A. E. Kingston, M. J. Murray, and F. J. Smith, “Recommended Data on the
  Electron Impact Ionization of Atoms and Ions: Fluorine to Nickel,”
  *Journal of Physical and Chemical Reference Data* 17, 1285--1363 (1988).
- **DOI:** `10.1063/1.555809`
- **Official NIST full text:**
  `https://srd.nist.gov/jpcrdreprint/1.555809.pdf`
- **Status:** PRIMARY NIST FULL TEXT READ + EQ. 6/TABLE 2 AUDITED
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

## Use decision

This source closes one real transcription ambiguity before the chlorine solver
is built. The executable law is base-10, not natural-log, and rejects
temperatures outside the published Eq.-6 domain. Lee's rounded coefficients
remain source-faithful; a later sensitivity can quantify whether replacing
them with Lennon's higher-precision row is material.
