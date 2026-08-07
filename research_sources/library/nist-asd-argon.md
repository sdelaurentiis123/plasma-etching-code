# nist-asd-argon

**NIST Atomic Spectra Database, neutral argon energy authority**

- **Citation:** A. Kramida, Yu. Ralchenko, J. Reader, and NIST ASD Team,
  *NIST Atomic Spectra Database*, version 5.12 (2024), National Institute of
  Standards and Technology.
- **Database DOI:** `10.18434/T4W30F`
- **Ionization-energy query:**
  `https://physics.nist.gov/cgi-bin/ASD/ie.pl?at_num_out=1&biblio=1&e_out=0&el_name_out=1&ion_charge_out=1&spectra=argon&unc_out=1&units=1`
- **Ar I level query:**
  `https://physics.nist.gov/cgi-bin/ASD/energy1.pl?biblio=on&conf_out=on&format=0&g_out=on&j_out=on&lande_out=on&level_out=on&output=0&page_size=15&perc_out=on&spectrum=Ar+I&submit=Retrieve+Data&term_out=on&unc_out=1&units=1`
- **Status:** PRIMARY NIST DATABASE QUERIES
- **Topic:** physical event energies for the argon reactor power ledger

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | Neutral argon’s first ionization energy is `15.7596119 ± 0.0000005 eV`. | Physical energy removed by ground-state ionization; not the `18.68 eV` exponent in Lee–Lieberman’s rate fit. |
| Q2 | The Ar I `3s2 3p5(2P°3/2)4s`, `J=2` level lies at `11.54835442 eV`. | Representative 4s metastable energy for excitation/superelastic bookkeeping in the Lee–Lieberman lumped `Ar*` deck. |
| Q3 | The corresponding metastable-to-ionization-limit energy is `4.21125748 eV`. | Step-ionization physical threshold by energy difference; not the `4.95 eV` rate-fit exponent. |

## Use decision

The database values prevent a category error in the reactor power balance:
fit exponents control the temperature dependence of rate coefficients, while
level/ionization energies determine energy removed per event. The lumped
`Ar*` state remains a reduced approximation because the 4s manifold contains
multiple levels; this value follows the metastable state represented by the
Lee–Lieberman deck and is not a full collisional-radiative model.
