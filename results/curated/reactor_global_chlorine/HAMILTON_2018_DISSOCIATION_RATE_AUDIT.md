# Hamilton 2018 Cl2 dissociation-rate audit

**Verdict: PASS**

The analytic Maxwellian integrator consumed all 50,000 source points for each
of eight dissociative excited states and was compared against the authors'
independently supplied Figure-5 `x=1` total-rate array at
235 temperatures from
0.3 to 5.0 eV.

- maximum absolute relative error:
  `0.4442%`
- mean absolute relative error:
  `0.1855%`
- engineering reproduction limit: `1.0000%`
- worst temperature:
  `0.310150 eV`

This is numerical source reproduction, not experimental validation and not a
rate fit. No reactor density, ion flux, etch rate, or feature depth selected a
coefficient.

## Physics gained

Each state retains its Table-2 vertical excitation energy, so particle
production and the electron-energy sink can be summed consistently instead
of treating Lee's `3.824 eV` Arrhenius exponent as a physical event energy.
All eight retained states dissociate to two ground-state Cl atoms.

## Boundary retained

The source uses fixed-nuclei R-matrix calculations for Cl2(v=0), then
transition-specific high-energy scaling. Hamilton et al. explicitly say
Cosby's vibrationally distributed experiment is not directly comparable and
publish no scalar uncertainty. The resulting rates are therefore
`semi_empirical`, Maxwellian-only, and fail outside the paper's stated
industrial 0.3--5 eV domain.
