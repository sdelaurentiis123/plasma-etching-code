# chantry-1987-wall-diffusion

**Partial-reflection diffusion boundary and low-density limit**

- **Citation:** P. J. Chantry, “A simple formula for diffusion calculations
  involving wall reflection and low density,” *Journal of Applied Physics*
  **62**, 1141–1148 (1987).
- **DOI:** `10.1063/1.339662`
- **Publisher record:**
  `https://pubs.aip.org/aip/jap/article/62/4/1141/18611878`
- **Status:** PRIMARY PUBLISHER ABSTRACT AND EQUATIONS READ; FULL PDF REQUESTED
- **Local extraction:**
  `research_sources/thesis_extracts/chantry_1987_wall_diffusion.txt`

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| C1 | A zero-density wall boundary can seriously misstate radical loss when the diffusion mean free path is comparable to the chamber or particles reflect. | Forbids using `D/Lambda0^2` as a universal neutral wall loss. |
| C2 | With reflection coefficient `R`, the linear extrapolation length is `lambda = 2(1+R) lambda_m / [3(1-R)]`. | With loss probability `gamma = 1-R` and `D = lambda_m vbar/3`, this independently reduces to `lambda = 2(2-gamma)D/(gamma vbar)`. |
| C3 | The exact diffusion lengths come from transcendental eigenvalue equations for each container shape. | The implementation solves the cylindrical radial and axial Robin roots directly. |
| C4 | Chantry's approximation `Lambda^2 = Lambda0^2 + (V/A) lambda` has the correct large- and small-`lambda` limits; its worst discrepancy across the tested shapes is 11%. | Retained as a diagnostic and a regression envelope, not used as the returned exact cylinder loss. |
| C5 | In the collisionless limit the predicted loss frequency agrees with the mean-chord result. | Supplies the surface-reaction asymptote `h A/V`, with `h = gamma vbar/[2(2-gamma)]`. |

## Executable decision

`neutral_transport.solve_cylindrical_neutral_wall_loss` solves the fundamental
cylindrical Robin eigenmode and separately reports the Chantry approximation,
absorbing-wall limit, and surface-reaction limit. The exact radial Bessel zero
is retained; rounded `2.405` is not used for certification.

The publisher abstract is sufficient to verify the extrapolation-length law,
limits, and approximation error. The exact cylinder root equations are an
independent separation-of-variables derivation and are tested by substitution
into both boundary equations. Retrieve and visually audit the full paper before
promoting any further shape-specific formula from this source.
