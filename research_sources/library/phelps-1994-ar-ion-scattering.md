# phelps-1994-ar-ion-scattering

**Phelps, consistent Ar+-Ar scattering cross sections for discharge sheaths**

- **Citation:** A. V. Phelps, “The application of scattering cross sections
  to ion flux models in discharge sheaths,” *Journal of Applied Physics* 76,
  747–753 (1994).
- **DOI:** `10.1063/1.357820`
- **Primary record:** `https://doi.org/10.1063/1.357820`
- **Status:** PRIMARY PUBLISHER RECORD + ABSTRACT READ; EQUATION
  CROSS-CHECKED AGAINST OPEN TECHNICAL REPRODUCTIONS
- **Topic:** Ar+ elastic/backward scattering, momentum-transfer cross section,
  ion-neutral sheath transport

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | Phelps constructs a mutually consistent Ar+-Ar scattering set rather than treating elastic scattering and symmetric charge transfer as unrelated collisions. | Prevents double counting charge transfer and elastic momentum loss. |
| Q2 | The suggested analytic momentum-transfer law is `Qm(E) = 1.15e-18 E^-0.1 (1 + 0.015/E)^0.6 m2`, with `E` the center-of-mass collision energy in eV. | Landed in `petch.reactor_global.transport`; the energy frame is explicit and tested. |
| Q3 | The scattering set is intended to improve calculated angular, energy, and time distributions of argon ions crossing discharge sheaths. | Supports the ion-neutral transport closure, not a surface-yield or reactor-density claim. |

## Use decision

The analytic momentum-transfer law is imported without digitizing or
redistributing the LXCat Phelps dataset. The original publisher article is
paywalled; its equation was independently cross-checked in open technical
implementations and in the public LXCat query interface. The closure remains
`published_model` until its reactor-level outputs pass an independent board.
