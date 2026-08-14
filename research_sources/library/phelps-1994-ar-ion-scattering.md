# phelps-1994-ar-ion-scattering

**Phelps, consistent Ar+-Ar scattering cross sections for discharge sheaths**

- **Citation:** A. V. Phelps, “The application of scattering cross sections
  to ion flux models in discharge sheaths,” *Journal of Applied Physics* 76,
  747–753 (1994).
- **DOI:** `10.1063/1.357820`
- **Primary record:** `https://doi.org/10.1063/1.357820`
- **Status:** PRIMARY PUBLISHER RECORD + PUBLIC LXCAT PHELPS PROCESS RECORDS
  3075/3076 READ 2026-08-13
- **Topic:** Ar+ elastic/backward scattering, momentum-transfer cross section,
  ion-neutral sheath transport

## Claims table

| # | verified source claim | use and boundary |
|---|---|---|
| Q1 | Phelps constructs a mutually consistent Ar+-Ar scattering set rather than treating elastic scattering and symmetric charge transfer as unrelated collisions. | Prevents double counting charge transfer and elastic momentum loss. |
| Q2 | The suggested analytic momentum-transfer law is `Qm(E_lab) = 1.15e-18 E_lab^-0.1 (1 + 0.015/E_lab)^0.6 m2`. The public LXCat record writes the energy coordinate as center-of-mass energy and evaluates this law at `2 E_cm`, which is the projectile laboratory energy for equal masses. | Landed in `petch.reactor_global.transport`; the equal-mass frame conversion is explicit and tested. |
| Q3 | The scattering set is intended to improve calculated angular, energy, and time distributions of argon ions crossing discharge sheaths. | Supports the ion-neutral transport closure, not a surface-yield or reactor-density claim. |
| Q4 | LXCat's Phelps Ar+/Ar records split the cross section into an isotropic component `Qi = 2e-19 / (sqrt(E_lab)(1+E_lab)) + 3e-19 E_lab / (1+E_lab/3)^2.3` and a backscatter component `Qb = (Qm-Qi)/2`, so `Qm = Qi + 2Qb` while the physical collision-event cross section is `Qi + Qb`. | Below 400 eV the deterministic sheath operator samples the isotropic COM branch and an exact charge-label-swap/backscatter branch with these energy-dependent weights. |

## Use decision

The analytic laws are transcribed with attribution from the public interactive
LXCat Phelps records, not redistributed as an LXCat data file. The original
publisher article is paywalled; the equation and component identity were
cross-checked in open technical reproductions and the public LXCat query
interface. The closure remains `published_model` until reactor outputs pass an
independent measured IEAD board.
