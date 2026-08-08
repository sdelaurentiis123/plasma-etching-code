# Malyshev Hamilton EEPF grid convergence

## Verdict

The held-out 500 W, 30%-absorbed-power sensitivity is **numerically converged**
from 415 to 813 actual threshold-aligned cells. The largest physical-output
change is `electron_density_m3` at
`0.0348%`, below the
preregistered `0.05%` receipt threshold.

| metric | 415 cells | 813 cells | relative change |
|---|---:|---:|---:|
| reduced_electric_field_Td | 880.54694 | 880.679829 | +0.0151% |
| mean_electron_energy_eV | 3.72110899 | 3.72130079 | +0.0052% |
| electron_density_m3 | 4.38987661e+16 | 4.38834686e+16 | -0.0348% |
| electronegativity | 1.6262401 | 1.6267436 | +0.0310% |
| cl_to_cl2_ratio | 0.934488625 | 0.934290256 | -0.0212% |
| modeled_relative_cl2_density_percent_eq11 | 68.1549754 | 68.1595829 | +0.0068% |
| total_positive_ion_axial_flux_m2_s | 4.79244352e+19 | 4.79104639e+19 | -0.0292% |
| clplus_ion_fraction | 0.35052793 | 0.35045289 | -0.0214% |

This is a numerical receipt only. It does not validate the reactor state,
wafer flux, or feature depth, and no feature observable selected either grid.
The raw collision decks remain user-supplied and uncommitted; their identities
are hash-gated in the JSON receipt.
