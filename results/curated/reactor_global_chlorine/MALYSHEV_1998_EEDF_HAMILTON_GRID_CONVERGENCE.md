# Malyshev Hamilton EEPF grid convergence

## Verdict

The held-out 500 W, 30%-absorbed-power sensitivity is **numerically converged**
from 415 to 813 actual threshold-aligned cells. The largest physical-output
change is `electron_density_m3` at
`0.0328%`, below the
preregistered `0.05%` receipt threshold.

| metric | 415 cells | 813 cells | relative change |
|---|---:|---:|---:|
| reduced_electric_field_Td | 1091.18039 | 1091.28396 | +0.0095% |
| mean_electron_energy_eV | 3.81860433 | 3.81876983 | +0.0043% |
| electron_density_m3 | 4.43078414e+16 | 4.42933244e+16 | -0.0328% |
| electronegativity | 1.53062211 | 1.53106298 | +0.0288% |
| cl_to_cl2_ratio | 0.78148106 | 0.781399101 | -0.0105% |
| modeled_relative_cl2_density_percent_proxy | 56.1330694 | 56.135652 | +0.0046% |
| total_positive_ion_axial_flux_m2_s | 4.87972348e+19 | 4.87835347e+19 | -0.0281% |
| clplus_ion_fraction | 0.322192437 | 0.322131642 | -0.0189% |

This is a numerical receipt only. It does not validate the reactor state,
wafer flux, or feature depth, and no feature observable selected either grid.
The raw collision decks remain user-supplied and uncommitted; their identities
are hash-gated in the JSON receipt.
