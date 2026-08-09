# Tian 2017 zonal partial-redistribution audit

## Verdict

Coupled spatial and partial-frequency transport closes the base mechanism without fitting atomic or depth parameters. Because the source target was visible during model selection and one line's emitter field is unpublished, this is a high-accuracy mechanism reproduction, not experimental validation. The boundary ledger is numerically usable, but absolute wafer photon flux additionally requires a line-specific upper-state population field.

| line | observed HPEM | deterministic zonal PFR | error % | complete-redistribution result |
|---|---:|---:|---:|---:|
| Ar_104.8_nm_trapping_factor | 214.463 | 215.277 | 0.38 | 438.289 |
| Ar_106.7_nm_trapping_factor | 360.331 | 336.796 | -6.53 | 264.337 |

Combined two-line MAPE: 3.46%.

## Numerical gate

| line | order 10→12 trapping change | order 10→12 wafer change | frequency refinement change | transition conservation error | terminal conservation error | GMRES residual |
|---|---:|---:|---:|---:|---:|---:|
| Ar_104.8_nm_trapping_factor | 0.345% | 0.065% | 0.00000% | 1.066e-14 | 1.805e-09 | 4.631e-09 |
| Ar_106.7_nm_trapping_factor | 0.917% | 0.523% | 0.00000% | 1.110e-14 | 4.725e-07 | 4.866e-09 |

## Exact blocker

Tian does not publish the Ar(1s2) emitter field or its numerical mesh. A raw source-state field export is required for a blind base-case grade and mixture-dependent spatial fields are required for the 22-point sweep.
