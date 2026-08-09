# Tian 2017 partial-frequency-redistribution audit

## Verdict

Partial frequency redistribution is necessary and materially improves the held-out two-line board, but a one-temperature homogeneous reactor is still insufficient. The base-case split moment shows that the next closure is nonlocal spatial coupling, not a fitted photon multiplier. No feature-depth parameter was touched.

All 22 mixture-sweep trapping markers stayed held out. No target selected an atomic or collision parameter.

| uniform gas temperature K | combined MAPE % | 104.8 MAPE % | 106.7 MAPE % | opposite line-order points |
|---:|---:|---:|---:|---:|
| 300 | 59.53 | 20.75 | 98.31 | 3 |
| 400 | 41.85 | 16.76 | 66.93 | 2 |
| 600 | 25.37 | 12.11 | 38.64 | 1 |
| 833 | 26.93 | 15.16 | 38.70 | 3 |

Best homogeneous PFR MAPE: 25.37% versus 38.35% for complete redistribution.

## Split spatial-moment mechanism diagnostic

At the printed 20% Cl2 base condition, a 300 K cold-absorber density moment plus 600 K source-zone Doppler moment gives 0.84% MAPE across the two lines. The source target was visible during this model-selection step, so this is mechanism reproduction, not validation.

## Numerical receipt

| line | grid change | half-range change | GMRES residual |
|---|---:|---:|---:|
| Ar_104.8_nm_trapping_factor | 0.00001% | 0.00117% | 4.821e-09 |
| Ar_106.7_nm_trapping_factor | 0.00000% | 0.00197% | 3.945e-09 |

## Required next state

- source-derived axisymmetric gas-temperature and Ar-density field
- line-specific Ar(1s2) and Ar(1s4) emitter spatial moments
- source-grade velocity-changing and quenching collision rates
