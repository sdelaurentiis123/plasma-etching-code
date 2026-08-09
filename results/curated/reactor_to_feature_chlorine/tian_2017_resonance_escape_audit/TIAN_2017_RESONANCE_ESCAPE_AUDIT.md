# Tian 2017 deterministic resonance-escape audit

## Verdict

A homogeneous complete-redistribution cylinder is not the next predictive rung. It misses the two-line trapping board and, over part of the mixture sweep, predicts the opposite ordering of 104.8- and 106.7-nm trapping. Tian's own fields explain why: the gas spans 300--833 K and the two radiating states have different, mixture-dependent spatial support. The deterministic transport kernel is viable, but it must consume spatial moments from a two-zone or axisymmetric reactor rather than one homogeneous 0-D state. No depth parameter was touched.

All 22 digitized trapping markers were held out; no target selected a parameter. The best of the four declared homogeneous-temperature sensitivities was 833 K at 38.35% combined MAPE. This is a failed reduced-model gate, not experimental validation.

| uniform gas temperature K | combined MAPE % | 104.8 MAPE % | 106.7 MAPE % | opposite line-order points |
|---:|---:|---:|---:|---:|
| 300 | 93.79 | 120.52 | 67.05 | 7 |
| 400 | 60.94 | 85.17 | 36.72 | 7 |
| 600 | 40.05 | 41.31 | 38.79 | 7 |
| 833 | 38.35 | 24.56 | 52.13 | 7 |

## Required next state

- spatial gas-temperature and ground-state-density field
- line-specific Ar(1s2) and Ar(1s4) emitter spatial moments
- partial frequency redistribution after resonant absorption
- state-specific collisional broadening and quenching

The appropriate deterministic replacement for Tian's photon Monte Carlo is fixed-quadrature ray/frequency transport over those spatial moments. A scalar escape factor or fitted broadband photon yield is not supported.
