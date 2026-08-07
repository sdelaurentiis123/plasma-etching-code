# Lee--Lieberman Figure 3 argon reproduction

**Verdict: PASS**

This is a no-fit reproduction of the published global-model curve, not an
independent reactor validation. The 18 reference points were digitized before
running this grade. Neither a target electron temperature, reactor density,
ion flux, etch rate, nor feature depth selected any coefficient.

| ion wall energy member | MAPE | maximum APE | monotonic | maximum balance residual | verdict |
|---:|---:|---:|:---:|---:|:---:|
| 5.0 | 9.206% | 14.669% | yes | 1.107e-14 | PASS |
| 8.0 | 9.210% | 14.669% | yes | 1.437e-14 | PASS |

Frozen limits were MAPE <= 10%, maximum
APE <= 20%, normalized particle/power
residual <= 1e-08, strictly decreasing electron temperature
with pressure, and positive finite densities/fluxes. Both endpoints of the
source's published 5--8 Te ion-wall-energy range had to pass.

## Claim boundary

The transport closure is still labeled `published_model`: Phelps' Ar+-Ar law
is source-backed and energy-averaged, but the NIST Ar-in-Ar self-diffusion
correlation is extrapolated from its stated 418 K ceiling to Lee and
Lieberman's 600 K condition. Independent, condition-specific reactor
measurements are the next gate. No Krueger depth result is fitted or changed
by this board.
