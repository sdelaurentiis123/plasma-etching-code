# Mahoney 1994 independent argon ICP board

**Verdict: FAIL**

The unchanged five-reaction argon closure was run on the chamber geometry and
all five 100 W operating rows transcribed from Mahoney et al. Table I. No
Mahoney temperature, density, flux, or profile selected a coefficient.

| gas T (K) | ion wall energy (Te) | Te MAPE | Te max APE | model/measured density | density-shape log RMSE | max residual | verdict |
|---:|---:|---:|---:|---:|---:|---:|:---:|
| 300 | 5 | 36.96% | 44.38% | 8.00–15.85 | 0.333 | 6.31e-15 | FAIL |
| 300 | 8 | 36.96% | 44.38% | 6.73–13.64 | 0.351 | 3.05e-15 | FAIL |
| 600 | 5 | 29.07% | 36.96% | 4.35–8.61 | 0.325 | 3.66e-15 | FAIL |
| 600 | 8 | 29.07% | 36.96% | 3.62–7.35 | 0.345 | 2.65e-15 | FAIL |

Frozen limits were Te MAPE <= 30%, Te maximum
APE <= 50%, model/measured peak-density
ratio inside [1,
5], normalized density-shape log RMSE
<= ln(2), the source trends, positive finite state variables, and normalized
balance residual <= 1e-08. Both neutral-temperature sensitivity
endpoints and both published ion-wall-energy endpoints had to pass.

## Claim boundary

Mahoney reports net generator power after reflected-power subtraction, not
calorimetric absorbed plasma power, and explicitly notes coil/matching-network
losses. This board therefore runs 100 W as an all-net-power-absorbed
upper-bound scenario. It cannot by itself validate absolute plasma density or
a net-to-absorbed transfer closure.

The paper also states that its electron-density diagnostic can read two to five
times below ion-density determinations. That source-declared interval is the
absolute-density comparison boundary; it was frozen before the run.

This is an independent plasma-state test of the argon chemistry/transport
closure. It does not validate C4F6 chemistry, a sheath IEAD, or Krueger depth.
