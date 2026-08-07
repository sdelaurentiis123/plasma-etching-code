# Mahoney 1994 independent argon ICP board

**Verdict: FAIL**

The unchanged five-reaction argon closure was run on the chamber geometry and
all five 100 W operating rows transcribed from Mahoney et al. Table I. No
Mahoney temperature, density, flux, or profile selected a coefficient.

| gas T (K) | ion wall energy (Te) | Te MAPE | Te max APE | model/measured density | density-shape log RMSE | max residual | verdict |
|---:|---:|---:|---:|---:|---:|---:|:---:|
| 300 | 5 | 37.78% | 45.10% | 9.21–19.09 | 0.368 | 3.37e-15 | FAIL |
| 300 | 8 | 37.78% | 45.10% | 7.67–16.45 | 0.387 | 4.88e-15 | FAIL |
| 600 | 5 | 29.86% | 37.67% | 4.86–9.92 | 0.350 | 3.46e-15 | FAIL |
| 600 | 8 | 29.86% | 37.67% | 4.05–8.48 | 0.371 | 2.85e-15 | FAIL |

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
