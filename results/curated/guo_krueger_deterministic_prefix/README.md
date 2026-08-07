# Guo/Kwon–Krüger deterministic-extruded prefix refinement

Status: **10 nm time refinement passed on the 0.015625/0.0078125 s
twofold pair. Spatial refinement is now authorized under the committed
addendum.**

The candidate moving-profile neutral-exchange authority is
`deterministic_extruded_2d`, appropriate only because Krüger's published
trench is translationally invariant. The authority and all thresholds were
frozen in `PREREGISTRATION.md` before the 0.03125 s rung completed.

## Four-rung time ladder

Both trajectories start from the same flat 10 nm geometry, Guo/Kwon
translating-mixed-layer state, published Krüger boundary with ion
normalization 1, exact projective occlusion, and common-refinement state
transfer.

| nominal step | depth at 0.5 s | mouth at 0.5 s |
|---:|---:|---:|
| 0.0625 s | 6.126733 nm | 85.591999 nm |
| 0.03125 s | 6.216545 nm | 85.595661 nm |
| 0.015625 s | 6.254472 nm | 85.580895 nm |
| 0.0078125 s | 6.267566 nm | 85.583387 nm |

Every adjacent pair is a fresh `t=0` trajectory. The original
0.0625/0.03125 s pair passes at the endpoint but fails the maximum matched
depth criterion at the first scored point. The 0.03125/0.015625 s pair
narrows that early discrepancy but still fails the same frozen criterion.
The 0.015625/0.0078125 s pair passes every gate:

| twofold pair | max matched depth difference | terminal depth difference | terminal mouth difference | result |
|---:|---:|---:|---:|---|
| 0.0625 / 0.03125 s | 19.8965% | 1.4447% | 0.003662 nm | FAIL |
| 0.03125 / 0.015625 s | 6.4026% | 0.6064% | 0.014766 nm | FAIL |
| 0.015625 / 0.0078125 s | 1.9944% | 0.2089% | 0.002492 nm | **PASS** |

At the first scored instant, 0.0625 s, the final pair gives 0.648163 and
0.635488 nm. Both values remain sub-cell, but the difference is inside the
frozen 5% trajectory gate without waiving or changing it. All solver-health
gates pass:

- exact material ledgers;
- maximum radiosity balance residual `2.20e-12`;
- maximum extrusion deviation `6.94e-18` mesh units;
- raw-face/resolved-grid speed ratio exactly 1;
- zero rejected trials, topology events, and asymmetric cells.

`time_gate.json` and `time_gate_refined.json` are deliberate failing
receipts; `time_gate_final.json` is the passing receipt. Native-resolution
profiles for the final pair were visually inspected: the exterior gas,
floor, and mouth remain connected with no isolated one-cell feature.

`PREREGISTRATION_ADDENDUM_TIME_REFINEMENT.md` freezes the ensuing 10/5 nm
spatial pair at 0.015625 s before the 5 nm trajectory is run. It changes no
spatial threshold and introduces no experimental observable.

This entire board is numerical. No experimental depth, yield, flux scale, or
target-dependent constant enters the gate, and a future pass cannot by itself
identify Krüger's unpublished reactor boundary.
