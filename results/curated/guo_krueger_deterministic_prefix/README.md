# Guo/Kwon–Krüger deterministic-extruded prefix refinement

Status: **10 nm time gate failed exactly one preregistered trajectory
criterion; finer rung required.**

The candidate moving-profile neutral-exchange authority is
`deterministic_extruded_2d`, appropriate only because Krüger's published
trench is translationally invariant. The authority and all thresholds were
frozen in `PREREGISTRATION.md` before the 0.03125 s rung completed.

## First two rungs

Both trajectories start from the same flat 10 nm geometry, Guo/Kwon
translating-mixed-layer state, published Krüger boundary with ion
normalization 1, exact projective occlusion, and common-refinement state
transfer.

| nominal step | depth at 0.5 s | mouth at 0.5 s |
|---:|---:|---:|
| 0.0625 s | 6.126733 nm | 85.591999 nm |
| 0.03125 s | 6.216545 nm | 85.595661 nm |

The terminal depth difference is 1.445%, terminal mouth difference is
0.0037 nm, and all solver-health gates pass:

- exact material ledgers;
- maximum radiosity balance residual `2.10e-12`;
- maximum extrusion deviation `6.94e-18` mesh units;
- raw-face/resolved-grid speed ratio exactly 1;
- zero rejected trials, topology events, and asymmetric cells.

The frozen maximum-trajectory depth gate nevertheless fails. At 0.0625 s,
the depths are only 0.827 and 0.690 nm—both far below one 10 nm cell—and
their relative difference is 19.90%. The discrepancy falls to 0.72% by
0.125 s and stays near 2% or below afterward. That pattern suggests an
initial interface/metric transient, but the source of a failure cannot waive
a preregistered gate.

`time_gate.json` is therefore a deliberate failing receipt, not a green
result. The next authorized action is a fresh 0.015625 s rung at 10 nm.
Spatial refinement and the 60 s sensitivity forecast remain prohibited until
a twofold time pair passes or the startup interface is repaired and
re-preregistered.

This entire board is numerical. No experimental depth, yield, flux scale, or
target-dependent constant enters the gate, and a future pass cannot by itself
identify Krüger's unpublished reactor boundary.
