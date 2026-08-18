# Zhu NPG80 self-bias/global-state ensemble v1

This board propagates every printed Oxford-80 self-bias witness through the
wall-resolved conserved reactor and the deterministic axisymmetric wafer
transfer. Each 200, 250, 276, 300, 360, 387, and 400 V state closes its
Maxwellian plasma-potential/sheath fixed point below 0.01 V, retains its full
numerical continuation chain, and uses no SEM, etch depth, or surface yield.

At fixed 90 W absorbed power, ion and F flux both decrease monotonically with
self-bias. The 200 V state has 1.141x the ion flux and 2.608x the F flux of the
400 V state: increasing sheath voltage consumes more of the fixed power at
the wall and reduces the bulk plasma sustained by the remainder. This is why
voltage cannot be varied independently of reactor density in the depth model.

The exact NGP80 source only states censored conditioning thresholds (>300 V
initially and <about 200 V finally). The receipt integrates the explicitly
linear 300-to-200 V witness with deterministic Simpson nodes at 300, 250, and
200 V. It remains a sensitivity history, not a reconstructed measurement.

Every state also includes a 48x16 axisymmetric central-optic flux. The smooth
radial correction stays small; no target serial-tool radial diagnostic or
TiO2/Cr surface law is inferred. Absolute profile-depth support remains false.

```bash
python scripts/audit_zhu_npg80_self_bias_global_ensemble.py --check
```
