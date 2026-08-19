# Zhu NPG80 absorbed-power ensemble v1

The Oxford screenshot supplies 150 W forward power, not absorbed plasma
power. This target-free ensemble holds the 276 V machine-family self-bias
sensitivity fixed and solves 60, 90, 105, and 120 W absorbed-power states.
Every accepted state closes the conserved open-reactor equations, the
powered/grounded sheath fixed point, and the deterministic axisymmetric lift.

Central-optic ion flux rises from about 1.13e19 to 1.74e19 m^-2 s^-1 over the
40--80% transfer envelope. Consequently, the exact 700 nm atom-count
requirement spans roughly 0.82--1.63 TiO2 formula units per incident positive
ion across density and power. The frozen earlier binary clearance call remains
part of the record, but it is not a robust current depth verdict after the
wall-power correction and this transfer envelope. Unique clearance/profile
depth remains unidentified without target TiO2/Cr surface response or a
measured absorbed-power/wafer-flux boundary.

The 120 W solve also found a model-form issue. It was rejected at the original
600 Td bound with a 0.45% CHF3 balance residual. Extending the represented-feed
coordinate to 900 Td closes conservation at about 690 Td, while the E/N formed
with the total neutral inventory remains below 100 Td. The disparity occurs
because only about 13% of the neutral inventory is represented in the electron
collision basis at that dissociated state. The 120 W row is therefore a field-
domain sensitivity and a direct pointer to the next reactor task: electron
collisions with daughter gases. It is not promoted as a stronger machine
prediction.

Rebuild the audit after reproducing the state continuation ladder:

```bash
python scripts/audit_zhu_npg80_absorbed_power_ensemble.py --check
pytest -q tests/test_zhu_npg80_absorbed_power_ensemble.py
```
