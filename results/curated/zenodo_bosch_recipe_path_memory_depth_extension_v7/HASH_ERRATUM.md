# Hash erratum for the frozen v7 preregistration

The `target_firewall.calibration_measurements_sha256` value written in
`preregistration.json` is a transcription error. It says

```text
c46b82fd7d15903dc61098d72823ba3992cf8ac7b14c69d4c4b419361ccbc317b05
```

The actual SHA-256 of
`data/experimental/zenodo_17122442/calibration_Si_Oxide_etch_89_points.csv`
at the v7 freeze commit `f784852` and at the current audit is

```text
c46b82fd7d15903dc61098d72823ba3992cf8ac7b14c69d71e3e90da4caa20eb
```

The correct value is independently present in the v2, v3, v5, and v6
preregistrations and in the v5 and v6 calibration receipts. The calibration
asset itself did not change. The v7 input loader used the repository asset and
the frozen calibration-key allowlist; it did not use the mistyped hash as an
asset selector.

The frozen v7 preregistration is deliberately not rewritten after execution.
This erratum records the literal metadata defect, the verification command,
and the correct value for all later receipts:

```text
git show f784852:data/experimental/zenodo_17122442/calibration_Si_Oxide_etch_89_points.csv | sha256sum
```

This correction changes no outcome, coefficient, split, gate, or heldout
firewall.
