# Deterministic electron-kinetics receipts

`two_term_bolos_oracle_v1.json` is the stdout receipt from:

```text
python scripts/audit_two_term_bolos_oracle.py \
  --bolos-path /private/tmp/petch-bolos-ref
```

The comparison uses a repository-defined manufactured elastic-plus-excitation
deck, not LXCat data. BOLOS 0.2 is an LGPL local oracle and is not a petch
dependency or redistributed artifact. At 2400 fixed energy cells, petch and
BOLOS agree to 0.458% in mean energy and 0.634% in excitation rate; the
flux-reduced-mobility residual is 0.634%, the scalar-reduced-diffusion
residual is 0.403%, and the weighted EEPF L1 difference is 0.122%. The
monotone refinement trend is part of the gate. The remaining discrepancy is
consistent with petch's current
piecewise-constant inelastic source reconstruction versus BOLOS's
piecewise-exponential reconstruction.

This is a numerical operator gate only. It does not validate a physical
collision deck, measurement-specific swarm coefficients, a reactor state,
wafer flux, or feature depth.
