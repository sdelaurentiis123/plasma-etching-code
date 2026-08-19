# Oxford TiO2 surface-topology contract

This target-free audit binds three independent TiO2 response axes before the
withheld Oxford SEM is used: Choi's DC-bias response, Ji's RF morphology
response, and Ji's spacing response. Together they reject a single
energy-independent removal scalar as a complete target surface law.

The existing conditional Oxford profiles remain physically meaningful inside
their declared removal-only scope. In particular, a wider lower section does
not require material growth: a shadowed lower wall can simply recede more
slowly than the upper wall. But numerical convergence of that geometry cannot
identify Freddie's actual mechanism because the present sentinel omits the
experimentally required TiO2 fluorination, passivation volume, Cr-mask motion,
and energy-dependent desorption.

The minimum next solver remains deterministic and differentiable. Each
surface element carries bounded fluorination and activation fractions, a
nonnegative passivation inventory, and a removed-material ledger. Deterministic
angular quadrature supplies species/energy/angle-resolved flux; analytic
bounded reaction updates and conservative level-set remapping evolve TiO2,
Cr, and the overlayer. No Monte Carlo operator is required.

These cross-process boards identify model topology and response signs, not the
Oxford coefficients. Absolute profile prediction still needs same-condition
self-bias/waveform, blanket TiO2 loss, Cr loss, target GDS/radius, and the SEM
answer key.

```bash
python scripts/audit_zhu_npg80_tio2_surface_topology.py --check
pytest -q tests/test_audit_zhu_npg80_tio2_surface_topology.py
```
