# Karahashi reactive-ion event audit

`audit.json` separates source-backed CF3+ impact facts from a hypothetical
Figure-4/Figure-10 stoichiometric join:

```bash
python3 scripts/audit_karahashi_reactive_ion_event.py --check
```

The audit is intentionally non-production. Figure 4 gives normal-incidence
total SiO2 removal yield, whereas Figure 10 does not report ion incidence
angle and gives only fractions among detected SiFx products. The source also
shows that product composition changes with angle. The conditional join is
therefore useful only to expose atom-balance requirements, including exchange
with a resident fluorocarbon surface layer.

The production conclusion is negative but constraining: mass-partitioned
fragment energy can enter a transport closure, but the measured beam yield
must not be combined with an additional live-film attenuation law until a
finite C/F/O/Si surface inventory and delayed-product kinetics are independently
identified. Doing so would risk counting the beam-conditioned surface twice.
