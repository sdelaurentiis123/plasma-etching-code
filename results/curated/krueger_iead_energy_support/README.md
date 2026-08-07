# Krüger IEAD surface-physics energy support

This provenance board asks how much of Krüger's published Figure 4(a)
positive-ion distribution lies inside each independently established
surface-physics energy domain. It performs no yield fit and never reads the
825 nm feature endpoint.

Run:

```bash
python scripts/audit_krueger_iead_energy_support.py --check
```

The result is restrictive:

| surface-physics domain | upper energy | IEAD at or below | IEAD above |
|---|---:|---:|---:|
| Guo/Yin beam-regressed TML coefficients | 370 eV | 0.000% | 100.000% |
| An et al. released DFT-trained NNP/ZBL outputs | 1000 eV | 5.102% | 94.898% |
| Karahashi direct mass-selected CFx+/SiO2 beams | 2000 eV | 13.183% | 86.817% |
| Tachi Si-target literature lead, not SiO2 support | 3000 eV | 27.142% | 72.858% |

The input distribution itself is not a measurement: it is the source's HPEM
output, digitized from a normalized, clipped color plot. It combines all
positive ions and does not publish their composition. Energy coverage
therefore cannot solve the projectile-identity problem, identify the C4F6
radical boundary, or transfer a C4F8 beam regression to C4F6 plasma.

The Tachi row is deliberately quarantined as a target-mismatched lead. The
accessible record describes mass-selected fluorocarbon bombardment of
elemental Si through 3 keV; it is neither a landed numerical SiO2 board nor
authority for a SiO2 yield law.

Consequently, the current finite-fluence feature trajectory is a valuable
numerical and mechanism sensitivity, but not an atomic-level depth
prediction. Closing that claim requires a validated high-energy
species-resolved reactive-event kernel, persistent mixed-layer/product
physics, and a species-resolved reactor boundary. Matching 825 nm may not be
used to choose any of those missing quantities.
