# Oxford NPG80 species-resolved feature boundary v1

This receipt closes the deterministic software join from the strongest
conserved Oxford reactor state to the microscopic feature boundary. All 20
positive ions and 37 thermal neutrals retain their names, absolute fluxes,
masses, and charges. Singly and doubly charged ions receive distinct
charge-resolved impact energies, and the two-component angular law is exported
as fixed quadrature rather than sampled particles.

This interface is not TiO2-specific. The adapter accepts arbitrary reactor
species, species-specific energy measures and IADFs, and arbitrary neutral
inventories. Oxford/TiO2 is the first audited configuration passed through it.

The result remains conditional. The 276 V self-bias is a same-family transfer,
not a diagnostic from Freddie's exact run. The Kim angular widths are direct
beam measurements but not from this Oxford chemistry, the 0.65 tail share is a
sensitivity, and neutral radial transfer is not yet solved. Most importantly,
the target TiO2/Cr reaction probabilities are not measured. Consequently the
boundary is executable for feature transport, but it does not certify a unique
absolute SEM profile.

Rebuild or check the checksum-bound receipt with:

```bash
python scripts/audit_zhu_npg80_species_resolved_feature_boundary.py
python scripts/audit_zhu_npg80_species_resolved_feature_boundary.py --check
```
