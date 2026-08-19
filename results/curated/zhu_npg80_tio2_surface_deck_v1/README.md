# Oxford TiO2 fail-closed surface deck

This receipt binds the executable TiO2-specific deck to the Oxford condition
without inventing a surface coefficient. The deck labels ALD TiO2 and TiO2
formula-unit inventory explicitly, converts an explicitly supplied film mass
density into formula-unit density, and requires provenance for every sticking,
passivation, oxygen-removal, and energetic-yield input.

The reduced common kernel may be executed only through an explicit sensitivity
override. It now carries a separate, bounded competitive oxygen site fraction:
oxygen adsorption and energetic cleanup are integrated analytically and
conservatively with fluorination and passivation. It still cannot certify an
absolute Oxford profile because the target coefficients are not identified and
chemistry-dependent roughness evolution remains unresolved. Even a full set of
evidence flags cannot override that model-form gap.

This is useful progress rather than a placeholder: future beam, blanket, or
same-tool measurements now land in named dimensional slots, and no SiO2 number
can become a TiO2 prediction through relabelling.

```bash
python scripts/audit_zhu_npg80_tio2_surface_deck.py --check
pytest -q tests/test_tio2_surface_deck.py \
  tests/test_audit_zhu_npg80_tio2_surface_deck.py
```
