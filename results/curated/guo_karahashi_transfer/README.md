# Guo/Kwon transfer to mass-selected reactive-ion beams

This board is a no-fit test of the frozen Guo/Kwon C4F8-Ar/SiO2 translating
mixed-layer closure against all 21 digitized Karahashi F+, CF+, CF2+, and CF3+
beam-yield points. The experiment is deliberately simple: one mass-selected
ion, known energy, normal incidence, and no neutral radicals.

Run:

```bash
python scripts/audit_guo_karahashi_transfer.py --check
```

The atom ledgers and steady-state solves close, but the physical transfer does
not. The model overpredicts F+ and CF+ removal and predicts deposition at the
measured 250 and 500 eV CF3+ etch points. CF2+ and CF3+ are closer at some
higher energies, but those points were not a preregistered validation subset.

Consequently the CF2+/CF3+ Krueger depth envelope remains a useful ion-identity
sensitivity, not a validated species-resolved prediction. The next closure
must distinguish projectile stopping, implantation/fragmentation, and
chemical removal by ion species and energy. The measured Karahashi table
remains the highest-evidence direct closure inside its exact support.
