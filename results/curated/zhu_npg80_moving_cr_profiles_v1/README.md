# Oxford TiO2 profiles with moving Cr mask

This board upgrades the pinned-mask conditional atlas to the common
multi-material feature engine. ALD TiO2 and the 45 nm Cr mask evolve on
separate material level sets under the same deterministic energy-angle
transport. The preregistered Janissen TiO2-rate interval and TiO2:Cr
selectivity interval remain cross-machine sensitivity inputs; Freddie's SEM
and depth are not used.

Production uses a 10 nm mesh so the initial mask spans 4.5 vertical cells. It
reports mask thickness and stops a controlled-profile claim if the Cr center
is exhausted. The board is not a microscopic Cr prediction: Nguyen et al.
show that real F/O Cr etching cycles through CrOx and CrFx with neutral
conversion and ion-assisted inhibitor removal. Those state coefficients are
not published for the Oxford condition.

```bash
python scripts/audit_zhu_npg80_moving_cr_profiles.py --check
pytest -q tests/test_audit_zhu_npg80_moving_cr_profiles.py
```
