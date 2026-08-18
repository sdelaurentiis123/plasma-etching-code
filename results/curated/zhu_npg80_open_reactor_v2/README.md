# Zhu NPG80 source-geometry correction v2

This rung replaces the unsourced `170 mm` v1 development diameter with the
manufacturer's `240 mm` PlasmaPro 80 electrode specification. The CUNY ASRC
inventory independently identifies the exact facility tool as an Oxford
PlasmaPro NPG80 RIE with a 300 W, 13.56 MHz RF source. The evidence and its
remaining gaps are pinned in
`data/experimental/zhu_2026_tio2_npg80/machine_geometry_evidence.json`.

Both states remain target-free conserved global sensitivities. They use the
same recipe, 90 W absorbed-power sensitivity, 30 mm active-height sensitivity,
wall closures, and v1 continuation state. `source_geometry_central.json` uses
the Voloshin 350 K CHF3+F rate; the alternate uses Lim's published 700 K rate,
which is 11.5 times lower at its stated temperature.

The geometry correction is not small. The central positive-ion flux is
`1.251e19 m^-2 s^-1`, 54.64% of v1, while the F flux is 38.26% of v1. At the
source-correct geometry, changing between the two published CHF3+F branches
moves total ion flux only 1.19% but F flux 37.33%. Thus the binary ion-dose
ledger is robust to this neutral-rate conflict and sensitive to the machine
volume; polymerization, Cr loss, and detailed profile remain sensitive to the
radical boundary.

For a 700 nm film and 1200 s exposure, the corrected central state requires
`1.143--1.459` TiO2 formula units per incident positive ion before feature
attenuation. That requirement is physically reachable, but a target-tool
surface yield, local wafer flux, feature transmission, and exact self-bias are
not measured. The reactor-only clearance call therefore fails closed as
unresolved. The preregistered v1 binary forecast remains immutable and is not
retrofitted.

Rebuild the two states, then the receipt:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python scripts/run_zhu_open_reactor.py \
  --source-workbook /path/to/o2_song_2026_supplement.xlsx \
  --initial-state-json results/curated/zhu_npg80_open_reactor_v1/central.json \
  --output results/curated/zhu_npg80_open_reactor_v2/source_geometry_central.json \
  --maximum-evaluations 300

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python scripts/run_zhu_open_reactor.py \
  --source-workbook /path/to/o2_song_2026_supplement.xlsx \
  --initial-state-json results/curated/zhu_npg80_open_reactor_v1/central.json \
  --output results/curated/zhu_npg80_open_reactor_v2/source_geometry_chf3_rate_alternate.json \
  --chf3-f-rate-branch lim_700K \
  --maximum-evaluations 300

python scripts/audit_zhu_npg80_geometry_correction.py --check
```
