# Zhu NPG80 open-reactor checkpoint v1

`central.json` is the first conservation-gated 66-species state for the
withheld-SEM Oxford NPG80 recipe.  It uses no feature depth or TiO2 response.
Its declared machine sensitivities are 90 W absorbed power, 350 K gas,
170 mm powered diameter, 30 mm active height, 100 um ion momentum mean free
path, and 250 eV mean all-wall positive-ion energy.

The nonlinear solve is initialized from
`hydrogen_closed_continuation.json`, the preceding 59-species tier before the
Lim CHF3/O2 daughter block was enabled.  That file is a numerical continuation
state, not an experimental calibration.  Its SHA-256 is pinned in
`central.json`.

The central state passes particle, charge, pressure, and power balances at a
maximum normalized residual of `3.90e-7`.  It is not a certified wafer-flux or
depth prediction: only `37.6%` of its neutral inventory is represented in the
current CHF3/SF6/O2 electron-collision basis, and all declared machine closure
inputs remain sensitivities rather than target-tool diagnostics.

Rebuild with the hash-locked Song O2 workbook:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python scripts/run_zhu_open_reactor.py \
  --source-workbook /path/to/o2_song_2026_supplement.xlsx \
  --initial-state-json \
    results/curated/zhu_npg80_open_reactor_v1/hydrogen_closed_continuation.json \
  --electrode-diameter-mm 170 \
  --output results/curated/zhu_npg80_open_reactor_v1/central.json \
  --maximum-evaluations 150
```

This explicit `170 mm` flag is now required to replay the immutable v1
development state. It was an unsourced sensitivity. Oxford's current
PlasmaPro 80 specification gives a `240 mm` electrode; the exact CUNY
inventory confirms that the target tool is an NPG80. The source-correct
follow-up is archived separately in `zhu_npg80_open_reactor_v2` so the blind
v1 forecast is never rewritten after the fact.
