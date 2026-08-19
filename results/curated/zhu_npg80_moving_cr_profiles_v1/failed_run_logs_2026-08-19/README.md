# Production board failure logs — 2026-08-19

Both launches of the 56-trajectory moving-Cr board died before any full
trajectory completed (the two cached w200 files were 12 s smoke probes).

- `zhu_production.log` — first launch: the sub-cell Cr island failure,
  fixed by the fixed-point cleanup in 4265e70.
- `zhu_production_fixed.log` — relaunch at the fixed tree: died at
  w80/low/tail0/sel14 with "marching-cubes surface contains 3 unmatched
  interior edges". Root cause: the Cr mask thins non-uniformly, so
  footprint-edge columns degenerate to sub-cell slivers while the v2
  centre-column exhaustion check still sees >1 vertical cell. Punch-through
  arithmetic (45 nm Cr clears at 630/811 normalized-nm against the 869 the
  fast rate endpoint requires) makes mask exhaustion the generic outcome of
  every cell, not a w80 quirk.

Fix: model revision v3 (commit 2ad547d) ends each trajectory as a declared
`cr_mask_below_vertical_resolution_in_footprint` event; the topology
certifier is untouched.
