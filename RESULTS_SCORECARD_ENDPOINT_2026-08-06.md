# The 60 s endpoints, measured (2026-08-06)

The scorecard's published criteria are endpoint quantities. Every previous
grading of the final mechanism was an early-time or extrapolated proxy:

- `ml23` (the Gray-laws confirmation) was cut at **t = 2.80 s** and the 60 s
  depth was a linear extrapolation, quoted as ~418 nm (−49 %);
- the 2026-08-06 scorecard reached only **t = 1.948 s** because the box drew a
  warp build with no CUDA device, so its grader explicitly refused clog,
  necking and absolute depth.

This pass ran all six conditions to **t = 60.000 s** on a GPU-verified box
(`nvidia-smi` + `torch.cuda.is_available()` + `warp.get_devices()` all checked
before deploy, `warp-lang` pinned to 1.15.0). Nothing below is extrapolated.

## Endpoint state

| condition | depth (nm) | opening (nm) | top CD (nm) | mask (nm) | closure/etch | x Krüger |
|---|---|---|---|---|---|---|
| base (6 kW, O2 1.0) | 346.8 | 39.82 | 64.26 | 850.21 | 0.0723 | 2.33 |
| O2 0.5 | 201.6 | 13.08 | 49.04 | 850.40 | 0.1908 | 6.15 |
| O2 1.5 | 322.4 | 33.60 | 44.17 | 850.27 | 0.0875 | 2.82 |
| O2 2.5 | 366.6 | 49.60 | 53.10 | 850.24 | 0.0551 | 1.78 |
| 4 kW | 315.1 | 34.86 | 52.68 | 850.25 | 0.0875 | 2.82 |
| 8 kW | 359.2 | 37.12 | 64.48 | 850.20 | 0.0736 | 2.37 |

## Graded criteria

| criterion | target | measured | verdict |
|---|---|---|---|
| **r(4/6)** | [0.84, 0.94] | **0.909** | **PASS** |
| **r(8/6)** | [0.97, 1.06] | **1.036** | **PASS** |
| base mask remaining | 850 +/- 2 % | 850.21 | **PASS** |
| necking absent at O2 2.5 | open | 49.60 | **PASS** |
| base etch depth | 825 +/- 5 % | 346.8 | **MISS (−58.0 %)** |
| base mask opening | 45 +/- 5 | 39.82 | **MISS** (0.18 nm below band) |
| clog at O2 0.5 | sealed (0.000) | 13.08 | **MISS** |
| O2 depth rank max at 1.5 | 1.5 highest | 2.5 highest | **MISS** |

Grader: `scripts/grade_scorecard_endpoint.py`, which marks any run that did not
reach 60 s INCOMPLETE rather than grading it early.

## What the endpoints settle

**The power transfer is certified, not indicated.** Scorecard-1 (ml9a
mechanism) read r(4/6) = 0.672 (MISS) and r(8/6) = 1.085 (near-miss). The
matched-time pass read 0.902 / 1.017 and was careful to call that an early
indication because the published bands are endpoint quantities. Measured at the
endpoint they are **0.909 and 1.036** — both in band, and within 0.02 of the
matched-time values, so the early indication was sound. This is the scorecard's
strongest result and it is now an endpoint result.

**Mask survival is exact across all six conditions** (850.20–850.40 nm against
850). The armored a-C mask reproduces at every oxygen ratio and both powers,
which is a consistency statement the base case alone could not make.

**The depth miss is worse than the extrapolation, and that is the honest
number.** `ml23`'s linear projection gave ~418 nm; the measured endpoint is
**346.8 nm**, −58.0 % against 825 +/- 5 %. The projection was optimistic
because it extrapolated a rate that was still decelerating: the run's own
aspect-ratio dependence — the physical behaviour the Gray laws restored — bends
the trajectory below any straight line drawn through its early span. Every
condition is consistent with this: all six sit at 202–367 nm where the
published etch is 825.

**The mouth no longer equilibrates at 50.9 nm.** Pre-Gray it did (drift
−17 pm/s, `ml19`). Under the final mechanism the base aperture passes *through*
45 nm at t ≈ 44 s and ends at 39.82 — 0.18 nm below the band. This is not an
independent defect: with the etch slowed, more lip closure accumulates per nm of
depth, and the run-average closure/etch is 0.0723 against Krüger's 0.0310
(2.33x). Depth and mouth are one statement — the etch is too slow relative to
lip closure — which is what the decomposed channel bound already says.

**Clog at O2 0.5 is directionally preserved but does not seal within 60 s.**
The starved case is by far the narrowest (13.08 nm, 3x narrower than base) and
was still closing at the endpoint (16.3 nm at t = 50.7 s), so the mechanism
produces the clogging *branch*; it does not reach zero in the published time
because the etch is slow. Recorded as a MISS on the published criterion rather
than reinterpreted.

**The oxygen ordering inverts at 1.0/1.5 at the endpoint too.** Aperture
ordering is 0.5 (13.1) < 1.5 (33.6) < 1.0 (39.8) < 2.5 (49.6): correct at the
extremes, inverted in the middle pair — the same inversion the matched-time
pass recorded at 0.60 nm, now 6.2 nm at the endpoint, so it is a real feature of
the mechanism and not sub-cell noise. The depth ranking has its maximum at
O2 2.5 rather than 1.5, which is a MISS on the published criterion.

## Cost

One RTX 3090, GPU verified before deploy. Base run solo at ~37 s/step; the five
conditions resumed from checkpoints and ran concurrently at 20–180 s/step.
~7 box-hours, **~$1.20**. Box 46992097 destroyed; `vastai show instances`
returns NONE.

## Provenance note

Six conditions, one commit, one mechanism. The runs carry the full final stack:
Gray's measured `sqrt(E)` yield laws with the co-regressed `s0`, deposition-driven
crosslinking with published per-material bond multiplicity, the measurement-bounded
oxide angular class, the two-component beam, and the sliver dissolution that
makes 60 s completions reachable at all — every 60 s run in this pass completed
without a topology stall.
