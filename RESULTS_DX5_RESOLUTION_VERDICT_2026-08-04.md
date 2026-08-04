# dx = 5 nm resolution verdict — the front-loaded taper is not a grid artifact

Run `ml17a-dx5-short` (box 46830758, 167 GB / 32-core, RTX 3090, **destroyed**),
HEAD `bd35d86` deployed clean via `git archive` (no live-tree patching — the
stale-module lesson of 2026-08-04, which cost three false relaunches).
Archived: `results/curated/ml17_dx5_resolution/`.

## Purpose

Campaign 5 + the neck regrade (`fcc98eb`) localised the whole remaining
validation defect to the top ~150 nm of mask, where petch closes **10.8×**
faster than Krüger. `RESULTS_LIP_CHANNELS_AND_CRITICAL_ANGLE_2026-08-04.md`
reduced that to a taper-propagation question: our simulated wall is tilted
17.3° in the top 50 nm and exhausted by 100 nm, where his sustains ~7°. The
standing suspect was discretisation — Krüger runs **1 nm voxels**, petch runs
10 nm, and one cell of lip film per side is 20 nm of aperture.

## Deployment note (the guard was never broken)

The three prior 5 nm launches died on `geometry variation along the declared
extrusion axis exceeds the projection guard; deviation=6.65e-05, guard=5e-05`
*after* the guard was raised to 0.02 cells (`3848197`, `bd35d86`). Reproducing
locally showed step 1 clean. The failures were a **stale module on the box** —
files patched onto a live tree while an old interpreter held the pre-fix code.
Clean `git archive` deployment passes the guard on the first attempt
(`grep -c "projection guard"` = 0). Guard values verified 0.02 in the extracted
tree *before* launch.

## Cost reality that forced the method

At dx = 5 nm the first step took **1583 s** wall (one-time exchange build) and
steady steps ~500 s at dt ≈ 0.071 s. Reaching t = 12 s needs ~120 steps ≈
**17 h**; t = 5 s ≈ 7 h. Rather than wait, the verdict was taken from
**matched-simulated-time comparison** against the archived 10 nm run
(`ml16a-verbatim-lift`), which is the same measurement the 12 s endpoint would
have supplied, evaluated over the interval both runs share.

## Result — matched simulated time, 5 nm vs 10 nm

| t (s) | depth 5 nm | depth 10 nm | Δ | cumulative closure 5 nm | 10 nm | Δ |
|---|---|---|---|---|---|---|
| 0.200 | 4.404 | 4.406 | −0.05 % | 1.910 | 1.969 | −3.0 % |
| 0.250 | 5.505 | 5.508 | −0.05 % | 2.287 | 2.461 | −7.1 % |
| 0.300 | 6.610 | 6.609 | +0.02 % | 2.643 | 2.953 | −10.5 % |
| 0.350 | 7.712 | 7.707 | +0.07 % | 3.014 | 3.397 | −11.3 % |

- **Floor etch rate is grid-converged to 0.05 %** (22.0 nm/s both). The
  substrate channel carries no discretisation error worth naming.
- **Mouth closure is ~8 % slower at 5 nm**, stable and converging (3.0 → 11.3 %
  as the film establishes). Richardson extrapolation to dx → 0 gives
  **−15 % (2nd order) to −22.5 % (1st order)** relative to the 10 nm closure.
- Throat position is identical at both resolutions: `throat_z` pins at
  2.645 µm — 5 nm below the 2.650 µm mask top — from the first step onward.
  The pinch starts at the mask top on **both** grids.

## Verdict

**Resolution is not the driver.** Halving dx buys 8–11 % less mouth closure and
converges; the discrepancy to be explained is **10.8× (+980 %)**. Even the
infinitely-refined limit recovers at most ~22 % of it — under 3 % of the gap.
The front-loaded taper is a property of the modelled physics/transport, not of
the mesh, and the 1 nm-voxel difference against Krüger is **acquitted as the
explanation** for the top-band excess.

Correspondingly, the ml16 depth undershoot (590 vs 825 nm) is also not a
resolution effect: the floor rate agrees to 0.05 % at both grids, so the depth
deficit is inherited from the throttled aperture, exactly as the neck regrade
implied.

## What this leaves

The top-band defect is now bounded on every side that has been tested:

| candidate | status |
|---|---|
| surface chemistry constants | matched 200–270 nm band to 1–11 % (`fcc98eb`) |
| O-radical channel | normalisation bug fixed, exact to 4 digits (`63cfefa`) |
| grazing ion removal law | faithful transcription; would need f = 61.6 vs 4.17 max (`aca2aeb`) |
| lip deposition delivery | faithful: visibility, sticking, isotropy all verified (`cbbd2d6`) |
| wall-slope / critical angle | falsified as sole mechanism (`3b7f104`); his neck holds at 1.86° |
| hot-neutral delivery | 0.1 % of top-band ion flux, already counted (`350c971`) |
| **grid resolution** | **falsified here — ≤22 % of a 980 % gap** |

What remains untested is the **time-dependence of taper propagation** — how the
tilt profile evolves over 0–30 s, where the 10 nm reference shows the throat
descending normally to 180 nm by t = 30 s and then *reversing* to 130 nm while
the aperture collapses 23.5 → 11.1 nm. That reversal (a second, higher
constriction overtaking the descending neck) is the specific unexplained event,
and it is a profile-evolution question, not a per-face budget question.

### 10 nm reference trajectory (the reversal)

| t (s) | opening (nm) | throat depth below mask top |
|---|---|---|
| 2 | 73.1 | 30 nm |
| 6 | 50.7 | 70 nm |
| 12 | **38.4** | 120 nm |
| 20 | 29.0 | 160 nm |
| 30 | 23.5 | **180 nm** |
| 40 | 16.8 | 130 nm ← reversal |
| 60 | 11.1 | 130 nm |

At t = 12 s petch reads 38.4 nm against the experimental 39.0 nm neck: the early
trajectory is *correct*. The defect is entirely a late-time takeover.

## Method note for reuse

Matched-simulated-time comparison against an archived trajectory answered in
~2.5 h of box time (≈ $0.50) what the specified 12 s endpoint would have cost
17 h. Any future grid or parameter sensitivity on this problem should be graded
this way rather than at a fixed endpoint.
