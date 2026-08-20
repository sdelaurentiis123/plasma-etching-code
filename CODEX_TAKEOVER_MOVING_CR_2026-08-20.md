# Codex takeover report — Oxford moving-Cr board

Status timestamp: 2026-08-20 13:44 EDT  
Repository: `plasma-etching-code`  
Branch: `codex/validation-first-multiphysics`  
Authoritative fix commit: `d8a8d02` (pushed)  

This report supersedes `HANDOFF_MOVING_CR_BOARD_2026-08-20.md` wherever the
two disagree.  The earlier handoff preserved useful cache and campaign
context, but its live-process diagnosis and recovery instructions were not
true on inspection.

## Executive finding

The missing `w320 / selectivity 14 / ion_low_tail_0p0` cell was not stuck in
an unidentified non-converging physics loop.  It hit a deterministic numerical
singularity in periodic strip symmetrization:

1. A sparse energetic event carries a dimensional rate equal to its flux
   density times the hit triangle's area.
2. The v3 symmetrizer divided that event rate equally among all triangles in
   the equivalent extruded strip.
3. It then divided each equal rate share by the recipient triangle area.
4. When marching cubes produced unequal-area strip triangles, a near-zero-area
   sliver received a near-infinite flux density.
5. On accepted outer step 39, the grid-resolved recession speed jumped from
   about `7.245e-4` to `3310.3929595` mesh units/s.  The CFL controller then
   requested `5,331,155` internal advection substeps instead of 2.  That is why
   one CPU core ran indefinitely while the GPU was idle.

The corrected v4 operator projects **uniform flux density**, which is the
actual translational-symmetry invariant.  Event rate is distributed in
proportion to recipient triangle area, so total rate is conserved without a
small-area singularity.  The implementation is also vectorized while
preserving deterministic event-major/member-major ordering.

The exact formerly failing trajectory completed under v4 in 381.65 wall
seconds, with two CFL substeps and speed approximately `7.25e-4` throughout.
This is a direct same-trajectory reproduction of both the defect and the
correction.

## What was wrong with the prior handoff

The prior handoff claimed:

- an instrumented `/root/run_w320_debug.py` existed;
- `/root/w320_debug.log` could be monitored;
- the live worker was a known single stuck cell inside an otherwise live
  campaign;
- the box cost `$0.176/hr`;
- 55 remote caches implied only the target was absent.

The observed state was:

- neither the claimed debug driver, log, nor pid file existed;
- the campaign parent had exited;
- eight `multiprocessing.spawn` workers plus one resource tracker were
  orphaned under PID 1 for more than 12 hours;
- one orphan used 100% CPU and about 3.87 GB RSS; the other workers waited on
  futexes or pipes; GPU utilization was 0%;
- the authoritative Vast price is `$0.2011111111/hr` (`$0.1733333` GPU plus
  `$0.0277778` disk), not `$0.176/hr`;
- the remote directory held 57 JSON files: 55 eligible v3 production caches
  plus obsolete v1 and v2 smoke caches.  The target v3 production cache was
  genuinely absent.

The nine orphan PIDs were individually identity-checked and terminated.  No
file, cache, repository state, or unrelated process was deleted.

## Evidence chain

### 1. Reproducible diagnostic, not an ad-hoc remote script

Commit `f1f6620` adds
`scripts/debug_zhu_npg80_moving_cr_trajectory.py`.  It selects the exact board
cell, prints its full content-addressed job specification, runs the same
`_execute` path as production, registers `SIGUSR1` all-thread stack dumps, and
automatically dumps stacks on a timer.  Its job selection is regression
tested.

Commits `ae41c8c` and `94f623d` add diagnostic-only tracing for:

- sparse and projected event counts in strip symmetrization;
- grid shape, step duration, CFL substep count, and extended-velocity extrema.

These wrappers run only in the diagnostic process and do not alter production
physics.

### 2. Stack capture

An early two-minute snapshot landed in
`boundary_transport_3d._symmetrize_face_resolved`; later repeated snapshots
showed normal progression through remap and advection.  At step 39, repeated
snapshots stayed in `threed.advect_3d`, proving the process was spending its
time in millions of ordinary advection iterations rather than one unbounded
Python `while` loop.

### 3. Cardinality and CFL capture

Normal steps showed roughly:

- 115k–180k sparse energetic events;
- 0.56M–2.88M projected strip events;
- 2 CFL substeps;
- maximum resolved speed around `0.0007245` mesh units/s.

The failing step printed exactly:

```text
ADVECTION_SHAPE grid=(41, 41, 106)
duration_s=4.8312883435582821 dx=0.01
substeps=5331155
speed_abs_max=3310.3929595316654
speed_min=0 speed_max=3310.3929595316654
```

The immediately preceding step used 2 substeps at
`0.00072451052354539122`.  A factor of roughly 4.57 million in one accepted
step is incompatible with the rate-normalized removal law and identifies
numerical amplification, not a surface-physics transition.

### 4. Physics correction

Commit `d8a8d02` changes `_symmetrize_face_resolved` from equal rate per
triangle to equal flux density over the strip.  For group areas `A_j` and an
incoming event rate `R`, each member now receives flux density

```text
Phi_j = R / sum_j(A_j)
```

and hence rate `R_j = Phi_j A_j`.  Therefore:

- every member has the same flux density;
- `sum_j R_j = R` exactly up to floating-point summation;
- a vanishing-area triangle receives vanishing rate rather than divergent
  density;
- the extrusion invariant agrees with the neutral-flux projection already in
  the same function.

A new unequal-area regression test fails under the v3 formula and passes under
v4.  The relevant focused suite is 99 tests green.

### 5. Same-path confirmation

The v4 replay passed the exact former step 39 and completed after 82 accepted
steps.  All printed CFL counts were 2 and all maximum speeds stayed near
`7.25e-4`.  There was no stack timeout or GPU/context error at the former
failure point.  The resulting cache is
`w320_s14.000_ion_low_tail_0p0_d037283bbcd7364e.json`.

The declared endpoint was local Cr-mask loss: the footprint minimum was zero
while 26.77 nm remained at the centre.  At that point the conditional TiO2
profile depth was 278.16 nm, particle balance error was exactly zero, and the
maximum conservative-remap residual was `6.59e-16`.

## Critical scope finding: this is not yet a 20-minute profile board

Every one of the 55 committed v3 trajectories terminated on
`cr_mask_below_vertical_resolution_in_footprint`; none reached the requested
1200-second process endpoint.  Their accepted process-equivalent durations
ranged from 391.33 to 713.85 seconds and their endpoint depths from 274.90 to
406.09 nm.  The corrected v4 target likewise stopped at local mask loss after
about 396.17 reference seconds.

That is not automatically a bad physical result.  A 45 nm Cr mask combined
with the preregistered low TiO2:Cr selectivity cannot survive a 20-minute
blanket exposure, and feature corners can fail earlier than the centre.  But
it changes the deliverable meaning:

- a completed 56-cell cache board certifies conditional evolution **up to
  first unresolved local Cr loss**;
- it does not yet predict Freddie's final 20-minute SEM;
- copying the terminal geometry forward to 1200 seconds would be a false
  prediction;
- the next feature-engine requirement is a certified continuation through
  local hard-mask perforation and eventual material-component extinction.

That continuation must remove sub-resolution Cr locally while preserving the
remaining resolved mask, rebuild/redistance a valid material level set, and
keep topology, mass/remap, and particle-balance certifiers active.  It must not
silently reinterpret mask loss as complete process success.

## Cache consequence

The old 55 production caches are valid receipts for the v3 operator but are
**not eligible for a v4 certification board**.  The area weighting can change
surface flux wherever strip triangles have unequal area, not only in the one
catastrophic cell.  `MODEL_REVISION` is therefore bumped from v3 to v4, giving
new content-addressed cache names.  Reusing v3 results would hide a changed
transport operator and is forbidden.

Do not delete the committed v3 caches; they are forensic receipts.  A clean v4
campaign must create a separate set of v4 cache files and assemble a v4
`audit.json`.

## Live compute state at report time

Vast instance:

- ID `48177892`
- label `petch-zhu-moving-cr-bigram`
- SSH `ssh -p 17892 root@ssh6.vast.ai`
- RTX 3090, 28 effective vCPUs, about 179 GB RAM limit
- actual billing `$0.2011111111/hr`
- tree `/root/petch-4b656fd`
- venv `/root/petch-venv`

The exact v4 diagnostic completed.  Its log is
`/root/w320_v4_debug.log`; its pid file is `/root/w320_debug.pid` but that PID
is no longer live.  A clean eight-worker v4 board was subsequently launched
under parent PID 5718 with log `/root/zhu_v4_board.log` and pid file
`/root/zhu_v4_board.pid`.  Always re-read the pid and verify the command line
before sending a signal.  Never use `pgrep -f` with the pattern in the probe
command because it can match the probe itself.

The remote tree began from commit `4b656fd` and has the committed v4 versions
of these files copied into place:

- `src/petch/boundary_transport_3d.py`
- `scripts/audit_zhu_npg80_moving_cr_profiles.py`
- `scripts/debug_zhu_npg80_moving_cr_trajectory.py`

Run with both the repository and `src` on `PYTHONPATH`:

```bash
PYTHONPATH=/root/petch-4b656fd:/root/petch-4b656fd/src
```

## Required landing sequence

1. Monitor the active clean v4 production board.  The exact target cache is
   already complete, so 55 cells remain.  Because all v3 caches are
   ineligible, the other cells must recompute.  Eight workers are active;
   monitor RSS because 12 workers did not fit a prior 70 GB machine.
2. Assemble v4 `audit.json`, pull all v4 caches and audit, then run local
   `--check` against the identical committed code.
3. Run the full local suite, render the conditional endpoint atlas, and commit
   the v4 receipts.
4. Destroy Vast instance `48177892` immediately after all artifacts are pulled
   and verified.  Confirm it disappears from `vastai show instances --raw`.
5. Treat post-mask continuation as the next implementation block before
   calling this a final 20-minute profile predictor.

## Scientific interpretation

This fix removes a numerical false singularity and makes the deterministic
feature solver faster and more physical.  It does **not** turn the Oxford
board into an atomic-accuracy prediction.

The Freddie/Oxford board remains deliberately conditional because its TiO2
and Cr surface response is rate-normalized from cross-machine witness data,
and the target tool's achieved self-bias/species-resolved surface law has not
been independently measured.  The target SEM was not used.  Correct claims
after v4 are:

- the supplied geometry and conditional transport axes can be evolved without
  the v3 numerical singularity;
- the board is a blind, preregistered conditional profile envelope;
- it can be scored against Freddie's SEM when received;
- it is not yet an absolute Oxford knobs-to-SEM certification.

Keep those boundaries explicit in any website or external writeup.

## Repository hygiene

The following pre-existing untracked user artifacts were observed and left
untouched:

- `results/curated/mixed_layer_feature_v1/ml20-bonds-12s.log`
- `results/curated/mouth_equilibrium_probe_dx/`
- `scratch_ignore_calc.py`

No reset, checkout, cache deletion, or broad process kill was used.
