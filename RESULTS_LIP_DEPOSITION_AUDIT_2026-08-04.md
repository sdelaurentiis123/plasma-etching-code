# The lip deposition side: no over-delivery, but the probe's transport was wrong (2026-08-04)

`RESULTS_LIP_REMOVAL_AUDIT_2026-08-04.md` (`aca2aeb`) falsified the angular
ion-removal law and concluded the residual ~0.8 x deposition must live on the
*deposition* side.  This pass measured that side directly -- delivered depositor
flux against an analytic arrival reference, and effective sticking against the
published rows -- per depth band, with the 200-270 nm band (which already matches
Krüger to 0.88-1.11x) as the control.

**The deposition side has no over-delivery bug.**  What the audit found instead
is that the probe itself was running the wrong transport configuration, and the
error is *angle selective*, so it silently distorted every removal/deposition
ratio the probe has ever reported -- including the ones the previous two passes
reasoned from.

## The defect: periodic profile, non-periodic transport

`scripts/mouth_equilibrium_probe.py` called `advance_feature_step_3d` with
`profile_periodic_lateral=True` but never set `ballistic_periodic_lateral=True`.
In `feature_step_3d.py` the ballistic flag is *not* inherited from the profile
flag:

```
periodic_ballistic = (periodic_neutral if ballistic_periodic_lateral is None
                      else bool(ballistic_periodic_lateral))
```

and `periodic_neutral` is false unless a radiosity backend is selected.  The
gather therefore ran its finite-source branch, which keeps a ray only if its
back-projection lands inside the source rectangle:

```
source_point = points + travel * reverse            # boundary_transport_3d
visible &= (source_point in source_bounds)
```

The probe's cell is 0.13 x 0.02 um with the source plane ~0.15 um above the mask
top.  A ray leaving at 45 deg is displaced ~0.15 um laterally -- far outside the
20 nm y-extent -- so essentially every off-vertical direction was discarded.

That cut is **angle selective**, which is what makes it dangerous rather than
merely wrong in scale:

| population | angular spread | fate under the finite-source cut |
|---|---|---|
| thermal neutrals (depositors, O) | cosine law, full hemisphere | ~97 % discarded |
| ions | sigma ~0.83 deg (corrected lift) | ~intact |

## Measured, before and after

Fully exposed mask-top face, delivered / source (must be 1.0), and the global
opaque-cell flux balance (must be 1.0), at 45 nm neck, dx = 0.01:

| quantity | before | after | required |
|---|---|---|---|
| mask-top delivered / source | **0.0281** | **1.0021** | 1.0 |
| global sum(delivered x area) / (source x source area) | **0.0206** | **1.0000** | 1.0 |

A 36x under-delivery of thermal neutrals against a near-intact ion beam.

Per-band audit with the corrected transport (`scripts/lip_deposition_audit.py`,
lateral mask faces, area weighted):

| band (nm) | wall tilt | delivered/source | A(n) | visibility | O/dep isotropy | s_eff | x_xl | meas/analytic slot | removal/dep | net (nm/s) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0-50 | 0.47 | 0.3719 | 0.5016 | 0.742 | 1.0000 | 0.0586 | 0.163 | 1.151 | 0.2163 | -1.768 |
| 50-100 | 2.24 | 0.1938 | 0.5180 | 0.375 | 1.0000 | 0.0487 | 0.427 | 1.074 | 0.5130 | -0.478 |
| 100-150 | 6.78 | 0.1247 | 0.5585 | 0.224 | 1.0000 | 0.0654 | 0.407 | 1.123 | 0.9931 | -0.006 |
| 200-270 | 4.79 | 0.0313 | 0.5281 | 0.057 | 1.0000 | 0.0786 | 0.228 | 2.083 | 0.8836 | -0.030 |

`A(n) = sum_s w_s max(-d_s . n, 0) / |d_s,z|` is the unobstructed arrival factor
built from the boundary's own quadrature nodes: 1 for an exposed horizontal face,
1/2 for an exposed vertical wall.  `visibility = (delivered/source) / A(n)` is
therefore a pure occlusion fraction in [0, 1].

Three things the table settles:

1. **No over-delivery.**  Visibility is 0.74 at the top band falling to 0.06 at
   the neck -- monotone, bounded by 1, and within 7-15 % of an independent
   analytic parallel-slot reference at the three mask bands (the 2.08 at
   200-270 nm is where the parallel-wall approximation is worst, not where the
   gather is worst).
2. **Sticking fires the published rows.**  `s_eff` = 0.049-0.079 sits exactly on
   the fresh/crosslinked blend of 0.1 and 0.02 at the measured crosslinked
   fractions (0.16-0.43); nothing anomalous, nothing imported.
3. **Isotropy is exact.**  O and the depositors deliver in a face-by-face ratio
   of 1.0000, so the geometry-free 0.1953 O share
   (`RESULTS_O_CHANNEL_2026-08-04.md`) survives untouched.

## What this invalidates

Everything the probe reported before this fix is quantitatively unsafe, because
the distortion is in the ion:neutral *ratio*, not a common scale factor:

- "removal is 85.5 % of deposition at the 39 nm neck"
  (`RESULTS_MOUTH_EQUILIBRIUM_PROBE_2026-08-02.md`)
- "no equilibrium aperture exists; the neck closes at every width" (same)
- the wall-angle sweep and its "no sign change over 350x modulation"
  (`RESULTS_WALL_SLOPE_FALSIFICATION_2026-08-04.md`)

The previously recorded methodology caveat -- probe absolute rates ~10x below
what evolution requires -- is now explained and superseded: the mechanism is the
36x neutral truncation, and it moves ratios, not just magnitudes.  With correct
transport the top-band net closure is -1.77 nm/s, the right order for the
~0.33 nm/s the evolution run needs, instead of the -0.034 nm/s the broken probe
reported.

**Not invalidated:** the 10.8x top-band closure excess and the top-150 nm
localisation, which come from the *evolution* runs regraded against Krüger's
digitised profile (`RESULTS_NECK_REGRADE_2026-08-04.md`, `fcc98eb`).  The
production pilot never had this defect -- `--radiosity-backend
deterministic_extruded_2d` sets periodic neutral transport itself -- so no box
result is affected.

Independently, the corrected probe *reproduces* that localisation from a
different instrument: with the transport fixed, the 100-150 nm band sits at
removal/deposition = 0.993 and the 200-270 nm band at 0.884 (both near balance),
while the top band sits at 0.216 and closes hard.

## Gates

`tests/test_periodic_neutral_delivery.py` (4 tests) pins the transport contract
from both sides so the two paths can never be silently swapped again:

- an exposed face under periodic transport receives exactly the source flux
  (rtol 1e-6), and the opaque-cell flux balance closes;
- delivered species ratios equal source ratios face by face (the invariant the
  geometry-free O share depends on);
- the non-periodic finite-source path is *confirmed* to truncate thermal
  neutrals while leaving the ion beam comparatively intact -- the angle
  selectivity is asserted, not merely described.

The directive's `<= 1.3x` top-band closure gate is **not met**: no change was
made to the physics in this pass, because the measurement showed the deposition
side to be faithful.  The 200-270 nm control band was likewise not disturbed.

## Next

The deposition and removal sides are now both audited and both faithful to the
published mechanism, and the transport that feeds them is verified against
analytic references.  What remains unexplained is a *geometric* difference: at
the top band our evolved profile forms a high constriction (tilt 17.3 deg at
0-50 nm falling to 2.2 deg by 100-150 nm, `RESULTS_WALL_SLOPE_FALSIFICATION`)
where Krüger's narrows monotonically.  The wall-angle sweep must be re-derived
with the corrected transport before any further mechanism hypothesis, since the
version that falsified the angle theory ran on the truncated neutral field.
