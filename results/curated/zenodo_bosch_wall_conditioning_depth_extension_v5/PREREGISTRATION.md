# Bosch wall-conditioning depth extension v5

This freezes the only new Bosch model form allowed after the calibration-only
v4 audit.  Higher positive source harmonics did not beat the leave-one-lot-out
empirical map baseline, while the official metadata declares systematic
carbon-cycle, Si, and SiO2 chamber preparation classes.

One coefficient triplet is shared across all wafers.  It maps the declared lot
type to a bounded neutral wall-loss multiplier.  The multiplier acts inside the
reactor: it shortens the effective F and film-precursor lifetimes and increases
their upper/sidewall loss velocities.  It leaves positive ions, lower-wafer
collection, and both surface mechanisms unchanged.  No depth correction or
per-lot offset is permitted.

The factor-four multiplier interval matches the neutral pressure-lifetime
interval already frozen in the reduced reactor.  The 3C state is the exact
unit-multiplier reference.  Model selection remains leave-one-entire-lot-out,
and the heldout outcome firewall remains closed until the physics model beats
both empirical baselines, passes grid refinement, and is hash-sealed.

Calibration outcomes and the prior model-form failures were visible before
this freeze.  No conditioning coefficient was fitted and no heldout numeric
outcome was examined before this file was written.
