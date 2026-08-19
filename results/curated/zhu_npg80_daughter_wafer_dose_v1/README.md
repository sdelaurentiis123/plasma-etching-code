# Zhu NPG80 daughter-reclosed radial wafer/dose board

This receipt lifts every accepted daughter-reclosed reactor state through the
fixed-topology axisymmetric transport operator and converts the centered
positive-ion dose into exact TiO2 formula-unit arithmetic.  It was generated
before receipt of the target SEM and does not select a surface coefficient
from a depth outcome.

## What the reactor now says

| absorbed power | centered 3 mm ion flux | smooth center enhancement | F/ion incident-flux ratio | 700 nm blanket yield requirement |
|---:|---:|---:|---:|---:|
| 60 W | `1.0262e19 m-2 s-1` | 0.73% | 104.5 | 1.393--1.779 formula units/ion |
| 90 W | `1.0637e19 m-2 s-1` | 1.15% | 396.8 | 1.344--1.716 formula units/ion |
| 105 W | `1.0340e19 m-2 s-1` | 1.60% | 646.9 | 1.383--1.765 formula units/ion |
| 120 W | `1.0475e19 m-2 s-1` | 2.34% | 880.3 | 1.365--1.743 formula units/ion |

The yield interval is only the published `3.25--4.15 g cm-3` ALD-TiO2 density
sensitivity.  It assumes blanket delivery.  A feature-floor transmission
`T < 1` raises the microscopic yield requirement by `1/T`.  For example, the
high-density endpoint needs an effective `yield * transmission` of roughly
`1.72--1.78` to clear the 700 nm film.  At a yield of 1.5 and full
transmission, the high-density endpoint reaches only about 590--612 nm; at a
yield of 2.0 and full transmission both density endpoints clear.

Those are conditional predictions, not fitted uncertainty bars.  The receipt
contains the full `yield x transmission` depth board so the future SEM can
grade the physics without rewriting the preregistration.

## What this changes about the SEM forecast

The original frozen binary call, "full film clearance expected," remains in
its historical receipt.  The corrected daughter-reclosed board does not
support that call across the whole admissible surface/transport space.  It
supports a narrower statement: full clearance is physically reachable, but it
depends decisively on the unmeasured TiO2/Cr surface response and feature-floor
delivery.  The independent Janissen and Hong TiO2 boards keep clearance
plausible; they do not identify this Oxford/ALD-TiO2 coefficient.

The smooth axisymmetric model predicts only a 0.7--2.3% center enhancement.
Spatially clustered fallen pillars in a micrograph are therefore not, by
themselves, evidence for a large reactor-scale radial flux gradient.  Layout
microloading, local CD/lithography variation, Cr-mask loss, adhesion, and
post-etch drying/capillary collapse remain competing explanations.  Mapping
the SEM field to its exact GDS geometry and wafer position is required to
separate them.

## Remaining decisive measurements

The smallest high-value target-run additions are:

1. achieved DC self-bias/readback during the 20 minute etch;
2. exact GDS/CD, pitch, layout family, and chip/wafer position for the SEM;
3. scale-bearing cross-section plus top-down SEMs, with whether Cr was stripped;
4. ideally, blanket TiO2 loss and residual Cr on the same run.

Rebuild/check with:

```bash
python scripts/audit_zhu_npg80_daughter_wafer_dose_board.py
python scripts/audit_zhu_npg80_daughter_wafer_dose_board.py --check
```
