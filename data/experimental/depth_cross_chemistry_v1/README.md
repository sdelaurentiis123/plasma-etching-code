# Value-blind absolute-depth campaign

This directory freezes the targets for the next depth campaign before new
Yoshie Figure 4--6 numbers are digitized or used by the model.  Run:

```bash
python scripts/validate_depth_cross_chemistry_preregistration.py
```

The protocol commits source checksums, panels, process conditions,
calibration/held-out roles, chemistry families, boundary-evidence tiers, and
admissible parameter bounds.  Numerical observations are absent by design.
The later value reveal must contain every and only these target IDs.

## What “predict depth” means here

The target is a fixed-duration feature depth predicted without fitting that
feature endpoint.  Agreement earns only the evidence tier of its incident
boundary:

- **A — measured species, energy, and angle:** the incident beam is measured
  with enough resolution to exercise surface physics directly.
- **B — facility-conditioned:** a boundary scale may be inferred from an
  independent blanket observable in the same reactor, then tested on held-out
  features.  This is useful process prediction, but it is not a
  first-principles knobs-to-flux result.
- **C — underidentified:** process knobs or an aggregate simulated flux are
  printed, but species-resolved wafer flux/IEAD information needed by the
  surface law is absent.  A numerical match can be reported only as
  retrospective agreement; it cannot pass formal predictive validation.

Krüger and Jeong are Tier C at preregistration.  No flux scale may be fit to
their feature depths.  Yoshie is Tier B only if its Figure-4 blanket data are
used to condition the same-reactor boundary before Figures 5--6 are revealed;
otherwise it is Tier C.  The Chang beam profiles have a Tier-A incident
boundary, but the source does not print exact exposure time.  They therefore
grade geometry transfer, not an absolute time-rate prediction.

## Why three chemistries

The held-out board deliberately spans:

1. C4F6/C4F8 fluorocarbon etching of SiO2;
2. alternating C4F8 passivation and SF6 removal of Si;
3. Cl/Cl2 ion-assisted etching of Si.

This prevents a correction that merely memorizes the Krüger reactor.  The
unchanged transport/evolution core must carry all three; only chemistry and
independently evidenced boundaries may vary.

## Source and visibility record

- Karahashi, Takada, Krüger, Jeong, and Chang numerical values were already
  visible in the repository before this protocol.  They are retrospective
  calibration/transfer targets, not newly blinded observations.
- Yoshie Figures 4--6 were visually inspected to establish axes, series, and
  process conditions before this commit, but no pixel-to-value dataset was
  created and no Yoshie value was supplied to the simulator.  Their complete
  width/timing grids are committed here before digitization to prevent
  favorable-point selection.
- Official Yoshie figure URLs and exact raster SHA-256 values are committed;
  the copyrighted raster files are not redistributed.

The physics standard is species- and energy-resolved atom accounting:
adsorption, surface activation/passivation, stoichiometric substrate removal,
volatile-product routing, feature transport, and moving geometry.  A reactor
normalization cannot compensate for a wrong surface law, and a surface yield
cannot compensate for missing wafer species fluxes.
