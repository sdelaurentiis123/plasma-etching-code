# Zhu / Oxford NPG80 conditional square-pillar profiles

This directory freezes the first evolving-profile board for the target-free
Oxford PlasmaPro NPG80 condition supplied by Freddie Zhu. It is a deterministic
conditional envelope, **not** a fitted SEM and **not** a certified Oxford
absolute-depth prediction.

## Frozen condition

- 700 nm ALD TiO2 on fused silica with a 45 nm Cr hard mask
- 55 / 5 / 1 sccm CHF3 / SF6 / O2
- 30 mTorr, 150 W forward RF, 20 C, 20 minutes
- square-pillar prior: 400 nm pitch and widths 80--320 nm in 40 nm steps
- four target-free ion transport sensitivities: low/high inferred sheath energy
  crossed with 0/0.65 two-component angular-tail fraction
- two TiO2 rate endpoints, 34.125 and 43.4667 nm/min, from the audited
  Janissen cross-machine process analog

The width/pitch board is a same-group public geometry prior, not a claim about
the unreleased target GDS. The TiO2 rates are cross-machine evidence, not fitted
target coefficients. The target SEM, target depth, and target-selected
coefficient are absent by construction.

## Result

The 28 independent rate-normalized trajectories produce 56 reported endpoints.
Because the declared surface law is autonomous and linear in blanket removal
rate, the two rate endpoints share one geometry trajectory at exact equal dose;
no profile interpolation is used.

| nominal width (nm) | cleared endpoints / 8 | 20 min depth range (nm) | interpretation |
|---:|---:|---:|---|
| 80 | 8 / 8 | 700 | clearance in every conditional case; lateral CD unresolved at 20 nm mesh |
| 120 | 8 / 8 | 700 | clearance in every conditional case |
| 160 | 8 / 8 | 700 | clearance in every conditional case |
| 200 | 8 / 8 | 700 | clearance in every conditional case |
| 240 | 8 / 8 | 700 | clearance in every conditional case |
| 280 | 6 / 8 | 622.6--700 | low-rate, broad-tail cases do not clear |
| 320 | 4 / 8 | 567.4--700 | high-rate cases clear; all low-rate cases do not |

For widths through 240 nm, conditional clearance occurs between approximately
920 and 1182 s. At 320 nm the high-rate cases clear between approximately 982
and 1006 s, while the four low-rate cases retain 107--133 nm of TiO2 after the
full 1200 s.

When clearance is detected, `profile.etched_depth_nm` is capped at the known
700 nm film thickness, the time is stored as a bracket, and the CD/sidewall
geometry is explicitly the **last pre-clearance** surface. The fused-silica
post-clear profile is not identified by this one-material law. When clearance
is not detected, the receipt contains the 1200 s endpoint geometry.

## What is physical and what is not yet certified

The solver preserves real geometric recession, corner rounding, bowing,
shadowing, and possible failure to clear. Particle accounting closes below
`2e-12`; the largest conservative surface-state remap residual is about
`3.2e-15`.

The production mesh is 20 nm and is only the first resolution rung. A CD is
marked resolved only when both the vertical relief is at least two cells and
the narrower of pillar width and neighbor gap spans at least six cells. Thus
the 80 nm pillar and the reciprocal 80 nm gap at 320 nm remain CD-unresolved.
The 80 nm high-tail case develops a grid-scale x/y difference even though the
continuum problem is symmetric; that raw result is retained as a numerical
warning rather than smoothed into a physical claim. A finer-grid sentinel is
required before quantitative CD certification.

The Cr mask is pinned in this evolving-profile rung. Mask survival is handled
separately by the blind selectivity board because no target-condition
species/energy-resolved Cr erosion law exists. A profile after conditional mask
exhaustion must not be interpreted as the target's final shape.

## Why the surface coefficient is still conditional

The open literature supplies tool-level TiO2 etch rates and selectivities, but
not the species-, energy-, angle-, coverage-, and material-state-resolved law
for ALD TiO2 and Cr under this exact CHF3/SF6/O2 boundary. This board therefore
transfers only an explicitly labeled rate interval. It omits neutral-radical
response, fluorocarbon polymer/passivation kinetics, TiO2 fluorination state,
species-dependent sputter yields, charging feedback, redeposition, and Cr
shape evolution. Those omissions are physical uncertainty, not optimizer
degrees of freedom.

The highest-value same-run measurements remain achieved DC self-bias or the
electrode waveform, blanket TiO2 loss, and residual Cr thickness. The target
SEM is an answer key for the frozen prediction, not a simulator input.

## Reproduction

```bash
python scripts/audit_zhu_npg80_conditional_profiles.py --check
pytest -q tests/test_zhu_npg80_conditional_profiles.py
```

`audit.json` is the aggregate receipt. `trajectories/` contains all 28
specification-bound, checksum-pinned trajectory caches so the board can be
audited or resumed without rerunning completed cases.
