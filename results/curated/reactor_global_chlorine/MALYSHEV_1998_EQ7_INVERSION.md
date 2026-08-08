# Malyshev 1998 Lam measured-state Eq.-7 inversion

**Verdict: diagnostic closure, not a wall fit or depth prediction**

The audit joins three independently digitized boards from the commercial Lam
Alliance chlorine reactor: Figure 3 electron temperature, Figure 11
volume-average electron density, and Figures 7--8 relative Cl2 density. Of 38
audited dissociation markers, 23 have physical dissociation and supported Te
and ne at the same gap, pressure, and power.

For each supported marker, the audit evaluates Hamilton's eight-state
Maxwellian neutral-dissociation rate, adds the retained Lee--Lieberman
dissociative-attachment rate, and inverts Malyshev Eq. 7:

`relative Cl2 = 1 / (1 + kd ne / (2 kr))`.

Every row reproduces its measured relative-Cl2 marker algebraically. The
required wall-return frequency has an all-row median of `172.5 s^-1`. Restricting
only to the 16 rows whose upper bound remains finite under the source's
non-statistical +/-25% absolute-density statement gives a median of
`148.6 s^-1` and a range of `65.1--304.8 s^-1`.

This does not identify a wall recombination probability. That mapping still
requires validated Cl-in-Cl2 diffusion at the Lam gas temperature and a
surface-state law applicable to the conditioned anodized-Al chamber. Te and
volume-average ne lack complete uncertainties, Hamilton publishes no scalar
physical uncertainty, and the attachment rate remains a compilation. No
formal pass/fail interval is therefore defensible yet.

## Source inconsistency caught

Footnote 14 prints `ne = 1e11 cm^-3`, `kd = 7e-9 cm3/s`, and `kd ne = 700
s^-1`. The article's own printed `kdis` law, evaluated at the highest measured
11 cm/10 mTorr Te and augmented by the article's maximum stated 1/7 attachment
contribution, gives at most `1.13e-9 cm3/s` or `113 s^-1`. The nominal footnote
rate is higher by `6.18x`; it is quarantined from calibration.

Neither the inversion nor this source audit supports local sheath density,
species-resolved wafer flux, etch rate, or feature depth.
