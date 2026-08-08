# Lam chlorine attachment-energy support audit

**Claim class: measured-Te support diagnostic; not a power or depth result**

The 500-dpi-audited NIST Table 16 supplies 42 Cl2 dissociative-attachment
cross sections over `0.05--11.8 eV`. All 62 Lam Figure-3 electron-temperature
markers (`1.2801--3.7783 eV`)
were evaluated without extrapolating that table.

At the inherited `1e-6` support tolerance, particle-rate support is complete
for `0/62` markers and incident-energy support is complete
for `0/62`. The missing constant-cross-section kernel
above 11.8 eV spans `0.00101384--0.181503`
for `<sigma v>` and
`0.00522943--0.396187`
for `<sigma v E>`. These are EEDF-kernel exposure diagnostics, not bounds on
the unknown cross-section-weighted error.

On printed support alone, the collision-conditioned incident energy spans
`2.18299--4.8142 eV`. That number is
reported as a partial moment, not substituted for a complete attachment
energy loss.

## Verdict

Table 16 advances the chlorine ledger because it separates the attachment
particle moment from the electron-removal energy moment. It does **not** close
the Lam electron-power balance: both tails remain exposed on every measured
temperature marker, and elastic, vibrational, non-dissociative electronic,
detachment, ion-pair, and molecular branching channels remain open.

No coefficient was selected against a reactor observable or feature depth.
No absorbed-power, wafer-flux, or etched-depth prediction is supported by this
audit.
