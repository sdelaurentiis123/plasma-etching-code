# SF6 evaluated electron-collision audit

The exact NIST recommended/suggested/deduced tables are now replayable in SI
units. The aggregate deck avoids total-scattering double counting and carries
explicit vibration, neutral-dissociation, ionization, and attachment closures.

The source-derived attachment-rate replay has a median absolute residual of
`15.46%` (maximum
`39.04%`). This is a
source-consistency check, not independent validation.

The raw flux-drift comparison has a median residual of
`14.32%`. It is deliberately
not graded as measurement-equivalent: in strongly attaching SF6 the solver's
local flux drift/temporal growth and the source's swarm drift/spatial Townsend
observables diverge. The corresponding critical fields are
`441.1 Td` temporal versus
`359.3 +/- 3.0 Td`
spatial. Closing that distinction requires a spatial-growth/bulk-transport
solver, not parameter fitting.

Representative-grid maximum changes are
`0.000%`
for flux drift and
`0.000%`
for attachment rate.

This authorizes SF6 as a bounded component of the mixed-gas EEDF. It does not
by itself authorize a unique Oxford reactor state, wafer flux, or feature depth.
