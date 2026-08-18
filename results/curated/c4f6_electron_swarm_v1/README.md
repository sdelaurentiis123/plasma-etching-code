# C4F6 electron-collision source-replay audit

Lan--Jeon's printed C4F6 tables are now checksum-locked and consumed directly
by the deterministic two-term solver. The 18 visually audited pure-C4F6
Figure-7 markers span 121--1197 Td.

Against that board, local flux drift has a mean absolute residual of
`7.82%` and a maximum of
`14.87%`. This is not
graded as a failed cross-section fit: the paper reports pulsed-Townsend average
drift `Wv`, while the reactor solver returns flux drift. Those observables
diverge in a gas with attachment and ionization.

Doubling the energy-grid density changes flux drift by at most
`0.078%`
and aggregate dissociation rate by at most
`0.794%`.
The remaining gap is therefore physical/observational, not grid error.

This authorizes the deck as a bounded C4F6 component for local EEDF and
aggregate-rate calculations. It does not resolve product branching, a Krueger
reactor state, wafer flux, or feature depth. The next gates are deterministic
density-gradient/PT transport and mass-resolved reactor validation against
Benck's ion-current board.
