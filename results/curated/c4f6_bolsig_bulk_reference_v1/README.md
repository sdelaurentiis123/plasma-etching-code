# C4F6 BOLSIG+ bulk-transport reference

The exact Lan--Jeon collision tables were exported without retuning and run in
the official deterministic BOLSIG+ `03/2016` density-gradient mode at all 18
Figure-7 fields. This is an independent comparator, not a production
dependency.

BOLSIG+ flux drift agrees with petch flux drift to at most
`0.234%` across the board.
That independently corroborates the local two-term implementation.

Changing the comparison quantity from flux to BOLSIG+ bulk drift does **not**
reproduce the plotted legacy PT average drift `Wv`: the mean absolute residual
is `8.18%`
and the maximum is
`29.08%`,
versus `7.77%`
for the independent flux replay. Bulk drift overcorrects the low-field,
attachment-dominated points and improves much of the high-field board.

This is a physical observable-definition result, not a numerical failure.
`casey-2021-pt-foundations` shows that the older PT transport property transforms as
`W_B = W_tilde + alpha_tilde D_L_tilde`. Lan--Jeon reports `Wv` but not the
same-study Townsend property and longitudinal diffusion needed for that
transformation. The source cross sections must not be retuned to force either
flux or modern bulk drift through a legacy quantity with missing co-observables.

The result validates the collision solver's local flux calculation. It does
not identify a C4F6 reactor state, wafer flux, or Krueger depth.
