# Benck C4F6 differential light-ion loss bound

The previous common-loss inverse produced an unphysical negative CF3 neutral
density in all three Ar-containing feed conditions. This audit fixes
`n(CF3)/n(CF) = 0`, identifies a nonnegative parent source only from the
independent measured `CF2+/CF+` ratio, and asks how much additional CF3+ loss
is then required by `CF3+/CF+`.

All 15 Ar-mixture temperature rows close with nonnegative sources. The required
CF3+-selective first-order loss is
`0.103`--
`0.614`
times the common light-ion wall/exhaust loss over 2--6 eV. The CF2 ratio replay
error is below `1.11e-16`.

This repairs the *sign* of the Ar-mixture inverse without clipping a density or
using depth. It does not identify a unique reaction. The pure-C4F6 condition
changes the required operator sign across the temperature grid, so one fixed
CF3+ loss coefficient is still rejected. The next forward reactor must evolve
condition-dependent ion-neutral conversion, heavy-fragment cascades, and
surface-product return, then pass both Benck feed and pressure boards.

No absolute current, Krueger depth, or feature result was used. This receipt
does not provide a Krueger boundary or wafer flux.
