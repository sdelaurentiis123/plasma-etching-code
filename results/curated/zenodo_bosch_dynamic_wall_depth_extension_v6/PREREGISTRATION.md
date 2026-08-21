# Bosch dynamic chamber-wall depth extension v6

This freezes the next Bosch model form before its operator or coefficients are
implemented. Calibration outcomes and the v5 residual audit are visible; the
chronological heldout outcomes remain unopened.

## Source correction

The dataset authors define `C` as the bare system **chuck**, not carbon. Their
conditioning sequence is two minutes of O2 plasma followed by two minutes of
O2/SF6 plasma, repeated one, three, or nine times on the chuck, a blank silicon
wafer, or a blank silicon-dioxide wafer. Ten Bosch wafers then run sequentially
without intermediate cleaning.

The v5 variable name `log_carbon_cycle_coefficient` was therefore semantically
wrong. Its arithmetic used the repeat count and remains numerically intact.
Version 6 calls that feature the conditioning-repeat count.

## Missing physics being tested

Every calibration lot loses measured silicon depth as wafer sequence advances.
The static v5 reactor predicts flat or rising depth. That coherent signed
residual, together with the authors' no-intermediate-cleaning protocol, supports
a between-wafer chamber state.

Version 6 introduces one bounded fluorocarbon wall occupancy relative to the
declared post-conditioning state. C4F8 dose fills unoccupied wall state; SF6
dose removes occupied state. The exact first-order update carries the end state
from one wafer into the next. The state is reset only by a declared conditioning
sequence, never from wafer number or a fitted lot offset.

The state changes the same neutral upper/sidewall loss closure used by v5. It
does not change positive ions, the lower-wafer collection boundary, or the
Belen and La Magna/Garozzo surface laws. The three static preparation
coefficients and three dynamic coefficients are shared by every lot.

## Why the model remains falsifiable

The SF6 and C4F8 doses vary only weakly across this dataset, so deposition,
cleaning, and wall-response coefficients may be poorly identified. A full-rank
Jacobian, bounded condition number, parameter-correlation limit, and bound
contact rule are frozen. Failure of that gate refuses a prediction seal rather
than manufacturing a unique wall history.

Whole-lot leave-one-out refits are mandatory. The physics path must beat the
same global-depth and mean-map baselines, improve the signed within-lot drift,
pass refinement, and be replayed through the exact reactor after any
interpolation-accelerated search. Only then may the chronological heldout
outcomes be opened.

Primary topology support is the Sayyed et al. target-tool protocol
(DOI `10.5281/zenodo.17122442`) and Schaepkens et al.'s independent
fluorocarbon-ICP wall study (DOI `10.1116/1.581316`). No numerical wall rate is
imported from the latter.
