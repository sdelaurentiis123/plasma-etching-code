# Bosch identifiable recipe-path memory depth extension v7

This freezes version 7 before its operator or coefficients are implemented.
Calibration outcomes and the failed v6 identifiability audit are visible. The
chronological heldout outcomes remain unopened.

## Why v6 is not the answer

Version 6 demonstrated that a carried chamber state improves mean absolute
depth and the within-lot depth drift. It did not identify separate C4F8
deposition, SF6 cleaning, and wall-response coefficients: the fit contacted two
dynamic bounds, its Jacobian condition number was 6.76e6, and its maximum
parameter correlation was 0.999997. The process traces explain why. Across the
calibration set, normalized C4F8 and SF6 doses and their ratio vary by only a few
parts per thousand. This experiment cannot resolve three independent rates.

Version 7 therefore represents only the combination the data can identify: net
wall loading accumulated along the observed, nearly fixed Bosch recipe path. It
does not claim that the state is literal fluorocarbon coverage or that separate
deposition and cleaning rates have been measured.

## Frozen state and response

Let `d_i` be the nonnegative integrated C4F8 flow of process trace `i`, divided
by the frozen calibration median dose. Each declared conditioning sequence sets
`H_start=0`. During a production wafer,

```
H_mean = H_start + 0.5*d_i
H_end  = H_start + d_i
```

`H_end` is carried into the next production wafer in the same lot, including
the one processed wafer that lacks an outcome. The recorded one-minute
plasma-off breaks hold `H` constant. There is no fitted initial state and no
wafer-number, lot-number, date, or outcome input.

The applied neutral wall-loss multiplier is

```
clip(exp(b_repeat*log(repeat_count/3)
       + b_Si*I[conditioning on blank Si]
       + b_SiO2*I[conditioning on blank SiO2]
       + b_H*H_mean), 0.25, 4.0)
```

The multiplier acts only on the existing zero-dimensional and cylindrical
neutral upper/sidewall loss closure. Positive ions, lower-wafer collection, and
the frozen Belen and La Magna/Garozzo surface laws remain unchanged.

## Scope and falsification

The four shared coefficients are the three existing static preparation
coefficients plus `b_H`. The frozen bound is `|b_H| <= ln(4)/10`, so ten
reference-dose production wafers cannot introduce more than the existing
factor-four response before clipping.

This is a recipe-path model valid only over the measured calibration ratio
domain `5.635257777059961 <= D_SF6/D_C4F8 <= 5.737394214801463`. It must refuse extrapolation
outside that domain. A varying-ratio experiment or direct wall diagnostic is
required to replace it with independently identified deposition and cleaning
physics.

Whole-lot leave-one-out refits, the unchanged empirical baselines, the
four-parameter identifiability gates, midpoint interpolation validation, exact
selected-parameter replay, and certification-grid refinement remain mandatory.
The response table uses thirteen exact Chebyshev nodes because the v6 nine-node
table missed the frozen interpolation tolerance. Failure of any gate leaves the
heldout outcomes sealed.

Primary topology support remains the Sayyed et al. protocol
(DOI `10.5281/zenodo.17122442`) and Schaepkens et al.'s independent
fluorocarbon-ICP wall study (DOI `10.1116/1.581316`). No numerical wall rate is
imported from the latter.
