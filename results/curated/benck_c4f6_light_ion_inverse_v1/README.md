# Benck C4F6 light-ion inverse

This target-free audit asks whether the smallest direct-plus-secondary C4F6
ion-source model can be mapped to Benck's measured CF+, CF2+, and CF3+ wall
currents using one common loss factor. It uses the paper's independent neutral
Figure-14(a) neutral ratio at each feed condition and solves for the remaining
parent source and `n(CF3)/n(CF)` at 2--6 eV.

The algebra replays the two measured ion ratios to within
`1.11e-16`. All three
Ar-containing feed entries require a negative inferred CF3 density throughout
the declared temperature grid. This remains true if every line-of-sight
CF2/CF ratio is halved, addressing the paper's warning that the local plasma
ratio is probably lower. No negative density was clipped or accepted. This is
a model-form failure: one common source-to-wall-current loss operator cannot
explain the co-conditioned neutral and ion feed boards together.

The next reactor rung must evolve species-dependent Bohm/wall/exhaust losses,
ion-neutral conversion, heavy-fragment cascades, and surface-product return.
This audit does not decide which one dominates, does not use absolute current
scale or feature depth, and does not supply a Krueger boundary or wafer flux.
