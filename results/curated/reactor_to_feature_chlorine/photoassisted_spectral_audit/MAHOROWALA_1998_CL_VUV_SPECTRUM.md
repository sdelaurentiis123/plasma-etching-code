# Mahorowala 1998 line-resolved chlorine VUV audit

## Verdict

The atomic-Cl source does not populate the band where the large 90--244 Si/photon yield was measured, so that measurement cannot close Mahorowala depth. The source instead predicts strong 109--111 and 118.88 nm lines, but their wavelength-resolved chlorinated-Si yield is unpublished. Depth is therefore experimentally blocked at a much narrower interface: those line yields, their RF anti-synergy, and feature-floor photon delivery.

The conserved non-photon full-feature board remains at 15.265% MAPE. This is the deterministic 40 nm feature calculation, not the planar surface-rate projection. Passing only the calculated 104.82--106.67 nm Cl I source through Du's measured yield gives a MAPE interval of 15.254--15.236%. It does not close the board.

| band nm | calculated direct-coronal rate coefficient cm^3/s | yield at unit delivery | residual rate / target-depth direct-floor flux | status |
|---|---:|---:|---:|---|
| 104.82--106.67 | 4.23e-15--7.17e-15 | measured 90--244 | not a reconstruction | supported but negligible source |
| 106.67--112 | 7.97e-12--1.43e-11 | 3.5--110 | 10.8--589 Si/photon | source sensitivity; surface response unmeasured |
| 112--120 | 1.53e-11--2.52e-11 | 1.89--59.8 | 5.86--321 Si/photon | source sensitivity; surface response unmeasured |

The dominant calculated shortwave lines are 118.88, 109.74, and 110.75 nm. OPEN-ADAS raw data are not redistributed; the exact physical-record hashes and license boundary are retained in the JSON.

The dominant 118.88-nm line is now propagated with conservative partial-frequency redistribution, ground-fine-structure absorption, alternate radiative branches, and a finite 200-mm wafer. Its wafer flux spans 4.08e+13--1.58e+14 cm^-2 s^-1. Closing only the positive depth residual would require 2.15--67.8 Si/photon at the wafer, or 6.64--365 Si/photon after the target-depth direct-floor sensitivity. Those are measurement targets, not fitted yields; velocity-changing and nonradiative collision frequencies are still unset.

The released 118.88-nm transport is quadrature checked: the 12x12 production surface/direction rule changes by only 0.160% against the 16x16 rule (1% gate: `True`). The finite-wafer radial discontinuity is split exactly rather than crossed by an unsplit Gaussian panel.

The exact absorbing-ray cylinder-to-line integral delivers only 0.169--0.323 of wafer-plane photons to the target-depth floors. This is a useful deterministic transport sensitivity, not the final photon boundary: the 310 nm opening is only 2.6--2.8 wavelengths across at the dominant lines, so wave-optical validation remains mandatory.

A strictly additive photon channel also cannot be the complete board repair: runs 6, 9, 13 are already at or above their measured depths. Any successful spectrum-resolved mechanism must reproduce the condition dependence and RF ion/photon anti-synergy, not add one global depth offset.

## Experiment that decides depth

- absolute line-resolved 109.2--110.9 and 118.88 nm wafer flux under each target condition
- chlorinated-poly-Si photoetch yield at 109.7, 110.8, and 118.9 nm at 60 C
- the same yields under the target 13.56 MHz RF ion waveform to resolve PAE/IAE anti-synergy
- electromagnetic validation of the direct geometrical floor sensitivity for the 310 nm opening (2.6--2.8 wavelengths wide)

No target depth was used to select a reactor, atomic, radiative, or surface parameter. The reported unmeasured-band yields are measurement targets, not fitted constants. Both a unity-delivery lower bound and a target-depth absorbing-ray sensitivity are shown; neither substitutes for wavelength-resolved electromagnetic feature validation.
