# Bosch current-conserving wafer-boundary map extension v8

This freezes version 8 before its operator, response table, coefficients, or
selection audit are implemented. Calibration outcomes and the v7 residual
decomposition are visible. The chronological heldout outcome remains unopened.

## Scientific hypothesis

Version 7 identifies a shared recipe-path wall-memory law and transfers wafer
mean depth across left-out lots. It still loses to the empirical mean-map
baseline on point and normalized-shape scores. The calibration-only
decomposition shows that one repeatable spatial map contains 86.2018% of the
squared v7 map residual. Its leading process-dependent component is
edge-positive and center-negative and tracks measured C4F8-phase platen Vpp.

Version 8 tests one specific physical interpretation: the remaining term is a
machine-specific positive-ion transmission fingerprint produced by the
electrode, chamber, sheath, and focus-ring boundary, plus a small Vpp-dependent
edge-focusing perturbation. It is not a surface-yield correction and it is not
an additive depth map.

## Frozen operator

On the official circular wafer, let `rho=r/R` and `phi=atan2(y,x)`. Let
`Z_k(rho,phi)` be the complete real-Zernike basis through a candidate maximum
order, excluding the piston. For wafer `i`,

```text
z_i = (Vpp_RMS_i - 637.4409584828442 V) / 3.816305957878358 V
g_i(rho,phi) = sum_k a_k Z_k + z_i sum_l b_l Z_l
u_i(rho,phi) = exp(g_i)
T_i = u_i * integral(Jion_base dA) / integral(Jion_base*u_i dA)
Jion_v8 = Jion_base*T_i
```

The static basis candidates are complete orders 1 through 10. The dynamic
basis candidates are absent or complete through order 2, excluding the piston.
The same coefficients apply to every calibration and heldout wafer. The
normalization is evaluated on the deterministic cylindrical finite-volume
wafer grid using the baseline ion-current density, not the 89 measurement
points. It must conserve total positive-ion current to relative error at most
`2e-14`.

Only the positive-ion spatial channel changes. Neutral F and film-precursor
fluxes, wafer-average positive-ion current, ion energy, v7 wall memory, Belen
silicon kinetics, and La Magna/Garozzo oxide/film kinetics remain unchanged.
The exponential guarantees positivity. There is no clipping in the forward
operator: a coefficient or field outside the frozen domain raises.

## Frozen domains and complexity

- Static and dynamic coefficients are each bounded to `[-0.1, 0.1]` log units.
- Across a dense disk certification grid and the complete frozen Vpp domain,
  `max(abs(g_i)) <= ln(2)` is mandatory.
- The measured C4F8 platen-Vpp-RMS domain is
  `[626.9533265149638, 643.534317555529] V`; application outside it is refused.
- The static tool map is valid only for this SPTS tool, wafer coordinate frame,
  electrode/focus-ring configuration, and measurement convention.
- Candidate selection chooses the lowest total coefficient count that passes
  every gate. Ties choose the lower static order and then omit the dynamic
  field. If no candidate passes, v8 fails and the heldout remains sealed.

The operator and conservation law are portable. The fitted coefficients are a
tool calibration, not universal plasma constants. Complete order 10 contains
65 static non-piston coefficients; its availability is a falsification ceiling,
not permission to relabel an arbitrary 89-point correction as physics.

## Calibration and whole-lot validation

The four v7 wall coefficients are not jointly traded against the spatial map.
The full-calibration v8 fit uses the already identified full v7 coefficients.
For each whole-lot leave-one-out fold, the v7 coefficients are first refit on
that fold's training lots using the frozen v7 procedure; they are then held
fixed while only v8 map coefficients are fit on the same training lots. No
test-lot outcome enters either stage.

Each candidate must satisfy all of the following:

1. All unchanged absolute silicon-mean, point, normalized-shape, oxide, and
   selectivity gates pass.
2. Aggregate whole-lot silicon-mean MAE beats `0.338486 um`.
3. Aggregate whole-lot silicon-point RMSE beats `0.486585 um`.
4. Aggregate whole-lot normalized-shape RMSE beats `0.636619%`.
5. Within-lot slope MAE does not exceed the v7 value `0.082903 um/wafer`.
6. The selected map Jacobian has full column rank, condition number at most
   `1e6`, maximum pairwise parameter correlation below `0.995`, and no
   coefficient within 0.1% of a bound.
7. Across folds, median coefficient standard deviation is at most `0.005`,
   maximum coefficient range is at most `0.05`, and maximum Zernike design
   condition number is at most `1000`.
8. Independent response-table midpoint interpolation error is below 0.05 of
   every frozen absolute gate.
9. The selected full fit and every reported seal score are replayed through the
   exact reactor, cylindrical transmission, and surface recurrence.
10. Cylindrical and Zernike-grid refinement change every observable by less
    than 0.25 of its frozen gate.

The prior all-calibration residual search over 131 machine features is declared
selection exposure. V8 freezes only C4F8 platen Vpp RMS as the dynamic input;
no second feature may be substituted after seeing v8 scores.

## Deterministic acceleration

Calibration search may use a tensor table over 13 Chebyshev log wall-loss
nodes in `[ln(0.25), ln(4)]` and 13 Chebyshev log local-ion-factor nodes in the
same interval. At each node, the unchanged surface recurrence is evaluated
from a physical reactor boundary. Pointwise interpolation is valid because the
surface states at the 89 wafer coordinates are independent given the boundary.
The table is an acceleration only; it cannot supply the final replay or heldout
seal.

## Target firewall and failure action

The heldout numeric outcome fields remain forbidden until every calibration,
whole-lot, identifiability, stability, interpolation, exact-replay, refinement,
hash, and runtime-audit gate passes. Failure leaves
`heldout_outcomes_read=false`, `heldout_prediction_written=false`, and
`eligible_for_prediction_seal=false`.

Forbidden shortcuts include direct depth corrections, changing surface yields,
changing total ion current, fitting per-lot/per-wafer maps, using wafer number
or outcome as an input, swapping the Vpp feature after scoring, accepting a
nonconservative 89-point normalization, or opening the heldout because an
output-space proxy happened to beat a baseline.

Primary topology support is Babaeva and Kushner's focus-ring sheath analysis
(DOI `10.1088/0022-3727/41/6/062004`) and Seong et al.'s experimental wafer-edge
focus-ring study (DOI `10.3390/nano12223963`). No numerical map coefficient is
imported from either source.
