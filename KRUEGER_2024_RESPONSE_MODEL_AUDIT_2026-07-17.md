# Krüger 2024 base response-model audit

Date: 2026-07-17
Status: completed bounded read-only numerical audit; no new simulation run
Scope: Krüger base-calibration opening/depth responses only

## Evidence boundary

This audit used only completed base-calibration endpoints and their base-only proposal artifacts.
The sealed oxygen, power, and transfer-profile observations were not opened or used. No simulation,
GPU job, parameter evaluation, held-out reveal, commit, or push was performed for this audit.

The calibration target remains:

- minimum mask opening: `45 nm`;
- etch depth: `825 nm`.

The response order below is always `[mask_opening_nm, etch_depth_nm]`. The parameter order is always
`[effective_mask_crosslinked_growth_fraction, oxide_etch_yield_scale]`, abbreviated `[f, s]`.

## Unique completed 10 nm base endpoints

Seven unique completed 60 s endpoints are present. The corresponding artifacts in
`results/vast_instance_45125477_archive/` are byte-identical copies of the primary files, not
independent replicates.

| Revision | Transport closure | `f` | `s` | Opening (nm) | Depth (nm) | Primary audit SHA-256 |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| R11-f0 | single published plane | 0.000000000 | 1.000000000 | 13.166934 | 980.933485 | `8091daeebb027eef72a774cca3b73ba42ecc10ecb649b57107bc671da18fbc65` |
| R11-f0.693 | single published plane | 0.693139320 | 1.000000000 | 14.686674 | 1191.470019 | `efe63f087f4989c8217dcc27a104df916746c204b2775c8d79ca980dec9c913f` |
| R11-f1 | single published plane | 1.000000000 | 1.000000000 | 59.092861 | 1420.756575 | `669bca08157462172cd4041eb815866e1c09e6c6482adfddc4097dfd5a82c497` |
| R12 | single published plane | 0.902613908 | 0.612022523 | 39.950256 | 900.404822 | `0542c57bf96430631bf05fbf304a54a44fba07d2aff7d07ef51fd98511b8314d` |
| R15 | axisymmetric uniform, order 16 | 0.934040082 | 0.519249144 | 50.307100 | 750.205922 | `de9d5251b6b655827f4ab6fc9c8f9812070f9ec5f6493798b0b932b9e66828c5` |
| R16 | axisymmetric uniform, order 16 | 0.921017052 | 0.574877578 | 49.890804 | 853.253317 | `aef745dcbba1daaec0266e22b1003be12aa1fa834e99f400481d18902bb9b267` |
| R17 | axisymmetric uniform, order 16 | 0.893405974 | 0.566763272 | 42.846799 | 837.905130 | `034045b79cda1fa55ec87ffcc57c74cd5647e476a64f4b2c83a17034fab0f160` |

R13/R14 added azimuthal-closure diagnostics at an existing endpoint; they did not add a new
parameter-to-60-s-response observation. The archived zero-time benchmark and incomplete `running`
artifacts are not calibration endpoints. The R17 5 nm interrupted/continued trajectory is also not
part of this 10 nm response fit because it is a mixed-operator development trajectory.

## Prospective test of the old 2×2 response

The response used to propose R17 from R16 was

```text
J_old = [[ 167.787546,   31.796806],
         [ 448.047344, 1957.314277]]
```

For the actual R16-to-R17 parameter move,

```text
delta_parameter = [-0.027611078, -0.008114306]
```

the model and engine gave:

| Response change | Predicted (nm) | Observed (nm) | Observed / predicted |
| --- | ---: | ---: | ---: |
| Opening | -4.890804 | -7.044005 | 1.440255 |
| Depth | -28.253317 | -15.348187 | 0.543235 |

The prediction error `observed - predicted` was `[-2.153201, +12.905130] nm`. The direction was
useful (cosine similarity `0.96668`), but the vector error was `77.5%` of the observed response
norm. The old model understated the local opening sensitivity and overstated the local depth
sensitivity.

Using one 5 nm fine-cell width as the declared scale for both base residuals, the R16 merit was
`5.734701` and the observed R17 merit was `2.616705`. The proposal therefore achieved `54.37%` of
the predicted reduction. This is an accept-and-update trust-region result, not evidence for an
expanded or full Newton step.

## Same-operator local affine response

R15, R16, and R17 are the three completed points that share the current 10 nm axisymmetric
transport closure. The unique affine interpolation through them has Jacobian

```text
J_local = [[240.750948,   48.878080],
           [ 10.743900, 1854.938292]]
```

The two displacement directions used to identify this matrix have raw condition number `4.13`
and trust-coordinate-scaled condition number `3.57`. Response conditioning is:

| Response matrix | Raw condition number | Previous-step/trust-coordinate-scaled condition number |
| --- | ---: | ---: |
| `J_old` | 12.8484 | 5.8942 |
| `J_local` | 7.7199 | 2.2713 |

Under the same scaling, the old response columns have cosine `0.9421`, while the local response
columns have cosine `0.0709`. Thus the observed local data suggest that `f` primarily controls the
opening and `s` primarily controls depth. This is more identifiable and physically interpretable
than the old response, but it is based on exactly three points and consequently has no residual
degrees of freedom.

Solving the local affine response for the target gives

```text
f = 0.90377432
s = 0.55974604
```

This is interpolation, not extrapolation. Its barycentric weights in the R15/R16/R17 parameter
triangle are approximately

```text
[0.169274, 0.126401, 0.704325]
```

which are all positive.

## Alternative-model and leave-one-out audit

Operator-aware global regressions were examined only as uncertainty diagnostics. They mix the
historical single-plane and current axisymmetric closures, so none can supersede the local
same-operator interpolation.

| Model | Scaled design condition | Axisymmetric leave-one-out RMSE, opening/depth (nm) | Target candidate `[f,s]` | Finding |
| --- | ---: | ---: | --- | --- |
| Affine in `f,s` + operator indicator | 3.46 | 4.29 / 34.19 | `[0.84228, 0.59037]` | Candidate lies far outside the local triangle; reject |
| Quadratic in `f`, linear in `s` + operator indicator | 11.56 | 1.56 / 42.47 | `[0.89921, 0.57486]` | Depth prediction is high-leverage and unstable |
| Full quadratic in `f,s` + operator indicator | 312.57 | 3.56 / 26.25 | `[0.90406, 0.55933]` | Seven coefficients fit seven points exactly; saturated and not predictive |
| R15/R16/R17 local affine | response condition 7.72 raw / 2.27 scaled | Not internally estimable | `[0.90377, 0.55975]` | Best local model, but zero residual degrees of freedom |

A same-operator quadratic is not identifiable: a quadratic in two parameters requires six
independent coefficients per response, while only three same-operator observations exist. Adding
historical closure data can make the algebra square, but cannot make the resulting model
scientifically local or predictive.

The cleanest empirical model-error estimate is therefore the genuinely prospective R17 miss:
about `2.15 nm` opening and `12.91 nm` depth. Propagating that error box through `J_local` gives
rough parameter half-widths of approximately `0.0104` in `f` and `0.0070` in `s`. Parameter
uncertainty is not yet closed.

## Candidate comparison

| Candidate | `[f,s]` | Step L2 / previous R16-to-R17 step L2 | Trust-scaled L-infinity | Barycentric status | Local 10 nm prediction, opening/depth (nm) | Decision |
| --- | --- | ---: | ---: | --- | --- | --- |
| Full old-J correction from the mixed 5 nm endpoint | `[0.935559, 0.518358]` | 2.230 | 5.965 | Outside; one negative weight | `50.63 / 748.57` | Reject |
| Existing safeguarded proposal | `[0.900472, 0.558649]` | 0.374 | 1.000 | Inside; weights `[0.1714,0.0037,0.8249]` | `44.15 / 822.93` | Evaluate once at 10 nm |
| Local-interpolation sensitivity candidate | `[0.903774, 0.559746]` | 0.435 | 0.865 | Inside; weights `[0.1693,0.1264,0.7043]` | `45.00 / 825.00` | Record as sensitivity only; do not launch a second run |

For the existing safeguarded proposal, `J_old` predicts `43.774 / 825.189 nm` and `J_local`
predicts `44.151 / 822.930 nm`. The model spread is only `0.38 nm` opening and `2.26 nm` depth.
The safeguarded and local-interpolation parameter pairs differ by only `0.00330` in `f` and
`0.00110` in `s`, corresponding under `J_local` to about `0.85 nm` opening and `2.07 nm` depth.
Creating and evaluating both would not be an informative use of an hour-scale simulation.

The full old-J candidate is not sanctioned. It is almost six trust radii in its limiting
coordinate and the updated same-operator response predicts a large overshoot.

## Exact next-run recommendation and gates

Evaluate exactly one run: the already checksum-bound safeguarded proposal at **10 nm proposal
fidelity only**.

Protocol compatibility must be explicit: the two extra 10 nm endpoints allowed by R1.6 were R16
and R17, so that numerical run allowance is exhausted. R1.8 reclassified 10 nm as proposal fidelity
but did not silently reset the count. Consequently this audit is a recommendation, not run
authorization, until a pre-run protocol amendment names this one fixed model-check point and closes
the 10 nm sequence afterward. Treating it as an informal diagnostic would evade the evidence
contract.

```text
effective_mask_crosslinked_growth_fraction = 0.9004722559883319
oxide_etch_yield_scale                    = 0.5586489665864749
ion_azimuthal_closure                     = axisymmetric_uniform
ion_azimuthal_order                       = 16
topology_change_policy                    = continue_gas_cavity
boundary_case                             = base
duration_s                                = 60
dx_um                                     = 0.01
```

No fine/AMR evaluation and no second candidate should be launched from this audit. The 10 nm run
tests the response model; it does not confer fine-grid authority.

### Pre-run gates

1. A written protocol amendment must authorize exactly this fixed 10 nm model check and forbid
   candidate churn afterward.
2. The continuous minimum-opening observable must be the reviewed/certified definition that will
   score both prediction and result.
3. The proposal checksum, source endpoint checksums, base target checksum, executable source
   checksums, numerical controls, and inactive legacy-quadrature provenance must be in the run
   manifest.
4. Only the base calibration table may be read. Oxygen, power, and transfer held-outs remain sealed.
5. Start at `t=0`; do not resume or splice an earlier parameter trajectory.
6. Use the current hard topology-continuation policy and conservative remap. No terminal-clog
   endpoint may be accepted as a 60 s response.
7. Impose an explicit wall/step budget. A budget checkpoint is resumable evidence, not a result.

### Completion gates

1. The engine must reach exactly `60 s` with status `complete`.
2. Material and surface ledgers, radiosity balance, topology events, unresolved-volume bounds, and
   all existing numerical diagnostics must be reported and pass their unchanged contracts.
3. Score the actual opening/depth with no post-run parameter or observable-definition change.
4. Report predictions from both committed diagnostics:
   - old response: `43.7744 / 825.1889 nm`;
   - local response: `44.1514 / 822.9295 nm`.
5. Report actual-versus-predicted response, scaled merit, and trust ratio. With R17 merit
   `m0=2.616705`, use the committed local prediction `m_pred=0.447528` and

   ```text
   rho = (m0 - m_actual) / (m0 - m_pred)
   ```

### Decision gates

- **Strong response-model pass:** actual opening and depth are each within `5 nm` of target and
  within the empirical R17 model-error envelope around the local prediction (`2.153 nm` opening,
  `12.905 nm` depth). Stop 10 nm calibration; do not evaluate the nearby interpolation candidate.
- **Accept/hold radius:** `rho >= 0.25` and actual merit is lower than R17. Update the local model
  with the actual result, but do not expand the radius or launch a fine run until the late
  coarse/fine discrepancy and authoritative numerical path are resolved.
- **Marginal/shrink:** `0.10 <= rho < 0.25`. Preserve the result and stop 10 nm calibration. Any
  further proposal would require a new explicit protocol decision; this audit does not authorize it.
- **Reject:** `rho < 0.10`, actual merit does not improve, the run is incomplete, or any numerical
  contract fails. Diagnose response transfer/numerics before another parameter evaluation.
- **Fine authority remains separate:** regardless of the 10 nm outcome, a final parameter pair
  belongs only to a clean `t=0` uniform-5-nm or certified-AMR operator. This proposal run cannot
  authorize a held-out reveal.

## Conclusion

The existing safeguarded proposal is appropriately sized. It should be evaluated once at 10 nm,
not at 5 nm, and it should not be replaced by a numerically near-duplicate proposal. The historical
response was directionally useful but quantitatively inaccurate; the same-operator local affine
response is better conditioned and places the target inside the observed response triangle. One
additional cheap, clean base response is needed to test that local model before any expensive
authority run.
