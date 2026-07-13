# Charging co-evolution C3 integration audit

Date: 2026-07-13. Contract revision: `CCA-2026-07-13-R2` (signed and in force).

## Outcome

- **The unified C3 engine path is implemented and passes its manufactured integration gates.** One
  code path now couples authoritative face sheet charge, compatible Q1 Poisson, exact hard-visibility
  charged transport, the conservative charged-response/re-impact cascade, surface chemistry, level-set
  motion, and C1 signed charge remap. Chemistry consumes the charging solve's final transport object;
  it does not retrace a second kinetic operator.
- **Quasi-static and waveform-resolved operation share the same operators.** Quasi-static mode must
  pass B1 and B2 before profile motion. Waveform mode advances the same charge ODE once per explicitly
  declared physical segment and does not pretend the charge saturated. Supplying a waveform in
  quasi-static mode is refused.
- **This is an engine-integration pass, not C3 scientific closure.** The designated real-trench
  timestep/grid/sample refinement, cold-versus-warm branch test, observable invariance, and independent
  high-sample exact-operator B5 audit remain pending. C4 is therefore not authorized by this report.

Replay artifact: `results/charging_coevolution_c3/audit_summary.json`.

## Unified transaction

For each quasi-static profile increment:

1. Extract the current material surface and rebuild its Q1 Poisson system.
2. Project authoritative face sheet charge to Q1 nodes.
3. Integrate `d sigma / dt = J_i - J_e` by fixed physical time, then optionally safeguarded SER once
   the residual enters its activation region.
4. Require potential-rate saturation and patch-current balance at no fewer than two physical scales.
5. Reuse that final exact charged/reflected transport for surface chemistry.
6. Advance the common level-set/material engine.
7. Remap positive and negative retained surface-charge inventories separately and ledger charge removed
   with the etched material.
8. Warm-start the rebuilt field on the next geometry.

The SER rule follows residual-ratio pseudo-transient continuation with bounded growth and a
residual-growth rejection safeguard. It changes only the explicit step size of the same conservative
charge ODE. No frozen-map root solver, smoothed visibility, volume Boltzmann term, or alternate final
operator is present. See [Kelley and Keyes](https://doi.org/10.1137/S0036142996304796) for the
pseudo-transient/SER basis.

Waveform-resolved mode replaces step 3's saturation loop with one physical update for each declared
`ResolvedBiasSegment3D`. The corresponding endpoint transport drives that segment's profile motion,
then the same remap and field rebuild occur. Segment durations must sum to the declared run duration;
species identity, charge, mass, ordering, and reference plane must remain consistent.

## Contract enforcement

| Signed clause | Engine enforcement | Manufactured result |
| --- | --- | ---: |
| B1 potential saturation | Maximum `abs(dV/dt)` reported; required in quasi-static mode | 0 V/s |
| B2 physical patch balance | `abs(Ji-Je)/Ji`; at least two distinct fixed physical scales required | 0 at 0.25 and 1.0 micrometers |
| B2 feature-scale rider | A feature claim is refused unless one patch scale is no larger than its declared extent | Refusal gate passes |
| B3 tolerance rider | Claim tolerance must not exceed combined experimental and digitization uncertainty | Refusal gate passes |
| B4 deposition conservation | Every evaluated step is checked at roundoff; failed evaluation history is retained | 0 relative error (quasi-static) |
| Signed C1 inventories | Retained/removed positive and negative charge reported independently | 0 ledger error |
| Exact operator | Hard visibility and caller-supplied response used for final transport | Exact object reused |
| Pulsed-bias no-go | Quasi-static mode refuses any supplied waveform | Refusal gate passes |
| Parameter provenance | Response parameters require identical parameter/bounds/source sets | Reflection manifest complete |
| Per-node diagnostic retention | RMS and worst-node current imbalance recorded in every charging history item | 0 / 0 in planar gate |

The B3 API uses `ExperimentalObservableTolerance3D`. Profile-feature claims additionally declare
`feature_extent_m`; the driver checks it against the run's physical patch scales. Runs with no
experimental claim do not invent an observable tolerance.

The historical symmetric `(Ji-Je)/(Ji+Je)` node and patch diagnostics are also retained, but they do
not gate B2. The audit caught and removed an initially too-permissive use of that symmetric statistic
before any real-trench C3 claim was run.

## Manufactured results

The deterministic planar gate uses equal 100 eV Ar+ and electron fluxes on a dielectric plane,
two physical patch scales, the bounded C2 reflection response, and physical sputter chemistry.

| Quantity | Quasi-static result |
| --- | ---: |
| B1 potential rate | 0 V/s |
| Retained node RMS / worst imbalance | 0 / 0 |
| Patch maximum imbalance, 0.25 / 1.0 micrometers | 0 / 0 |
| Deposition conservation relative error | 0 |
| Surface-response conservation relative error | 0 |
| C1 remap conservation relative error | 0 |
| Final charged/reflected transport reused by chemistry | yes |
| Wall clock on recorded CPU | about 0.03 s |

The waveform smoke gate uses a 1 ns ion-rich segment followed by a 1 ns electron-rich segment.
Both segments complete exactly one physical charge update; neither is required or reported to be
saturated. Deposition conservation closes to `2.43e-16` and `1.21e-16`; signed charge removed with
the two etched surface increments is itemized separately. This earns a tested co-simulation path,
not pulsed-process validation.

The full local regression suite after integration is **364 passed, 1 skipped**. The skip is the
existing unavailable-CUDA condition on this CPU-only build.

## Legacy checkpoint migration refusal

The refined pre-C3 checkpoint stores nodal charge but not authoritative face sheet charge. Reusing
it would require an inverse of the compatible Q1 face-to-node projection. On its identical archived
176-face mesh, that projection has shape 765 by 176, numerical rank 121, and condition number
`1.41e18`. The least-squares face field misses the archived nodal state by relative L2 `0.584` and
relative Linf `0.681`.

Migration is therefore refused. This is not a new physical failure: the legacy refined trajectory
was initialized and evolved as a nodal state and contains components outside the face-sheet image.
There is no unique conservative surface charge that can be remapped from it. The real C3 campaign
must start from zero face charge or a checkpoint written by this C3 face-authoritative path. No
minimum-norm, regularized, or guessed inverse is admitted.

## Evidence and provenance

Audit config hash: `d3b5485aff03a950c82f6fb4a0161e76532b5120dc7f5a075bb340a7a4c444fc`.

The JSON artifact records the base git revision, SHA-256 for the C3 driver, shared feature step, C1
remap, charged response, and audit runner, plus hardware, Python version, seeds, trajectory horizon,
patch scales, hard-visibility statement, full response parameter values/bounds/sources, both signed
charge ledgers, wall clock, and all retained node/patch diagnostics.

## Remaining C3 gates

Before scientific closure or C4:

1. Run a designated real-trench mid-etch and final-profile checkpoint from both zero charge and the
   remapped warm start generated by C3 itself; require the same stationary branch.
2. Compare fixed physical time and safeguarded SER on the same exact operator; perform timestep
   halving and schedule refinement.
3. Report B2 on at least two physical patch scales and demonstrate tightening under one grid
   refinement, with a scale no larger than each claimed feature.
4. Demonstrate sample-level, timestep, and grid invariance of the declared observables. For any
   experimental claim, cap its tolerance at the benchmark uncertainty including digitization.
5. Run the independent high-sample B5 audit with exact hard visibility and report retained per-node
   RMS/worst diagnostics even if the amended patch/observable contract passes.

If explicit timestep refinement and a stable PTC/implicit treatment do not approach the same
stationary state, the signed kill condition applies: stop and open the discrete-equilibrium audit.
