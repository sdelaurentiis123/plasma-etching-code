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
- **The production trajectory step has been corrected and centrally certified.** A later charged
  checkpoint exposed that the former `0.005` fixed trajectory step could cross a material surface
  and report a back-face ion hit. The shared field transport now derives impact direction and cosine
  from terminal velocity plus the declared gas normal and refuses every solid-side hit. A subsequent
  refinement exposed a rarer float32 shared-edge miss: only a hit that fails that certificate is now
  replayed from its original state with the same fixed-step Verlet scheme and edge-inclusive float64
  hard-triangle visibility. Every replay is counted. All earlier `0.005` trench trajectories are
  retained only as controller-mechanics evidence, not physical-time refinement evidence.

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
saturated. Deposition conservation closes to `8.09e-17` and `4.05e-17`; signed charge removed with
the two etched surface increments is itemized separately. This earns a tested co-simulation path,
not pulsed-process validation.

The full local regression suite after the trajectory-lineage integration is **372 passed, 1
skipped**. The skip is the
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

## Bounded real-trench pilot

A two-update coarse real-trench run now exercises C3 from zero authoritative face charge with the
existing separately selected 40-face ion estimator map, exact hard visibility, bounded C2 grazing
reflection, and the 0.25/0.50 micrometer B2 scales. No legacy nodal charge is reused. Config
`f188a7eb1eb6a7476313ffa44af810c47432fe71c3f940a8a8313f3478c7b96e` writes the first replayable
C3 face checkpoint to `results/charging_coevolution_c3_trench_pilot/`.

The first strict run correctly refused the response cascade: deterministic weighted reflection can
leave a positive but geometrically vanishing tail forever, so merely raising the 16-bounce cap cannot
make the population literally empty. A declared conservative tail closure was added with default
tolerance zero. After at least one reflected flight, if the remaining absolute charge rate is below
the declared fraction of the primary absolute rate, it is absorbed on its current impact faces. The
global charge ledger remains exact. The normalized L1 error of the spatial current distribution is
rigorously bounded by twice the closed tail fraction and is reported on every evaluation. A
nondecaying perfect-specular cavity still refuses; the tolerance cannot hide a trapped finite tail.

This is a deterministic bounded-error closure, not stochastic roulette. Classical transport codes
normally use weight cutoffs with Russian roulette so the mean carried weight stays unbiased; see the
[Los Alamos Monte Carlo error analysis](https://www.osti.gov/servlets/purl/6286976-DagIii/). C3 keeps
the deterministic replay and exact per-run charge ledger instead, and requires tail-tolerance
refinement as a separate numerical error dimension.

| Quantity | Initial evaluation | After 1 update | After 2 updates |
| --- | ---: | ---: | ---: |
| Physical time (microseconds) | 0 | 0.125 | 0.250 |
| Node RMS / worst diagnostic | 0.760 / 0.995 | 0.533 / 0.959 | 0.513 / 0.930 |
| B2 max, 0.25 micrometers | 1111.15 | 0.9785 | 0.9628 |
| B2 max, 0.50 micrometers | 662.94 | 0.9784 | 0.9626 |
| Maximum potential-rate magnitude (V/s) | 2.217e8 | 2.158e7 | 2.115e7 |
| Tail spatial-current L1 bound | 2.40e-12 | 4.62e-11 | 6.54e-11 |

The huge initial B2 value is real under the signed contract normalization: some active patches
receive far more electron than ion current, so division by local ion current is severe. It is not the
historical symmetric metric. The pilot remains far from the 0.08 B2 gate and makes no convergence
claim.

Tightening the tail tolerance from `1e-10` to `1e-12` changes final face sigma by relative L2
`1.44e-14` and potential by relative L2 `4.53e-15`; the tighter run's largest reported L1 current
bound is `8.52e-14`. This passes the bounded tail-refinement check at the pilot state. Both runs close
charge deposition below `1.78e-16` relative and surface-transfer charge below `4.68e-15`. The
machine-readable comparison, including hashes of both source summaries and face checkpoints, is
`results/charging_coevolution_c3_trench_pilot/tail_refinement_comparison.json`. Each two-update run
takes about 1.14 seconds on the recorded CPU; both exact wall-clock values are in their manifests.

## Charging-timestep refinement and safeguarded SER/PTC

The coarse trench was advanced from zero face charge to the same 2.5 microsecond physical time with
three fixed steps. Every run used the same source checksums, hard-visibility transport, frozen
estimator-method map, and `1e-12` response-tail tolerance.

| Step | Updates | Node RMS / worst | B2 max, 0.25 / 0.50 micrometers | Potential rate (V/s) |
| ---: | ---: | ---: | ---: | ---: |
| 125 ns | 20 | 0.3963 / 0.8919 | 19.813 / 17.074 | 1.508e6 |
| 62.5 ns | 40 | 0.3956 / 0.8887 | 19.685 / 17.335 | 1.554e6 |
| 31.25 ns | 80 | 0.3952 / 0.8887 | 19.674 / 17.327 | 1.549e6 |

Successive halving changes face sigma by `3.76%` then `1.32%`, nodal charge by `2.12%` then
`0.843%`, and potential by `0.261%` then `0.164%`. The sequence tightens and its face-charge
successive-difference order is about 1.50, but the finest face state is not yet invariant enough to
close the timestep gate. The B2 rebound after the initial two updates survives refinement. It is not
an explicit-step blow-up; its local ion-normalized denominator changes as accessibility evolves.

These runs all used the former `0.005` **particle-trajectory** step. The later trajectory audit below
shows that setting becomes invalid as the field strengthens. The table therefore remains evidence
that the charge-ODE timestep controller and conservation ledgers behave consistently on the early
trajectory, but it no longer supports a physical transient-refinement claim.

The first real SER run exposed three engine issues that are now corrected without changing the
kinetic operator or any acceptance tolerance:

1. The old safeguard rejected a step when the signed B2 ratio grew even if the dimensional current
   residual fell. B2 remains the final contract gate, but only the absolute ODE current residual now
   safeguards the PTC trajectory.
2. The old rejection was detected at the candidate state without restoring the prior charge. SER
   now retains the prior authoritative face/nodal state and clocks, rolls a rejected trial back,
   halves its timestep, and retries. Both face and independent nodal candidates use that same trial
   timestep; their maximum relative mismatch in the real run is `2.46e-16`.
3. Charge-ledger roundoff had been normalized by near-cancelled signed charge. Each step now retains
   positive, negative, absolute-throughput, and signed-net inventories separately; conservation is
   normalized by positive-plus-negative throughput. The prior false failure was a `9.12e-34` C
   residual divided by a `9.42e-22` C signed net. The corrected 80-step SER run's worst relative
   deposition-ledger error is `9.53e-17`.

At approximately equal elapsed time, SER reaches 4.3549 microseconds with 80 accepted steps and two
rolled-back trials; the 31.25 ns reference reaches 4.3750 microseconds with 140 steps. SER differs
from fixed time by relative L2 `0.339%` in face sigma and `0.176%` in potential, uses 42.9% fewer
accepted steps, and is 1.65x faster on the recorded CPU. This passes the bounded same-operator PTC
mechanics check. It does not pass the stationary-state or schedule-refinement gates: the SER endpoint
still has node RMS/worst `0.350/0.790`, B2 maxima `11.50/10.59`, and potential rate `1.46e6` V/s.

The complete paired evidence is
`results/charging_coevolution_c3_trench_refinement/comparison.json`; it hashes every source, summary,
and face checkpoint. The runner also writes a replayable face checkpoint, full diagnostics, source
hashes, and wall clock before exiting nonzero on any `SurfaceChargingSaturationError`, so rejected
campaigns no longer disappear into tracebacks.

## Particle-trajectory resolution and back-face refusal

At the first later checkpoint where fixed-time continuation failed, the shared field kernel reported
an Ar+ event on top mask face 3 with stored incidence cosine `+0.998886`, but its terminal velocity
pointed along the declared `+z` gas normal. The authoritative geometric cosine was therefore
`-0.998886`: the fixed-step trajectory had crossed the solid and intersected the surface from behind.
The old response-level consistency check refused this event before it could affect reflection or
charge. That check has now been moved into the shared field transport, so primary boundary particles,
neutral reuse, bidirectional forward estimates, and surface-emitted re-impact cascades all use the
same terminal-velocity/gas-normal certification. The response check remains as defense in depth.

The exact face checkpoint was replayed without a charging update over a fixed phase-space sample set.
The `0.005` level refuses. All finer levels shown below finish without truncation and close deposition
and surface-transfer charge at roundoff.

| Particle step | Status | Node RMS / worst | B2 max, 0.25 / 0.50 micrometers | Potential rate (V/s) |
| ---: | --- | ---: | ---: | ---: |
| 0.005 | refused: solid-side hit | -- | -- | -- |
| 0.0025 | accepted | 0.34376 / 0.78032 | 10.6601 / 10.3011 | 9.635e5 |
| 0.00125 | accepted | 0.34463 / 0.78228 | 10.7789 / 10.4142 | 7.289e5 |
| 0.000625 | accepted | 0.34515 / 0.78328 | 10.8407 / 10.4722 | 8.238e5 |
| 0.0003125 | accepted | 0.34574 / 0.78376 | 10.8721 / 10.5021 | 8.791e5 |
| 0.00015625 | accepted | 0.34577 / 0.78402 | 10.8878 / 10.5169 | 9.030e5 |
| 0.000078125 | accepted | 0.34584 / 0.78414 | 10.8954 / 10.5244 | 9.162e5 |

The final halving changes node RMS/worst by `0.0216% / 0.0159%`, the two exact B2 maxima by
`0.0697% / 0.0708%`, and maximum potential rate by `1.44%`. This first audit established local
timestep convergence at one fixed checkpoint, but it did not prove that a raw float32 mesh query
would remain certified everywhere along a changing charge trajectory. The initially selected
`0.0003125` level is therefore a bounded campaign setting, not a globally certified visibility
scale or a waiver of final B1 refinement.

Machine-readable evidence, source hashes, all per-node and patch diagnostics including the refused
run, and the durable input checkpoint are in
`results/charging_coevolution_c3_trajectory_refinement/audit.json`. This audit changes the next step:
long continuation must restart with the resolved transport, not continue a coarse-trajectory state.

## Flight-horizon separation and exact hard-visibility replay

Restarting from zero charge separated two numerical effects that the preceding fixed-checkpoint
audit could not distinguish:

1. At `dt=0.0003125`, a 50,000-step horizon (dimensionless flight time `15.625`) refuses because
   slow adjoint electrons remain unresolved. A 128,000-step horizon (`40`) passes, and doubling it
   to `80` produces identical reported currents. The bounded runner default is therefore 128,000
   steps; horizon and timestep are recorded independently.
2. At 2.375 microseconds, the raw float32 Warp mesh query misses a wall intersection near a shared
   triangle edge and later reports a solid-side hit on top-mask face 38. Halving the particle step
   once removes the failure, but halving it again produces a different solid-side hit. That
   non-monotone `refuse / pass / refuse` pattern is a floating intersection degeneracy, not evidence
   of a physical smoothing length or ordinary ODE truncation error.

The engine response is deliberately narrow. The Warp float32 fixed-step Verlet path remains the fast
operator. Its terminal velocity and declared level-set gas normal form an independent lineage
certificate. Only a hit that fails that certificate is replayed from the original incident or
surface-emitted state using the same fixed timestep, the same trilinear Q1 electric field, and an
edge-inclusive float64 hard-triangle intersection. The replay may land on another exact face or
escape; it may not be softened, dropped, or accepted incomplete. Primary and reflected/re-emitted
charged paths share this repair, and `transport_lineage_replay_count` is retained in every charging
evaluation.

The exact failure checkpoint was then evaluated at three particle steps:

| Particle step | Replays | Node RMS / worst | B2 max, 0.25 / 0.50 micrometers | Potential rate (V/s) |
| ---: | ---: | ---: | ---: | ---: |
| 0.0003125 | 1 | 0.401044 / 0.892723 | 20.9533 / 17.2026 | 1.427e6 |
| 0.00015625 | 0 | 0.401382 / 0.892854 | 20.9813 / 17.2261 | 1.401e6 |
| 0.000078125 | 1 | 0.401436 / 0.892920 | 20.9943 / 17.2376 | 1.387e6 |

Across both halvings, node RMS changes `0.0975%`, worst node `0.0220%`, the two B2 scales
`0.196% / 0.203%`, and potential rate `2.89%`. Replay incidence is non-monotone (`1 / 0 / 1`), while
the physical diagnostics refine smoothly. A manufactured ray aimed exactly at the diagonal shared
edge of two coplanar triangles also lands with unit incidence and closes its lineage.

Finally, the repaired zero-charge trajectory completes all 20 requested 125 ns charging steps to
2.5 microseconds. Exactly one lineage is replayed, at evaluation 20 (physical time 2.375
microseconds), where the unrepaired operator refused after 18 accepted steps. The repaired endpoint
has node RMS/worst `0.396962 / 0.886989`, B2 `19.8109 / 17.4492`, potential rate `1.498e6` V/s,
deposition conservation `1.51e-16`, and surface-transfer balance `4.51e-15`. These values remain far
from B1/B2; this is transport certification and bounded progress, not charging convergence.

Machine-readable hashes, the refusal/repair pair, horizon pair, failure-state refinement, replay
counts, and decisions are in
`results/charging_coevolution_c3_lineage_replay/audit.json`.

Historical residual baselines (`0.788`, `0.627`, and the later approximately `0.30` values) are not
used as quantitative comparators for this repaired C3 operator: face-authoritative charge, response
physics, particle timestep, flight horizon, and hit certification differ. The shared-edge defect makes
uncertified float32 histories potentially vulnerable to rare hit leakage, but this audit does not
claim that every historical sample contained one. Likewise, the earlier stuck-map horizon audit found
zero unresolved trajectories and no current change under 4x/8x horizon extension; the new horizon
failure arose only when timestep was reduced without preserving total flight time. The post-repair
baseline begins with the artifacts in this section.

## Repaired-operator decay, restart, and PTC schedule audit

The exact repaired operator was advanced from zero charge for 60 fixed 125 ns updates. A second path
ran 20 updates, wrote the authoritative face checkpoint, reloaded it, and ran 40 more. At the common
7.5 microsecond endpoint, every checkpoint array is bitwise equal and both NPZ files have SHA-256
`91122190f3e306f119ffff62ba7b2140559a7ebf45555883a2219de0ddb2954e`. This certifies restart
invariance. It is not yet the signed cold-versus-remapped-warm stationary-branch gate, because neither
path is stationary and both originate from the same zero-charge trajectory.

| Physical time (microseconds) | Node RMS / worst | B2 max, 0.25 / 0.50 micrometers | Maximum `abs(dV/dt)` (V/s) |
| ---: | ---: | ---: | ---: |
| 0 | 0.7600 / 0.9951 | 1111.05 / 662.87 | 2.217e8 |
| 2.5 | 0.3970 / 0.8870 | 19.8109 / 17.4492 | 1.498e6 |
| 5.0 | 0.3398 / 0.7712 | 10.1449 / 9.9566 | 1.410e6 |
| 7.5 | 0.2999 / 0.7229 | 8.4359 / 7.5235 | 1.027e6 |

B2 is a dimensionless ion-normalized ratio, not a percentage-valued column: `8.4359` means about
`843.6%`, while the contract gate is `0.08` or `8%`. The endpoint is therefore about 105 times the
strict B2 gate, and B1 is about 1027 times its `1000 V/s` pilot tolerance.

Nonnegative fits of `c + A exp(-(t-2.5 microseconds)/tau)` over 2.5--7.5 microseconds describe node
RMS well (`R^2=0.997`, `tau=5.71 microseconds`) but place its fitted floor at `0.240`. The two B2 fits
place their floors at `9.66` and `8.58`. Most importantly, maximum potential rate is not described by
a single exponential (`R^2=0.300`). These short-window fits are diagnostics, not asymptotic proofs;
they do not support a saturation-time projection or a fit-and-jump proposal.

Two safeguarded SER schedules started from the identical repaired 2.5 microsecond checkpoint:

| Schedule | Allowed residual growth | Accepted / rejected | Pseudo-time advanced | Outcome |
| --- | ---: | ---: | ---: | --- |
| A | 0.5% | 30 / 4 | 2.188 microseconds | Repeated halvings collapse the step to 8.77 ns; slower than fixed time |
| B | 2% | 30 / 0 | 4.206 microseconds | Modest same-operator acceleration; about 10.8% fewer accepted steps |

At schedule B's 6.706 microsecond total equivalent time, the nearest fixed checkpoint is 6.75
microseconds. RMS differs `0.059%`, worst node `0.379%`, B2 `0.054% / 0.128%`, and potential rate
`3.27%`. Schedule B earns continued bounded use, but not final promotion: potential rate remains the
least-refined diagnostic and a matched endpoint/state comparison is still required. Schedule A is
rejected as an accelerator, not mistaken for a physical failure.

The requested replay canary now records its denominator throughout primary and charged re-impact
transport. At the exact failure state, one of 8653 eligible field lineages replays (`0.0116%`); at
7.5 microseconds, zero of 9079 replay. The fixed trajectory has isolated one-lineage events at 2.375,
3.0, and 7.125 microseconds, not monotone growth with charge.

The proposed potential-space rewrite is not promoted. The physical engine already computes
`dV/dt` by applying the exact Poisson charge-to-voltage response to conservative `dQ/dt`. Making
voltage authoritative would require inverting the rank-deficient face-charge projection and would
lose the unique signed surface ledger. PTC theory permits a nonsingular scaling/mass matrix, but a
future matrix pseudo-time preconditioner must therefore be defined in authoritative face-charge
space, remain easily invertible, and reproduce physical-time stationary states under schedule
refinement; none is assumed here. See [Coffey, Kelley, and Keyes](https://doi.org/10.1137/S106482750241044X).

Machine-readable fits, hashes, schedule comparisons, unit interpretation, replay fractions, and
accelerator decisions are in `results/charging_coevolution_c3_decay_audit/audit.json`.

## Independent-scramble and sample-level audit at 7.5 microseconds

The repaired 7.5 microsecond face checkpoint was evaluated without another charge update using the
same hard-visibility operator and frozen estimator-method map. The audit uses nested same-scramble
sample doubling at forward/adjoint levels `10/8`, `11/9`, and `12/10`; eight independent scrambles
at `11/9`; and a second nested `11/9 -> 12/10` pair. A runner defect discovered during setup is now
fixed: `--seed` previously changed forward phase-space samples but left the adjoint proposals frozen
at 79/83. It now controls the Ar+ adjoint proposal as well, with the electron proposal declared as
`seed+4`. The default seed-79 evaluation remains identical, while independent runs now refresh every
sample family.

| Forward / adjoint level, seed | Node RMS / worst | B2 max, 0.25 / 0.50 micrometers | Maximum `abs(dV/dt)` (V/s) |
| --- | ---: | ---: | ---: |
| 10 / 8, 79 | 0.299887 / 0.722924 | 8.43589 / 7.52354 | 1.027e6 |
| 11 / 9, 79 | 0.297124 / 0.719829 | 7.93320 / 7.17026 | 1.110e6 |
| 12 / 10, 79 | 0.291805 / 0.698030 | 7.26521 / 6.77372 | 8.741e5 |
| 11 / 9, 179 | 0.291437 / 0.697093 | 7.12022 / 6.76942 | 1.164e6 |
| 12 / 10, 179 | 0.292649 / 0.687451 | 6.59268 / 6.41140 | 1.120e6 |

The eight-scramble `11/9` ensemble gives node RMS `0.29145 +/- 0.00247`, worst node
`0.69934 +/- 0.01126`, B2 `7.228 +/- 0.511` and `6.481 +/- 0.527`, and maximum potential rate
`(1.126 +/- 0.133)e6 V/s` (95% Student-t intervals on the mean). Every confidence interval is far
outside B1/B2, so the non-convergence conclusion is robust. The numerical value of the apparent
late-time state is not sample-refined: the two same-seed level doublings reduce B2 by 5.3--8.4%, and
the potential-rate change ranges from 3.8% to 21.3%. A longer transient may still be informative,
but these data do not earn a precise saturation projection, plateau claim, or fit-and-jump.

The runner now stores `current_audit.npz` with positive, negative, and net face currents, nodal
currents, physical face area, and patch membership at both scales. Localization shows that the B2
maximum is the lower mask sidewall (material 2, approximately z=1.33--1.42 micrometers). Its
0.5-micrometer patch carries about 3% of total current throughput; it is not an inactive numerical
speck. Electron collection is approximately 8--9 times ion collection there, and the relatively
small ion denominator controls much of the remaining B2 sampling variation. This is an operator
diagnostic, not grounds to exclude, merge, or loosen that patch.

Across the eight `11/9` audits, exact replay is used once among 145,269 eligible lineages
(`0.000688%`). Maximum deposition and surface-transfer ledger errors are `2.89e-17` and `5.03e-15`.
Machine-readable protocol, intervals, paired changes, patch localization, and decisions are in
`results/charging_coevolution_c3_sample_audit/audit.json`.

## Fixed versus fresh-scramble physical time

The fixed level-10/8 reference was extended from 7.5 to 15 microseconds without a replay or ledger
failure. Node RMS/worst fell to `0.2215 / 0.5697`, and B2 fell to `4.135 / 3.710`; the earlier fitted
late-time floor was therefore a short-window artifact. Maximum potential rate did not fall
monotonically and ended at `1.53e6 V/s`. A further 1.25-microsecond fixed continuation at 125 ns and
62.5 ns agrees within `0.031%` in face sigma and `0.057%` in potential, so the sign-changing rate is
not a global charging-timestep instability.

Independent level-11/9 current audits reveal the more important issue. At 15 microseconds the
ensemble-mean floor-center voltage rate is `-3.30e6 V/s`; at 16.25 microseconds it is still
`-2.74e6 V/s`. The low-sample frozen trajectory instead drives floor potential from `9.65` to
`11.28 V`. A finite frozen sample operator is therefore steering this interval opposite to the
ensemble-expected current, even though every individual update is exactly conservative. Conservation
prevents particle/charge loss; it does not make a finite-sample drift unbiased.

The unified C3 integrator now supports an explicit `frozen` or `fresh` scramble mode. Fresh mode
regenerates forward, Ar+ adjoint, and electron adjoint samples at every accepted physical update;
records the seed epoch on every evaluation; and gives the final state a new, unused audit epoch.
Fresh mode is restricted to fixed physical time. It refuses SER because residual changes between
independent scrambles cannot safely control deterministic trial rollback. Frozen behavior and its
default seed sequence are unchanged.

Eight fresh-scramble realizations were advanced from the identical 15-microsecond checkpoint to
16.25 microseconds:

| Charge timestep | Mean node RMS | Mean B2, 0.25 / 0.50 micrometers | Mean floor potential (V) | Mean maximum rate (V/s) |
| ---: | ---: | ---: | ---: | ---: |
| 125 ns | 0.1892 +/- 0.0076 | 3.320 +/- 0.267 / 2.930 +/- 0.245 | 9.208 +/- 0.331 | `(1.686 +/- 0.767)e6` |
| 62.5 ns | 0.1912 +/- 0.0045 | 3.404 +/- 0.450 / 2.927 +/- 0.278 | 9.386 +/- 0.309 | `(2.022 +/- 1.058)e6` |

Intervals are two-sided 95% Student-t intervals on eight realization means. The timestep-halved
ensemble means differ by `0.218%` in face sigma and `1.86%` in potential L2; the largest absolute
potential difference is `0.185 V`. This bounded interval passes a statistical timestep check for the
state trajectory. It remains far outside B1/B2, and maximum rate remains the least precise measure.

Two proposed diagnostics were narrowed after code and literature review. First, `1/sqrt(N)`
Richardson extrapolation is not a valid scrambled-Sobol law: scrambled-net variance is integrand- and
resolution-dependent and can decrease faster than Monte Carlo, while B2 is a nonlinear maximum of
ratios. Independent scrambles and nested powers of two remain the evidence contract; see
[Owen 1997](https://epubs.siam.org/doi/10.1137/S0036142994277468). Second, the controlling mask ion
faces use the separately frozen forward estimator, so selective *adjoint* retracing would not refine
their B2 denominator without a new independent estimator-certification pilot.

A diagonal pseudo-time scaling remains mathematically plausible but held. Kelley and Keyes permit a
nonsingular, often diagonal scaling matrix to equilibrate local CFL numbers, but their convergence
analysis assumes stable sufficiently regular dynamics; it does not supply a relaxation matrix for
this noisy nonlocal Poisson/current operator or justify a claimed 10--50x speedup. Any future local
face scaling must be built from ensemble-mean response, safeguarded, and reproduce the physical-time
stationary state under refinement; see [Kelley and Keyes](https://repository.lib.ncsu.edu/bitstreams/339a4ac5-c0a2-48ad-b4a3-c4ec0085f0aa/download)
and [Coffey, Kelley, and Keyes](https://doi.org/10.1137/S106482750241044X).

Machine-readable trajectories, ensemble intervals, timestep differences, and method decisions are in
`results/charging_coevolution_c3_stochastic_transient_audit/audit.json`.

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
2. Restart fixed physical time and safeguarded SER with the resolved particle trajectory; compare the
   same exact operator under charge-timestep halving and SER schedule refinement.
3. Report B2 on at least two physical patch scales and demonstrate tightening under one grid
   refinement, with a scale no larger than each claimed feature.
4. Demonstrate sample-level, timestep, and grid invariance of the declared observables. For any
   experimental claim, cap its tolerance at the benchmark uncertainty including digitization.
5. Run the independent high-sample B5 audit with exact hard visibility and report retained per-node
   RMS/worst diagnostics even if the amended patch/observable contract passes.

If explicit timestep refinement and a stable PTC/implicit treatment do not approach the same
stationary state, the signed kill condition applies: stop and open the discrete-equilibrium audit.
