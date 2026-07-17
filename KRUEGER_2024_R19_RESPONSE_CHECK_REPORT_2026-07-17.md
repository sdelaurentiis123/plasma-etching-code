# Krüger 2024 R1.9 response-check report

Date: 2026-07-17
Scope: one precommitted 10 nm base-only development run; no held-out profile was opened; no
validation or fine-grid authority claim is made.

## Executive result

The repaired engine completed the full 60 s process normally. The **physical run passed** its
execution contracts; the **old endpoint response model failed** its precommitted prediction gate.

| Quantity | Experiment | R19 result | Error | Status |
| --- | ---: | ---: | ---: | --- |
| Minimum mask opening | 45.000 nm | 45.085 nm | +0.085 nm (+0.19%) | inside one 5 nm fine cell |
| Etch depth | 825.000 nm | 853.219 nm | +28.219 nm (+3.42%) | outside one 5 nm fine cell |

The very close opening match is encouraging but not authority. At 10 nm the throat is represented by
one interior extrusion row and the late trajectory contains grid-scale throat-selection excursions.
The final number must therefore retain the opening-refinement qualification from
`KRUEGER_2024_MASK_OPENING_AUDIT_2026-07-17.md`.

The run took `467.55 s` wall time on an RTX 4090, accepted `436` adaptive steps, reached exactly
`60.0 s`, recorded no topology event or terminal refusal, and closed the material ledger exactly.
Maximum reported neutral-radiosity balance error was `2.75e-12`.

## What was committed before execution

R1.9 authorized one and only one additional 10 nm response check at

```text
effective_mask_crosslinked_growth_fraction = 0.9004722559883319
oxide_etch_yield_scale                      = 0.5586489665864749
```

The same-operator historical affine response predicted `44.151 / 822.930 nm`; the older response
predicted `43.774 / 825.189 nm`. The empirical model-error envelope was fixed at `2.153 nm` opening
and `12.905 nm` depth. The exact launch manifest is
`results/krueger_2024_r19_response_check/remote_artifacts/launch_manifest.json`, SHA-256
`821301496f133a16498e5051fe491827e5bdf2fd430f562cbf8e0ecf2d9e7aa5`.

The actual response missed the local prediction by `+0.934 nm` opening and `+30.289 nm` depth.
Scaled merit worsened from the historical R17 value `2.6167` to `5.6438`; the committed trust ratio
is `rho = -1.3955`. The mechanical verdict is therefore `reject_response_model`. R1.9 closes the
10 nm candidate sequence regardless of this outcome; there is no second informal parameter run.

The evaluation receipt is
`results/krueger_2024_r19_response_check/evaluation.json`, file SHA-256
`15047f554ca9d96a6efab4973e1433f552913a3656850ba08c26244b2c87dd77`.

## Why the response model failed

Two effects are now separated.

### 1. The response is path-nonlinear

The lower oxide-yield scale initially behaves as expected. Relative to R17, the R19 candidate is
`0.098 nm` shallower at `0.5 s` and reaches a maximum deficit of `4.952 nm` at `11.77 s`. The depth
trajectories cross at `14.89 s`. The R19 path is then deeper and ends `15.314 nm` deeper than R17,
despite the lower nominal oxide yield.

During the same interval the mask/access geometry separates. At `30 s` the R19 throat is
`4.445 nm` wider than R17. Transport, surface coverage, and geometry therefore feed back on each
other strongly enough to overwhelm the local endpoint derivative. Near a necking/access threshold,
two scalar endpoint sensitivities are not globally linear even when the parameter point lies inside
the historical triangle.

### 2. The historical response belonged to a pre-repair operator epoch

The stopped R17 worker preserved its executable sources, so this is no longer guesswork.

| Runtime component | R17 versus R19 |
| --- | --- |
| Amorphous-carbon mask chemistry | byte-identical |
| Material mechanism | byte-identical |
| Neutral radiosity | byte-identical |
| Boundary transport | byte-identical |
| Published boundary tables and compressed IEAD | checksum-identical |
| Feature evolution/remap | **changed:** R17 `a1e8d02f...`, R19 `61dbef70...` |

The feature-step repair makes material-local state remapping periodic across the 20 nm extrusion
seam, validates distance against the represented triangles rather than arbitrary triangle
centroids, and adds bounded subcell-ownership handling. Periodic neighbor selection participates in
every accepted state remap, even when no topology event occurs. An endpoint Jacobian learned under
the old nonperiodic remap cannot be assumed to describe the repaired operator.

At the same R17 parameter pair, old and current operators agree at `0.5 s` to `0.0045 nm` in depth.
That excludes a gross boundary or initial-rate regression. It does not make the old late response
portable: cumulative remap/state feedback becomes visible only after the profile evolves.

### 3. Frozen-checkpoint decomposition localizes access feedback, but `dt=0` is oxide-blind

A bounded current-operator screen evaluated the R17 and R19 completed checkpoint geometries/states
with both R17 and R19 mechanism parameter pairs. All four cells used the same base-only reduced
quadrature, `duration_s=0`, and no profile evolution. The matrix completed in `49.37 s` on local CPU;
every cell had an exact material ledger and neutral-radiosity error below `2.62e-12`.

The accumulated checkpoint effect dominates the instantaneous mask-parameter effect: at fixed
parameters the magnitude of the net mask rate falls by about `55--59%` from the R17 checkpoint to
the R19 checkpoint, whereas changing the parameter pair moves the mask rate by `7.45%` at the R17
state and `15.70%` at the R19 state. The R19 geometry/state also directs about `5.98%` more ion rate
to SiO2, `2.87%` less ion rate to the mask, and modestly less O rate. This supports the path-feedback
diagnosis without another moving-profile run.

The screen cannot decide the depth response. Both the one-point screen and a selected three-point
quadrature returned zero instantaneous SiO2 velocity even with nonzero SiO2 ion arrival. This is not
a missing trajectory: every saved SiO2 face carries finite polymer inventory, and the exact surface
kernel defines `dt=0` removal velocity as zero. In a real positive step, the conservative chemistry
first evolves/depletes that film and only then exposes oxide removal. The earned diagnostic is
therefore a frozen-geometry, positive-horizon chemistry micro-step with fixed transport—not a longer
profile or another parameter fit.

The low-fidelity screen receipt is
`results/krueger_2024_r19_response_check/frozen_checkpoint_2x2/audit.json`, SHA-256
`9eba876f153c9acc1881edc6d4894f6626682f8b251bac9b8bca975ced4e4237`. The selected q3 receipt is
`results/krueger_2024_r19_response_check/frozen_checkpoint_q3_r19_benchmark/audit.json`, SHA-256
`c352ceb5354c6e55fc2d82ab09e0603fb3cd2de5e50b1b4ef1a958ec0d15b8f6`.

The R17 source subset and complete source manifest are preserved under
`results/krueger_2024_r19_response_check/r17_source_snapshot/`. The complete current source tree
executed on the GPU matches the local 76-file Python tree after excluding inert macOS AppleDouble
metadata files.

### 4. Coupled surface-state/radiosity microstep restores the expected oxide-yield direction

The first positive-time chemistry audit showed why `dt=0` was blind: at only
`0.00779536 s`, polymer state changed enough to move a local neutral reaction probability by about
`5.96%`. The production feature step had been holding its diffuse-neutral radiosity solution fixed
over a much longer outer step. Because polymer state controls sticking and reflected neutrals alter
the next surface flux, that split was not numerically resolved.

A follow-up diagnostic kept the exact R1.9 geometry fixed, captured one production periodic direct
transport and one production form-factor operator, and then re-solved radiosity while subcycling the
exact material mechanisms. Each candidate chemistry step was graded by embedded step-doubling: one
full step versus two half steps with a midpoint radiosity solve; only the conservative fine path was
accepted. The comparison covered every state increment, all exchange inventories, and per-face
integrated recession/growth.

The first two horizons passed every declared gate:

| Horizon | R17 oxide removal | R19 oxide removal | R19/R17 | Max gross motion | Tight/nominal difference |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `7.795 ms` | `2794.390` units | `2753.804` units | `0.985476` | `0.02104 dx` | `<0.0065%` |
| `15.591 ms` | `5405.351` units | `5325.748` units | `0.985273` | `0.03991 dx` | `<0.0094%` |

Thus the R1.9 parameter pair produces about `1.45--1.47%` less oxide removal than R17 on the same
late geometry/state, the physically expected direction from its lower oxide-yield scale. Material
ledgers are exact, accepted embedded errors stay below `0.76%` (`0.20%` on the tighter path), and
radiosity balance stays near `2e-12`. The dominant coupling error is the evolving SiO2 polymer
inventory, not oxide removal itself.

The predeclared `480 s` CPU budget expired while evaluating the `31.181 ms` horizon. That timeout is
implementation/controller evidence only; it is not a physics failure. The largest fully completed
frozen-geometry horizon remains `15.591 ms`. Audit SHA-256:
`0f337c8516542f4603926cc8b1d28d156abba19859d68ac2d6b6b6231a97706e`.

One additional API hazard was found in the diagnostic: periodic ballistic first-hit wrapping had
been inferred from the neutral-radiosity option. Disabling radiosity to obtain a
"direct" result also disabled periodic transport and lost `25--76%` of integrated incident rate.
The accepted diagnostic instead captures the direct transport *inside* the exact production periodic
call, then proves the cached replay identical face-by-face and hash-by-hash. The old seed-unbound q3
receipt was regenerated with explicit transport seed `241` and radiosity seed `10241`; the corrected
receipt is under `frozen_checkpoint_q3_seed241_current/`. The common engine now exposes an explicit
`ballistic_periodic_lateral` control, retaining the old radiosity-inherited default for replay while
allowing periodic first-hit transport without enabling diffuse re-emission.

This closes the causal direction, not the calibration. A reusable cached surface/radiosity integrator,
manufactured DAE tests, nested form-factor refinement, and a short current-operator profile comparison
must pass before another endpoint is purchased. Held-out oxygen/power observations remain sealed.

## Visual evidence

![Bounded checkpoint and coupled-chemistry evidence](results/krueger_2024_r19_response_check/krueger_bounded_endpoint_diagnostics.png)

The first panel separates the accumulated geometry/state effect from the parameter-only effect. The
second shows that the R19 direction is stable across both completed positive horizons. The third
puts every numerical/scope diagnostic on its own declared-limit scale: both coupled horizons remain
left of the pass line, the old one-shot frozen-flux approximation fails, and the unresolved longer
horizon is labeled as controller/runtime evidence only. The plot is regenerated from the three
checksummed base-only audit JSON files by `scripts/plot_krueger_2024_bounded_diagnostics.py`; its
provenance sidecar sits beside the PNG.

![R17/R19 response comparison](results/krueger_2024_r19_response_check/r17_r19_response_comparison.png)

The upper panels show stable full trajectories. The lower-left panel shows that both precommitted
models expected the target neighborhood while the actual repaired-operator depth landed higher.
The lower-right panel shows the sign reversal in depth response at `14.89 s`, the signature of a
trajectory-level nonlinear feedback rather than a constant scale error.

![R19 final profile](results/krueger_2024_r19_response_check/remote_artifacts/profile.png)

The final profile is a physically plausible deep, narrow trench with no sealed cavity. Its smooth
completion, exact conservation, and bounded runtime are evidence for engine operation—not yet for
predictive accuracy.

## What this result does and does not earn

It earns:

- confidence that the repaired 10 nm engine runs the complete Krüger base case unattended in under
  eight minutes on one 4090;
- an excellent coarse-grid mask-opening development result;
- a hard rejection of stale response/Jacobian reuse across operator epochs;
- closure of the 10 nm candidate sequence before a costly fine run;
- promotion of source-epoch binding and the AMR/fine numerical-authority problem.

It does not earn:

- a calibrated parameter freeze;
- a held-out oxygen or power reveal;
- a formal Krüger validation claim;
- permission to tune another 10 nm point;
- evidence that the remaining 3.42% depth error is chemistry rather than late grid discrepancy.

The last point is decisive. The existing late 10/5 evidence is not a clean same-operator pair, and
its observed difference is larger than the current `28.2 nm` depth miss. Adjusting chemistry before
late numerical authority would fit a numerical discrepancy.

## Earned next branch

No further long profile is launched immediately. The bounded order is:

1. bind every calibration response and held-out run to an executable operator-epoch checksum;
2. complete the bounded positive-horizon frozen-chemistry micro-step needed to expose the oxide-yield
   direction that `dt=0` cannot observe;
3. preserve the completed frozen R17/R19 access diagnostic and do not promote its reduced quadrature
   to endpoint authority;
4. certify sparse narrow-band/AMR geometry, surface-state transfer, and the paper-defined opening on
   manufactured translation/recession/pinch-off cases;
5. reproduce the existing uniform 10/5 nm `0.5 s` operator agreement with AMR;
6. only then choose one current-operator authoritative base confirmation at uniform 5 nm or
   certified AMR;
7. keep all held-out profile outcomes sealed until that base authority gate passes.

This is the adaptive branch defined by the large execution program: a failed cheap response model
promotes numerical-authority work, not solver roulette or parameter churn.
