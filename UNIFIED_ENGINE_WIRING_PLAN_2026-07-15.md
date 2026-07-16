# Unified engine wiring plan and status

Date: 2026-07-15

This record starts from executable common-engine code. It does not authorize a validation claim,
alter `CCA-2026-07-13-R2`, reopen frozen-map root solvers, or redirect the legacy `Process` facade.
All work is local; no remote push is part of this increment.

## One engine, two scientifically different operating modes

```text
 PlasmaBoundaryState + material geometry + declared mechanisms
                              |
                              v
        hard-visibility neutral/charged transport + certified replay
                              |
                  reflection / SEE response cascade
                              |
                conservative face-resolved arrivals
                              |
          +-------------------+-------------------+
          |                                       |
  mean-field physical time                    finite arrivals
  fresh ensemble currents                 Poisson(flux*area*time)
          |                                       |
          +---------- Q1 charge / field ----------+
                              |
                    MaterialMechanismRouter3D
                     /                      \
             substrate mechanism       mask mechanism
                     \                      /
                      emitted product ledger
                              |
              optional bounded diffuse redeposition
                 (same material only in v1)
                              |
             signed level-set motion + C1 charge/state remap
                              |
                        next geometry
                              |
          notch/bow observables or twist ensemble statistics
```

Mean-field mode predicts the expected charged profile. Finite-arrival mode predicts a distribution of
profiles and is the only mode eligible for a stochastic-twist claim. They share transport, Poisson,
surface mechanisms, interface evolution, conservation ledgers, manifests, and checkpoints.

## Wiring status

| Work package | Status | Executable result | Claim still withheld |
| --- | --- | --- | --- |
| W1 continuation/replay | complete | In-memory continuation and versioned `allow_pickle=False` checkpoint preserve geometry, face charge, material-namespaced surface state, mesh fingerprint, and source-manifest hash | Fresh-environment migration remains a release gate |
| W2 profile observables | complete | Geometry-native notch depth, left/right asymmetry, bow width, hole centerline, displacement, equivalent diameter, and onset AR have manufactured refinement gates and ensemble confidence intervals | No experimental profile has been scored |
| W3 C5 statistical runner | complete infrastructure | Nested N/2N and paired sample-doubling runner plus isotropy/systematic-direction gates; public ensemble requires physical Poisson arrivals | N>=30 AR campaign has not run, so twist probability is not validated |
| W4 C4 packaging | operational replay complete | Checksum/pixel replay ingest, commit/reveal split, maximum-two-parameter calibration, exact held-out coverage, uncertainty ladder, charging-off causality contract, and the installed `petch-nozawa` open-area replay all execute through the common engine | The paper does not report measurement uncertainty; signed charging saturation, numerical ladders, and untouched held-out solves have not passed, so validation remains withheld |
| W5 materials/mass loop | partial but executable | Material-ID router independently advances mask and substrate laws; product ledgers remain face-exact; opt-in diffuse product return converts explicitly bounded same-material sticking to signed growth | Production mask chemistry, reactive SiO2 product branching, and cross-material film creation still require evidence and a new material layer |
| W6 hardening | substantially complete | Versioned run/checkpoint schemas, geometry/boundary/operator hashes, machine-readable mechanism/redeposition provenance, device declaration, inline-recovery/error budget, durable campaign supervisor, and heartbeats | Fresh-install replay and an accuracy-matched published CPU/GPU manifest remain release gates |
| W7 ordinary energetic-response/profile bridge | complete | Ordinary and charged profile modes now share the certified material response/reimpact cascade. Response-enabled field-free ions use the zero-field certified tracer; periodic transport now drives periodic Godunov motion and wrapped redistancing, preventing moving seam artifacts | Reflection remains a bounded sensitivity until material-specific differential data and independent morphology validation exist |

## ARDE migration increment closed on 2026-07-15

The git-history reconciliation found that evolving profiles, Belen kinetics, Knudsen reductions,
reflection, and legacy redeposition had already been built in the compatibility engine. Only the proven
pieces were migrated: the common Belen mechanism and radiosity were retained; legacy calibrated Knudsen
and AR-shaped passivation closures were not promoted; the existing certified response cascade was wired
into ordinary feature evolution.

The source-correct reflected-ion Figure-9 development diagnostic improved the twelve-point RMSE from
3.545 to 2.557 micrometres but did not validate the model. A nested three-level moving-profile gate then
completed with roundoff-level charge/state ledgers and exact energy closure. Its latest sampling delta
contracted to 0.000220 micrometres, and the moving profile finished 0.001184 micrometres shallower than
the identical frozen-geometry rate counterfactual. These are mechanism-development results; all Figure-9
markers remain exposed development data.

## The deliberately narrow redeposition closure

The engine now closes the loop only when all of these are explicit:

1. the surface mechanism names a population and closes it against removed material face-by-face;
2. its energy and angular launch law is declared;
3. every target material has a bounded, sourced sticking probability;
4. the deposited material density is bounded and sourced;
5. reacted material lands on an already evolving face of the same material.

The returned material is split exactly into deposited plus escaped inventory. Deposited thickness is
converted to a negative recession velocity, so the existing level-set and charge-remap path sees real
growth rather than a diagnostic-only flux. A positive cross-material sticking event refuses because a
mask-derived coating on substrate is a new material state/level set, not “more substrate.” This is the
smallest honest implementation that can support physical resputtering and redeposition studies.

## Execution plan from here

### P1 — Close the C4 scientific gates, then score it

- Preserve the ingested Nozawa/Hwang source images, checksums, axis/pixel transforms, marker pixels,
  and digitization uncertainty; the exact installed CSVs are now checksum-bound.
- Preserve the already committed calibration family and at-most-two-parameter protocol.
- Resolve the absent experimental measurement uncertainty without inventing one; until then the
  numerical engine may be developed but the headline validation scorer must continue to refuse.
- Run charging-off/on, grid, timestep, and sample ladders through `PhysicalChargingProcess`.
- Score the two untouched families once. Outside combined uncertainty produces a decomposed miss, never
  held-out retuning.

### P2 — Run the C5 ensemble campaign

- Use `physical_time_resolved` plus dimensional Poisson arrivals on a 3-D hole AR sweep.
- Use base N>=30, nested N doubling, paired `n_position` doubling, and a symmetric isotropy control.
- Report probability and confidence interval for onset AR/displacement. A single realization remains a
  visualization, not a forecast.

### P3 — Promote additional material physics only when P1/P2 demand it

- Add a real mask law from material-specific yield/chemistry data.
- Add reactive SiO2 branching only from sourced product identities and launch distributions.
- Add cross-material coatings as a new level-set material with conservative state transfer; do not weaken
  the v1 refusal.
- Run the already specified SEE/conduction parameter box last, and promote neither channel without a
  charging-off/on morphology causality result.

### P4 — Release hardening

- Replay one signed manifest in a clean environment.
- Publish accuracy-matched CPU/CUDA event, conservation, and profile comparisons.
- Exercise supervisor recovery at each recoverable refusal and corruption refusal at each hard-stop path.
- Freeze the public schemas only after migration/version-refusal tests pass.

## Current claim boundary

The code can execute a first-principles mean-field charged notch/bow mechanism, material-specific mask
and substrate motion, conservative same-material sputter redeposition, and a physically dimensioned
finite-arrival 3-D ensemble through one engine. It has the numerical packaging needed to test Nozawa
notching and statistical twisting, and the Nozawa open-area condition now has a packaged user replay
with an operational exact-operator smoke and scientific refusal checkpoint. It cannot yet claim validated notch depth, predictive twist
probability, universal mask chemistry, reactive-product redeposition, or cross-material film growth.
