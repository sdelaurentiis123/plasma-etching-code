# Codex takeover report — reactor-to-depth campaign, exact state at 2026-08-21 16:05 UTC

Repository: `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code`

Branch: `codex/validation-first-multiphysics`

Pre-report committed HEAD and remote: `e3bc275` (`Implement post-seal Bosch
heldout scorer`)

This is the authoritative continuation report. It supersedes the live-state
claims in `CODEX_TAKEOVER_REPORT_2026-08-21_1515UTC.md` and every earlier
Oxford moving-Cr handoff. Older files remain forensic records, but their
process counts, Bosch status, and Oxford diagnosis are stale.

## Executive answer: what is actually happening

There are three different experiments, not one tangled failure.

1. **Bosch Si/SF6-C4F8 reactor-to-wafer depth has passed a genuine sealed
   chronological heldout test.** The prediction was written, hashed, committed,
   and pushed before numeric outcome reveal. On 13 unseen wafers with 89 radial
   measurements each, the model reaches 0.239 um wafer-mean Si-depth MAE
   (0.541% MAPE), 0.330 um point RMSE, 0.262% normalized radial-shape RMSE,
   0.024 um oxide-loss MAE, and 3.16% selectivity MAPE. Every frozen absolute
   gate and empirical-baseline gate passes. This is the strongest absolute-
   depth result in the repository.
2. **Oxford/Freddie is stopped on a newly localized numerical bookkeeping
   defect, not unknown physics and not a hang.** One sparse remap coefficient
   becomes `-2.2204460492503131e-16`. That is one binary64 ulp of negative
   roundoff created by closing a normalized row onto its final coefficient.
   The fail-closed validator correctly rejects it. The implementation must
   close the row through a safely positive coefficient while retaining exact
   sum and conservation; the validator must not simply be weakened. Twenty-
   three of 56 clean v6 trajectories are already quarantined and committed.
3. **Three Krueger C4F6/Ar/O2 no-fit feature forecasts are healthy and still
   running.** At the last inspection they were at 25.4-27.0 of 60 physical
   seconds, with depths of about 226-259 nm. They are nonlinear evolving-
   geometry calculations and cannot be extrapolated honestly. The published
   825 nm target is not matched yet. The old website's 790-811 nm match remains
   retracted because it came from canceling bugs.

The mission is therefore neither a failure nor finished. One chemistry/tool
now has a real heldout reactor-to-depth win. Freddie still lacks a frozen,
complete feature board and target SEM score. Krueger still lacks a supported
species-resolved reactor boundary and its final independent forecasts.

## The claim boundary

The current evidence supports this statement:

> The stack can deterministically map measured machine time series through a
> conserved reactor/wafer boundary to absolute wafer-scale Si and oxide depths,
> selectivity, radial shape, and within-lot drift on one heldout Bosch tool
> dataset. The common feature engine separately has strong analytic transport,
> conservation, charging, moving-material, and direct-beam receipts.

It does **not** yet support these stronger statements:

- arbitrary tool knobs produce atomically accurate feature profiles;
- Freddie's TiO2 pillars have been matched;
- Krueger's 825 nm trench has been reproduced from published inputs;
- the Bosch heldout result validates feature charging, scallops, sidewalls,
  ARDE, or a different machine;
- any chemistry can be added by naming gas species without an independently
  evidenced gas/surface reaction deck.

Atomic grid spacing is not atomic predictive accuracy. A defensible prediction
must carry boundary uncertainty, surface-law uncertainty, numerical
convergence, and experimental metrology uncertainty. The architecture can
approach atomic-layer accuracy for an ALE process with the required data; this
campaign has not demonstrated that accuracy for these plasma etches.

## Repository and safety state

Use only the repository and branch named at the top. Do not switch to older
branches or treat `reactorlab`, `petch-torchsim`, or website-era HTML as the
authoritative implementation.

At the start of this report the branch matched its remote at `e3bc275`.
Campaign-owned, intentionally uncommitted post-reveal Bosch artifacts were:

- `data/experimental/zenodo_17122442/revealed_heldout_Si_Oxide_etch_89_points.csv`;
- `data/experimental/zenodo_17122442/revealed_heldout_Si_Oxide_etch_89_points_manifest.json`;
- `results/curated/zenodo_bosch_wafer_boundary_map_depth_extension_v8/heldout_score.json`;
- the post-reveal README and artifact-test updates made with this report.

They must be committed with exact paths after their tests pass.

The following unrelated user-owned paths must never be edited, staged,
deleted, reset, cleaned, or absorbed into campaign commits:

- `results/curated/mixed_layer_feature_v1/ml20-bonds-12s.log`;
- `results/curated/mouth_equilibrium_probe_dx/`;
- `scratch_ignore_calc.py`.

Do not run destructive Git commands. Stage exact files only.

## Paid compute

Only Vast instance `48177892` belongs to this campaign.

- SSH: `ssh -p 17892 root@ssh6.vast.ai`
- Oxford tree: `/root/petch-4b656fd`
- Krueger tree: `/root/petch-d852a1f`
- environment: `/root/petch-venv`
- hardware: RTX 3090, 24 GiB GPU, 28 vCPUs
- last recorded price: approximately `$0.2011/hour`

Do not touch any other instance. Do not destroy `48177892` until the finished
Oxford and Krueger artifacts are copied locally, hash-verified, audited,
committed, and pushed. Inspect exact PID files or `ps -p`; never use `pgrep -f`
as proof of liveness because the probe can match itself.

## Track A — Bosch: successful heldout reactor-to-wafer depth

### Physical task

The public Sayyed et al. Zenodo 17122442 record contains SPTS Omega i2L DSi
Rapier machine traces and independent wafer maps for a Bosch Si etch:

- alternating SF6 etch and C4F8 passivation phases;
- source/platen RF, reflected power, platen Vpp/DC, pressure, gas, current,
  thermal, and backside-helium channels;
- 89 spatial Si-depth and oxide-mask measurements per measured wafer;
- chronological process records spanning calibration and heldout dates.

The model uses the measured machine trajectory, deterministic reduced
reactor/wall memory, frozen Belen silicon and oxide mechanisms, and a positive
wafer ion-transmission field. The field is a complete real-Zernike expansion
whose finite-volume normalization conserves total positive-ion current. It
changes ion delivery only: it does not apply a depth correction, change
neutral flux, alter ion energy, or refit the surface law.

### Preregistered selection and exact replay

Twenty calibration-only candidate bases were evaluated. The selected lowest-
capacity passing model has static order 9, dynamic order 2 driven by measured
C4F8 platen Vpp, and 59 coefficients. The exact Jacobian is full rank 59/59,
condition number 9.954, maximum parameter correlation 0.714, no bound contact,
and maximum log field 0.14744 below the frozen `ln(2)` ceiling.

It was replayed through the exact reactor/surface recurrence. Whole-lot
leave-one-lot-out results were:

| Metric | Physics | Frozen empirical baseline | Result |
|---|---:|---:|---|
| Si wafer-mean MAE | 0.270179 um | 0.338486 um | pass |
| Si point RMSE | 0.431199 um | 0.486585 um | pass |
| normalized shape RMSE | 0.635267% | 0.636619% | pass, narrow |
| within-lot slope MAE | 0.077810 um/wafer | 0.082903 gate | pass |

The independent response-tensor midpoint audit and denser reactor/Zernike
refinement gates also pass. Maximum refinement movement is 0.043559 of any
frozen gate.

### Seal and heldout result

The prediction SHA-256 is:

```text
56ed2429832fe77280762fbca86cb6ffa4de3fd9687aa84f3b5cfd4ca99a3b1a
```

It was committed and pushed in `fd39277` before a separate extractor parsed
heldout numeric outcomes. Thirteen of 20 chronological heldout process records
have 89-point outcomes; seven are absent from the source and are explicitly
reported missing, never imputed.

Heldout score:

| Metric | Physics | Frozen gate |
|---|---:|---:|
| Si wafer-mean MAE | 0.238906 um | <=1.0 um |
| Si wafer-mean MAPE | 0.541245% | <=3.0% |
| Si 89-point RMSE | 0.329751 um | <=1.5 um |
| normalized radial-shape RMSE | 0.261855% | <=2.0% |
| oxide wafer-mean MAE | 0.023908 um | <=0.08 um |
| selectivity MAPE | 3.161040% | <=12.0% |
| Si wafer-mean correlation | 0.873410 | descriptive |

The physics result also beats:

- the calibration-global mean-depth baseline: 0.238906 vs 0.400599 um MAE;
- the calibration mean-map baseline: 0.329751 vs 0.477852 um point RMSE;
- the calibration mean-map shape baseline: 0.261855% vs 0.294344%;
- the zero-slope within-lot baseline: 0.083713 vs 0.130151 um/wafer slope MAE.

Twenty-thousand wafer-level bootstrap replicates preserve positive 95%
improvement intervals for mean depth and pointwise depth. The shape-improvement
interval is `[-0.0118, 0.0714]` percentage points, so the shape baseline is
beaten by the full-sample score but the advantage is not statistically
resolved at 95%. Do not suppress that qualification.

Post-reveal scoring never changes the sealed prediction. The local score and
reveal checks pass. The scope is one tool and one Bosch process family.

### Why this is a real breakthrough

This is not an optimizer fitting a target depth after seeing it. The chronological
answer was sealed before reveal, exact replay and refinement were binding, and
the model beats both absolute tolerances and empirical baselines. It establishes
that a deterministic, differentiable, conserved machine-to-wafer boundary can
carry enough reactor physics and wall history to predict absolute depth and
radial transfer on heldout wafers.

It is the reactor-side proof needed by the wider program. The next scientific
extension is to feed this validated boundary into feature-scale geometries and
test ARDE/scallop/profile outcomes without changing the core.

## Track B — Oxford/Freddie TiO2/Cr square-pillar board

### Frozen supplied process

- Oxford PlasmaPro NPG80 RIE;
- 55/5/1 sccm CHF3/SF6/O2;
- 30 mTorr;
- 150 W forward table RF;
- 20 C table temperature;
- 1,200 s process;
- 700 nm ALD TiO2 on fused silica;
- 45 nm Cr hard mask;
- square pillars, approximately 400 nm pitch;
- nominal widths 80, 120, 160, 200, 240, 280, and 320 nm.

The target SEM has not been used for tuning. The current board propagates two
independently sourced TiO2:Cr selectivities, low/high self-bias energy cases,
and core/core-plus-tail angular cases: 7 widths x 2 selectivities x 2 energies
x 2 angular cases = 56 trajectories.

The reactor/wafer boundary is much more than the old forward-power heuristic:
it carries 67 plasma species, electron kinetics, sheath closure, parent and
daughter collision chemistry, and radial axisymmetric transport into a
species-resolved feature boundary. The feature solver then evolves TiO2 and Cr
interfaces, surface state, charging, sputtering/chemical removal,
polymerization, and material extinction deterministically in 3-D.

### What was fixed before the current stop

The historical 55/56 v3 CPU spin was real but has been solved. Its chain
included a marching-cubes sliver, unequal-area flux projection, undefined mask
extinction behavior, unsupported material-owner nodes, and non-idempotent
regional cleanup. The corrected v6 engine has explicit material-extinction
lifecycle, exact owner projection, conservative ledgers, and regression
receipts. The formerly failing acceptance trajectory completed.

One committed acceptance case, width 200 nm/selectivity 18.016664/high energy/
zero tail, reached:

- depth 679.407775 nm;
- middle CD 196.542208 nm;
- bottom CD 203.951939 nm;
- sidewall angle 81.4639 degrees;
- bow 94.566238 nm;
- particle residual zero;
- remap residual `9.37e-16`.

That is a conditional pre-SEM result, not an experimental match.

### Exact current failure

The clean v6 campaign stopped after 23/56 unique caches on:

```text
width_nm    = 120
selectivity = 18.016664028610727
ion_energy  = 296.20650777233976 eV
tail        = 0.0
duration    = 1200 s
```

The original aggregate error was instrumented without changing acceptance
rules. The exact reproduction now fails with:

```text
invalid sparse surface-transfer weights:
negative_weights count=1, first_entry=24516,
minimum=-2.2204460492503131e-16
```

The root mechanism is visible in `build_surface_transfer_3d`: inverse-distance
weights are normalized, then the final coefficient is assigned
`1 - sum(previous coefficients)` to close the row. In this row the preceding
floating-point coefficients round to slightly more than one, so the final
coefficient becomes one ulp negative. A probability/interpolation weight that
is truly negative would be unphysical; this magnitude and construction prove
the observed value is numerical closure roundoff.

The correct repair is narrow:

1. apply the normalization residual to a safely positive coefficient, normally
   the largest coefficient, instead of forcing it into an arbitrarily tiny
   final coefficient;
2. prove all weights nonnegative and the row sum closed under the exact stored
   summation order;
3. pin the production checkpoint as a deterministic regression;
4. prove intensive identity/manufactured transfer, extensive per-material
   conservation, material locality, CPU/CUDA agreement, and fingerprint
   determinism;
5. compare the corrected trajectory to a high-precision reference;
6. bump the model revision and recompute all 56 if any stored trajectory can
   change beyond diagnostics. Never mix revisions.

Do **not** merely clip negatives without re-closing the row, and do not relax
the validator tolerance. The guard found a real implementation flaw even
though the triggering magnitude is roundoff.

### Current Oxford artifacts

The 23 clean v6 caches, failure log, inventory, and hashes are quarantined in:

`results/quarantine/zhu_npg80_v6_failure_20260821/`

They are committed in `0dbbd08`. Whether they remain cache-compatible depends
on exact-equivalence evidence after the repair. The remote diagnostic log is
`/root/zhu_v6_w120_diag.log`. The standalone worker has terminated after
producing the named failure; it is not silently billing a GPU loop.

### What must happen before an Oxford SEM claim

1. Land the remap correction with the production regression and focused suite.
2. Decide revision/cache compatibility honestly.
3. Complete all 56 frozen trajectories and assemble the audit/atlas.
4. Commit and push the frozen board before opening the target SEM.
5. Digitize target depth, top/middle/bottom CDs, sidewall angle, bow, Cr
   survival, collapsed/missing-pillar rate, and spatial position.
6. Score the frozen ensemble; do not select whichever uncertainty corner looks
   closest after reveal.

### What to ask Freddie

The supplied recipe is enough to run the blind uncertainty board. These items
turn it into a uniquely identified machine/profile test:

- cross-section and top-down SEMs for this exact recipe, with scale bars;
- the width/pitch/GDS mapping for each imaged region;
- sample radius and orientation on the powered electrode;
- recorded DC self-bias or electrode voltage/current waveform;
- blanket TiO2 thickness loss under the same recipe;
- remaining Cr thickness after the same recipe;
- whether Cr was stripped before imaging and the rinse/dry history.

The SEM is an answer key, not a fit input. Self-bias constrains ion energy;
blanket TiO2/Cr losses isolate the surface law from feature transport. Those
three measurements are the highest-value additions.

Clustered fallen pillars do not by themselves prove an uneven top-down plasma
flux. Smooth reactor nonuniformity, local pattern loading, mask adhesion or
undercut, feature weakening, micromasking, and rinse/dry collapse can all make
spatial clusters. The simulator should test those hypotheses against sample
coordinates and the intact profile instead of assuming the cause.

## Track C — Krueger C4F6/Ar/O2 SiO2 trench

### Honest historical verdict

Krueger reports approximately 825 nm after 60 s. The old website-era 790-811
nm agreement is retracted: it came from two bugs canceling. Under the corrected
published-input mechanism the measured 60 s result is approximately 346.833
nm. Direct beam evidence implies that reaching 825 nm at the paper's printed
aggregate ion flux requires a removal yield above the supported physical
ceiling.

That does not prove the experiment is wrong. It proves the published boundary
is insufficient or internally inconsistent for absolute depth. Plausible
missing quantities include the true ion-flux normalization, species-resolved
positive-ion mixture, stable C4F6/reactive neutral channels, and ion-energy/
angular distributions at the wafer.

### Current independent forecast

The active forecasts transfer Guo planar/mixed-layer surface laws without
reading the Krueger depth target during execution. Three frozen ion-identity
hypotheses are run because the paper reports an aggregate positive-ion channel:

- nominal unresolved aggregate ion;
- all aggregate ion declared CF2+;
- all aggregate ion declared CF3+.

At 2026-08-21 16:02 UTC all three supervisors and one child each were healthy:

| Case | Accepted time | Depth | Mask opening |
|---|---:|---:|---:|
| nominal unresolved | 26.969 s | 226.444 nm | 18.342 nm |
| all CF2+ | 25.703 s | 254.390 nm | 22.547 nm |
| all CF3+ | 25.391 s | 259.062 nm | 21.809 nm |

The children checkpoint every 30 wall minutes and resume under supervisors
10284-10286. All must reach a physical terminal result. Do not extrapolate
linearly and do not choose the target-closest ion case after completion.

### Why the current Krueger runs are not yet certified predictions

Their own receipts set:

```text
parameter_evidence_supports_prediction = false
within_declared_scope = false
```

Reasons include neutrals outside Guo's printed source list, energy beyond the
<=370 eV regression board, a declared angular-law typesetting repair, and
unresolved aggregate ion composition. The Karahashi mass-selected-beam transfer
also shows that the frozen Guo law overpredicts F+/CF+ and predicts deposition
at measured 250/500 eV CF3+ points. Therefore the CF2+/CF3+ cases are valuable
sensitivity/falsification runs, not validated species-resolved closures.

The path to a defensible Krueger match is one of:

1. species-resolved wafer flux/IEDF/IADF data from that reactor;
2. a validated C4F6 global/spatial reactor model constrained by independent
   diagnostics and blanket rates, then tested on the trench;
3. same-tool blanket SiO2 and mask-loss measurements that identify the boundary
   normalization without fitting the feature depth.

An optimizer-adjusted yield chosen to land at 825 nm is not a prediction.

## Other validation and generalization state

The common core is generalizable by mechanism/deck, not automatically
chemistry universal.

What transfers across chemistries:

- deterministic particle transport and analytic Clausing limits through
  aspect ratio 200:1;
- conserved species-resolved wafer boundary objects;
- differentiable reduced reactor, sheath, and radial transport operators;
- self-consistent feature charging and ion deflection;
- moving multi-material 3-D level sets;
- intensive/extensive surface-state remapping with per-material ledgers;
- polymer deposition/crosslink/breakage state machinery;
- mask erosion/extinction and topology handling;
- exact JVP/VJP-compatible smooth closures where implemented;
- target firewalls, preregistration, hashing, replay, refinement, and
  heldout scoring.

What must be supplied or independently inferred per chemistry/tool:

- gas-phase species and rate coefficients over the actual EEDF/temperature;
- electron-impact cross sections and ionization/dissociation branching;
- wall recombination/sticking/desorption and wall-history laws;
- sheath/electrode transfer or measured waveforms;
- ion- and neutral-resolved surface reaction probabilities/yields versus
  energy, angle, coverage, temperature, and material;
- mask and substrate material response;
- independently measured tool geometry and boundary/blanket observables.

Existing cross-chemistry receipts remain useful but have different strength:

- analytic molecular transport and conservation: strong numerical/physics
  validation;
- Gray/Kim/Karahashi direct-beam boards: direct surface/beam constraints with
  named passes and misses;
- SF6/O2 silicon/de Boer: transfer validation, with source constants declared
  profile-fitted rather than blind;
- Bosch SF6/C4F8 Si: now a sealed heldout absolute reactor/wafer-depth pass;
- Oxford CHF3/SF6/O2 TiO2/Cr: blind conditional feature board unfinished;
- Krueger C4F6/Ar/O2 SiO2: endpoint and boundary identifiability unresolved.

The library should remain the provenance spine. Every new reaction coefficient
must enter with source, exact support, unit conversion, evidence grade, and
fit/validation exposure; it should never become an unlabeled JSON knob.

## Immediate continuation sequence

1. **Land Bosch post-reveal evidence.** Run extractor, scorer, prediction seal,
   and focused tests in check mode; stage only the five campaign files plus
   this report; commit and push. The model and sealed prediction are immutable
   after reveal.
2. **Fix Oxford remap closure.** Add the captured production regression, apply
   the residual to a safely positive weight, prove nonnegativity/conservation/
   determinism/high-precision agreement, and decide revision compatibility.
3. **Restart/complete Oxford.** Deploy only the reviewed commit to the existing
   box, finish all 56 cells, retrieve/hash/check locally, render the blind
   atlas, and commit before SEM reveal.
4. **Continue monitoring Krueger.** Leave the three healthy supervisors alone;
   retrieve terminal results for all cases and grade the preregistered envelope
   without target-based selection.
5. **Run focused then broad certification.** At minimum include Bosch artifact
   gates, sparse-transfer and Oxford moving-mask regressions, Guo/Krueger
   forecast checks, then the full suite. Record exact counts and commit hashes.
6. **Score Freddie only after freeze.** The target SEM becomes a heldout
   measurement. If the board misses, decompose reactor boundary, blanket
   surface rate, mask survival, feature transport, charging, and post-etch
   mechanics rather than fitting the final contour wholesale.

## Terminal definition for the overall goal

The active goal is not complete until all are true:

- the Bosch heldout artifacts are landed and reproducible;
- the Oxford blind board completes numerically and is scored against the exact
  target SEM without post-reveal parameter selection;
- the Krueger forecasts terminate and their evidence-domain failures are
  reported;
- either Krueger depth is reproduced using independently constrained reactor/
  surface data, or a concrete minimal experiment is demonstrated necessary;
- at least two additional chemistries show heldout absolute depth/profile
  transfer with the unchanged common core;
- numerical convergence, conservation, provenance, uncertainty, and target
  firewalls remain green.

The scientifically correct current headline is:

> A machine-to-depth breakthrough has occurred on heldout Bosch wafers. The
> Freddie feature-scale prediction is one narrow numerical remap repair and a
> completed blind board away from scoring. Krueger remains an unresolved
> published-boundary problem, with independent no-fit feature forecasts still
> running. The architecture is general; atomic-level and arbitrary-chemistry
> accuracy are not yet demonstrated.
