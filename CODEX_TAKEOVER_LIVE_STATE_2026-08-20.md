# Codex takeover — live reactor-to-depth state

Snapshot: 2026-08-20 22:39 EDT / 2026-08-21 02:39 UTC

Repository: `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code`

Branch: `codex/validation-first-multiphysics`

Pushed HEAD before this report: `910f981`

This is the current operational handoff.  It supersedes the live-state and
next-action sections of `CODEX_TAKEOVER_FULL_STACK_2026-08-20.md` and
`HANDOFF_MOVING_CR_BOARD_2026-08-20.md`.  Those older files remain useful
for forensic history, but their Oxford cache counts, Bosch implementation
status, and recovery instructions are stale.

## Executive verdict

The campaign is not lost, and the old “55 of 56 trajectories plus an unknown
CPU spin” account is no longer the active state.

Three distinct scientific tracks are underway:

1. **Oxford TiO2/Cr blind profiles:** two real numerical defects were found,
   fixed, tested, and pushed.  The corrected v5 board reused its valid caches
   and reached 39 of 56 cells.  It then stopped safely on a new deterministic
   sub-grid material-ownership cleanup failure in one cell.  This is an engine
   topology defect to reproduce and fix; it is not an experimental mismatch.
2. **Krueger C4F6/Ar/O2 SiO2 depth:** three frozen, no-target-fit feature
   forecasts are healthy at about 9 simulated seconds of 60.  Their current
   depths are 108.746, 132.258, and 140.104 nm for the nominal unresolved-ion,
   all-CF2+, and all-CF3+ cases.  These prefixes cannot be extrapolated linearly
   or reported as final answers.  The verified published-input common model
   remains about 346.833 nm versus the paper's 825 nm.  No Krueger match is
   currently claimed.
3. **SPTS Bosch Si/C4F8 reactor-to-wafer depth:** the deterministic
   waveform-to-species-to-2-D-wafer-map-to-surface path is implemented.  On all
   75 calibration wafers, one shared bounded equipment parameterization reaches
   0.893% Si mean-depth MAPE, 0.840 um pointwise RMSE, 1.647% normalized-map
   RMSE, 6.269% oxide-loss MAPE, and 6.112% selectivity MAPE without changing
   the Belen/La Magna surface laws.  It passes the absolute calibration gates.
   It does **not** yet beat leave-one-lot-out empirical baselines, so the heldout
   outcomes remain sealed and no predictive pass is claimed.

The honest overall state is therefore: the software now contains a serious,
conservative, deterministic reactor-to-depth architecture and calibration
absolute scale is reachable without a depth multiplier.  Absolute heldout
depth/profile accuracy is not yet closed across the mission's three required
chemistries.

## Repository safety and ownership

Use only:

```text
/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code
codex/validation-first-multiphysics
```

Do not switch to the older branches or use another repository as the campaign
source of truth.

The only non-clean paths at this snapshot are unrelated, pre-existing user
work and must remain untouched:

- `results/curated/mixed_layer_feature_v1/ml20-bonds-12s.log`;
- `results/curated/mouth_equilibrium_probe_dx/`;
- `scratch_ignore_calc.py`.

Do not reset, clean, delete, stage, or absorb them.  Before this report was
written, no campaign file was uncommitted; the latest Bosch calibration
searches existed only in the console and are recorded below.

## What changed after the earlier full-stack handoff

The pushed sequence after the earlier `e5aed0f` checkpoint is:

- `24c423c` — refresh the full-stack reactor/depth handoff;
- `32a6613` — couple cylindrical point fluxes to the exact fused surface law;
- `cb16a9d` — freeze the Bosch cylindrical v2 extension and record the
  axisymmetry model-form falsification;
- `58b41ff` — reuse one cylindrical sparse factorization for source fitting;
- `02ec3da` — vary the shared radial source without refactorization;
- `8b5d797` — preregister the species-resolved radial-source v3 extension;
- `787672f` — implement species-resolved radial source moments;
- `167bafc` — preregister the wafer-edge sheath v4 extension;
- `910f981` — implement a current-conserving sub-grid wafer-edge ion focus.

The reusable cylindrical solve reduced calibration cost from about 1.65 s to
about 0.012 s per source evaluation.  This is an acceleration of the same
finite-volume operator, not a surrogate replacement.

## Live paid compute

Only Vast instance `48177892` belongs to this task:

- label `petch-zhu-moving-cr-bigram`;
- SSH `ssh -p 17892 root@ssh6.vast.ai`;
- Oxford tree `/root/petch-4b656fd`;
- Krueger tree `/root/petch-d852a1f`;
- environment `/root/petch-venv`;
- RTX 3090, 28 vCPU, approximately 251 GiB host memory;
- last recorded price about `$0.2011/hour`.

Do not touch any other Vast instance.  Do not destroy `48177892` until all
useful Oxford caches and all Krueger outputs/checkpoints have been copied,
hash-verified locally, and committed as appropriate.

Always inspect the exact PID from its PID file.  Do not use `pgrep -f`; an
earlier liveness probe matched its own command line and produced a false
positive.

## Oxford/Freddie blind-board state

### Frozen experimental condition

- Oxford PlasmaPro NPG80 RIE;
- 55 / 5 / 1 sccm CHF3 / SF6 / O2;
- 30 mTorr;
- 150 W forward table RF;
- 20 C table;
- 1200 s;
- 700 nm ALD TiO2 on fused silica;
- 45 nm Cr hard mask;
- square-pillar prior, approximately 400 nm pitch;
- widths 80, 120, 160, 200, 240, 280, and 320 nm.

The current conditional boundary contains 67 plasma species, including 20
positive ions and 37 thermal neutrals.  Its central positive-ion flux is about
`1.457e19 m^-2 s^-1`.  The conditional sheath solution uses approximately
296 V across the powered-electrode sheath, about 299 eV for singly charged
ions.  The uncertainty board explicitly carries ion energy, angular-tail, and
TiO2:Cr selectivity rather than selecting one answer from Freddie's SEM.

### Bugs already fixed

The original alleged spin was captured and explained.  Periodic strip
symmetrization preserved integrated particle rate by splitting it equally
among triangles and then dividing by recipient triangle area.  A tiny
marching-cubes sliver therefore received an enormous flux density, created a
local recession speed near 3310 mesh units/s, and requested 5,331,155 CFL
substeps.  `d8a8d02` changed the projection to preserve uniform flux density
and distribute integrated rate in proportion to triangle area.  The exact hard
cell completed afterward.

The corrected board then exposed complete Cr loss.  `84c2c01` added a
certified material-extinction lifecycle, and `44956d0` allowed this declared
Oxford campaign to continue conservatively through Cr extinction while
preserving the remaining material, topology, and ledgers.  The exact hard v5
acceptance trajectory was committed in `bdabd8a`.

### New failure at this snapshot

The v5 parent PID file still contains PID `8980`, but that process is no longer
alive.  Exactly 39 trajectory JSONs have a
`job_spec.model_revision` ending in `v5`; do not count all JSONs because the
directory intentionally contains v1-v4 forensic caches.

The parent stopped on this exact cell:

```text
width_nm=200
selectivity=18.016664028610727
impact_energy_eV=296.20650777233976
impact_energy_label=high
tail_fraction=0.0
elapsed=719.862 s
accepted_steps=149
Cr centre remaining=19.850306886455172 nm
Cr interior minimum=0.000 nm
```

The actual guardrail was:

```text
subcell material-island cleanup did not reach a fixed point;
passes=12;
coordinates=((14,12,80),(15,12,80),(15,28,80),(28,15,80),(28,27,80))
```

The path is `_restore_unresolved_material_ownership` inside
`advance_feature_step_3d`, around `src/petch/feature_step_3d.py:2989-3019`.
Redistancing after one subcell ownership repair can expose another subcell
component.  The code iterates that repair and deliberately refuses to continue
if it has not reached a fixed point after 12 passes.

Do **not** merely raise the pass cap and call the problem solved.  First make a
single-cell deterministic reproduction from this exact job spec; record the
unresolved-node set and owner before/after each pass; determine whether the
operator converges slowly, oscillates between ownership states, or creates a
new component through redistancing.  The acceptable fix is an idempotent or
monotone ownership/topology update with a proved postcondition and particle/
volume ledger preservation.  Add a regression for the exact cell, then rerun
the v5 driver.  The 39 content-addressed v5 caches should be reused unchanged.

The Oxford board is still a blind physical envelope, not an SEM match.  The
answer-key SEM has not been used for tuning.  The most valuable measurements
remain exact-run DC self-bias or electrode waveform, blanket TiO2 loss,
remaining Cr, exact GDS/width/pitch, sample radius/orientation, and SEMs with
scale/view metadata.

## Krueger state

Three independent, frozen Guo-to-Krueger transfers are running.  All target
outcomes remain excluded from parameter selection.  Exact supervisor PID files
and PIDs are:

| Case | PID file | PID | Snapshot prefix |
| --- | --- | ---: | ---: |
| nominal unresolved aggregate ion | `/root/krueger_guo_60s_nominal.pid` | 10284 | 9.516 s, 108.746 nm |
| all aggregate ions treated as CF2+ | `/root/krueger_guo_60s_cf2.pid` | 10285 | 9.125 s, 132.258 nm |
| all aggregate ions treated as CF3+ | `/root/krueger_guo_60s_cf3.pid` | 10286 | 8.953 s, 140.104 nm |

Each supervisor is alive and has a healthy feature worker.  Checkpoints and
`audit.json` files were updating at the snapshot.  Each worker advances about
0.015625 simulated seconds per accepted step at roughly 37-42 wall seconds per
step.  This is why the run is long.

Do not linearly extrapolate these values.  Mask opening, surface state,
transport, and topology evolve.  Do not choose the composition case closest to
825 nm after reveal.  Grade all three frozen cases and preserve aggregate-ion
composition as uncertainty.

The current science remains:

- published-input common-engine result about 346.833 nm versus 825 nm;
- old website 790-811 nm result retracted because two defects canceled;
- Guo/Kwon planar transfer effective yield 2.613 versus the 2.521 implied by
  Krueger depth, a promising but not sufficient approximately 5% agreement;
- independent full-feature forecasts still needed because early prefixes and
  planar arithmetic are not a profile prediction.

## Bosch machine-waveform depth gate

### Data firewall and frozen gates

The official Sayyed et al. SPTS Omega i2L DSi Rapier dataset contains 96 process
records with measured 5 Hz machine waveforms and 89-point wafer outcomes.  The
target broker created a calibration-only asset with 75 measured wafers and
6,675 rows.  Thirteen measured heldout wafers and 1,157 outcome rows are absent
from the fitting asset.  The heldout CSV has not been opened by the calibration
path.

Frozen heldout gates are:

- Si wafer-mean MAE `<=1 um`;
- Si wafer-mean MAPE `<=3%`;
- pointwise Si RMSE `<=1.5 um`;
- normalized spatial-shape RMSE `<=2%`;
- oxide mean-loss MAE `<=0.08 um`;
- selectivity MAPE `<=12%`;
- beat a leave-one-lot-out global-mean depth baseline;
- beat a leave-one-lot-out mean 2-D map baseline.

No heldout reveal is allowed until a versioned model receipt includes hashes,
parameter bounds/values, leave-one-lot-out scores, refinement, and a frozen
prediction file.

### Implemented unchanged-core path

The path is:

```text
measured 5 Hz power/pressure/gas/bias waveforms
  -> exact reduced F / film-precursor / ion population dynamics
  -> conservative positive cylindrical (r,phi,z) wafer transport
  -> species-resolved point fluxes at all 89 measured coordinates
  -> exact fused Belen-Si / La-Magna-film-and-oxide recurrence
  -> Si depth, oxide loss, selectivity, and spatial map
```

The finite-volume lift has periodic azimuth, Robin lower/upper/side losses, an
M-matrix solve, exact inventory ledgers, positive source moments, analytic JVPs,
and exact axisymmetric recovery at zero harmonics.  The surface evaluator has
machine-precision parity against the rich canonical mechanisms.  The sub-grid
wafer-edge focus conserves area-integrated ion current exactly.

### Model-form discoveries

The old shared axisymmetric model has a mathematical lower bound of 2.128943%
normalized map error.  Even a per-wafer axisymmetric oracle has median 2.053148%
error, both above the frozen 2% gate.  A shared unconstrained 2-D map reaches
0.621753%.  Axisymmetry is therefore falsified by the calibration data; the
cylindrical extension was required rather than optional.

Species-specific radial source moments did not materially improve the result.
An extreme-diffusion/edge-focus basin lowered map error but collapsed mean
depth to about 15 um and was rejected as a nonphysical trade.  At the
depth-valid diffusion coefficients `(2,1,10) m2/s`, a current-conserving ion
edge focus plus first azimuthal harmonic was sufficient to pass every absolute
calibration gate.

The current strongest depth-valid shared candidate uses:

```text
reference lifetimes: F=0.412333985 s
                     film=2.99382399e-5 s
                     ion=1.53294139e-4 s
diffusion:            F/film/ion=(2,1,10) m2/s
source ring:          radius=0.12 m, width=0.02 m
ion central fraction: 0.07954333
edge focus:           amplitude=0.5
                      onset=0.0907157698 m
                      width=0.00705256433 m
first harmonic, F:    cos=-0.0356084449, sin=-0.405495811
first harmonic, film: cos=-0.273223495,  sin=-1.49999895
first harmonic, ion:  cos=0.000805789873, sin=-0.0191321177
```

Across all 75 calibration wafers, exact scores are:

| Metric | Physics candidate | Frozen gate |
| --- | ---: | ---: |
| Si mean-depth MAE | 0.389181 um | <=1 um |
| Si mean-depth MAPE | 0.89312% | <=3% |
| pointwise Si RMSE | 0.839958 um | <=1.5 um |
| normalized map RMSE | 1.647094% | <=2% |
| oxide mean-loss MAE | 0.039568 um | <=0.08 um |
| oxide mean-loss MAPE | 6.269% | diagnostic |
| selectivity MAPE | 6.1117% | <=12% |

But the leave-one-lot-out empirical baselines are stronger:

| Baseline | Error |
| --- | ---: |
| global mean depth | MAE 0.338486 um; MAPE 0.7754% |
| shared mean 2-D map | point RMSE 0.486585 um; normalized shape 0.636619% |

The physics candidate does not yet beat either baseline.  Its Si wafer-mean
correlation is about 0.041, showing that it mostly reaches the correct scale
and shared shape rather than explaining lot-to-lot drift.

A final calibration-only search added positive source harmonics through order
four.  Third-order normalized map error was 1.372983%; fourth-order was
1.372168%.  The negligible change falsifies “add more harmonics” as the next
closure.  Those exploratory coefficients were not committed and should not be
promoted.  The preregistered next mechanism is one shared wall-conditioning
law using at most three coefficients: log carbon-count plus Si and SiO2
condition indicators, acting as a bounded multiplier on effective neutral
wall-loss velocity.  It must be implemented physically through reactor state,
not as a per-lot depth offset.

The conditioning law should be evaluated with true leave-one-lot-out refits.
If it cannot beat both baselines while satisfying absolute gates and grid
refinement, retain the result as a model-form failure and do not reveal
heldout outcomes.

## Validation and generality

Existing unchanged-core evidence includes direct SF5+ to Si/SiO2 beam points
at about 5.88% MAPE, Cl2/Ar+ to Si ALE points at about 12.88% maximum error,
a Mahorowala Cl2 feature board near 19% MAPE without a formal pass gate, and
direct fluorocarbon beam agreement at feature energy with an existing 4.7%
grade.

These are real mechanism checks, but they do not yet satisfy the mission's
requirement for Krueger plus at least two additional formal heldout
reactor-to-depth/profile gates using the unchanged core.  Bosch can become one
such gate after baseline-beating calibration and sealed reveal.  Oxford can
become another when the board is complete and Freddie's exact answer key is
available.  A second independent feature chemistry should still be promoted
in case Oxford data remain incomplete.

The engine is general by deck and conserved operator, not omniscient by JSON.
A new chemistry still requires evidence for gas reactions, electron-impact
rates, wall losses, sheath/ion-energy transfer, material-specific sticking,
desorption, sputter/removal yields, polymer formation/removal, geometry, and
heldout measurements.  DFT/MD can close selected barriers and yields; it cannot
infer missing machine boundary conditions from forward power alone.

“Atomic accuracy” currently means atom/formula-unit ledgers and physically
resolved reaction accounting where data exist.  The Oxford production grid is
10 nm and the surface laws are reduced kinetics.  No atom-by-atom final-SEM
claim is warranted.  The defensible goal is deterministic, differentiable,
uncertainty-propagating knobs-to-profile prediction at the experimental error
floor.

## Exact takeover sequence

1. **Reproduce the Oxford v5 cleanup failure locally or on the task-owned box.**
   Use the exact job spec above and instrument ownership sets per cleanup pass.
   Do not change physics parameters or consume the SEM.
2. **Fix the topology operator, not the assertion.**  Require idempotence or a
   monotone termination argument, preserve ledgers, add the exact regression,
   and run focused plus full tests.
3. **Resume the Oxford board content-addressably.**  Reuse all 39 valid v5
   caches.  Require exactly 56 v5 caches, remote native `--check`, hash-verified
   retrieval, local `--check`, and then commit/push the final blind envelope.
4. **Formalize the Bosch calibration candidate.**  Put the current parameter
   set, exact scores, baselines, and hashes in a version-controlled audit
   script and receipt; do not rely on console heredocs.
5. **Implement only the preregistered conditioning closure.**  Use one shared
   physical wall-state law, true leave-one-lot-out fitting, bounded parameters,
   and no per-wafer/per-lot depth offsets.
6. **Run refinement before reveal.**  Compare the fitting grid with the frozen
   certification grid and require observable changes below the preregistered
   tolerance.  Seal code/data/parameter/prediction hashes.
7. **Reveal Bosch only after it beats both baselines.**  Score every frozen
   gate once.  Do not tune after reveal.
8. **Leave all three Krueger cases frozen and running.**  Retrieve and
   hash-verify checkpoints/results only after completion or a declared physical
   terminal event.  Grade every case, including domain violations.
9. **Promote another heldout feature chemistry.**  Surface points alone do not
   complete the mission.
10. **Destroy only instance `48177892` after verified retrieval.**  Confirm the
    task artifacts exist locally first.  Never touch other account instances.

## Validation status

Focused tests for the latest cylindrical, wafer-depth, target-firewall, and
edge-focus paths passed.  The latest edge-focus checkpoint had 15 focused tests
passing.  The last repository-wide run before the newest Bosch commits was
2161 passed and 7 skipped; a new full-suite run is still required before a
heldout seal.

## Claims that must not be made

- Oxford has not matched Freddie's SEM; it is still a blind conditional board.
- Krueger 825 nm has not been reproduced by the verified common model.
- The old website 790-811 nm match is retracted and must not be restored.
- Bosch has not passed heldout; only calibration outcomes have been used.
- Higher Fourier order did not solve Bosch generalization.
- The model is conservative and atom-balanced where mechanisms are resolved,
  but it is not an atom-by-atom exact-geometry oracle.

## One-sentence handoff

Fix one reproducible Oxford sub-grid ownership cycle and reuse 39 valid caches;
keep the three healthy Krueger no-fit forecasts frozen; and close Bosch
lot-conditioning with one preregistered physical wall-state law before sealing
and revealing the heldout machine-to-depth test.
