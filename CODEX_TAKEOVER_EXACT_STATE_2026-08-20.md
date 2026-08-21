# Codex takeover — exact absolute-depth campaign state

Snapshot: 2026-08-20 23:25 EDT / 2026-08-21 03:25 UTC

Repository: `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code`

Branch: `codex/validation-first-multiphysics`

Pushed scientific HEAD before this report: `2a03b3a`

Read this file first. It supersedes the live-state and next-action sections of
all earlier 2026-08-20 takeover and moving-Cr handoffs. Those files remain
for forensic history only.

## The short, honest answer

The campaign is alive, version controlled, and scientifically productive, but
it has **not yet matched Krueger depth, completed Freddie's blind profile board,
or earned the sealed Bosch heldout reveal**.

What looked like one mysterious Oxford CPU spin was actually a sequence of
separable numerical defects exposed as the board became more physical. Three
were fixed with tests: a mesh-area flux singularity, the Cr-extinction
lifecycle, and invalid recombination of independently redistanced material
fields. The current fourth failure is deterministic, bounded, and checkpointed:
the subcell material-owner cleanup does not reach a fixed point in one exact
trajectory. It is an engine-topology problem, not an experimental mismatch.

Meanwhile, the Bosch reactor-to-depth stack achieved all frozen **absolute
calibration** tolerances with unchanged literature surface laws, but it still
loses to simple leave-one-lot-out empirical baselines. The first static wall
conditioning closure was therefore refused rather than used to open heldout
data. A new calibration-only residual audit found a very strong within-lot
wafer-sequence signature, pointing to a missing dynamic chamber/tool memory
state rather than another static lot label.

Three frozen no-target-fit Krueger feature forecasts are still advancing on
the paid box. Their prefixes are healthy but must not be extrapolated or called
final depths.

## Mission and acceptance condition

The active goal is:

> Achieve defensible absolute-depth prediction rather than a one-target fit:
> identify and implement the missing reactor/surface closures, validate them
> on independent measurements, and demonstrate heldout depth/profile accuracy
> across Krueger C4F6/Ar/O2 SiO2 plus at least two additional chemistries using
> the unchanged multiphysics core.

The goal is not complete. A calibration match, a conditional uncertainty
envelope, or one target-specific fit is not sufficient.

## Repository safety

Use only the repository and branch at the top of this file. Do not switch to
the older multiphysics repository/branch mentioned in historic conversations.

The following untracked paths are pre-existing user work. Do not stage, edit,
delete, clean, reset, or absorb them:

- `results/curated/mixed_layer_feature_v1/ml20-bonds-12s.log`;
- `results/curated/mouth_equilibrium_probe_dx/`;
- `scratch_ignore_calc.py`.

At this snapshot, all task-owned work through the Bosch residual discovery is
committed and pushed. The protected paths above are the only worktree dirt.

## What was actually completed after the previous handoff

Pushed commits after the earlier live report:

- `f10ecdf` preserves an authoritative material union and owner map instead of
  independently redistancing the combined union after owner cleanup;
- `01abf02` adds an exact Oxford single-cell reproducer;
- `921e62a` freezes the first Bosch wall-conditioning law before fitting;
- `852df0b` implements that bounded wall-state law in the reduced and
  cylindrical reactors without changing ions or surface laws;
- `8ccd278` adds Oxford pre-failure checkpoints and pass-by-pass cleanup
  histories;
- `9c2c009` records the Bosch conditioning calibration and explicitly refuses
  a prediction seal;
- `2a03b3a` records the calibration-only Bosch residual-feature audit and
  preserves the heldout firewall.

Focused checks at the snapshot:

- 93 Oxford/feature tests passed after `f10ecdf` and `8ccd278`;
- 23 Bosch data/reduced-reactor tests passed after `2a03b3a`;
- the Bosch residual receipt was regenerated exactly from its inputs;
- the last recorded broad suite earlier in the Oxford campaign was 2,125
  passing tests. A new full-suite run is still required after the current
  Oxford topology fix is completed.

## Live paid compute

Only Vast instance `48177892` belongs to this task:

- SSH: `ssh -p 17892 root@ssh6.vast.ai`;
- GPU: RTX 3090;
- Oxford source tree: `/root/petch-4b656fd`;
- Krueger source tree: `/root/petch-d852a1f`;
- environment: `/root/petch-venv`;
- last recorded price: approximately `$0.2011/hour`.

Do not touch another Vast instance. Do not destroy this one until the Oxford
v5 caches, Oxford diagnostic checkpoint/log, and all Krueger terminal outputs
have been copied locally and hash-verified.

Use exact PID files for liveness. Do not use `pgrep -f`; a previous probe
matched its own command line and created a false-positive status report.

## Track A — Oxford NPG80 TiO2/Cr blind board

### Frozen physical condition

- Oxford PlasmaPro NPG80 RIE;
- 55/5/1 sccm CHF3/SF6/O2;
- 30 mTorr;
- 150 W forward table RF;
- 20 C table temperature;
- 1,200 s process time;
- 700 nm ALD TiO2 on fused silica;
- 45 nm Cr hard mask;
- blind square-pillar board, approximately 400 nm pitch;
- widths 80, 120, 160, 200, 240, 280, and 320 nm.

The reactor boundary is species-resolved: 67 plasma species, including 20
positive ions and 37 thermal neutrals. The current conditional solution has a
central positive-ion flux of about `1.457e19 m^-2 s^-1` and a powered-electrode
sheath near 296 V, giving about 299 eV for singly charged ions. Because Freddie
provided forward RF power but not the achieved electrode waveform/self-bias,
the board propagates ion-energy and angular-tail uncertainty rather than
pretending that one exact IEDF is known.

The board also propagates the experimentally unresolved TiO2:Cr selectivity.
These are named physical uncertainties, not knobs selected from the target SEM.

### Cache truth

The local trajectory directory contains:

- 55 v3 trajectories;
- 56 complete v4 trajectories;
- 39 v5 trajectories copied from the paid box at this snapshot.

The same 39 of 56 v5 trajectories remain on the paid box and are now committed
locally with this report. The v3/v4 caches are valuable forensic evidence, but
**must not be graded as the current v5 board**.

All 39 local v5 files were checked against the current preregistration-derived
job specs and contain nonempty profile lists. Seventeen current job specs are
still missing; the failed `w200/s18.017/ion_high/tail0` cell is the first one.

### Bugs fixed before the current blocker

1. Periodic strip symmetrization divided equal integrated rate by unequal
   marching-cubes triangle areas. A microscopic sliver acquired an enormous
   flux density and requested millions of CFL substeps. The corrected operator
   preserves uniform flux density and distributes rate by triangle area.
2. Complete Cr-mask loss was previously an undefined topology event. The engine
   now has a certified material-extinction lifecycle and a declared
   campaign-specific continuation rule.
3. Independently redistancing the repaired authoritative material union could
   manufacture unsupported owner nodes. The union/owner invariant is now
   explicit and regression-tested.

The third fix is a valid invariant repair, but the exact failing cell reproduced
unchanged. It was therefore not the root cause of the present failure.

### Current deterministic failure

Exact cell:

```text
width_nm = 200
selectivity = 18.016664028610727
scenario = ion_high_tail_0p0
impact_energy = 296.20650777233976 eV
tail_fraction = 0
failure elapsed = 719.862 s
accepted steps = 149
Cr centre remaining = 19.850306886455172 nm
Cr interior minimum = 0 nm
```

Guardrail:

```text
subcell material-island cleanup did not reach a fixed point;
passes=12;
coordinates=((14,12,80),(15,12,80),(15,28,80),(28,15,80),(28,27,80))
```

This is a bounded certification failure in
`_restore_unresolved_material_ownership`, not an unbounded spin. Do not raise
the 12-pass limit. The acceptable fix must make the ownership update monotone
or idempotent, preserve the material union, and conserve the existing particle
and volume ledgers.

### Captured diagnostic

Exact PID `18157` terminated on the expected guardrail and wrote the
pre-failure checkpoint. Both artifacts were copied off the paid box and are
version controlled:

```text
tests/data/zhu_v5_prefailure_8ccd278.npz
results/curated/zhu_npg80_moving_cr_profiles_v1/zhu_v5_exact_trace_8ccd278.txt
```

Their SHA-256 hashes are:

```text
95d43ec8df8372ef14df57cd28df93ac90fe5c2e72da84cceaf66965f7c50998  checkpoint
20addfad8d0c815bbed05e7d868a015baa5af4a02a54bb01ab1ba531b8e5d8fa  trace
```

The checkpoint is immediately before the failed step at elapsed
`719.8619631901854 s`, after 149 accepted steps, with next step duration
`4.831288343558282 s`.

The committed checkpoint was replayed locally on CPU after retrieval. It
reproduced the same coordinates and complete twelve-pass history in 8.8 wall
seconds, despite the different CPU/Warp backend. The failure is therefore a
small, deterministic, cross-platform regression rather than a cloud-only race
or GPU artifact.

The twelve-pass history rules out a simple two-state oscillation. Every
reported island has owner 2 (the Cr layer). Node counts evolve as:

```text
393 -> 172 -> 22 -> 1 -> 22 -> 9 -> 11 -> 4 -> 4 -> 2 -> 3 -> 5
```

The first four passes nearly erase the original unresolved set, but global
redistancing then exposes new one-cell-wide Cr fragments around the mask
perimeter. This is a non-idempotent cleanup wave: each local sign repair plus
global redistance can create the next subcell component. Merely increasing the
pass limit might let this one grid converge but would leave an unbounded,
mesh-dependent erosion operator in the physics path.

Next actions:

1. resume the committed checkpoint for a one-step local reproduction;
2. replace the repeated global-redistance cleanup with a monotone/idempotent
   ownership repair that cannot create new positive Cr nodes outside the
   authoritative pre-repair support;
3. prove the material-union, resolved-volume, and ledger postconditions in an
   exact checkpoint regression;
4. replay the checkpoint, then the full exact cell;
5. sync the fix to the box and resume the v5 parent, reusing all 39 caches;
6. copy all 56 v5 results locally, run the audit in check mode, render the
   frozen atlas, and commit the prediction package before SEM reveal.

### What can and cannot be claimed

Can claim: a deterministic knobs-to-species-to-moving-mask profile envelope is
implemented, serious numerical singularities have been found and fixed, and
the target SEM has not been used for tuning.

Cannot claim: Freddie's SEM has been predicted or matched. The current v5 board
is incomplete, and the result is conditional on unmeasured self-bias/IEDF and
TiO2/Cr surface-response uncertainty.

Freddie's eventual SEM is an answer key, not a simulator input. Highest-value
additional measurements are exact-run self-bias/electrode waveform, blanket
TiO2 loss, remaining Cr thickness, exact GDS/width/pitch, sample radius and
orientation, and SEM scale/view metadata.

## Track B — Krueger C4F6/Ar/O2 SiO2 depth

Three frozen, no-target-fit Guo-to-Krueger feature forecasts are running on
CPU. Their supervisor PID files are:

```text
/root/krueger_guo_60s_nominal.pid
/root/krueger_guo_60s_cf2.pid
/root/krueger_guo_60s_cf3.pid
```

Latest snapshot prefixes:

| Case | Simulated time | Depth | Mask opening |
| --- | ---: | ---: | ---: |
| unresolved aggregate ion | 10.531 s | 119.250 nm | 42.334 nm |
| aggregate treated as CF2+ | 10.094 s | 142.597 nm | 43.315 nm |
| aggregate treated as CF3+ | 9.875 s | 150.185 nm | 42.521 nm |

All three target 60 simulated seconds and checkpoint every 30 wall minutes.
Each accepted 0.015625 s feature step costs roughly 38–41 wall seconds. They
are slow but healthy.

Do not linearly extrapolate these prefixes. Transport, mask opening, surface
state, and topology evolve. Do not select whichever ion-composition case ends
closest to 825 nm. All three were frozen before the final result and must be
graded together as composition uncertainty.

Current defensible Krueger statements:

- the published-input common-engine result is about 346.833 nm versus the
  paper's 825 nm;
- the old website's 790–811 nm result is retracted because two implementation
  errors canceled;
- direct planar Guo/Kwon transfer gives effective yield 2.613 versus the 2.521
  implied by Krueger's depth, an encouraging approximately 5% arithmetic
  agreement;
- only the terminal full-feature forecasts can say whether this independent
  surface transfer predicts the depth and profile under feature transport.

Therefore: **Krueger depth is not matched yet.**

## Track C — SPTS Bosch machine-waveform-to-depth validation

### What exists

The official Sayyed/SPTS dataset supplies 5 Hz machine waveforms and 89-point
wafer maps. The implemented deterministic path is:

```text
measured power/pressure/gas/Vpp waveforms
  -> exact reduced F, film-precursor, and positive-ion population dynamics
  -> conservative positive cylindrical wafer transport
  -> species-resolved flux at all 89 measurement points
  -> unchanged Belen Si and La Magna film/oxide surface recurrence
  -> Si depth, oxide loss, selectivity, and spatial profile
```

The calibration broker exposes 75 measured wafers/6,675 outcome rows. Thirteen
heldout wafers/1,157 rows are absent from the fit asset. Heldout outcomes have
not been opened.

### Absolute calibration result

The static wall-conditioned candidate reaches:

| Metric | Candidate | Frozen gate | Pass |
| --- | ---: | ---: | --- |
| Si mean MAE | 0.3756 um | <=1.0 um | yes |
| Si mean MAPE | 0.8607% | <=3% | yes |
| Si point RMSE | 0.8412 um | <=1.5 um | yes |
| normalized shape RMSE | 1.6546% | <=2% | yes |
| oxide mean MAE | 0.03957 um | <=0.08 um | yes |
| selectivity MAPE | 5.994% | <=12% | yes |

The surface laws were unchanged and no per-wafer/per-lot depth offsets were
used.

### Why the seal was refused

The candidate's Si mean-depth correlation is only 0.061 and oxide correlation
is only 0.073. It does not beat the leave-one-lot-out global mean or mean-map
baselines:

| Metric | Physics | Empirical baseline | Physics wins |
| --- | ---: | ---: | --- |
| Si mean MAE | 0.3756 um | 0.3385 um | no |
| Si point RMSE | 0.8412 um | 0.4866 um | no |
| normalized shape RMSE | 1.6546% | 0.6366% | no |

The fitted wall multipliers stayed between about 0.9995 and 1.016. That means
the preregistered static conditioning categories do not explain the missing
variation. `eligible_for_prediction_seal` is false; no heldout prediction was
written; no heldout outcome was read.

### New scientific discovery

The calibration-only residual audit tested outcome-free machine summaries
against residuals while retaining a hard no-reveal flag. The Si mean residual
has:

- raw Pearson correlation 0.935 with within-lot wafer number;
- within-lot correlation 0.956 with wafer number;
- within-lot correlations of magnitude 0.85–0.89 with measured reflected
  power, pressure, and Vpp summaries.

Observed Si depth falls by roughly 0.11–0.16 um per sequential wafer in every
lot, about 1 um over a run, while the present physics predicts a nearly flat or
slightly increasing depth. An exploratory one-feature leave-one-lot-out
correction using wafer sequence reduces Si mean MAE to about 0.130 um, but this
is discovery only: the feature was selected after looking at calibration
residuals and cannot authorize a heldout reveal.

The next defensible closure is a dynamic chamber/tool state carried across
wafers—most plausibly fluorocarbon wall inventory/seasoning, with possible
thermal and matching-network observables—not a direct wafer-number depth
offset. The recurrence form and parameter bounds must be frozen before a new
fit. The next audit should first partial out wafer sequence and rank remaining
machine signals, so a measured transfer drift is not misidentified as wall
chemistry.

## Cross-chemistry status

The unchanged core has substantial component validation: exact conservative
transport tests, direct-beam surface comparisons at feature-relevant energy,
charging and moving-material ledgers, and multiple chemistry decks. But the
mission requires terminal heldout depth/profile evidence, not component count.

Current terminal board:

| Chemistry / tool | Current evidence | Terminal heldout depth/profile pass |
| --- | --- | --- |
| Krueger C4F6/Ar/O2 on SiO2 | full frozen feature forecasts running | no |
| Oxford CHF3/SF6/O2 on TiO2/Cr | blind v5 board 39/56; topology blocker | no |
| SPTS SF6/C4F8 Bosch Si/SiO2 | all absolute calibration gates pass; baselines fail | no |

No language implying “atomic accuracy” is justified yet. The correct target is
physics-constrained accuracy with explicit uncertainty and progressively
tighter intervals as independent diagnostics close the boundary and surface
state.

## Exact next sequence

1. Finish the Oxford diagnostic replay and repair the owner-cleanup operator
   without weakening its certifier.
2. Complete, retrieve, and freeze all 56 Oxford v5 trajectories before the
   SEM answer key is used.
3. Let all three Krueger forecasts reach their declared terminal state;
   retrieve and grade all cases without selecting by closeness to 825 nm.
4. Extend the Bosch residual audit to partial correlations after within-lot
   sequence detrending.
5. Freeze one dynamic chamber-memory recurrence from physical evidence, fit it
   only on calibration lots, run leave-one-lot-out physics refits and grid
   refinement, and reveal heldout only if every seal prerequisite passes.
6. Run the broad suite after the Oxford engine fix, audit all versioned
   receipts, and update the public HTML only from the final certification
   table. The stale 790–811 nm Krueger claim must not reappear.

## Takeover verdict

Nothing has been deleted, no heldout firewall has been crossed, and the current
failures are informative rather than hidden. The most important new physics is
that static recipe labels are not sufficient for Bosch: sequential tool state
dominates the remaining depth error. The most important engine fact is that
Oxford is blocked by one reproducible topology-cleanup invariant, with the
instrumentation needed to resolve it. The most important Krueger fact is that
the first honest independent-transfer feature runs are still in flight.

Continue from those facts. Do not revive any superseded depth match, fit a
target-specific scale factor, open sealed outcomes early, or grade stale model
revisions as current results.
