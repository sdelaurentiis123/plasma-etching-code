# Codex takeover report — exact live state at 2026-08-21 15:15 UTC

Repository: `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code`

Branch: `codex/validation-first-multiphysics`

Committed HEAD and remote at inspection: `07daf1d` (`Tabulate exact Bosch wall
ion response`)

This is the authoritative continuation report for the current session. It
supersedes the live-state sections of:

- `HANDOFF_MOVING_CR_BOARD_2026-08-20.md`;
- `CODEX_TAKEOVER_REPORT_2026-08-21.md`;
- earlier Oxford v3-v5 progress reports.

Those files remain useful forensic records. Their process counts and statements
that the Oxford board is running are no longer current.

## Executive answer: what the fuck is happening

There are three separate scientific tracks, and two different Oxford failures
have been conflated in the older narration.

1. **The alarming 55/56 Oxford spin is historical and fixed.** It occurred in
   an older revision. It was deterministically reproduced, localized to
   non-idempotent TiO2/Cr material-owner cleanup, fixed in `980800a`, and pinned
   by a full successful acceptance trajectory in `fefd970`. That old spin is
   not consuming compute now.
2. **The clean Oxford v6 board later encountered a new failure after 23/56
   unique v6 caches.** The parent has exited. The new failure is in sparse
   surface-state remapping for a different cell; it is not the old infinite
   loop and is not an out-of-memory event. No Oxford workers are alive, so the
   machine is not secretly burning compute on this board.
3. **The three Krueger no-fit forecasts are healthy and still running.** Their
   supervisors sleep while one exact child per case uses a CPU core. They are
   around 24-26 of 60 physical seconds and must not be extrapolated linearly.
4. **Bosch v8 has made real progress locally.** The conservative differentiable
   ion-transmission operator and exact wall/ion response tensor are committed.
   An uncommitted calibration audit has selected the first preregistered basis
   that beats all pre-replay calibration and whole-lot baselines, and the
   independent interpolation audit passes. Exact reactor replay, refinement,
   prediction hashing, and the heldout reveal have not happened.

Nothing has been deleted or reset. The unfinished Bosch files and the user's
pre-existing paths are still present. The scientific mission is not complete:
Freddie's SEM has not been matched, Krueger's 825 nm has not been matched, and
Bosch's chronological heldout result remains sealed.

## Active mission and non-negotiable standard

The active goal is:

> Achieve defensible absolute-depth prediction rather than a one-target fit:
> identify and implement missing reactor and surface closures, validate them
> against independent measurements, and demonstrate heldout depth/profile
> accuracy for Krueger C4F6/Ar/O2 SiO2 plus at least two additional chemistries
> using the unchanged multiphysics core.

Completion requires numerical convergence and conservation, independently
constrained reactor boundaries and surface laws, a frozen prediction before
target reveal, heldout profile/depth metrics that beat preregistered empirical
baselines, and cross-chemistry transfer. Atomic mesh resolution is not atomic
predictive accuracy.

## Repository state and safety boundaries

At inspection, the committed branch and remote both point to `07daf1d`.

Campaign-owned work that is intentionally **uncommitted** because its exact
replay gate is unfinished:

- modified
  `scripts/audit_bosch_wafer_boundary_map_calibration.py`;
- new
  `results/curated/zenodo_bosch_wafer_boundary_map_depth_extension_v8/full_calibration_capacity.json`;
- new
  `results/curated/zenodo_bosch_wafer_boundary_map_depth_extension_v8/calibration_fit.json`;
- new
  `results/curated/zenodo_bosch_wafer_boundary_map_depth_extension_v8/interpolation_validation.json`.

Preserve these byte-for-byte until the v8 replay is checked or deliberately
recomputed. The script compiles. It now contains `build_exact_replay` and the
`--exact-replay` CLI, but no `exact_replay.json` exists and the path has not yet
been executed or scientifically certified.

Pre-existing user work that must never be edited, staged, deleted, cleaned,
reset, or absorbed into a campaign commit:

- `results/curated/mixed_layer_feature_v1/ml20-bonds-12s.log`;
- `results/curated/mouth_equilibrium_probe_dx/`;
- `scratch_ignore_calc.py`.

Do not use destructive Git commands. Stage exact paths only.

## Paid compute

Only Vast instance `48177892` belongs to this campaign.

- SSH: `ssh -p 17892 root@ssh6.vast.ai`
- Oxford tree: `/root/petch-4b656fd`
- Krueger tree: `/root/petch-d852a1f`
- environment: `/root/petch-venv`
- hardware: RTX 3090, 24 GiB
- last recorded rate: approximately `$0.2011/hour`

Do not touch other instances. Do not destroy `48177892` until the 23 Oxford v6
caches and failure log plus all terminal Krueger outputs have been copied
locally, hash-verified, audited, committed, and pushed.

Use PID files and `ps -p`; do not use `pgrep -f`, because the historical probe
matched its own command and falsely reported liveness.

## Track A — Oxford NPG80 TiO2/Cr blind board

### Frozen physical condition

- Oxford PlasmaPro NPG80 RIE;
- 55/5/1 sccm CHF3/SF6/O2;
- 30 mTorr;
- 150 W forward table RF;
- 20 C table temperature;
- 1,200 s process;
- 700 nm ALD TiO2 on fused silica;
- 45 nm Cr hard mask;
- square pillars, approximately 400 nm pitch;
- nominal width board: 80, 120, 160, 200, 240, 280, and 320 nm.

The prediction remains blind to Freddie's target SEM. Missing experimental
metadata includes exact GDS mapping, sample radius/orientation, achieved
self-bias/electrode waveform, blanket TiO2 loss, and remaining Cr thickness.
The solver propagates named uncertainty cases rather than choosing values to
match a profile.

### Resolved historical v5 failure

The old 55/56 story was real but stale. The failure chain included a
marching-cubes sliver, undefined mask extinction, unsupported material-owner
nodes, and finally non-idempotent Cr/TiO2 regional cleanup. The decisive defect
was fixed by freezing the locally repaired owner and projecting regional
fields against the exact authoritative union. The union remains bitwise
unchanged and the owner update is idempotent. The former failing full
trajectory completed under v6.

Committed v6 acceptance cell:

- width 200 nm;
- selectivity 18.016664;
- high-energy, zero-tail scenario;
- depth 679.407775 nm;
- middle CD 196.542208 nm;
- bottom CD 203.951939 nm;
- sidewall angle 81.4639 degrees;
- bow 94.566238 nm;
- particle transport residual zero;
- remap residual `9.37e-16`.

That is a conditional simulation result, not an SEM match.

### New live-state finding: v6 stopped at 23/56

At 2026-08-21 15:11 UTC:

- PID file still contains parent PID `27994`, but `ps -p 27994` returns no
  process;
- no Oxford child workers remain;
- exactly 23 unique v6 trajectory caches exist remotely;
- `/root/zhu_v6_board_980800a.log` ends in a deterministic worker exception;
- GPU use is zero because Oxford is stopped and Krueger is CPU-only.

Failing preregistered cell:

```text
width_nm       = 120
selectivity    = 18.016664028610727
ion energy     = 296.20650777233976 eV (high)
angular tail   = 0.0
duration       = 1200 s
revision       = two-material-moving-tio2-cr-owner-projection-v6
```

Exact exception path:

```text
advance_feature_step_3d
  -> _remap_surface_state_with_indexed_transfer
  -> build_surface_transfer_3d
  -> SurfaceTransferWeights3D.__post_init__
  -> ValueError: invalid sparse surface-transfer weights
```

This is a fail-closed invariant, not a silent corrupt result. The constructor
rejects one of several possible conditions: malformed CSR row offsets,
non-finite/negative weights or distances, invalid predecessor indices,
non-unit row sums at an `8e-16` tolerance, invalid scalar diagnostics, or
malformed metadata. The current aggregate message does not identify which
predicate failed. Therefore the cause is not yet proven and must not be called
physics or roundoff by assumption.

### Correct Oxford recovery sequence

1. Copy the remote v6 log and all 23 caches to a temporary local quarantine;
   hash them before any restart.
2. Reproduce only the failed 120 nm cell at the same commit and revision.
3. Instrument `SurfaceTransferWeights3D.__post_init__` so it reports the exact
   violated predicate and the step/checkpoint fingerprint. Do not relax all
   validation tolerances just to make the cell pass.
4. Add a deterministic regression fixture at the failing remap.
5. If the cause is nonphysical numerical bookkeeping, fix it narrowly and
   prove conservation, material locality, CPU/CUDA agreement, determinism, and
   mesh/refinement stability. If it is a real feature-topology limit, report
   the physical termination rather than forcing continuation.
6. Decide cache compatibility from exact equivalence evidence. If the fix can
   change trajectories, bump the model revision and recompute all 56; do not
   mix old and new revisions. If it is diagnostics-only or bitwise-equivalent,
   document why the 23 caches remain valid.
7. Finish all 56, verify the preregistered spec set, run check mode, render the
   blind atlas, run focused and broad suites, then commit and push the complete
   frozen board before looking at the target SEM.

Do not restart the old parent blindly. The v6 parent is dead and the exception
will recur unless diagnosed.

## Track B — Krueger C4F6/Ar/O2 SiO2 depth

Krueger's reported target is approximately 825 nm. The old website result of
roughly 790-811 nm is retracted because it came from cancelling bugs. The
pure-published-input feature result is approximately 346.833 nm. Direct beam
data still imply that reaching 825 nm with the paper's printed ion flux would
require a surface yield above the measured physical ceiling.

The current three runs transfer independently fixed Guo beam/planar surface
laws under three declared positive-ion composition hypotheses. They are not
post-target fits:

- `nominal_unresolved`;
- `all_cf2`;
- `all_cf3`.

At 2026-08-21 15:11 UTC, the three supervisor PIDs `10284`, `10285`, and
`10286` are alive. Their current child PIDs are alive at approximately 100% of
one CPU core each. Latest complete prefixes were:

| case | accepted step | process time | depth | mask opening |
|---|---:|---:|---:|---:|
| nominal unresolved | 1661 | 25.953125 s | 222.500871 nm | 19.024338 nm |
| all CF2 | 1584 | 24.750000 s | 249.345094 nm | 21.977417 nm |
| all CF3 | 1564 | 24.437500 s | 254.123894 nm | 23.466897 nm |

The target time is 60 s. These are nonlinear, geometry-evolving prefixes; do
not linearly extrapolate them or choose whichever terminal case lands closest
to 825 nm. All three must finish and be reported as the preregistered envelope.

The live receipts also explicitly say
`parameter_evidence_supports_prediction = false` and
`within_declared_scope = false`. Reasons include neutral species outside Guo's
printed source list, ion energy beyond the <=370 eV regression board, a
declared angular-law typesetting repair, and unresolved aggregate ion
composition. The runs are valuable transfer/falsification tests, not yet a
certified absolute prediction.

Krueger depth is **not matched** at this snapshot.

## Track C — Bosch SF6/C4F8 reactor-to-depth transfer

### What is committed

The chronological heldout outcome remains unopened. The following v8 work is
committed and pushed:

- `3b2ef76`: preregistered positive, smooth, complete-real-Zernike
  ion-transmission law with static orders 1-10 and optional order-2 Vpp mode;
- `2fed62d`: conservative current-normalized wafer-boundary implementation;
- `26931da`: differentiable accelerated factors and analytic Jacobian;
- `07daf1d`: exact wall/ion response tensor and equivalence tests.

The operator multiplies only the positive-ion wafer transmission, is strictly
positive, and is normalized against the baseline cylindrical ion-current grid
so total positive-ion current is unchanged. It does not alter neutrals,
energies, wall memory, or surface laws.

The exact response tensor covers 13 wall nodes x 13 local-ion nodes x 75 wafers
x 89 positions. Its SHA-256 is:

```text
d27104f3bfdf4e72279fffcb1ed78ecb3ca600b8ba95edd3a9b28e39da5d3b4e
```

At unity local-ion factor it matches the independently generated v7 response
table to about `4.74e-20 m` silicon and `6.35e-22 m` oxide.

### What the uncommitted v8 calibration found

Twenty preregistered static/dynamic basis candidates were fit on the 75
calibration wafers. The selection rule was frozen as the lowest coefficient
count passing every pre-replay gate. The selected candidate is:

- static complete-real-Zernike maximum order 9;
- dynamic complete-real-Zernike maximum order 2 driven only by standardized
  C4F8 platen Vpp RMS;
- 59 total coefficients;
- full-fit Jacobian rank 59/59;
- full-fit condition number 9.954;
- maximum pairwise parameter correlation 0.714;
- no coefficient bound contact;
- local ion factors remain near unity rather than hiding a scalar depth fit.

Full calibration metrics:

- silicon mean MAE: `0.227620 um`;
- silicon mean MAPE: `0.520465%`;
- silicon point RMSE: `0.389184 um`;
- normalized shape RMSE: `0.619833%`;
- oxide mean MAE: `0.039242 um`;
- selectivity MAPE: `6.017391%`.

Whole-lot leave-one-out metrics:

- silicon mean MAE: `0.269005 um`, beating the `0.338486 um` global-depth
  baseline;
- silicon point RMSE: `0.430195 um`, beating the `0.486585 um` mean-map
  baseline;
- normalized shape RMSE: `0.635023%`, narrowly beating the `0.636619%`
  mean-map baseline;
- within-lot slope MAE: `0.079842 um/wafer`, better than v7;
- median fold coefficient standard deviation: `0.000196`;
- maximum fold coefficient range: `0.003866`;
- maximum fold Jacobian condition number: `10.062`.

All preregistered pre-replay gates pass. The shape margin is only about
`0.001596` percentage points, so exact replay is decisive and cannot be
skipped.

Independent tensor-midpoint interpolation validation used 144 combinations
across the full factor-four wall and local-ion domain. It passes: the worst
error is `0.036988` of a frozen gate, below the frozen `0.05` limit.

### Bosch gates still closed

- No exact selected full-calibration replay receipt exists.
- No exact whole-lot replay receipt exists.
- No cylindrical-grid/refinement receipt exists.
- No final code/data/fold/prediction hash bundle exists.
- No chronological heldout prediction has been written.
- The heldout outcome has not been read.
- The candidate is not eligible for a prediction seal.

The immediate local step is to review and test the newly added exact-replay
path, run it without opening the heldout, and accept the candidate only if the
exact solver and refinement preserve every gate. If the narrow shape advantage
disappears, v8 fails honestly.

## Generalization status

The architecture is general at the solver level: chemistry decks, conserved
reactors, species/energy/angle-resolved wafer boundaries, charging, surface
state, and moving multi-material geometry are modular and deterministic.
Uncertainty and geometry boards are embarrassingly parallel.

That does not mean an arbitrary chemistry becomes predictive by editing a
JSON. Each new system still needs electron-impact and gas chemistry data,
wall/sticking data, IEDF/IAD closure, species- and energy-resolved surface
yields, passivation/polymer kinetics, mask response, and independent reactor
and surface validation.

The three current systems are:

1. Krueger C4F6/Ar/O2 on SiO2, with direct beam-transfer constraints but an
   unresolved published boundary;
2. Oxford CHF3/SF6/O2 on TiO2/Cr, with a reactor-to-feature conditional board
   that is currently numerically blocked before blind SEM scoring;
3. SPTS SF6/C4F8 Bosch silicon/oxide, with strong calibration/LOLO progress and
   an intact chronological heldout firewall.

This is meaningful cross-system science, but the requested three-chemistry
heldout accuracy claim has not yet been earned.

## Exact takeover order

1. Confirm this branch and preserve all seven uncommitted paths listed above.
2. Copy and hash the Oxford v6 log and 23 caches from instance `48177892`.
3. Reproduce and instrument the exact Oxford 120 nm remap failure; land a
   minimal regression and only then decide whether a revision bump/full rerun
   is required.
4. Continue monitoring the three Krueger PID files and exact children. Do not
   restart healthy runs or extrapolate prefixes.
5. Review the uncommitted Bosch exact-replay implementation, add focused tests,
   and run exact replay/refinement without reading heldout outcomes.
6. If Bosch replay passes, commit and push the complete v8 calibration receipt;
   then write and hash the chronological prediction before any outcome read.
7. Finish and certify the complete Oxford board before target SEM scoring.
8. Retrieve and certify all three terminal Krueger cases as one declared
   envelope.
9. Run focused tests and the full suite after the numerical work lands.
10. Destroy Vast instance `48177892` only after all remote-only evidence is
    local, verified, committed, and pushed.

## Claims allowed now

Allowed:

- the old Oxford 55/56 spin was diagnosed and fixed;
- the subsequent v6 board produced 23 caches and failed closed on a distinct
  sparse surface-transfer invariant;
- the three Krueger no-fit forecasts are alive and advancing;
- Bosch v8's conservative ion-boundary model is committed and its selected
  pre-replay candidate beats all frozen calibration and whole-lot baselines;
- Bosch interpolation validation passes and its heldout remains sealed;
- the core is deterministic, differentiable where intended, and modular
  across reactor, boundary, surface, and feature layers.

Not allowed:

- "Oxford is still running";
- "the Oxford board is complete";
- "we matched Freddie's SEM";
- "we matched Krueger's 825 nm";
- "Bosch heldout passed";
- "the simulator is atomically accurate";
- the retracted website claim of 790-811 nm Krueger depth;
- selecting a physical uncertainty case because it is closest to a target.

## Bottom line

The project is not destroyed and the old 55/56 failure is not the current
problem. Oxford v6 stopped cleanly on a new, diagnosable remap invariant after
23 cells. Krueger is healthy but incomplete. Bosch v8 is the strongest current
scientific advance: it turns the previously localized machine fingerprint into
a positive, current-conserving, differentiable wafer-boundary law that now
beats the whole-lot point and shape baselines before replay. But the exact
replay and heldout firewalls are doing their job, so no final depth/profile
breakthrough should be claimed yet.
