# Codex takeover report — exact reactor-to-depth state at 2026-08-21 16:50 UTC

Repository: `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code`

Branch: `codex/validation-first-multiphysics`

Committed local and remote HEAD at report start:
`53d35643e26563e1ca195b4e080fd00f2d9eb21a`

This is the authoritative live-state handoff. It supersedes the operational
state in `HANDOFF_MOVING_CR_BOARD_2026-08-20.md` and
`CODEX_TAKEOVER_REPORT_2026-08-21_1605UTC.md`. The older reports remain useful
for forensic detail, but their Oxford counts, diagnosis, and PIDs are stale.

## Executive answer: what the hell is happening

The program is not sitting on the old 55/56 mystery hang.

1. The old Oxford v3 55/56 CPU spin was real, was instrumented, and its chain
   of topology/remapping/mask-extinction defects was fixed.
2. A clean v6 recomputation then exposed a different fail-closed defect after
   23/56 cells: one sparse interpolation weight was
   `-2.2204460492503131e-16`, exactly one binary64 ulp below zero. The cause was
   deterministic row closure onto an arbitrarily tiny final coefficient, not
   plasma physics and not a random hang.
3. The closure algorithm is now repaired without clipping, weakening a
   validator, or hiding a conservation error. It closes through the largest
   safely positive coefficient. The exact failed production checkpoint passes
   on local arm64 CPU and remote x86_64 CUDA with the same transfer and mesh
   fingerprints.
4. Because even ulp-scale stored weights and fingerprints can change, the
   Oxford model revision was bumped to v7 and all 56 cells are being recomputed.
   No v6/v7 cache mixing is allowed.
5. The v7 run is healthy. At the 16:47 UTC inspection it had eight committed
   cells and a ninth completed remotely. Eight worker processes were each near
   100% CPU and held real CUDA allocations. This is slow deterministic feature
   evolution, not a silent deadlock.
6. Separately, the Bosch reactor-to-wafer depth track has achieved the campaign's
   first real sealed chronological heldout win: 0.239 um wafer-mean Si depth
   MAE and 0.541% MAPE across 13 unseen wafers, with 89 points per wafer.
7. Three Krueger no-fit forecasts are healthy and checkpointing. They are only
   around 26-28 of 60 physical seconds. None may be extrapolated or selected
   post hoc to match the paper's 825 nm target.

The correct headline is:

> A deterministic machine-to-depth model has passed a sealed heldout Bosch
> wafer test. Oxford's blind feature board is numerically repaired and
> recomputing from scratch. Krueger remains an unresolved published-boundary
> problem with three independent sensitivity forecasts still running.

## Safety and repository rules

Use only this repository and branch. Do not switch to the old website branch,
`reactorlab`, `petch-torchsim`, or any other stale implementation.

At report start, local HEAD and the branch remote are identical at `53d3564`.
The only untracked paths are unrelated user work and must never be edited,
staged, removed, reset, cleaned, or absorbed into a campaign commit:

- `results/curated/mixed_layer_feature_v1/ml20-bonds-12s.log`
- `results/curated/mouth_equilibrium_probe_dx/`
- `scratch_ignore_calc.py`

Do not use destructive Git commands. Stage exact campaign paths only. Commit
and push scientific checkpoints as they land.

## Paid compute

Only Vast instance `48177892` belongs to this campaign.

- SSH: `ssh -p 17892 root@ssh6.vast.ai`
- Oxford tree: `/root/petch-4b656fd`
- Krueger tree: `/root/petch-d852a1f`
- Python environment: `/root/petch-venv`
- accelerator: RTX 3090, 24 GiB
- last observed price: about `$0.2011/hour`

Do not touch any other instance. Do not destroy this one until completed
Oxford and Krueger artifacts have been retrieved, hash-verified, checked
locally, committed, and pushed.

Never use `pgrep -f` alone as liveness evidence: the probe can match itself.
Inspect the exact parent and children with `ps -p` and `ps --ppid`, and verify
recent artifact timestamps.

## Oxford/Freddie TiO2/Cr blind square-pillar board

### Frozen supplied condition

- Oxford PlasmaPro NPG80 RIE
- 55/5/1 sccm CHF3/SF6/O2
- 30 mTorr
- 150 W forward table RF
- 20 C table temperature
- 1,200 s process
- 700 nm ALD TiO2 on fused silica
- 45 nm Cr hard mask
- square-pillar board, approximately 400 nm pitch
- widths 80, 120, 160, 200, 240, 280, and 320 nm

The target SEM remains withheld from tuning. The board propagates two sourced
TiO2:Cr selectivity endpoints, low/high ion-energy cases, and core/core-plus-
tail angular cases: `7 x 2 x 2 x 2 = 56` trajectories.

This is a full conditional chain: a 67-species reactor chemistry, electron
kinetics, sheath closure, parent/daughter collision chemistry, radial wafer
transport, species-resolved boundary, and deterministic three-dimensional
moving TiO2/Cr feature evolution with charging, sputtering/chemical removal,
polymerization, surface-state remapping, mask erosion, and material extinction.

It is not yet a unique atomic-accuracy machine prediction because exact-run
self-bias and same-tool TiO2/Cr surface-rate measurements are not known. The
uncertain quantities are propagated as a frozen physics ensemble rather than
silently fitted to the SEM.

### Exact numerical failure and repair

The v6 failed checkpoint and diagnostic evidence are committed under:

`results/quarantine/zhu_npg80_v6_failure_20260821/`

Exact checkpoint SHA-256:

`1fbd52695a1df0379c67d22d9743d8e74bee29f2322142fd41737cafd7881dd6`

Exact diagnostic-log SHA-256:

`6683045c63d1efd0978a46d8940325043d758bb8c5a5aaac446c5bf092fc9a16`

Failure state:

- width 120 nm
- selectivity 18.016664028610727
- high ion energy, zero angular tail
- 133 accepted steps
- elapsed physical time 642.5613496932526 s
- next time step 4.831288343558282 s
- one negative CSR weight at entry 24516
- minimum weight `-2.2204460492503131e-16`

Root cause: inverse-distance weights were normalized and the final coefficient
was assigned `1 - sum(previous)`. In this row the true final physical tail was
about `1e-27`, while the preceding binary64 sum rounded to
`1 + 2.2204460492503131e-16`.

Commit `84a786e` introduced `_normalize_nonnegative_weights`, which applies the
closure residual through the largest positive coefficient. The validator was
not relaxed. Unit witnesses, dynamic-range fuzz tests, the exact production
checkpoint, locality, nonnegativity, intensive identity, extensive material
conservation, and CPU/CUDA behavior are pinned.

The former failed checkpoint now advances with:

- transfer fingerprint
  `64dafad5621689e785800294ec3408842939d24095c7c85d23c54e1abf498d96`
- next-mesh fingerprint
  `950c8f9a1176e4184020ef9f6a72d9e4ff6c61d201ba58efa4351a6478264bac`
- unresolved-node reassignment count: zero
- maximum remap residual: `1.371e-16` local CPU,
  `2.742e-16` remote CUDA

Focused tests passed 149 cases. The complete repository suite passed:

`2229 passed, 7 skipped in 1232.94 s`

The model revision is now:

`two-material-moving-tio2-cr-positive-row-closure-v7`

All v6 caches are invalid for the active board.

### Exact live v7 state

Launch implementation commit: `84a786e`

Source hashes deployed to the remote matched local byte for byte:

- `src/petch/surface_transfer_3d.py`:
  `c634011c047697615b0ec38cc351d409bd403ba0687ffa6e9ec399ade45f4276`
- `scripts/audit_zhu_npg80_moving_cr_profiles.py`:
  `7d7cb39fc854421a7b8a45223bbffe5216d922cd71a29479395ee7ec34436973`
- `scripts/reproduce_zhu_npg80_moving_cr_cell.py`:
  `02277fbc38e63b4267374ab740d873394ac9a728cbb2d0eb3ad3150bffe24281`

Processes observed at 16:47 UTC:

- wrapper shell: PID 33307
- board parent: PID 33309
- resource tracker: PID 33424
- eight workers: PIDs 33425, 33428, 33429, 33430, 33431, 33432, 33433,
  and 33434

All eight workers were runnable at about 100% CPU, used roughly 0.9-1.45 GiB
resident memory each, and held 322-354 MiB CUDA allocations. There was ample
host and device memory. The wrapper shell remains alive because the original
SSH launch command backgrounds the Python process inside a shell that also
contained follow-up monitoring commands. This is an operational oddity, not
the worker-liveness signal and not a reason to kill PID 33307 while the board
is healthy.

Remote log:

`/root/zhu_v7_board_84a786e.log`

The log only contains Warp initialization because workers write cache files
rather than per-cell progress lines. Liveness must be checked by worker CPU,
GPU allocations, and new cache timestamps.

At 16:43 UTC, all eight width-80 cells were locally committed in `53d3564`.
At 16:47:36 UTC, the ninth v7 cache landed remotely:

`w120_s14.000_ion_low_tail_0p0_9713fd1a430912a4.json`

Do not call the board 9/56 locally until that exact file is copied and hashed.
The local version-controlled state is 8/56; the remote execution state was
9/56 at the inspection instant.

The eight committed width-80 trajectories have zero exact particle-balance
error and maximum surface-state remap residuals of `9.21e-16` to `1.42e-15`.
Their final conditional etched depths span approximately 683.71-684.09 nm.
That is a pre-SEM physics prediction, not an observed match.

### Oxford continuation procedure

1. Leave the current parent and eight workers alone while CPU/GPU usage and
   output timestamps move.
2. Retrieve new v7 JSON files by exact filename. The trajectory directory also
   contains historical revisions; do not copy or count by an indiscriminate
   wildcard.
3. Confirm each retrieved JSON declares the v7 model revision and validate its
   embedded fingerprints, conservation ledgers, and cache key.
4. Commit and push recoverable batches, updating the live execution receipt.
5. When all 56 exist, allow the remote `--write` driver to assemble the final
   audit, copy every active artifact, then run the production audit locally in
   `--check` mode.
6. Render and inspect the full atlas. Run focused tests and the repository-wide
   suite.
7. Commit and push the frozen board before opening Freddie's exact target SEM.
8. Only then digitize and score depth, top/middle/bottom CD, sidewall angle,
   bow, Cr survival, missing/collapsed pillars, and spatial position. Never
   select the best-looking ensemble corner after reveal.

The highest-value additional measurements from Freddie remain exact-run DC
self-bias, blanket TiO2 loss, remaining Cr thickness, GDS/width mapping, sample
radius/orientation, and exact cross-section/top-down SEMs. These refine and
test identifiable submodels; they are not permission to fit the final contour.

## Krueger C4F6/Ar/O2 SiO2 trench

The historical paper target is about 825 nm at 60 s. The old website's
790-811 nm apparent match is retracted because it came from canceling bugs.
The corrected published-input endpoint is about 346.833 nm. Direct-beam
evidence indicates that 825 nm at the printed aggregate ion flux requires a
removal yield beyond the supported physical ceiling.

Three frozen no-fit Guo-transfer forecasts are running:

- nominal unresolved aggregate ion: supervisor PID 10284
- aggregate ion declared all CF2+: supervisor PID 10285
- aggregate ion declared all CF3+: supervisor PID 10286

At 16:47 UTC each supervisor had one live worker near 100% CPU and fresh
`checkpoint.npz`/`audit.json` files. Approximate latest accepted states were:

| Case | Physical time | Etch depth | Mask opening |
|---|---:|---:|---:|
| nominal unresolved | 27.88 s | 230.26 nm | 17.46 nm |
| all CF2+ | 26.56 s | 258.25 nm | 22.16 nm |
| all CF3+ | 26.27 s | 263.72 nm | 20.49 nm |

The final audit write was still advancing while these values were read, so
the table is a timestamped progress snapshot, not an endpoint. Do not linearly
extrapolate evolving topology to 60 s.

These receipts intentionally declare
`parameter_evidence_supports_prediction=false` and
`within_declared_scope=false`. The unresolved aggregate ion composition,
neutral channels outside Guo's printed source list, energy support, angular-law
repair, and direct mass-selected-beam misses make the cases sensitivity and
falsification calculations—not certified species-resolved predictions.

A defensible Krueger closure still requires either species-resolved wafer
flux/IEDF/IADF diagnostics, a separately validated C4F6 reactor model plus
blanket checks, or same-tool blanket SiO2 and mask-loss measurements. Adjusting
a yield until the trench lands at 825 nm is prohibited.

## Bosch heldout machine-to-depth breakthrough

This result is complete, committed, and pushed in `860fbe1`.

Prediction SHA-256 sealed before numeric reveal:

`56ed2429832fe77280762fbca86cb6ffa4de3fd9687aa84f3b5cfd4ca99a3b1a`

Thirteen unseen chronological wafers have 89 measured points each. Seven other
heldout process records lack source outcomes and are reported missing, never
imputed. Heldout results:

| Metric | Result |
|---|---:|
| Si wafer-mean depth MAE | 0.238906 um |
| Si wafer-mean depth MAPE | 0.541245% |
| Si 89-point RMSE | 0.329751 um |
| normalized radial-shape RMSE | 0.261855% |
| oxide wafer-mean loss MAE | 0.023908 um |
| selectivity MAPE | 3.161040% |
| Si wafer-mean correlation | 0.873410 |

Every preregistered absolute gate and empirical-baseline gate passes. Bootstrap
improvement intervals are positive for mean and pointwise depth. The full-
sample radial-shape score beats the mean-map baseline, but its 95% bootstrap
improvement interval crosses zero; retain that qualification.

The result validates deterministic reactor/wall-memory/wafer-boundary transfer
for absolute Si and oxide depth, radial shape, selectivity, and drift on one
Bosch tool/process family. It does not validate feature sidewalls, charging,
ARDE, scallops, or arbitrary tools and chemistries.

## Generalization and claim boundary

The common architecture is generalizable by evidenced mechanism/deck. It is
not automatically chemistry universal.

Reusable core:

- deterministic particle and molecular transport
- conserved species-resolved reactor/wafer boundaries
- differentiable reduced reactor, sheath, and radial transport operators
- self-consistent feature charging and ion deflection
- moving multi-material three-dimensional interfaces
- intensive/extensive surface-state remapping and material ledgers
- polymer deposition/crosslink/breakage machinery
- mask erosion/extinction and topology handling
- target firewalls, hashing, preregistration, replay, and heldout scoring

Per-chemistry/tool evidence still required:

- electron-impact cross sections and branching over the actual EEDF
- gas-phase and wall reaction rates
- wall sticking/recombination/desorption and chamber history
- sheath/electrode transfer or measured waveforms
- species/energy/angle/coverage/temperature-resolved surface laws
- mask and substrate response
- tool geometry and independently measured blanket/boundary observables

Atomic mesh spacing is not atomic predictive accuracy. The current system has
not demonstrated atomic-level accuracy for Oxford or Krueger. It has
demonstrated a heldout sub-1% absolute wafer-depth result on Bosch and strong
conservation/direct-beam/numerical receipts elsewhere.

## Exact next actions

1. Monitor Oxford by exact PIDs and output timestamps; retrieve and commit v7
   batches without mixing historical caches.
2. Complete all 56 Oxford cells, run local check mode, render the atlas, run
   focused and full suites, and freeze/push before SEM reveal.
3. Monitor all three Krueger supervisors to physical termination; retrieve and
   grade the complete artifacts without extrapolation or target-based case
   selection.
4. Score Freddie's exact SEM only after the blind freeze. Decompose any miss
   into reactor boundary, blanket TiO2/Cr law, feature transport/charging, mask
   survival, and post-etch mechanics.
5. Use the Bosch heldout result as the validated reactor-side anchor, then add
   independent feature-profile tests and another chemistry without changing
   the common core.
6. Keep the literature library as the provenance spine. Every new coefficient
   needs source support, unit conversion, evidence grade, and explicit
   calibration/validation exposure.

## Overall completion criterion

The campaign goal remains active. It is not complete until:

- Oxford's blind board is numerically complete and scored against the exact
  withheld SEM without post-reveal parameter selection;
- all Krueger forecasts terminate and their evidence-domain status is graded;
- Krueger is either closed using independently constrained boundary/surface
  data or reduced to a demonstrated minimal experiment;
- Bosch plus at least two additional chemistries show heldout absolute
  depth/profile transfer with the unchanged multiphysics core;
- conservation, convergence, provenance, uncertainty, and target firewalls
  remain green.

The takeover agent should continue the live runs, not restart the campaign,
not revive the old HTML number, and not turn a conditional envelope into a
claimed experimental match.
