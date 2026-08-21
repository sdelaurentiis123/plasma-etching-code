# Codex takeover report — reactor-to-depth campaign at 2026-08-21 17:00 UTC

Repository: `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code`

Branch: `codex/validation-first-multiphysics`

This report is the current operational handoff. It supersedes the live counts,
PIDs, and diagnosis in `HANDOFF_MOVING_CR_BOARD_2026-08-20.md` and the earlier
August 21 takeover reports. Those older files remain useful forensic records.

## The actual mission

Achieve defensible absolute-depth prediction rather than a one-target fit:
identify and implement missing reactor and surface closures, validate them on
independent measurements, and demonstrate heldout depth/profile accuracy for
Krueger C4F6/Ar/O2 SiO2 plus at least two additional chemistries using the
unchanged multiphysics core.

Do not mark the campaign complete merely because one trajectory board finishes.
The finish line is a validated machine-knobs-to-wafer-boundary-to-feature-depth
chain, or a rigorous statement of the minimum real experiment that remains
necessary to identify it.

## Executive answer: what is happening

The August 20 story that one Oxford cell was stuck at 55/56 is stale.

1. The old v3 CPU spin was real and reproducible. It exposed a chain of
   topology, surface-remapping, and mask-extinction defects. Those were fixed
   rather than timed out or hidden.
2. A clean v6 rerun then failed closed after 23/56 cells for a different reason:
   one sparse interpolation coefficient was
   `-2.2204460492503131e-16`, one binary64 ulp below zero. This was a numerical
   row-closure defect, not plasma physics and not a stochastic hang.
3. Commit `84a786e` repairs row closure through the largest safely positive
   coefficient. No coefficient clipping was added and no validator or
   conservation criterion was relaxed. The exact failed checkpoint advances on
   arm64 CPU and x86_64 CUDA with pinned fingerprints.
4. Because even roundoff-scale stored weights changed, the active Oxford model
   revision is v7 and all 56 cells are being recomputed. v6 and v7 caches must
   never be mixed.
5. At the 16:58 UTC snapshot, Oxford v7 had produced 15/56 trajectories. All 15
   are now retrieved locally, checked, and included in this checkpoint. Eight
   workers remain healthy at about 100% CPU each with real CUDA allocations.
6. Bosch has produced the first sealed chronological machine-to-wafer-depth win:
   0.239 um wafer-mean depth MAE and 0.541% MAPE over 13 unseen wafers, each with
   89 radial points. This is a real heldout result, but it is wafer depth rather
   than a feature-profile validation.
7. Three Krueger no-fit sensitivity forecasts are healthy and checkpointing at
   about 26.5-28.1 of 60 physical seconds. They are deliberately outside the
   declared predictive evidence scope because Krueger published an unresolved
   aggregate ion boundary and the transferred surface response is not fully
   supported at his energies and chemistry.

The honest headline is:

> The old Oxford hang is fixed and the clean v7 blind board is progressing.
> Bosch has passed a sealed heldout wafer-depth test. Krueger is still an
> unresolved reactor-boundary and surface-transfer problem, not a matched
> 825 nm result.

## Repository safety

Use only this repository and branch. Do not switch to the old website branch,
`reactorlab`, `petch-torchsim`, or another stale implementation.

At the start of this checkpoint, local and remote were both at `438f536`.
The only unrelated untracked paths are user work. Never edit, stage, remove,
reset, clean, or absorb them:

- `results/curated/mixed_layer_feature_v1/ml20-bonds-12s.log`
- `results/curated/mouth_equilibrium_probe_dx/`
- `scratch_ignore_calc.py`

Use exact-path staging. Do not use destructive Git commands. Commit and push
recoverable scientific checkpoints as they land.

## Paid compute

Only Vast instance `48177892` belongs to this campaign.

- SSH: `ssh -p 17892 root@ssh6.vast.ai`
- Oxford tree: `/root/petch-4b656fd`
- Krueger tree: `/root/petch-d852a1f`
- Python environment: `/root/petch-venv`
- accelerator: RTX 3090, 24 GiB
- last observed price: about `$0.2011/hour`

Do not touch another instance. Do not destroy this instance until both Oxford
and Krueger campaigns are terminal, all active artifacts have been retrieved
and hash-verified, local production checks pass, and the results are committed
and pushed.

Never use `pgrep -f` alone as liveness evidence; the probe can match itself.
Use exact known PIDs, `ps --ppid`, CUDA allocation, and recent artifact mtimes.

## Oxford/Freddie blind TiO2/Cr square-pillar board

### Frozen experimental condition

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

Freddie's exact target SEM is withheld from tuning. The frozen uncertainty board
propagates two sourced TiO2:Cr selectivity endpoints, two ion-energy cases, and
two angular cases: `7 x 2 x 2 x 2 = 56` trajectories. Each trajectory carries
two dose-normalized TiO2 rate endpoints.

The chain includes 67-species reactor chemistry, electron kinetics, sheath
closure, parent/daughter collision chemistry, radial wafer transport,
species-resolved wafer boundary, and deterministic 3D moving TiO2/Cr feature
evolution with charging, chemical/physical removal, polymerization, material
extinction, mask erosion, and conservative surface-state remapping.

It is a conditional blind prediction, not yet a unique atomic-accuracy tool
prediction. The exact-run self-bias and same-tool TiO2/Cr blanket response are
not known. Those uncertainties are exposed as a preregistered ensemble rather
than silently fit to the target SEM.

### Numerical repair certification

The captured v6 failure is quarantined under
`results/quarantine/zhu_npg80_v6_failure_20260821/`.

- failed checkpoint SHA-256:
  `1fbd52695a1df0379c67d22d9743d8e74bee29f2322142fd41737cafd7881dd6`
- diagnostic log SHA-256:
  `6683045c63d1efd0978a46d8940325043d758bb8c5a5aaac446c5bf092fc9a16`
- failed state: width 120 nm, selectivity 18.016664, high energy, no angular
  tail, 642.5613497 s accepted physical time
- active model revision:
  `two-material-moving-tio2-cr-positive-row-closure-v7`
- repaired transfer fingerprint:
  `64dafad5621689e785800294ec3408842939d24095c7c85d23c54e1abf498d96`
- repaired next-mesh fingerprint:
  `950c8f9a1176e4184020ef9f6a72d9e4ff6c61d201ba58efa4351a6478264bac`
- unresolved-node reassignments: zero
- local/remote maximum remap residual on the captured step:
  `1.371e-16` / `2.742e-16`
- focused suite: 149 passed
- repository suite: `2229 passed, 7 skipped`

### Exact v7 state at 16:58 UTC

Processes:

- wrapper shell: PID 33307
- board parent: PID 33309
- resource tracker: PID 33424
- workers: PIDs 33425, 33428, 33429, 33430, 33431, 33432, 33433, 33434

All eight workers were runnable at about 100% CPU. Every worker held a 322-354
MiB CUDA allocation. Host and device memory were healthy. The remote log only
contains Warp initialization because cell completion is expressed by atomic
cache writes; use exact processes and mtimes for liveness.

Fifteen active v7 caches have been retrieved. The authoritative hash ledger is:

`results/curated/zhu_npg80_moving_cr_profiles_v1/V7_LIVE_EXECUTION_2026-08-21.md`

Validation of all 15:

- exact v7 revision and frozen preregistration hash
- exactly two rate endpoints per trajectory
- exact zero transport particle-balance error
- maximum surface-state remap residual below `1.47e-15`
- width-80 conditional etched depth: approximately 683.71-684.09 nm
- seven completed width-120 cases: approximately 684.05-684.32 nm

Some high-rate endpoint profiles terminate by `domain_gas_breakthrough` after
the TiO2 clearing front reaches the modeled domain; the paired trajectory still
preserves its conservative ledgers. This must be interpreted explicitly in the
final atlas, not mislabeled as an ordinary requested-duration endpoint.

### Oxford continuation

1. Leave the current parent/workers alone while CPU/CUDA use and artifact mtimes
   move.
2. Retrieve only exact new v7 filenames; the directory contains historical
   revisions and must not be wildcard-counted or copied wholesale.
3. For every batch, verify revision, preregistration hash, both rate endpoints,
   terminal reason, particle ledger, remap ledger, and file hash.
4. Append the hash ledger and commit/push exact batches.
5. At 56/56, allow the remote writer to assemble the final audit; copy it and
   all active artifacts; run the production audit locally in `--check` mode.
6. Render and inspect the atlas, run focused tests and the full suite, and freeze
   the board in Git before revealing Freddie's SEM.
7. After reveal, digitize and score depth, top/middle/bottom CD, taper, bow,
   Cr survival, missing/collapsed pillars, and radial position. Never select an
   ensemble corner post hoc to manufacture agreement.

Highest-value future Freddie measurements: exact-run DC self-bias, blanket
TiO2 loss, remaining Cr thickness, GDS/width mapping, sample radius/orientation,
and matched top-down/cross-section SEMs. These identify submodels; they are not
permission to fit the final contour.

## Krueger C4F6/Ar/O2 SiO2 trench

The paper endpoint is about 825 nm after 60 s. The old website's 790-811 nm
match is retracted: it arose from canceling implementation errors. The corrected
published-input endpoint is about 346.833 nm. Direct-beam evidence indicates
that 825 nm at the printed aggregate ion flux requires removal above the
supported physical-yield ceiling.

Three frozen Guo-transfer sensitivity cases remain live:

- nominal unresolved aggregate ion: supervisor PID 10284
- aggregate ion treated as CF2+: supervisor PID 10285
- aggregate ion treated as CF3+: supervisor PID 10286

At 16:58 UTC each supervisor had one live child around 100% CPU and fresh
`audit.json` and `checkpoint.npz` files. The latest accepted records were:

| Case | Physical time | Step | Etch depth | Mask opening |
|---|---:|---:|---:|---:|
| nominal unresolved | 28.140625 s | 1801 | 230.967 nm | 17.240 nm |
| all CF2+ | 26.796875 s | 1715 | 259.231 nm | 21.684 nm |
| all CF3+ | 26.484375 s | 1695 | 264.862 nm | 20.281 nm |

Do not extrapolate these evolving-topology snapshots to 60 s and do not select
the closest branch after seeing the paper target. All three receipts correctly
declare `parameter_evidence_supports_prediction=false` and
`within_declared_scope=false`.

The unresolved physics is not merely runtime:

- Krueger reports an aggregate positive-ion population rather than a
  species-resolved ion flux and spectrum.
- neutral channels extend beyond the Guo source list.
- the trajectory samples energies above the validated <=370 eV response board.
- the angular-law typesetting repair is declared rather than independently
  measured for this reactor.
- mask and translating-layer coefficients cross chemistry and remain
  nonpredictive parameters.

The existing C4F6 evidence already rejects the shortcut of treating a direct
70 eV mass spectrum as the reactor ion mixture. Benck's reactor data show strong
pressure/feed-dependent secondary chemistry and heavy-ion current outside the
retained light ions. A common-loss inverse gives unphysical negative inferred
CF3 density; a fixed species-selective loss repairs only the Ar mixtures and
fails pure C4F6. The smallest defensible next global model must retain at least
C4F6, C3F3, CF, CF2, CF3 and their important ions, secondary ionization,
ion-neutral conversion, pressure-dependent residence, wall/Bohm loss, and
heavy-fragment current closure. It must be graded on the independent Benck
feed/pressure board before transfer to Krueger. It must not be tuned to 825 nm.

## Cross-chemistry validation state

### Bosch cyclic SF6/C4F8 silicon: sealed heldout win

The chronological heldout prediction was sealed before reading the 13 test
wafers. Across 89 radial points per wafer:

- wafer-mean Si depth MAE: 0.238906 um
- wafer-mean Si depth MAPE: 0.541245%
- point RMSE: 0.329751 um
- shape error: 0.261855%
- oxide depth MAE: 0.023908 um
- selectivity MAPE: 3.161%
- predicted/observed correlation: 0.87341

All preregistered absolute and baseline gates pass. The shape-improvement
bootstrap interval crosses zero, so do not claim a decisive shape-baseline win.
This is one tool and a wafer-depth board, not feature-scale universality.

### Additional chemistry evidence

- Tinacba SF5+ on Si/SiO2 and Vella-Hao Cl2/Ar+ silicon ALE provide direct
  retrospective depth-per-dose boards. They constrain surface response but are
  not yet formal heldout feature-profile passes.
- Mahorowala pure-Cl2 reactor-to-feature currently gives about 15.0-15.5%
  surface-plane MAPE and about 19.3% evolving-feature MAPE. It lacks a formal
  heldout gate and is still limited by species-resolved/radial ion flux and IEAD.
- Yoshie 2023 is the strongest next heldout feature dataset: 49 feature depths
  were preregistered before digitization, with same-reactor blanket rates and
  electron-density/OES/XPS diagnostics. Its timing-dependent feature/blanket
  ratio and rank changes prove a scalar blanket-rate rescaling is insufficient;
  the model must carry cycle and material history. Existing stratified C/F/Si
  state and the calibration-excluded Humbird-Graves response are the appropriate
  physics base, but a complete frozen reactor boundary and feature forecast have
  not yet been assembled.

The core is generalizable by species/material decks and conserved states. It is
not truthful to say arbitrary chemistry is predictive merely because a JSON
deck can be written. Each new chemistry still needs reaction provenance,
machine-boundary identifiability, surface-response support, and heldout grading.

## Immediate order of work

1. Continue retrieving, checking, committing, and pushing Oxford v7 batches.
2. Finish all 56; locally certify/render/freeze before opening the target SEM.
3. Let all three Krueger branches reach terminal 60 s; retrieve and grade them
   without extrapolation or branch selection.
4. Implement the smallest atom/charge-conserving C4F6/Ar global reactor closure
   authorized by the Lan-Jeon, NIST, and Benck evidence. Calibrate/validate model
   form on Benck, not on Krueger depth.
5. Transfer only the independently supported boundary uncertainty into Krueger
   and determine whether 825 nm becomes physically reachable. If not, state the
   missing measurement precisely.
6. Assemble the Yoshie cyclic SF6/C4F8 heldout feature forecast using the
   unchanged core and frozen preregistration. Do not equate OES intensity with
   absolute ground-state radical flux.
7. Only after these independent gates, update any public HTML claims.

## What must never be claimed yet

- Oxford matches Freddie's SEM: the exact target is still withheld.
- Oxford is atomically accurate: self-bias and same-tool surface response are
  still interval inputs, and a continuum profile cannot validate atom-by-atom
  dynamics from an SEM.
- Krueger's 825 nm depth is matched: it is not.
- The old website match remains valid: it does not.
- One Bosch tool proves universal knobs-to-profile prediction: it does not.
- Any chemistry is predictive merely by adding a deck: it is not.

The science is progressing because every failure is being made deterministic,
quarantined, repaired against unchanged invariants, and rerun from a clean
revision. The remaining failures are now being separated into numerical defects
that can be fixed in code and physical non-identifiability that requires either
independent literature support or a specific machine measurement.
