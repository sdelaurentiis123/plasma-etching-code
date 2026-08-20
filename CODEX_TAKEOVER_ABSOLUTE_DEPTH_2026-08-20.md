# Codex takeover — absolute depth and Oxford blind profiles

Status snapshot: 2026-08-20 15:42 EDT

Repository: `plasma-etching-code`

Branch: `codex/validation-first-multiphysics`

Baseline pushed checkpoint inspected: `d852a1f`

Read this file first when resuming the absolute-depth campaign.  Then read
`PROGRESS_STATE_2026-08-20.md` for the full scientific record.  This report
supersedes the live-process instructions in both
`HANDOFF_MOVING_CR_BOARD_2026-08-20.md` and
`CODEX_TAKEOVER_MOVING_CR_2026-08-20.md`.  Those two files are retained as
forensic history; do not follow their old v3/v4 recovery sequences.

## Executive state

The alarming “55 of 56 complete, one unexplained CPU spin” state is over.
There were two distinct problems, and both are now identified:

1. The original v3 trajectory suffered a deterministic unequal-area strip
   projection singularity.  A tiny marching-cubes triangle received an
   enormous flux density and requested 5,331,155 ordinary advection substeps.
   Commit `d8a8d02` fixes the projection while conserving total event rate.
2. The corrected v4 calculation then revealed that every board cell stopped
   at first local Cr-mask loss.  This was not another hang.  It exposed a real
   missing engine capability: certified continuation when one material
   disappears while another survives.  Commit `84c2c01` implements that
   lifecycle, and `44956d0` advances the Oxford driver through it.

The exact formerly hard cell has now passed under v5.  It crossed complete Cr
extinction once, conservatively retired the Cr surface state, initialized the
newly exposed TiO2 surface through common-refinement remap, and continued to
either the requested 1200 seconds or a declared physical domain breakthrough.
That result is committed in `bdabd8a`.

The full 56-cell v5 board is live on the one authorized Vast instance.  At
this snapshot it had nine completed v5 content-addressed caches—the exact
acceptance cell plus eight board cells—and eight healthy workers computing the
remaining cells.  No target SEM or target depth has been used to tune the
calculation.

## Live execution update at 15:42 EDT

The stale `55 of 56 / CPU spin` handoff is forensic history, not an active
failure.  A direct process inspection at this snapshot found:

- Oxford parent PID `8980` alive for 22 minutes with eight runnable workers;
- all eight Oxford workers near one fully used CPU each, stable resident
  memory, and no runaway-substep signature;
- exactly **9 of 56** v5 content-addressed trajectory caches present;
- no final v5 audit is accepted yet—the board parent is still running;
- host memory has 232 GiB available;
- the GPU has a valid 2.7 GiB context and was idle at the instant sampled,
  which is expected during CPU-heavy surface and level-set phases.

Three additional preregistered Krueger forecasts are running in parallel on
the same authorized box from the separately verified `/root/petch-d852a1f`
tree:

| Case | Supervisor PID | Child PID | Last sealed point at snapshot |
| --- | ---: | ---: | --- |
| published aggregate ion, identity unresolved | `10284` | `10288` | step 6, 0.09375 s, 0.9830 nm |
| all aggregate ions declared CF2+ | `10285` | `10289` | step 6, 0.09375 s, 0.8547 nm |
| all aggregate ions declared CF3+ | `10286` | `10290` | step 6, 0.09375 s, 0.8278 nm |

All three runtime audits say `experimental_outcomes_read=false`.  The CF2+
and CF3+ calculations are composition-envelope endpoints, not fitted reactor
mixtures.  Their early transient ordering must not be interpreted as the
60-second ordering and must not be changed to move toward the 825 nm target.

These are deliberately expensive converged feature integrations: 3,840
nominal steps per 60-second trajectory.  The first six accepted steps took
roughly 45 seconds each under the current shared-host load.  A naive linear
projection is therefore about 48 hours per full trajectory; a declared
physical clogging event could end a trajectory earlier.  The jobs checkpoint
and automatically resume every 30 wall-clock minutes.  Do not promise an
overnight Krueger answer, and do not replace the frozen resolution or time
step merely for speed without a preregistered convergence/performance study.

The archive used for those forecasts is `/private/tmp/petch-d852a1f.tgz`
locally and `/root/petch-d852a1f.tgz` remotely, with SHA-256:

`153b416c236ab1493389fd050181d5afdb1279bdb719d522fa8c03eefff524e0`

Logs and outputs:

- `/root/krueger_guo_60s_nominal.log` and
  `/root/krueger_guo_60s_nominal/`;
- `/root/krueger_guo_60s_cf2.log` and `/root/krueger_guo_60s_cf2/`;
- `/root/krueger_guo_60s_cf3.log` and `/root/krueger_guo_60s_cf3/`.

Use the exact PID files under `/root/krueger_guo_60s_*.pid`; never use
`pgrep -f` for liveness.  Do not destroy instance `48177892` until the Oxford
artifacts and every Krueger checkpoint/result worth retaining have been
retrieved and hash-verified.

## What the software now does

The working stack is:

```text
machine settings
  -> conserved reactor chemistry and electron kinetics
  -> sheath / ion-energy model
  -> radial wafer transfer
  -> species-, energy-, and angle-resolved wafer boundary
  -> material-routed surface kinetics and self-consistent charging
  -> conservative multi-material state evolution
  -> evolving three-dimensional profile
```

The common feature engine now supports a resolved material disappearing
without silently freezing the geometry, relabeling its state as another
material, or weakening topology checks.  Retirement is allowed only when:

- the material owns no geometry nodes;
- it contains no complete eight-corner positive volume cell;
- remaining positive level-set support is below a roundoff-scaled bound;
- that material's component count reaches exactly zero;
- surviving materials, solid components, gas cavities, and breakthrough
  state remain topologically consistent;
- the caller explicitly requests the material-extinction policy.

Default policies still refuse material loss.  Resolved hidden support and
above-roundoff residue still refuse.  Retired fields are integrated into an
explicit ledger and never transferred to the surviving material.

## Oxford/Freddie condition being predicted blindly

The frozen condition is:

- Oxford PlasmaPro NPG80 RIE;
- 55 / 5 / 1 sccm CHF3 / SF6 / O2;
- 30 mTorr;
- 150 W forward table RF;
- 20 C table temperature;
- 1200 s process duration;
- 700 nm ALD TiO2 on fused silica;
- 45 nm Cr hard mask;
- square-pillar prior, approximately 400 nm pitch;
- preregistered widths 80, 120, 160, 200, 240, 280, and 320 nm.

The current conditional reactor boundary contains 67 species, including 20
positive ions and 37 thermal neutrals.  Its central positive-ion flux is
approximately `1.457e19 m^-2 s^-1`; the conditional powered-electrode sheath
drop is approximately 296 V, giving about 299 eV for a singly charged ion.
The ion mixture remains species resolved rather than being collapsed into a
single generic ion.

This is a blind prediction envelope, not yet a unique Oxford-tool answer.  The
recipe screen gives forward generator power, but not the achieved electrode
waveform/self-bias or absorbed power.  The TiO2 and Cr absolute rate axes are
transferred from independent process witnesses because state-resolved
TiO2/Cr coefficients for this exact CHF3/SF6/O2 boundary are not published.

## Completed Oxford evidence

### Corrected v4 forensic board

Commit `30c3054` preserves all 56 corrected v4 trajectories and
`audit_v4.json`.  The remote native check passed, and the transfer archive
matched locally and remotely at:

`1fdad5f88f83c3e15704193968673de48358830c5c44483aab1c4648b03063ee`

All v4 paths stopped at first local Cr loss between 391.334 and 713.846
process-equivalent seconds.  Their conditional depth range was 274.896 to
406.087 nm.  Particle balance was exact and maximum remap residual was
`1.683e-15`.  These are valid pre-extinction receipts, not final 1200-second
profiles and not eligible as v5 caches.

### Exact v5 acceptance trajectory

The acceptance cell is:

- width 320 nm;
- TiO2:Cr selectivity 14;
- low ion-energy axis, 146.539 eV;
- zero angular-tail fraction;
- 10 nm production grid;
- model revision `two-material-moving-tio2-cr-dose-factorization-v5`.

Its local cache is:

`results/curated/zhu_npg80_moving_cr_profiles_v1/trajectories/w320_s14.000_ion_low_tail_0p0_593697ba778b8990.json`

Local and remote cache hash:

`4feeec492c52af860298495bd190e50106370306d208acd8fe24a7170cd1647f`

Cr extinction occurred at 560.429 s.  The final apparent positive Cr residue
was `1.937e-14` mesh units versus a `1.137e-13` numerical tolerance and had
zero resolved volume cells.  The event retired 2,284 Cr faces, preserved
particle balance exactly, and had maximum remap residual `7.166e-16`.

The two independent, no-target-fit surface-rate endpoints are:

- 34.125 nm/min: requested 1200 s reached; depth 646.856 nm; middle/bottom CD
  313.08/321.78 nm; sidewall angle 76.03 degrees from the wafer;
- 43.4667 nm/min: domain-gas breakthrough at 951.764 s; last accepted depth
  653.356 nm; middle/bottom CD 312.05/321.73 nm; sidewall angle 76.17 degrees.

The reported top CD of zero after mask loss is not a conventional zero-width
pillar.  It means the original top reference plane no longer intersects the
surviving profile.  Grade the full geometry and terminal class, not that
single scalar.

## Live v5 campaign

Only the following paid instance belongs to this task:

- Vast instance `48177892`;
- SSH `ssh -p 17892 root@ssh6.vast.ai`;
- repository `/root/petch-4b656fd`;
- venv `/root/petch-venv`;
- RTX 3090;
- billing `$0.2011111111/hour`;
- parent pid file `/root/zhu_v5_board.pid`;
- parent PID at this snapshot `8980`;
- log `/root/zhu_v5_board.log`;
- eight deterministic workers.

At inspection, all eight workers were runnable and each used approximately
100% of one CPU, total host memory usage was only 16 GiB of 251 GiB, and each
worker held a 322 MiB CUDA context.  Instantaneous GPU utilization was zero
during the CPU-heavy surface/level-set phase.  That is not evidence of the old
singularity.  The campaign writes a trajectory cache when a cell completes,
so the v5 count can remain unchanged for many minutes while all workers are
making progress.

Never infer liveness with `pgrep -f`, because the probe can match itself.
Read the pid file and inspect that exact PID and its children.  Count only
JSON files whose `job_spec.model_revision` ends in `v5`; the directory also
contains 55 v3 and 56 v4 forensic caches.

Do not stop or destroy the instance while the parent and healthy workers are
running.  Do not touch any other Vast instance.

## Required landing sequence for v5

When the parent exits normally:

1. Verify that exactly 56 v5 trajectory caches exist and that `audit.json`
   was written.
2. Run the remote native `--check` with the same tree and venv.  It should
   reuse all content-addressed caches rather than recompute them.
3. Package only the v5 caches plus the v5 audit, hash the archive remotely,
   retrieve it, and verify the local hash before accepting it.
4. Preserve v3 and v4 receipts; do not overwrite `audit_v4.json`.
5. Run the native local `--check` against the current committed code.
6. Summarize terminal classes, depth/CD/angle ranges, extinction times,
   particle balance, remap residuals, and any refusals.
7. Run the full repository test suite.
8. Commit and push the complete v5 board and updated result README.
9. Destroy only Vast instance `48177892`, then verify it no longer exists.

If the parent exits abnormally, preserve the log and every completed cache,
identify the exact failing content-addressed cell, and reproduce that cell
alone.  Do not relaunch the whole board blindly and do not weaken a
certification gate to make it finish.

## What remains before Oxford is experimentally matched

The eventual SEM is the held-out answer key, not a fitting input.  The most
valuable independent same-run information is:

1. achieved DC self-bias or electrode voltage/current waveform;
2. same-run blanket TiO2 thickness loss;
3. same-run remaining Cr thickness;
4. exact GDS or width/pitch map and sample radius/orientation;
5. rinse/dry and Cr-strip history;
6. cross-section and top-down SEMs with scale bars.

Self-bias constrains ion energy; blanket TiO2 loss constrains target-tool
surface removal without feature transport; remaining Cr constrains
selectivity; GDS and placement constrain geometry and radial boundary.  Those
measurements can collapse the current physical envelope without fitting the
target profile.

The axisymmetric reactor model predicts only a small smooth flux variation
across the central optic.  Spatially clustered fallen pillars therefore do
not, by themselves, prove a smooth top-down plasma-flux gradient.  Pattern
loading, mask undercut/adhesion, local micromasking, sample placement, and
rinse/dry mechanics remain viable causes.  The existing mechanical audit
finds intact pillars roughly 8–11 times stiffer than the simple capillary
collapse threshold, making etch-weakened bases or local process defects more
plausible than spontaneous collapse of ideal pillars.

## Krueger status

Krueger is not honestly matched under the paper's published aggregate wafer
boundary.  The current published-input result is approximately 346.833 nm
versus the reported 825 nm.  The old 790–811 nm website-era agreement was
retracted because two implementation errors canceled.

An independent surface-transfer route is promising and must be run to full
feature duration next:

- Guo/Kwon planar transfer predicts effective removal yield 2.613 versus the
  2.521 implied by the Krueger depth, about 5% apart;
- a finite-fluence planar sensitivity lands near 855 nm;
- a deterministic-extruded feature prefix predicts 11.93–12.27 nm/s, about
  11–13% below the 13.75 nm/s depth-average requirement;
- no Krueger depth was used to select the transferred surface parameters.

This is not yet an authorized match.  Krueger's IEAD extends beyond much of
the independent fit domain, and the paper omits species-resolved ion
composition, C4F6 parent flux, and enough neutral-boundary detail.  Run the
existing deterministic-extruded `GuoC4F8ArSiO2FeatureMechanism` forecast
without selecting parameters from the 825 nm target, then grade both its
prediction and its domain-of-validity violations.

The strongest eventual closure is a validated C4F6 reactor boundary or the
authors' species-resolved HPEM/PCMCM wafer output—not an optimizer-adjusted
yield chosen to hit 825 nm.

## Cross-chemistry status

Existing unchanged-core evidence includes:

- direct `SF5+ -> Si/SiO2` surface-dose points at about 5.88% MAPE;
- `Cl2/Ar+ -> Si` ALE points at about 12.88% maximum error;
- a Mahorowala Cl2 feature board near 19% MAPE without a formal pass gate;
- direct fluorocarbon beam agreement at feature energy with an existing 4.7%
  grade.

These checks are real, but surface points are not feature-profile validation.
The mission still requires at least two additional formal held-out
feature-depth/profile gates using the unchanged core.

## Generality and atomic-accuracy boundary

The numerical engine is general across chemistry and material decks.  New
chemistries can supply species/reaction data, machine coupling, sheath and
wafer-transfer evidence, surface laws, geometry, and held-out measurements
without changing transport, charging, conservative state evolution,
multi-material routing, topology checks, or level-set motion.

That does not mean arbitrary named chemistry is automatically predictive.
The evidence-bearing reactor and surface coefficients must exist or be
measured.  Today the production grid is 10 nm and the surface state is a
reduced kinetic description.  The system is atom-balanced, but it is not an
atom-by-atom final-geometry oracle.  DFT/MD can close selected barriers and
yields; it cannot infer missing self-bias, wall condition, radical flux, or
wafer placement from a recipe screen.

The defensible endpoint is a fast, deterministic, differentiable,
uncertainty-propagating knobs-to-profile predictor whose held-out error
approaches the measurement floor.

## Version-control and validation map

Key pushed checkpoints:

- `d8a8d02` — fix unequal-area strip flux amplification;
- `30c3054` — land all corrected v4 receipts;
- `84c2c01` — certified continuation through material extinction;
- `44956d0` — advance Oxford board past Cr loss;
- `242177b` — update the v5 content-addressed cache sentinel;
- `bdabd8a` — land the passed exact v5 continuation receipt.

The lifecycle-focused group is 138 passed.  The v5 driver/core group is 118
passed.  Before the final cache-sentinel update, the full repository run was
2131 passed, 7 skipped, and one stale expected v4 filename; the updated
debug/driver sentinel group is 9 passed.  A final full run remains mandatory
after the 56-cell v5 board lands.

The following pre-existing untracked paths belong to the user and remain
untouched:

- `results/curated/mixed_layer_feature_v1/ml20-bonds-12s.log`;
- `results/curated/mouth_equilibrium_probe_dx/`;
- `scratch_ignore_calc.py`.

Do not reset, delete, or absorb them into this campaign.

## Immediate priorities after Oxford v5 lands

1. Finish and certify the full deterministic Guo/Krueger feature forecast.
2. Turn two independent chemistry datasets into formal held-out feature-depth
   gates with the same common engine.
3. Replace Oxford's cross-machine rate axes with independent target-tool
   self-bias, blanket TiO2 loss, and Cr survival when Freddie supplies them.
4. Freeze the predicted board before opening the target SEM, then score depth,
   top/middle/bottom CD where defined, sidewall angle, bowing, mask survival,
   terminal class, and spatial variation.

The current accomplishment is substantial but bounded: the full stack can now
evolve a real two-material mask/film system through complete mask extinction
without a numerical singularity or silent freeze.  Oxford absolute rates are
still conditional; Krueger still needs its full independent transfer forecast;
and cross-chemistry feature-profile validation is not yet complete.
