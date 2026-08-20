# Codex takeover — reactor-to-depth and blind-profile campaign

Status snapshot: 2026-08-20 16:24 EDT

Repository: `plasma-etching-code`

Branch: `codex/validation-first-multiphysics`

Pushed implementation checkpoint: `aaeec38`

This is the canonical resume document.  It supersedes the live-state portions
of `CODEX_TAKEOVER_ABSOLUTE_DEPTH_2026-08-20.md`,
`CODEX_TAKEOVER_MOVING_CR_2026-08-20.md`, and
`HANDOFF_MOVING_CR_BOARD_2026-08-20.md`.  Keep those files as forensic history,
but do not follow their old recovery instructions without first checking this
report and the actual process IDs.

## Read this first: the honest state

There is no longer an unexplained Oxford CPU spin.  That failure was identified
as a deterministic unequal-area flux-projection singularity, fixed, and
regression-tested.  The exact formerly failing trajectory completed under the
corrected operator.  The feature engine was then extended to continue
conservatively through complete Cr-mask extinction instead of silently
freezing at first mask loss.

The current state is three separate experiments:

1. **Oxford TiO2/Cr blind board:** the corrected v5 board is actively running
   and healthy.  At the live inspection for this report, 24 of 56 v5 cells had
   completed.  One exact hard cell has already passed through Cr extinction.
   The board is a conditional surface-rate and ion-boundary envelope; it is not
   yet a match to Freddie's held-out SEM.
2. **Krueger SiO2 trench:** the published-input model still predicts about
   346.833 nm rather than 825 nm.  The old website-era 790--811 nm agreement is
   retracted because two bugs canceled.  Three independent Guo-to-Krueger,
   no-target-fit forecasts are actively running and had reached 1.000 s of a
   requested 60 s.  No Krueger match is claimed yet.
3. **SPTS Bosch Si/C4F8 depth gate:** the official machine dataset, chronological
   holdout, reduced deterministic reactor, radial wafer transfer, conserved
   C4F8 film memory, unchanged Belen silicon law, and oxide-mask law now exist
   in code.  The path runs in seconds and closes its ledgers, but its untouched
   defaults predict only 0.289 um mean Si depth for the first calibration wafer
   versus roughly 52 um in that wafer's allowed 89-point calibration record.
   The scaffold is therefore not calibrated and has not passed the holdout.

The bottom line is substantial but bounded: the software architecture now
spans measured machine waveforms to radial wafer depth and evolving
multi-material 3-D profiles.  Absolute depth is not yet closed across all
three experiments.

## Correct repository and protected working-tree state

Use only:

```text
/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code
codex/validation-first-multiphysics
```

Do not switch to the older branches or the other repositories for this
campaign.  At the checkpoint for this report, `aaeec38` is pushed to the
private remote.

The following unrelated, pre-existing untracked paths belong to the user and
must remain untouched:

- `results/curated/mixed_layer_feature_v1/ml20-bonds-12s.log`;
- `results/curated/mouth_equilibrium_probe_dx/`;
- `scratch_ignore_calc.py`.

Do not reset, clean, delete, stage, or absorb them into this campaign.

## Live compute: exact state and ownership

Only Vast instance `48177892` belongs to this task:

- label `petch-zhu-moving-cr-bigram`;
- SSH `ssh -p 17892 root@ssh6.vast.ai`;
- RTX 3090, 28 vCPUs, approximately 251 GiB host memory;
- current price `$0.2011/hour`;
- Oxford tree `/root/petch-4b656fd`;
- Krueger tree `/root/petch-d852a1f`;
- virtual environment `/root/petch-venv`.

Other Vast entries visible in the account are unrelated and must not be
touched.  Do not destroy `48177892` until useful Oxford artifacts and all
Krueger checkpoints/results have been retrieved and hash-verified.

### Oxford live process

At the 16:19 EDT inspection:

- parent PID file `/root/zhu_v5_board.pid` contained PID `8980`;
- PID `8980` was alive with eight runnable workers;
- workers used roughly one CPU each, with normal memory usage;
- host memory had about 234 GiB available;
- the GPU held a valid context and was lightly used at the instant sampled;
- 24 of 56 content-addressed v5 caches existed;
- the existing `audit.json` predates completion and must not be accepted as a
  final v5 audit merely because the file exists.

Count v5 cells by opening each trajectory JSON and checking that
`job_spec.model_revision` ends in `v5`.  The same directory intentionally
contains v1--v4 forensic caches.  Never use the raw JSON file count.

### Krueger live processes

At the same inspection:

| Case | Supervisor PID file | State at 1.000 s |
| --- | --- | ---: |
| published aggregate ion, identity unresolved | `/root/krueger_guo_60s_nominal.pid` | 12.438 nm |
| all aggregate ions declared CF2+ | `/root/krueger_guo_60s_cf2.pid` | 15.700 nm |
| all aggregate ions declared CF3+ | `/root/krueger_guo_60s_cf3.pid` | 17.296 nm |

The supervisor PIDs were `10284`, `10285`, and `10286`.  Each had one healthy
feature worker.  The workers checkpoint and resume every 30 wall-clock minutes.
The runtime audits say `experimental_outcomes_read=false`.

Do not linearly extrapolate the 1 s prefix to 60 s and call it a result.  Mask
opening, transport, surface state, and topology evolve.  At approximately
40 wall-seconds per accepted 0.015625 s step, a complete trajectory remains a
roughly 42--48 hour calculation unless a declared physical terminal event
ends it earlier.

Always read the exact pid files and inspect those PIDs.  Do not use `pgrep -f`
for liveness because the probe can match its own command line.

## What actually happened to the alleged Oxford hang

The stale handoff said 55 of 56 cells were complete and one worker was spinning
in an unknown loop.  That description applied to v3 and is no longer the live
problem.

The captured failure was:

1. A sparse energetic event carried an integrated rate equal to local flux
   density times the hit triangle area.
2. Periodic strip symmetrization divided that rate equally among equivalent
   triangles and then divided by each recipient's area.
3. A tiny marching-cubes sliver therefore received an enormous flux density.
4. The local recession speed jumped from approximately `7.25e-4` to
   `3310.39` mesh units/s.
5. The ordinary CFL integrator requested 5,331,155 advection substeps.

Commit `d8a8d02` projects uniform flux density and distributes integrated rate
in proportion to recipient area.  It preserves total event rate without the
small-area singularity.  The exact formerly hard cell then completed normally.

The corrected v4 board subsequently revealed a different, physical engine
gap: all cells reached local Cr-mask loss before the requested 1200 s.  Commit
`84c2c01` added a certified material-extinction lifecycle, and `44956d0`
enabled the Oxford driver to continue after Cr disappears.  The policy retires
a material only after its resolved geometry and volume support are gone,
preserves the surviving material and topology, and keeps particle/remap ledgers
active.  The default engine policy still refuses undeclared material loss.

## Oxford/Freddie blind board

The frozen supplied condition is:

- Oxford PlasmaPro NPG80 RIE;
- 55 / 5 / 1 sccm CHF3 / SF6 / O2;
- 30 mTorr;
- 150 W forward table RF;
- 20 C table temperature;
- 1200 s etch;
- 700 nm ALD TiO2 on fused silica;
- 45 nm Cr hard mask;
- square-pillar prior, approximately 400 nm pitch;
- widths 80, 120, 160, 200, 240, 280, and 320 nm.

The current boundary is species resolved: 67 plasma species, including 20
positive ions and 37 thermal neutrals.  The conditional central positive-ion
flux is approximately `1.457e19 m^-2 s^-1`.  The conditional sheath calculation
uses about a 296 V powered-electrode drop, or roughly 299 eV for a singly
charged ion.

The exact v5 acceptance cell is width 320 nm, TiO2:Cr selectivity 14, low ion
energy 146.539 eV, zero angular-tail fraction, and a 10 nm grid.  It passed Cr
extinction at 560.429 s with exact particle balance and a maximum remap
residual of `7.166e-16`.

Its two no-target-fit rate endpoints were:

- 34.125 nm/min: reached 1200 s, depth 646.856 nm, middle/bottom CD
  313.08/321.78 nm, sidewall angle 76.03 degrees from the wafer;
- 43.4667 nm/min: reached a domain-gas breakthrough at 951.764 s, last
  accepted depth 653.356 nm, middle/bottom CD 312.05/321.73 nm, sidewall angle
  76.17 degrees.

These numbers certify numerical continuation and define a conditional blind
envelope.  They do not establish experimental agreement.  The target SEM for
this exact run has not been used as an input or tuning target.

### Why Oxford is still conditional

Freddie's recipe screenshot supplies forward generator power, not the achieved
DC self-bias or electrode voltage/current waveform.  Forward power alone does
not uniquely determine ion energy because matching-network loss, plasma
impedance, electrode area asymmetry, and wall state matter.  The exact
CHF3/SF6/O2 TiO2 and Cr state-resolved surface coefficients are also not
published for this tool.  The current absolute rate axes are transferred from
independent witnesses and carried as an envelope rather than selected to match
the target SEM.

The most valuable same-run data are:

1. achieved DC self-bias or the electrode voltage/current waveform;
2. blanket TiO2 loss from the same recipe;
3. remaining Cr thickness from the same recipe;
4. exact GDS or nominal width/pitch map and sample radius/orientation;
5. cross-section and top-down SEMs with scale bars, viewing angle, etch time,
   and Cr-strip/rinse/dry history.

The SEM is the answer key, not an input needed to make a blind prediction.
Exact mask geometry is needed to predict one exact feature rather than a width
board.  Blanket TiO2 and Cr measurements isolate the surface law from feature
transport and are more diagnostic than fitting a final profile.

Clustered fallen pillars do not by themselves prove a smooth top-down plasma
flux gradient.  The current axisymmetric model predicts only a small smooth
variation across the central sample.  Pattern loading, mask undercut or
adhesion, micromasking, etch-weakened bases, and rinse/dry mechanics remain
viable.  The existing mechanical audit places ideal intact pillars about
8--11 times above a simple capillary-collapse threshold.

## Krueger absolute-depth state

The paper's published aggregate wafer boundary does not close its own reported
depth under the measured surface ceiling.  The honest common-engine result is
approximately 346.833 nm versus 825 nm.  The old HTML's 790--811 nm claim came
from two implementation errors canceling and must never be restored.

The independent transfer path is scientifically promising:

- Guo/Kwon planar transfer gives effective removal yield 2.613 versus the
  2.521 implied by Krueger's depth, about 5% apart;
- a finite-fluence planar sensitivity lands near 855 nm;
- a deterministic-extruded prefix gave 11.93--12.27 nm/s versus the required
  depth-average 13.75 nm/s;
- no Krueger depth was used to select those surface parameters.

That evidence motivated the three currently running full-feature forecasts.
The nominal, all-CF2+, and all-CF3+ cases are a composition envelope because
Krueger reports an aggregate ion channel without resolving the ion identity.
Do not select whichever case lands closest to 825 nm after reveal and call that
a prediction.  Grade all frozen cases, report their domain violations, and
retain the unresolved composition as uncertainty.

The strongest final closure is a validated C4F6 reactor boundary or the
authors' species-resolved HPEM/PCMCM wafer output.  An optimizer-selected yield
chosen from the target depth is not an acceptable closure.

## SPTS Bosch machine-waveform depth gate

The Bosch path is the cleanest current reactor-to-depth validation because the
public source includes actual 5 Hz machine waveforms and independent radial
wafer measurements.

### Source and frozen test

The official Sayyed et al. Zenodo 17122442 record contains:

- 96 SPTS Omega i2L DSi Rapier process records;
- approximately 100 cycles per wafer;
- roughly 4.5 s SF6 etch and 1.5 s C4F8 passivation phases;
- 89-point Si depth and oxide-mask maps;
- a separate 9-point outcome file;
- source/platen RF, reflected power, platen Vpp/DC, pressure, gas, current,
  backside helium, and thermal channels.

Commit `43e80ea` vendors the compact official assets with verified MD5s,
decodes the machine record without opening outcomes, and freezes an
execution-held-out chronological split:

- 76 calibration process records from July 2 through August 7;
- 20 held-out process records from August 21--22;
- one equipment-transfer parameter set shared by every wafer;
- leave-one-lot-out model selection on calibration lots only;
- held-out outcomes forbidden until a model receipt is hash-sealed.

This is explicitly not a pre-exposure blind test because the CSV files existed
in the repository before preregistration.  It is still a valid
execution-held-out test if the loader and receipt firewall are obeyed.

### Implemented physical path

Commit `638ac34` adds a conserved Bosch silicon surface closure:

- unchanged Belen SF6 bare-silicon removal;
- finite La Magna/Garozzo C4F8-derived film deposition and ion removal;
- exact within-step substrate exposure after film depletion;
- Belen silicon removal multiplied by the exposed-time fraction;
- independent silicon and film material ledgers;
- the La Magna SiO2 removal rate discarded from the Si law rather than
  accidentally double-counted;
- the same La Magna law retained for the oxide mask.

Commit `aaeec38` adds the measured-waveform reactor-to-depth scaffold:

- effective atomic F, C4F8 film precursor, and positive-ion populations;
- exact piecewise-constant integration of `dn/dt = S - n/tau` at 5 Hz;
- production capped by both absorbed RF power and inlet particle supply;
- pressure-dependent neutral lifetimes;
- ion energy conditioned on measured platen Vpp;
- deterministic axisymmetric finite-volume radial transfer from an annular
  ICP source;
- exact inventory ledgers and an analytic density-to-radial-flux JVP;
- all radial annuli advanced through the same surface mechanisms;
- no depth regression and no wafer-specific multiplier.

The current focused validation is 35 tests passed.  On the first official
calibration trace, the measured-waveform boundary solves in roughly 0.1 s and
the full radial surface evolution in roughly 4.6 s on the local CPU.  Its
ledger residual is zero at reported precision.

### Current Bosch result is a failed initialization, not a validation pass

For calibration wafer `2024-07-02_01`, the untouched defaults produce:

- mean Si depth `0.288922 um`;
- mean oxide loss `0.002426 um`;
- Si:oxide selectivity `119.07`;
- strong predicted radial falloff toward the wafer edge.

The allowed calibration record has Si depths near 52 um at the first measured
89-point locations.  The default scale is therefore about two orders of
magnitude too small.  This is not hidden or patched.  It means the physically
bounded initialization values for absorbed coupling, effective lifetimes,
wall loss, source geometry, and sheath transfer are not yet tool calibrated,
and possibly that the frozen surface law's declared domain must be tested more
carefully at Bosch conditions.

Do not repair this by multiplying every predicted depth by a fitted number.
The next step is a calibration-only physical fit of shared equipment-transfer
parameters, with leave-one-lot-out checks and explicit parameter bounds.  If
the unchanged surface law cannot reach the calibration scale within those
bounds, that is a scientific model-form failure to diagnose—not permission to
fit the held-out depths.

### Frozen Bosch acceptance gates

After sealing the calibration receipt, score the heldout board against:

- wafer-mean Si depth MAE `<= 1 um`;
- wafer-mean Si depth MAPE `<= 3%`;
- pointwise Si-depth RMSE `<= 1.5 um`;
- normalized radial-shape RMSE `<= 2%`;
- wafer-mean oxide-loss MAE `<= 0.08 um`;
- selectivity MAPE `<= 12%`;
- beat both the calibration-global-mean depth baseline and the
  calibration-mean radial-map baseline.

A Bosch pass validates absolute wafer depth, selectivity, radial transfer, and
drift on one tool.  It does not validate feature charging, ARDE, sidewall
angle, scallops, or a TiO2 chemistry.

## Cross-chemistry validation and generality

Existing unchanged-core evidence includes:

- direct SF5+ to Si/SiO2 surface-dose points at about 5.88% MAPE;
- Cl2/Ar+ to Si ALE points at about 12.88% maximum error;
- a Mahorowala Cl2 feature board near 19% MAPE without a formal pass gate;
- direct fluorocarbon beam agreement at feature energy with an existing 4.7%
  grade.

Those are useful mechanism and surface checks, but they do not yet satisfy the
mission's requirement for at least two additional formal held-out
feature-depth/profile gates with the unchanged core.  Bosch can become one
formal reactor-to-depth gate.  A second independent feature-profile chemistry
still has to be promoted and preregistered.

The core is general in the engineering sense: a new chemistry deck can supply
species/reaction data, reactor coupling, sheath and wafer-transfer evidence,
surface laws, geometry, and held-out measurements without changing the common
transport, charging, conservative remap, material routing, topology, and
level-set operators.  It is not predictive for an arbitrary material merely
because a JSON deck can name it.

“Atomic accuracy” must be stated carefully.  The engine is atom- and
formula-unit-balanced where reactions are resolved, but the production Oxford
grid is 10 nm and the surface laws are reduced kinetics.  It is not an
atom-by-atom final-geometry oracle.  DFT or MD can close selected barriers,
reaction probabilities, and energy-dependent yields; they cannot infer a
missing self-bias, radical flux, wall condition, or sample placement from a
recipe screenshot.  The defensible target is a deterministic,
differentiable, uncertainty-propagating knobs-to-profile predictor whose
held-out error approaches the experimental measurement floor.

## Exact next actions for the taking-over Codex

1. **Do not relaunch the live jobs.** Refresh the exact Oxford and Krueger PIDs,
   v5 cache count, logs, and available memory.  Preserve all checkpoints.
2. **Land Oxford only when complete.** Require exactly 56 v5 caches, run the
   remote native `--check`, package only v5 artifacts plus the v5 audit,
   hash-verify remote and local copies, run local `--check`, summarize terminal
   classes and conservation metrics, test, commit, and push.
3. **Let all three Krueger cases remain frozen.** On completion, retrieve and
   hash-verify all outputs.  Grade every case and its validity domain.  Do not
   choose a composition endpoint based on closeness to 825 nm.
4. **Build a calibration-only Bosch outcome loader.** It may materialize only
   the 76 preregistered calibration keys.  The mixed source CSV contains
   heldout rows, so check the experiment key before converting any numeric
   outcome fields.  Record in the receipt that heldout numeric fields were not
   parsed during fitting.
5. **Fit shared Bosch equipment transfer, not depth.** Candidate physical
   parameters are absorbed-power coupling, F/film/ion lifetimes, wall-loss
   velocities, source ring/width, and sheath transfer within explicit bounds.
   Use leave-one-lot-out calibration and retain one parameter set for every
   wafer.  Test whether the unchanged surface law can physically reach both
   absolute Si depth and oxide loss.
6. **Seal the Bosch model receipt before reveal.** Hash the code, source data,
   calibration keys, parameter bounds, fitted values, cross-validation scores,
   and prediction file.  Only then open and score the heldout outcomes.
7. **Promote a second independent heldout feature chemistry.** Surface points
   alone are not enough for the stated mission.
8. **When Freddie's answer key arrives, preserve blindness.** Freeze the
   Oxford predictions first, then digitize and score depth, top/middle/bottom
   CD where defined, sidewall angle, bowing, mask survival, terminal class,
   and spatial variation.
9. **Destroy only the task-owned instance after retrieval.** Verify the
   Oxford and Krueger artifacts locally first, then destroy only `48177892` and
   confirm it no longer exists.

## Version-control and validation map

Key pushed checkpoints:

- `d8a8d02` — fix unequal-area strip flux amplification;
- `30c3054` — preserve all corrected v4 trajectories;
- `84c2c01` — certified continuation through material extinction;
- `44956d0` — advance Oxford past Cr loss;
- `bdabd8a` — pass the exact hard v5 Oxford trajectory;
- `d852a1f` — preregister the no-target-fit Guo-to-Krueger forecasts;
- `43e80ea` — vendor and preregister the SPTS Bosch holdout;
- `638ac34` — conserved Bosch Si/C4F8 surface closure;
- `aaeec38` — deterministic measured-waveform Bosch reactor-to-depth scaffold.

Focused Bosch/common-surface validation at this snapshot: 35 passed.  A fresh
repository-wide run including the newest Bosch scaffold completed with
**2161 passed and 7 skipped in 1211.85 s**.  No full-suite failure is pending at
this checkpoint.

## One-sentence handoff

The numerical Oxford failure is fixed and the corrected blind board is running;
the honest Krueger match is still being tested through frozen independent
surface transfer; the first real machine-waveform Bosch full-stack path now
runs fast and conservatively but exposes a roughly two-order-of-magnitude
default depth-scale failure that must be resolved through calibration-only
physical equipment transfer before any heldout claim.
