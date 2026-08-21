# Codex takeover report — reactor-to-feature absolute-depth campaign

Snapshot: 2026-08-21 09:00 EDT / 13:00 UTC

Repository: `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code`

Branch: `codex/validation-first-multiphysics`

Scientific HEAD at snapshot: `fefd97066c039c86834ab55658ad2a2de29c78e0`

Remote branch at snapshot: the same commit.

This is the current takeover document. It supersedes the live-state and
next-action sections of `CODEX_TAKEOVER_EXACT_STATE_2026-08-20.md` and
`HANDOFF_MOVING_CR_BOARD_2026-08-20.md`. Those older files remain useful as a
forensic record of the Oxford v5 failure, but they do not describe the current
v6 engine or live campaign.

## Executive answer: what is actually happening

Nothing has been deleted, silently fitted to a target, or lost on the rented
machine. The campaign is alive, version controlled, and running.

The alarming handoff at `d9162fb` described a final Oxford cell spinning after
55 of 56 older trajectories. That was a real failure, but it is no longer the
current state. The failure was converted into a deterministic checkpoint,
reproduced locally, traced to non-idempotent ownership cleanup where the Cr and
TiO2 regional level sets meet, and fixed in `980800a`. The exact former failing
cell then completed cleanly under the new v6 revision. Its trace and cache are
committed in `fefd970`.

At this snapshot:

- the corrected 56-cell Oxford v6 blind-board audit is running from scratch on
  four workers;
- one exact v6 acceptance trajectory is already committed and pinned by a
  regression test;
- three target-free Krueger 60-second forecasts are still advancing through
  resumable checkpoints;
- the Bosch reactor-to-depth model passes its preregistered absolute
  calibration tolerances, but it has correctly refused the heldout reveal
  because it still loses to simple leave-one-lot-out baselines;
- the Bosch residual record identifies a missing dynamic chamber/tool-memory
  closure and measured impedance response rather than a need for another
  arbitrary depth multiplier;
- the goal is not complete: there is not yet a sealed heldout profile match for
  Freddie, a final Krueger depth match, or an independently validated third
  chemistry using the unchanged core.

The short scientific verdict is therefore:

> The computational path is substantially healthier than the old handoff
> suggests, and it now supports genuine blind conditional forecasts. It has
> not yet earned an atomic-accuracy or general absolute-depth claim.

## Mission and completion standard

The active goal is:

> Achieve defensible absolute-depth prediction rather than a one-target fit:
> identify and implement missing reactor and surface closures, validate them
> against independent measurements, and demonstrate heldout depth/profile
> accuracy for Krueger C4F6/Ar/O2 SiO2 plus at least two additional chemistries
> using the unchanged multiphysics core.

Completion requires all of the following:

1. Numerically certified feature evolution: convergence, conservation,
   topology invariants, material extinction, and mesh/refinement checks.
2. A wafer boundary derived from machine observables or independently measured
   flux/IEDF data, with uncertainty propagated rather than hidden.
3. Surface coefficients fixed from independent beam/blanket/material data, not
   chosen from the target SEM or target depth.
4. A frozen prediction written before target reveal.
5. Heldout depth and profile metrics that beat preregistered empirical
   baselines.
6. Transfer across at least three chemistry/material systems without changing
   the core solver merely to rescue a target.

Atomic-level numerical resolution is not the same as atomic-level predictive
accuracy. The latter also needs atomic-scale surface laws and exact-run chamber
boundary data. No report should conflate those two claims.

## Repository safety and version-control state

Use only the repository and branch named at the top. Do not switch to the old
multiphysics checkout or historic branch mentioned in older conversations.

At the snapshot, HEAD and origin are synchronized. The only worktree dirt is
pre-existing user work at these paths:

- `results/curated/mixed_layer_feature_v1/ml20-bonds-12s.log`;
- `results/curated/mouth_equilibrium_probe_dx/`;
- `scratch_ignore_calc.py`.

Do not edit, stage, delete, clean, reset, or absorb those paths. There is no
need for destructive Git commands. All campaign-owned work through the Oxford
v6 acceptance and Bosch v2 residual audit is committed and pushed.

Recent authoritative commits:

- `fefd970` records the clean Oxford v6 acceptance trajectory and the Bosch
  within-lot sequence-drift evidence;
- `980800a` makes material-owner cleanup idempotent while preserving the exact
  material union;
- `35e7966` preserves the Oxford v5 diagnostic checkpoint, trace, and caches;
- `2a03b3a` adds the Bosch calibration-only residual discovery audit;
- `9c2c009` records the static Bosch wall-conditioning fit and refuses a seal;
- `8ccd278` captures the exact Oxford pre-failure state;
- `852df0b` carries the first bounded Bosch wall state through the reactor;
- `921e62a` freezes that closure before calibration.

The latest relevant focused run has 112 passing Oxford/Bosch tests. The last
recorded broad suite before the final Oxford owner-projection fix had 2,125
passing tests. A new full-suite run is still required after the Oxford v6 board
is retrieved; do not substitute the older broad-suite count for that final
gate.

## Live paid compute

Only Vast instance `48177892` belongs to this campaign.

- SSH: `ssh -p 17892 root@ssh6.vast.ai`
- GPU: RTX 3090, 24 GiB
- approximate last recorded price: `$0.2011/hour`
- Oxford tree: `/root/petch-4b656fd`
- Krueger tree: `/root/petch-d852a1f`
- environment: `/root/petch-venv`

Do not touch any other Vast instance. Do not destroy this instance until all
Oxford v6 results and all terminal Krueger outputs have been copied locally,
hash-verified, audited, committed, and pushed.

Use PID files and `ps -p`. Do not use `pgrep -f`: the previous handoff's
liveness probe matched its own command and produced a false-positive report.

## Track A — Oxford NPG80 TiO2/Cr blind square-pillar board

### Frozen experiment condition

- Oxford PlasmaPro NPG80 RIE;
- 55/5/1 sccm CHF3/SF6/O2;
- 30 mTorr;
- 150 W forward table RF;
- 20 C table temperature;
- 1,200 s process duration;
- 700 nm ALD TiO2 on fused silica;
- 45 nm Cr hard mask;
- square-pillar board, approximately 400 nm pitch;
- nominal widths 80, 120, 160, 200, 240, 280, and 320 nm.

Freddie's target SEM has not been used to select the profile coefficients. The
board is intentionally blind. Exact GDS dimensions, sample radius, exact-run
self-bias/electrode waveform, blanket TiO2 loss, and remaining Cr thickness
remain missing experimental metadata.

### Reactor and boundary model

The current path is not a single fitted etch rate. It includes:

- a conserved 67-species CHF3/SF6/O2 global chemistry state;
- electron kinetics and power closure;
- species-resolved positive-ion and neutral wafer fluxes;
- a collisional sheath/IEDF boundary;
- axisymmetric radial wafer transport;
- feature-scale angular transport, surface chemistry, charging, moving TiO2,
  and moving Cr-mask geometry.

The conditional central solution is approximately:

- total positive-ion flux: `1.457e19 m^-2 s^-1`;
- powered-electrode sheath: approximately 296 V;
- singly charged ion impact energy: approximately 299 eV.

Those are conditional on an absorbed-power/self-bias closure. `150 W forward`
does not uniquely specify absorbed plasma power or the achieved electrode
waveform. The board therefore propagates high/low ion-energy, angular-tail,
and TiO2:Cr-selectivity scenarios. These are named physical uncertainties, not
values chosen to land on the SEM.

### What the apparent “spin” really was

The old failure was a chain of four numerical issues exposed by the increasingly
physical two-material board:

1. A marching-cubes sliver received a finite integrated rate divided by a tiny
   triangle area, requesting millions of CFL substeps.
2. Complete Cr-mask extinction was an undefined lifecycle event.
3. Independently redistanced material regions could manufacture unsupported
   owner nodes while preserving the apparent union.
4. Local repair followed by global redistance allowed one-cell Cr ownership to
   walk around the TiO2/Cr seam rather than reach a fixed point.

The fourth issue was captured at exactly `719.8619631901854 s`, after 149
accepted steps. Cleanup node counts evolved
`393 -> 172 -> 22 -> 1 -> 22 -> 9 -> 11 -> 4 -> 4 -> 2 -> 3 -> 5`, proving
that merely raising an iteration limit would hide a mesh-dependent erosion
operator.

The v6 fix freezes the locally repaired regional owner before redistance, then
projects redistanced regional fields so the desired owner retains the exact
authoritative union value and all competitors are strictly below it. The union
is bitwise unchanged and the owner update is exact and idempotent.

Regression evidence:

- `tests/data/zhu_v5_prefailure_8ccd278.npz` reproduces the old defect locally;
- the corrected checkpoint advances one exact step on CPU and CUDA;
- repaired owner nodes: 393;
- corrected fingerprint: `d53c...` in the committed test receipt;
- the former failing full trajectory completes cleanly under v6.

### Clean v6 acceptance result

Committed trajectory:

`results/curated/zhu_npg80_moving_cr_profiles_v1/trajectories/`
`w200_s18.017_ion_high_tail_0p0_c904f1b4f4500a0e.json`

Frozen cell:

- width 200 nm;
- TiO2:Cr selectivity 18.016664;
- high-energy, zero-tail ion scenario;
- requested duration 1,200 s;
- model revision `two-material-moving-tio2-cr-owner-projection-v6`.

The two surface-rate endpoints reach the same dose-limited final geometry by
different terminal paths:

- 34.125 nm/min endpoint: requested duration, 1,200 s accepted;
- 43.4667 nm/min endpoint: domain gas breakthrough after 942.101 s equivalent;
- etched depth: 679.407775 nm;
- middle CD: 196.542208 nm;
- bottom CD: 203.951939 nm;
- sidewall angle: 81.4639 degrees;
- bow: 94.566238 nm;
- reference Cr-mask extinction time: 724.693 s;
- particle transport balance residual: zero;
- remap residual: `9.37e-16`.

This is a successful conditional trajectory, not yet an experimental match.
Its receipt explicitly keeps
`parameter_evidence_supports_prediction = false` because the exact target-tool
surface coefficients and electrode waveform are not independently measured.

### Current v6 board run

At 2026-08-21 13:00 UTC:

- parent PID file: `/root/zhu_v6_board_980800a.pid`;
- parent PID: 27994;
- command: `scripts/audit_zhu_npg80_moving_cr_profiles.py --write
  --transport-device cuda:0 --workers 4`;
- log: `/root/zhu_v6_board_980800a.log`;
- parent state: alive;
- four child workers: alive, each near 100% CPU;
- one v6 cache exists: the previously accepted cell;
- the other 55 exact cells are being recomputed from scratch.

The workers spend most time in CPU feature evolution and call CUDA transport in
bursts. A single `nvidia-smi` snapshot can therefore show 0% GPU while all four
workers are healthy. That is not itself a stall. Judge progress by exact child
PIDs, cache creation, logs, and wall time, not instantaneous GPU utilization.

The local trajectory directory also contains 55 v3, 56 v4, and 39 v5 files.
They are forensic caches only. The v6 audit must not mix or grade them.

### Oxford landing sequence

1. Monitor PID 27994 and its four exact children without restarting them.
2. If the parent exits nonzero, preserve the log and any v6 caches before
   changing code.
3. On successful exit, copy `audit.json`, all 56 v6 trajectory caches, logs,
   and rendered artifacts to a temporary local path first.
4. Verify 56 unique preregistered specs, nonempty profiles, revision v6, and
   hashes before replacing or adding any curated files.
5. Run the audit locally in check mode.
6. Render and inspect the complete blind atlas.
7. Run focused tests and the full suite.
8. Commit and push the complete frozen prediction before target SEM reveal.
9. Only then score depth, top/mid/bottom CD, taper, bow, mask survival, bottom
   clearance, and spatial variation against Freddie's SEM.

## Track B — Krueger C4F6/Ar/O2 SiO2 depth

### What is known

Krueger reports approximately 825 nm feature depth. The earlier website story
claiming roughly 790–811 nm was retracted because it came from two cancelling
implementation errors. It must not be resurrected.

The pure-published-input feature result is approximately 346.833 nm. Direct
beam data constrain the physical removal yield below what is required to reach
825 nm using Krueger's published ion flux. That remains a real boundary-data
inconsistency: the paper's blanket rate and reported flux imply an effective
surface yield above the independently measured ceiling.

The current runs are a stricter transfer test. Their surface laws were fixed
from independent Guo beam/planar evidence before target scoring. Three bounded
species-allocation hypotheses are run without selecting the one closest to 825
nm:

- `nominal_unresolved`;
- `all_cf2`;
- `all_cf3`.

### Live Krueger state

At the snapshot, all supervisor PID files are alive and their logs continue to
advance through resumable half-hour segments:

| case | PID | elapsed process time | depth prefix | mask opening prefix |
|---|---:|---:|---:|---:|
| nominal unresolved | 10284 | 23.297 s | 210.779 nm | 20.299 nm |
| all CF2 | 10285 | 22.203 s | 235.923 nm | 25.940 nm |
| all CF3 | 10286 | 21.891 s | 240.383 nm | 23.148 nm |

Each target duration is 60 s. These prefixes must not be linearly extrapolated,
ranked by proximity to 825 nm, or reported as final depths. Geometry, mask
opening, transport, and surface state evolve nonlinearly.

Output directories:

- `/root/krueger_guo_60s_nominal`;
- `/root/krueger_guo_60s_cf2`;
- `/root/krueger_guo_60s_cf3`.

Each currently contains a resumable checkpoint, a large evolving `audit.json`,
and profile/metric images. Let all three finish. Retrieve and score all three as
a declared envelope. Do not select one after target comparison.

### Honest Krueger verdict

Krueger depth has not yet been matched. The current work may close the gap, or
it may strengthen the conclusion that the paper's missing/reported boundary is
incompatible with direct surface measurements. Either outcome is scientific.
If all no-fit transferred forecasts finish well below 825 nm, the next valid
closure is an independent reactor boundary or blanket-rate normalization with
uncertainty; the invalid move is increasing a surface yield until the feature
lands at 825 nm.

## Track C — Bosch SF6/C4F8 reactor-to-depth transfer

### What has succeeded

The reduced/cylindrical SPTS Bosch stack uses measured process traces, bounded
absorbed-power and transport closures, radial/azimuthal equipment response,
fused silicon/oxide surface recurrence, and unchanged surface laws.

The static wall-conditioning v5 calibration passes every preregistered absolute
calibration gate on 75 wafers:

- silicon mean MAE: 0.375596 um;
- silicon mean MAPE: 0.860673%;
- silicon point RMSE: 0.841164 um;
- normalized shape RMSE: 1.654573%;
- oxide mean MAE: 0.039568 um;
- selectivity MAPE: 5.993638%.

This is useful absolute-depth performance, but it is not enough for a prediction
seal.

### Why the seal was refused

Simple leave-one-lot-out empirical baselines are still better:

- global-mean silicon MAE: 0.338486 um;
- mean-map point RMSE: 0.486585 um;
- mean-map normalized shape RMSE: 0.636619%.

The fitted conditioning multipliers are only about 0.9995–1.016. They reproduce
static conditioning classes but do not explain the dominant sequential drift.
Therefore:

- `heldout_outcomes_read = false`;
- `heldout_prediction_written = false`;
- `eligible_for_prediction_seal = false`.

This refusal is correct. Passing absolute in-sample tolerances while losing to a
trivial heldout baseline is not evidence of a useful reactor predictor.

### New physical discovery

All eight calibration lots show declining measured silicon depth within the
lot:

`-0.1232, -0.1367, -0.1194, -0.1605, -0.1138, -0.1125, -0.1116,
-0.0558 um/wafer`.

The static physics model predicts nearly flat or rising sequences:

`+0.0217, +0.0204, +0.0072, -0.0107, +0.0467, +0.0280, +0.0638,
+0.0526 um/wafer`.

After removing lot identity and linear wafer sequence, the silicon residual is
still strongly related to measured platen/tool electrical channels:

- SF6 platen Vpp median: partial Pearson 0.58836;
- SF6 platen reflected power mean: 0.56150;
- SF6 platen Vpp RMS: 0.55251;
- C4F8 platen Vpp q90: 0.54494.

Gas-dose variation is tiny across calibration wafers, while Vpp/reflected-power
variation retains explanatory signal. This is evidence for two missing pieces:

1. a dynamic chamber wall/seasoning state carried from wafer to wafer; and
2. an impedance/sheath-transfer closure that uses the measured waveform and
   reflected-power response rather than forward power alone.

It is not evidence for a wafer-number depth offset. Wafer number is a proxy for
state evolution, not a physical input.

### Correct next Bosch iteration

Before fitting again, preregister a minimal exact chamber-memory recurrence.
A suitable first closure is a bounded scalar fluorocarbon wall coverage `q`:

```text
dq/ds = k_dep * D_C4F8 * (1 - q) - k_clean * D_SF6 * q
```

Use the analytic interval solution to carry `q_start`, `q_mean`, and `q_end`
through each wafer in sequence. Let `q_mean` modulate only a physically named
neutral wall-loss or sticking term in the first test. Do not change the ion or
surface laws simultaneously.

Because C4F8 and SF6 doses vary very little in this dataset, keep the parameter
count minimal and audit identifiability. A four-parameter coverage model can be
degenerate; prefer a reparameterized equilibrium coverage, relaxation rate,
initial state, and one bounded log-response coefficient, or freeze one from
literature/conditioning evidence.

Fit only calibration lots with whole-lot leave-one-out folds. Preserve the
heldout firewall. If a dynamic neutral wall state closes sequential drift but a
Vpp residual remains, preregister a separate measured impedance/sheath-transfer
closure in the next iteration. Do not introduce both at once and call the
mechanism identified.

## Track D — generalization and other chemistries

The architecture is generalizable at the solver level:

- chemistry mechanisms are data/deck driven;
- reactors produce species-resolved wafer boundaries;
- the sheath and feature transport operate on species/energy/angle
  distributions;
- surface state and moving material ownership are modular;
- deterministic ensembles are parallel over reactor settings, chemistry
  hypotheses, uncertainty corners, and geometries.

That does not mean arbitrary chemistry is predictive merely by changing a
JSON. Each new material/chemistry still needs:

- electron-impact cross sections and gas reaction rates;
- wall-recombination/sticking data;
- negative-ion and attachment physics where relevant;
- ion/neutral energy and angular distributions;
- energy- and species-resolved surface yields;
- passivation/polymer deposition and removal laws;
- mask response;
- independent reactor and surface validation data.

The current independent systems are strongest in:

1. Krueger fluorocarbon SiO2 feature transport and Guo beam-transfer surface
   evidence;
2. Oxford CHF3/SF6/O2 TiO2/Cr reactor-to-feature conditional forecasting;
3. SPTS SF6/C4F8 Bosch silicon/oxide reactor-to-depth calibration and sealed
   heldout protocol.

They share the core multiphysics architecture, but only the final two completed
heldout grades can establish the requested cross-chemistry claim. At present,
Oxford lacks the unrevealed SEM score and Bosch has not earned its reveal.

## What data would most improve predictive accuracy

For Oxford/Freddie, highest value in order:

1. exact-run DC self-bias or electrode voltage/current waveform;
2. blanket TiO2 loss from the same recipe and remaining Cr thickness;
3. exact GDS or nominal width/pitch mapping and sample radius/orientation;
4. cross-section and top-down SEMs with scale bars, tilt, and whether Cr was
   removed before imaging;
5. multiple wafer radii or coupons for radial response;
6. repeat runs to separate deterministic radial response from run-to-run drift.

For Krueger, highest value:

1. measured total ion flux for the reported reactor and condition;
2. species-resolved ion flux and IEDF/IAD;
3. blanket etch rate under the exact feature condition;
4. in-feature floor flux/energy or a second geometry sweep;
5. mask-loss history.

For Bosch, highest value:

1. exact wafer ordering and chamber conditioning/clean events;
2. matched V/I probe waveforms and reflected-power traces;
3. wall optical/emission or mass-spectrometry proxy for FC coverage;
4. a deliberate seasoning/cleaning perturbation sequence;
5. heldout lots kept sealed until the dynamic model beats baselines.

Those measurements are not admissions that the simulator “cannot do physics.”
They are the boundary and surface state needed to select one physical solution
from many solutions consistent with forward power, pressure, and flow alone.

## Immediate takeover checklist

1. Confirm the branch and protect the three user paths.
2. Monitor Oxford PID 27994 and its children; do not restart a healthy board.
3. Monitor the three exact Krueger PID files and logs; do not extrapolate
   partial depths.
4. While compute runs, write and hash the Bosch dynamic-wall v6
   preregistration before implementation or fitting.
5. Implement only the frozen recurrence and its conservation/sequence tests.
6. Run calibration-only whole-lot cross-validation; keep heldout sealed.
7. Retrieve and certify Oxford v6 immediately after completion.
8. Retrieve all three Krueger terminal audits after completion.
9. Run focused and full test suites.
10. Commit and push each complete scientific checkpoint.
11. Destroy instance 48177892 only after all required artifacts are local,
    hash-verified, committed, and pushed.

## Claims allowed now

Allowed:

- the engine implements a deterministic reactor-to-wafer-to-moving-feature
  multiphysics path;
- the old Oxford apparent spin was diagnosed and fixed as a certified
  material-topology defect;
- a clean v6 Oxford conditional trajectory now completes with conservative
  ledgers;
- the full blind Oxford board is running;
- Krueger transferred no-fit forecasts are running;
- Bosch reaches sub-percent mean-depth calibration error while correctly
  refusing heldout due to baseline failure;
- Bosch data expose a repeatable dynamic chamber/electrical-response signal;
- the target SEMs and heldout outcomes have not been used to tune the active
  prediction paths.

Not allowed:

- “we matched Freddie's SEM”;
- “we matched Krueger's 825 nm”;
- “the simulator is atomically accurate”;
- “any chemistry works by changing the deck”;
- “the Oxford board is complete”;
- “the Bosch heldout passed”;
- the retracted 790–811 nm website depth claim;
- choosing a Krueger species hypothesis or Oxford uncertainty corner because
  it is closest to the target.

## Bottom line

The project is not in a mysterious failure state. It is in the middle of three
well-separated scientific gates:

1. Oxford: the topology engine is repaired and the clean blind v6 board is
   running;
2. Krueger: the frozen surface-transfer hypotheses are only about one third of
   the way through their full nonlinear process time;
3. Bosch: absolute calibration works, but dynamic chamber memory and measured
   impedance response must be added and cross-validated before heldout reveal.

The right takeover behavior is to preserve those firewalls, finish the live
runs, and improve only closures supported by independent observables. The
wrong behavior is to treat a conditional 679 nm Oxford trajectory, a partial
240 nm Krueger prefix, or an in-sample Bosch pass as the requested general
absolute-depth breakthrough.
