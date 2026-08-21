# Codex takeover report — reactor-to-feature absolute-depth campaign

Snapshot: 2026-08-21 10:00 EDT / 14:00 UTC

Repository: `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code`

Branch: `codex/validation-first-multiphysics`

Scientific HEAD at snapshot: `cd83816`

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
- the Bosch reactor-to-depth model now contains an identifiable, preregistered
  recipe-path wall-memory closure; it passes all absolute calibration gates,
  beats the leave-one-lot-out global-depth baseline, and improves within-lot
  drift, but correctly refuses heldout reveal because it still loses the
  radial point and shape baselines;
- a subsequent calibration-only decomposition localizes `86.2018%` of the
  squared Bosch normalized-map residual to one repeatable spatial tool map;
  the leading dynamic residual is an edge-positive, center-negative mode tied
  to measured C4F8-phase platen Vpp, so the next closure is now materially more
  specific than a generic "radial/electrical response" guess;
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
v6 acceptance, Bosch v7 calibration audit, and Bosch spatial-residual
localization is committed and pushed.

Recent authoritative commits:

- `cd83816` records the Bosch spatial-boundary residual decomposition and
  output-space capacity diagnostics without changing the physical model or
  opening the heldout;
- `2cbec26` records the identifiable Bosch v7 calibration/LOLO audit and exact
  thirteen-node reactor-to-surface response table;
- `abb209b` implements the fixed-recipe wall-memory operator and tests;
- `f784852` freezes the reduced v7 law and its validity domain before code;
- `19c7db8` records the useful but structurally unidentifiable v6 result;
- `e64d5bc` implements the exact bounded v6 coverage recurrence;
- `2612a07` freezes the v6 physical wall-coverage hypothesis;
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

The Bosch-focused suite at this snapshot has 62 passing tests. The last
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

At 2026-08-21 14:00 UTC:

- parent PID file: `/root/zhu_v6_board_980800a.pid`;
- parent PID: 27994;
- command: `scripts/audit_zhu_npg80_moving_cr_profiles.py --write
  --transport-device cuda:0 --workers 4`;
- log: `/root/zhu_v6_board_980800a.log`;
- parent state: alive;
- four child workers: alive, each near 100% CPU;
- eighteen unique v6 caches are complete;
- the other 38 exact cells are being recomputed from scratch.

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

At 14:00 UTC, all supervisor PID files are alive and their logs continue to
advance through resumable half-hour segments:

| case | PID | elapsed process time | depth prefix | mask opening prefix |
|---|---:|---:|---:|---:|
| nominal unresolved | 10284 | 24.516 s | 216.178 nm | 19.804 nm |
| all CF2 | 10285 | 23.375 s | 242.399 nm | 23.367 nm |
| all CF3 | 10286 | 23.063 s | 246.620 nm | 24.297 nm |

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

### Stack and target firewall

The reduced/cylindrical SPTS Bosch stack uses all measured process traces,
bounded absorbed-power and neutral-loss closures, species-resolved radial wafer
transfer, and the unchanged Belen silicon plus La Magna/Garozzo oxide/film
surface laws. Seventy-five wafers have 89-point measured depth maps; a 76th
processed wafer without an outcome is still carried through chamber history.

The chronological heldout outcome remains unopened. Every v5, v6, and v7
receipt states:

- `heldout_outcomes_read = false`;
- `heldout_prediction_written = false`;
- `eligible_for_prediction_seal = false`.

### What v6 learned and why it was rejected

Version 6 implemented the physically attractive bounded occupancy law

```text
dq/ds = k_dep*D_C4F8*(1-q) - k_clean*D_SF6*q
```

and carried the exact analytic mean/end state between wafers. It improved
calibration mean-depth MAE from 0.375596 to 0.233331 um and reduced the
within-lot slope error from 0.145388 to 0.083769 um/wafer. That proved dynamic
memory matters.

It did not identify its physics coefficients:

- `k_clean` landed exactly on its lower bound;
- wall response landed exactly on its upper bound;
- Jacobian condition number: `6.762e6`, above the frozen `1e6` limit;
- maximum parameter correlation: `0.999997`, above `0.995`;
- the nine-node interpolation audit reached `0.06091` of a frozen gate,
  above the `0.05` limit.

This is structural, not optimizer bad luck. Across all 76 calibration traces,
C4F8 dose CV is 0.002336, SF6 dose CV is 0.000501, and the SF6/C4F8 dose-ratio
CV is 0.002572. One nearly fixed recipe cannot separate deposition, cleaning,
and response rates. Version 6 therefore correctly refused a seal.

### Version 7 identifiable recipe-path closure

Version 7 was frozen before implementation in `f784852`, implemented in
`abb209b`, and audited in `2cbec26`. It retains only the combination supported
by the experiment: cumulative normalized C4F8 production exposure `H` since a
declared conditioning sequence:

```text
H_mean = H_start + 0.5*d_C4F8
H_end  = H_start + d_C4F8
```

The four shared coefficients are three static conditioning terms and one
bounded `b_H` response. The law changes only neutral upper/sidewall loss. It
does not change positive ions, lower-wafer collection, surface yields, or depth
directly. It has no fitted lot state and never uses wafer number, date, lot
number, or target depth.

This is deliberately a **recipe-path** model, not a universal deposition and
cleaning law. It refuses any process trace outside the measured ratio domain
`5.635257777059961 <= D_SF6/D_C4F8 <= 5.737394214801463`. Independent
varying-ratio experiments or wall diagnostics are required before that scope
can be expanded.

### Version 7 result

The identification problem is closed:

- fitted `b_H`: `0.00661551` log wall-loss per reference wafer;
- Jacobian rank: `4/4`;
- condition number: `14.518`;
- maximum pairwise parameter correlation: `0.60234`;
- no parameter bound contact;
- thirteen-node/twelve-midpoint interpolation maximum: `0.01626` of a frozen
  gate, below `0.05`.

Absolute calibration performance on 75 wafers is:

- silicon mean MAE: `0.232214 um`;
- silicon mean MAPE: `0.530134%`;
- silicon point RMSE: `0.789018 um`;
- normalized shape RMSE: `1.664987%`;
- oxide mean MAE: `0.039568 um`;
- selectivity MAPE: `6.046622%`.

Whole-lot leave-one-out performance is:

- silicon mean MAE: `0.285331 um`, beating the `0.338486 um` global-depth
  baseline;
- silicon point RMSE: `0.814735 um`, losing to the `0.486585 um` mean-map
  baseline;
- normalized shape RMSE: `1.669293%`, losing to the `0.636619%` mean-map
  baseline;
- within-lot depth-slope MAE: `0.082903 um/wafer`, improved from the v5
  `0.145388 um/wafer`.

This is genuine progress: shared physical state now transfers mean depth across
whole left-out calibration lots and is numerically identifiable. It is not a
prediction seal because radial map accuracy still loses to a trivial empirical
mean map. Exact selected-parameter replay, grid refinement, and heldout hashing
are intentionally deferred rather than performed after a known prerequisite
failed.

### Spatial residual localization after v7

Commit `cd83816` adds a calibration-only model-form audit. It does not modify
the reactor, wafer boundary, or surface law and it does not read the heldout.
Its main result is unusually decisive:

- the shared mean spatial residual contains `86.2018%` of the squared v7
  normalized-map residual;
- raw v7 normalized shape RMSE is `1.664987%`;
- subtracting that all-calibration shared map leaves `0.618475%` RMSE;
- a whole-lot output-space shared-map proxy reaches `0.633474%` shape RMSE
  and `0.405513 um` point RMSE, beating the frozen mean-map baselines of
  `0.636619%` and `0.486585 um` while retaining the reactor-predicted wafer
  mean;
- adding one standardized C4F8 platen-Vpp-RMS slope improves the proxy to
  `0.620340%` shape and `0.399671 um` point RMSE;
- that slope map is edge-positive and center-negative, with `0.804859`
  Pearson correlation to the frozen edge basis;
- the measured C4F8 Vpp-RMS domain is `626.9533--643.5343 V`, with mean
  `637.44096 V` and standard deviation `3.81631 V`.

This is evidence that the dominant missing term is a stable, machine-specific
wafer ion-transmission fingerprint, plausibly produced by chamber, electrode,
and focus-ring geometry. The smaller dynamic term is consistent with
Vpp-driven sheath curvature or edge focusing. It is not evidence for another
scalar etch rate, wall-loss multiplier, or surface-yield adjustment.

The capacity audit also sets an important warning. Complete real-Zernike
output-space bases do not narrowly beat the strict shape baseline until order
10 (`65` non-piston coefficients, `0.635259%` shape RMSE). That is enough to
show representational capacity, but it is not yet a physical reactor-to-wafer
closure and must not be presented as one.

### Next Bosch physics: version 8, not yet frozen

The next iteration should be preregistered as a positive, smooth,
current-conserving wafer-boundary operator applied only to the positive-ion
channel before the unchanged feature and surface recurrence:

1. a static tool fingerprint in the official wafer coordinate frame;
2. a low-order Vpp-dependent edge mode using only the measured C4F8-phase
   waveform statistic;
3. exact area/current normalization so total positive-ion current is
   unchanged;
4. strict coefficient, field-amplitude, voltage-domain, identifiability, and
   interpolation gates;
5. whole-lot refits with no heldout leakage and the same point/shape baselines;
6. no changes to v7 wall memory, neutral chemistry, or surface yields.

The portable science is the operator and conservation law. The fitted static
map is a calibration of this specific SPTS tool and coordinate frame. If the
minimum basis capable of passing the baseline is too high-order or unstable
under whole-lot refits, v8 must fail rather than laundering an empirical depth
map into plasma physics.

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

They share the core multiphysics architecture, but only completed heldout grades
can establish the requested cross-chemistry claim. At present, Oxford lacks the
unrevealed SEM score, Bosch has not earned its chronological reveal, and
Krueger's no-fit terminal forecasts are still running.

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
4. Treat Bosch v7 as a completed, useful non-seal; do not reopen the heldout or
   refit its wall law.
5. Freeze the current-conserving static-tool-map plus low-order Vpp edge-mode
   v8 hypothesis before implementing it; keep wall and surface laws unchanged.
6. Retrieve and certify Oxford v6 immediately after completion.
7. Retrieve all three Krueger terminal audits after completion.
8. Run focused and full test suites after live artifacts land.
9. Commit and push each complete scientific checkpoint.
10. Destroy instance 48177892 only after all required artifacts are local,
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
- Bosch v7 reaches 0.53% calibration and 0.285 um whole-lot mean-depth error
  with an identifiable shared recipe-path memory law;
- Bosch v7 beats the global-depth baseline while correctly refusing heldout
  because radial point/shape baselines still win;
- Bosch data localize 86.2% of squared map residual to one repeatable tool
  fingerprint and expose a smaller Vpp-dependent edge mode;
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
3. Bosch: identifiable dynamic recipe-path memory now transfers mean depth;
   the dominant remaining map error is localized to a stable tool fingerprint,
   but that fingerprint still has to be implemented as a conservative physical
   wafer-boundary operator and beat point/shape baselines before heldout reveal.

The right takeover behavior is to preserve those firewalls, finish the live
runs, and improve only closures supported by independent observables. The
wrong behavior is to treat a conditional Oxford trajectory, a partial Krueger
prefix, or Bosch's improved calibration/LOLO mean depth as the requested
general absolute-depth breakthrough.

## Source-semantics correction preserved in v6/v7

The Bosch source README was visually audited. In the lot labels, `C` means
conditioning on the bare system
**chuck**, not a carbon cycle. The v5 static calculation used the numerical
repeat count correctly, but its `log_carbon_cycle` variable and this report's
earlier wording were semantically wrong. The authoritative v6 and v7 records
use conditioning-repeat semantics and record the actual O2 then O2/SF6
conditioning sequence.
