# Multiphysics / absolute-depth progress state — 2026-08-20

Status timestamp: 2026-08-20 15:23 EDT

Repository: `plasma-etching-code`

Branch: `codex/validation-first-multiphysics`

Baseline pushed checkpoint inspected: `bdabd8a`

This is the detailed scientific record for the Oxford/Freddie moving-mask work
and the broader reactor-to-feature absolute-depth mission.  For the concise
operational takeover sequence, read
`CODEX_TAKEOVER_ABSOLUTE_DEPTH_2026-08-20.md` first.  This file supersedes
`HANDOFF_MOVING_CR_BOARD_2026-08-20.md` and updates
`CODEX_TAKEOVER_MOVING_CR_2026-08-20.md`.  The older files remain useful as a
forensic history, but their live-process instructions must not be followed.

## Executive verdict

The project has a real deterministic full-stack architecture:

`machine condition -> conserved reactor chemistry -> sheath and radial wafer transfer -> species/energy/angle-resolved boundary -> material-routed surface kinetics -> evolving 3-D profile`.

The Oxford board exposed two numerical defects and one missing lifecycle
capability.  The numerical defects are now reproduced, fixed, regression
tested, committed, and pushed.  The material lifecycle is now implemented and
tested.  A complete corrected v4 board has landed as a forensic/conditional
receipt.  The exact v5 acceptance trajectory has now continued through Cr
extinction to its physical/requested endpoints, and the full v5 board is live.

No target SEM, target depth, or target-derived coefficient has been used.
Therefore the Oxford calculation remains a blind conditional prediction.  It
is not yet an experimentally graded final SEM and it is not an atomic-accuracy
claim.

The fundamental remaining science limitation is separate from these numerical
repairs: Oxford supplies generator power but not the achieved electrode
waveform/self-bias, and direct state-resolved TiO2/Cr response coefficients for
this CHF3/SF6/O2 boundary are not published.  The current absolute surface rate
and selectivity axes are transferred from independent cross-machine witnesses.
That supports a physically constrained envelope, not a unique target-tool
absolute answer.

## Frozen Oxford blind condition

The supplied target condition is:

- Oxford PlasmaPro NPG80 RIE;
- `55 / 5 / 1 sccm CHF3 / SF6 / O2`;
- `30 mTorr`;
- `150 W` forward table RF;
- `20 C` table temperature;
- `1200 s` process duration;
- `700 nm` ALD TiO2 on fused silica;
- `45 nm` Cr hard mask;
- square-pillar prior, approximately `400 nm` pitch;
- preregistered width board `80, 120, 160, 200, 240, 280, 320 nm`.

The exact GDS, sample radius/orientation, achieved DC self-bias, same-run
blanket loss, and target SEM remain outside the model.  This is useful: the
eventual SEM can be a genuine answer key rather than a fitting target.

## Current reactor-to-feature state

The strongest current Oxford boundary is a conditional 67-species conserved
reactor/sheath solution with deterministic radial wafer transfer.  Its central
optic prediction carries:

- 20 positive-ion species;
- 37 thermal-neutral species;
- absolute species-resolved particle fluxes;
- ion masses and charge states;
- charge-resolved impact energies;
- deterministic angular core/tail quadrature;
- central positive-ion flux approximately `1.457e19 m^-2 s^-1`;
- conditional powered-electrode sheath drop approximately `296 V`;
- singly charged mean impact energy approximately `299 eV`.

Dominant ions in that conditional state are approximately `CF3+ 75.3%`,
`CF+ 9.25%`, `H+ 5.54%`, `H2+ 4.06%`, `O2+ 2.50%`, and `CF2+ 2.12%`.
The interface does not collapse those populations into an anonymous ion.

The axisymmetric model predicts only a tiny smooth radial variation across the
central few-millimetre optic.  Therefore clustered fallen pillars in a local
SEM are evidence of spatial variation but not, by themselves, proof of a
smooth reactor-flux gradient.  Local pattern loading, mask adhesion/undercut,
sample placement, micro-masking, strip/rinse/dry mechanics, or unresolved
fine-scale tool nonuniformity remain plausible.

## What went wrong in the moving-Cr run

### 1. Subcell material cleanup was not a fixed point

One-pass cleanup removed a Cr island, redistancing exposed a new adjacent
one-node island, and the next mesh failed.  The same physically declared rule
now iterates to a fixed point.  It removes only components that own no complete
hexahedral volume cell.  Any resolved component still reaches the hard
topology gate.

### 2. Periodic strip projection amplified tiny triangles

The v3 strip symmetrizer divided an event rate equally among equivalent
triangles and then divided each share by recipient area.  A marching-cubes
sliver therefore received a near-infinite flux density.  On the exact failing
step, recession velocity jumped from about `7.25e-4` to `3310` mesh units/s and
the CFL controller requested `5,331,155` advection substeps instead of 2.

Commit `d8a8d02` projects uniform flux density instead.  Event rate is
distributed in proportion to recipient area, conserving total rate without a
small-area singularity.  The exact formerly failing trajectory then completed
with two CFL substeps throughout.

### 3. The driver confused first local mask loss with the final process

All old v3 and corrected v4 trajectories stopped when Cr became thinner than
one vertical cell somewhere in the original footprint.  That is a useful
physical/numerical event, but it is not the requested 20-minute endpoint.
Continuing required an explicit material lifecycle rather than copying the
terminal geometry forward.

## Complete corrected v4 board

Commit `30c3054` lands all 56 v4 trajectory receipts plus `audit_v4.json`.
The remote authority passed its native `--check`; the retrieved tarball hash
matched byte-for-byte between local and remote:

`1fdad5f88f83c3e15704193968673de48358830c5c44483aab1c4648b03063ee`.

The v4 board contains:

- 56 complete trajectory files;
- 112 rate-scaled endpoints;
- all 112 terminal reasons equal to
  `cr_mask_below_vertical_resolution_in_footprint`;
- accepted process-equivalent time range `391.334–713.846 s`;
- conditional etched-depth range `274.896–406.087 nm`;
- maximum particle-balance error exactly `0`;
- maximum surface-state remap relative residual `1.683e-15`.

Interpretation: the board is a complete, numerically corrected certificate of
the profile envelope up to first local Cr loss.  It is not a final 1200-second
profile board.  V3 and v4 caches remain immutable forensic receipts and are not
eligible as v5 results.

## New certified material-extinction lifecycle

Commit `84c2c01` adds the common-engine lifecycle used by v5.

A material may retire only when all of the following hold:

1. no geometry node is owned by that material;
2. its nonnegative level-set support contains no complete eight-corner volume
   cell;
3. any remaining positive level-set value is bounded by an explicit
   floating-point roundoff tolerance;
4. its periodic material-component count goes from positive to exactly zero;
5. solid components, gas cavities, domain breakthrough, and every surviving
   material-component count remain invariant;
6. the caller explicitly selects
   `continue_gas_cavity_and_material_extinction`.

The exact Oxford replay reported registered materials `[1, 2, 3]` but surviving
owners `[1, 3]`.  The disappearing Cr layer had:

- minimum level-set value `-0.05` mesh units;
- maximum `1.937e-14` mesh units;
- 111 roundoff-nonnegative nodes;
- no owned Cr nodes.

The maximum is approximately `2e-12` of one 10 nm cell.  Treating its sign as
physical chromium would be a floating-point artifact.  The v5 gate additionally
requires zero resolved volume cells, so a real hidden Cr island cannot pass.

At an accepted extinction event:

- retired Cr fields are integrated over the final Cr surface and written to an
  explicit retirement ledger;
- those fields are never transferred to TiO2;
- surviving TiO2 history remains material-local;
- newly exposed TiO2 receives the mechanism's declared fresh-surface state;
- common-refinement overlap, rather than nearest-neighbour borrowing, is used
  for the extinction step;
- the configured fast indexed remap resumes on later ordinary steps.

Default and gas-cavity-only policies still refuse material loss.  Above-roundoff
hidden support and any resolved hidden volume also refuse.  The focused
lifecycle/remap/feature group is `138 passed`; the v5 driver/core group is
`118 passed`.

## Passed v5 exact continuation

The exact hard cell completed successfully:

- width `320 nm`;
- selectivity `14`;
- low energy (`146.539 eV`);
- zero angular-tail fraction;
- requested duration `1200 s`;
- 10 nm grid;
- model revision `two-material-moving-tio2-cr-dose-factorization-v5`;
- log `/root/w320_v5.log`;
- wall time `1405.950 s`;
- cache
  `w320_s14.000_ion_low_tail_0p0_593697ba778b8990.json`.

Its retrieved cache hash matches the remote receipt:

`4feeec492c52af860298495bd190e50106370306d208acd8fe24a7170cd1647f`.

The Cr layer retired once at reference time `560.429 s`.  Its maximum residual
level-set value was `1.937e-14` mesh units, the roundoff gate was
`1.137e-13`, and the support contained zero resolved volume cells.  The event
used common refinement, retired 2,284 final Cr surface faces, and resumed
ordinary indexed remap afterward.  Particle balance remained exactly zero and
the maximum remap residual was `7.166e-16`.

The two no-target-fit conditional endpoints are:

- `34.125 nm/min` surface-rate axis: requested `1200 s` reached, predicted
  etched depth `646.856 nm`, middle/bottom CD `313.08/321.78 nm`, and sidewall
  angle `76.03 degrees` from the wafer;
- `43.4667 nm/min` axis: physical domain-gas breakthrough at
  `951.764 s`, preceding the requested 1200 s, with last accepted etched depth
  `653.356 nm`, middle/bottom CD `312.05/321.73 nm`, and sidewall angle
  `76.17 degrees`.

The zero top-CD metric after mask loss means the original top reference no
longer intersects a surviving pillar; it must not be described as a normal
trapezoidal top CD.  The eventual SEM grade must compare the full profile and
terminal class, not cherry-pick one CD.

## Live full v5 board

The full 56-cell v5 campaign is now running with eight deterministic workers:

- parent PID `8980`;
- log `/root/zhu_v5_board.log`;
- pid file `/root/zhu_v5_board.pid`;
- the passed exact cell is content-addressed and reused;
- the remaining 55 cells are new v5 computations.

The live box is Vast instance `48177892`, RTX 3090, at
`$0.2011111111/hour`.  Repository tree is `/root/petch-4b656fd` and the venv is
`/root/petch-venv`.  Only this project instance is in scope; unrelated Vast
instances must not be touched.

The first v5 trajectory proved:

- Cr extinction is accepted exactly once;
- newly exposed TiO2 remap balances;
- subsequent state contracts contain no retired Cr fields;
- the process advances beyond the old `~396 s` endpoint;
- the next physical terminal or the full requested duration is reached;
- no result is silently frozen at material loss.

After useful v5 artifacts are pulled and locally verified, instance `48177892`
must be destroyed and disappearance confirmed.  Do not destroy it while the
exact continuation or a deliberately launched v5 board is active.

## Surface coefficients: why the Oxford answer is still conditional

The software interface is general; the numerical coefficient evidence is not
yet target-complete.  The current TiO2/Cr board is rate-normalized from
independent process witnesses.  Direct literature supports the required model
topology—fluorination, oxygen blocking/cleanup, ion-assisted product desorption,
passivation growth/removal, Cr mask recession—but does not identify one unique
state-resolved Oxford deck from Freddie's recipe screen.

Still missing for a unique target-tool absolute prediction are independent
constraints on:

- achieved DC self-bias or electrode voltage/current waveform;
- absorbed plasma power rather than only forward generator power;
- same-run blanket TiO2 removal;
- same-run remaining Cr thickness/selectivity;
- species/energy-resolved bare and fluorinated TiO2 removal yields;
- fluorination and oxygen blocking/cleanup probabilities;
- fluorocarbon passivation sticking, density, and sputter yield;
- Cr response under the same mixed chemistry;
- actual mask/GDS dimensions and sample radius.

This is not a request for the target profile as a fitting input.  The best
minimal same-run calibration set is self-bias + blanket TiO2 loss + remaining
Cr, followed by a blind profile grade.

## Krueger C4F6/Ar/O2 depth

Krueger is not presently depth-matched under the paper's published aggregate
boundary.  The honest value remains approximately `346.833 nm` versus the
reported `825 nm`.  The old apparent `790–811 nm` match was retracted because
two implementation errors canceled.

There is now a promising independent surface-transfer route:

- the Guo/Kwon no-target-fit planar transfer predicts effective yield `2.613`
  versus the `2.521` needed by the Krueger depth, within about 5%;
- finite-fluence planar sensitivities land near `855 nm`;
- a deterministic-extruded feature prefix predicts `11.93–12.27 nm/s`, about
  `11–13%` below the `13.75 nm/s` feature-average requirement;
- no Krueger depth was used to solve those surface parameters.

It is not yet an authorized validation result because much of Krueger's IEAD
extends beyond the Guo fit domain, the source does not publish the required
species-resolved ion composition/C4F6 parent flux, and neutral uptake remains
uncertain.  The next high-value Krueger computation is a full deterministic-
extruded evolving feature forecast using the existing
`GuoC4F8ArSiO2FeatureMechanism`, explicitly labeled as a transfer sensitivity.
The scientifically strongest closure remains a validated C4F6 reactor boundary
or the authors' species-resolved HPEM/PCMCM wafer outputs.

## Other chemistry/depth validation

Cross-chemistry evidence currently includes:

- direct/planar `SF5+ -> Si/SiO2` surface-dose points at about `5.88% MAPE`;
- `Cl2/Ar+ -> Si` ALE points with about `12.88%` maximum error;
- a Mahorowala Cl2 feature board near `19% MAPE`, without a formal pass gate;
- direct beam agreement for the principal fluorocarbon surface law at feature
  energy, including the existing `4.7%` grade.

These are meaningful physics checks, but the repository still lacks two
additional completed formal held-out **feature-profile/depth** predictions with
the unchanged core.  That remains part of the absolute-depth mission.  Surface
points must not be relabeled as feature validation.

## What is genuinely general

The common engine is not hard-coded to TiO2 or Oxford.  A new process supplies:

1. gas species/reaction/collision data;
2. machine geometry, pressure/flow/power coupling, walls, and diagnostics;
3. sheath and wafer-transfer evidence;
4. species-resolved material surface laws and density;
5. mask/substrate laws and geometry;
6. untouched measurements for grading.

The deterministic transport, charging, conservative state update,
multi-material routing, topology gates, level-set motion, remap contracts, and
implicit/differentiable interfaces remain shared.  Arbitrary chemistry is
possible at the interface, but unsupported chemistry does not become
predictive merely because a JSON deck can name it.

## Atomic-accuracy boundary

No atomic-accuracy claim is supportable today.

The engine conserves atom/material-unit inventories and can use atomistically
derived barriers or yields, but the production feature grid is 10 nm, surface
states are reduced kinetic populations, and several target coefficients are
not measured.  DFT/MD can reduce selected surface-law uncertainty; they cannot
infer the exact chamber self-bias, wall condition, radical flux, or wafer
nonuniformity from generator power and flows alone.

The defensible ambition is a physics-closed, uncertainty-propagating,
experimentally graded knobs-to-profile predictor whose errors approach the
measurement floor.  “Atomic level” should mean atom-balanced mechanisms and
evidence-bearing coefficients, not a false claim that every final atom is
known from an incomplete recipe.

## Required next sequence

1. Monitor v5 board parent PID `8980` and count only v5 content-addressed
   caches; old v1-v4 files are not completion evidence.
2. Assemble/check/pull the v5 audit and every v5 trajectory, then run the local
   native check against identical code.
3. Rerun the full local suite, commit and push the complete board.
4. Destroy Vast instance `48177892` after every needed artifact is safe.
5. Run the full deterministic Guo/Krueger feature forecast and grade its domain
   of validity without selecting from the 825 nm target.
6. Promote at least two additional chemistries from surface evidence to formal
   held-out feature-depth/profile gates.
7. Ask Freddie for self-bias, blanket TiO2 loss, remaining Cr, exact GDS/sample
   position, rinse/dry history, and then the SEM answer key.

## Version-control and validation state

Pushed checkpoints relevant to this handoff:

- `44956d0` — advance Oxford board beyond Cr mask loss;
- `30c3054` — land complete corrected v4 moving-Cr receipts;
- `84c2c01` — continue through certified material extinction;
- `4ceda55` — expose material lifecycle mismatch diagnostics;
- `f65a445` — add no-cache post-mask exploration;
- `d8a8d02` — fix unequal-area strip flux amplification.

The full repository suite with the new lifecycle produced `2131 passed, 7
skipped` plus one expected stale v4 cache-filename assertion.  Commit `242177b`
updates that assertion to the v5 content hash; its debug/driver sentinel group
is `9 passed`.  A final full suite rerun remains mandatory after the complete
v5 board lands.

The only untracked local paths are pre-existing unrelated user artifacts and
remain untouched:

- `results/curated/mixed_layer_feature_v1/ml20-bonds-12s.log`;
- `results/curated/mouth_equilibrium_probe_dx/`;
- `scratch_ignore_calc.py`.
