# Multiphysics progress state — 2026-08-19

Branch: `codex/validation-first-multiphysics`

Current pushed checkpoint: `2b9bf62` (`Preserve species-specific surface response`)

This is the authoritative state after resuming the Oxford/Freddie and reactor-to-feature work. The unrelated untracked mixed-layer log and `mouth_equilibrium_probe_dx/` directory still predate this campaign and remain untouched.

## Executive verdict

The stack is now genuinely general at its principal software interfaces:

`machine condition -> conserved reactor state -> radial wafer transfer -> species-resolved deterministic feature boundary -> material-routed surface law -> evolving 3D profile`.

The Oxford/TiO2 work is one configuration and one validation target, not hard-coded physics. The boundary can carry arbitrary named ions and neutrals, species-specific masses, charge states, energy measures, angular distributions, and absolute fluxes. The oxide kernel can now carry a separate energetic-yield law for every incident ion species. The common level-set engine evolves arbitrary material IDs and the sub-grid topology repair applies to every moving-material campaign.

That does **not** mean the Oxford profile is already an absolute prediction. The target process still lacks direct measurements of the achieved self-bias/IEAD and the TiO2/Cr surface coefficients. The presently running profile board is a blind, cross-machine conditional ensemble. It is useful and physically constrained, but it is not atomic-accuracy certification and it is not fitted to Freddie's SEM.

Krueger is also not depth-matched. The honest published-boundary result remains about `346.833 nm` against the reported `825 nm`; the paper does not publish enough species-resolved wafer-boundary information to identify the missing factor without inventing it.

## Frozen Oxford condition

The supplied process is:

- Oxford PlasmaPro NPG80 RIE;
- `55/5/1 sccm CHF3/SF6/O2`;
- `30 mTorr`;
- `150 W` forward RF power;
- `20 C` table temperature;
- `1200 s` etch;
- `700 nm` ALD TiO2 on fused silica;
- `45 nm` Cr mask;
- blind geometry prior: square-pillar board, approximately `400 nm` pitch, widths `80--320 nm`.

The exact GDS and target SEM remain withheld/unavailable, which is appropriate for a blind prediction. The geometry prior should be replaced, not silently retained, when Freddie provides the actual dimensions.

## Reactor and wafer-boundary state

The strongest current Oxford reactor state is the wall-resolved sheath/global fixed point conditioned on:

- `90 W` absorbed-power sensitivity for the `150 W` forward setpoint;
- a same-Oxford-family `276 V` self-bias transfer;
- approximately `20.47 V` plasma potential;
- `296.47 V` powered-electrode sheath drop;
- conserved 67-species daughter chemistry;
- solved electron kinetics and particle/power balances;
- deterministic axisymmetric ion transfer to the wafer.

For the central 3 mm optic, that state predicts:

- positive-ion flux `1.457224881e19 m^-2 s^-1`;
- full-electrode/axisymmetric closure within about `0.36%`;
- final-grid central-flux change about `0.0025%`;
- smooth variation between the two modeled annuli intersecting the central optic only about `0.00031%`.

Therefore the smooth reactor-scale radial model does **not** explain strongly clustered pillar collapse inside a small field. Local loading, GDS/CD variation, mask adhesion/undercut, sample position, post-etch stripping/drying, or finer-scale plasma/tool nonuniformity remain viable. The existing SEM observation of spatially clustered collapse is evidence of spatial variation, but not by itself proof that incident ion flux is the cause.

The exact-run absolute boundary is still conditional because forward generator power is not absorbed plasma power and does not uniquely determine electrode voltage or self-bias. The current family transfer is far better than pretending `150 W` implies one ion energy, but it is not a measurement from Freddie's exact run.

## Species-resolved reactor-to-feature join

Commit `180a961` closes the software join without aggregating the plasma into one anonymous ion and one radical channel. The checksum-bound receipt is:

`results/curated/zhu_npg80_species_resolved_feature_boundary_v1/audit.json`.

It exports:

- 20 positive-ion species;
- 37 thermal-neutral species;
- 3,252 deterministic quadrature nodes per declared angular-tail case;
- ion mass and charge number;
- charge-resolved impact energy using `Z_i V_s + Te/2`;
- absolute particle flux and current density;
- two declared angular-tail sensitivities (`0` and `0.65`);
- no Monte Carlo sampling.

The dominant central ion fluxes in the conditional state are approximately:

- `CF3+`: `75.30%`;
- `CF+`: `9.25%`;
- `H+`: `5.54%`;
- `H2+`: `4.06%`;
- `O2+`: `2.50%`;
- `CF2+`: `2.12%`.

The neutral thermal-flux inventory is dominated by HF, followed by feed CHF3 and H2. This is scientifically important: a reactor may conserve dozens of species while only a subset has a known surface channel. The interface now preserves that fact rather than hiding it through aggregation.

Singly charged ions have conditional mean impact energy near `298.6 eV`; `SF2++` and `SF4++` receive about `595.1 eV`. Boundary particle flux closes to machine precision. The 2025 Kim core/tail widths are represented deterministically, but their tail fraction is not a target-tool measurement and is explicitly kept as a sensitivity.

## New common-engine repairs

### Fixed-point sub-grid material cleanup

The first production moving-Cr trajectory failed with an unmatched-edge/material-component error. Local reproduction showed the precise sequence:

1. Cr redistancing created a 14-node sub-cell component;
2. the existing physically defined cleanup removed it;
3. that redistance exposed a new adjacent one-node Cr island;
4. the old one-pass postcondition did not inspect the newly created neighbor.

Commit `4265e70` makes the same resolution rule iterate to a fixed point. It removes only components that do not own all eight corners of any physical hexahedral volume cell. Any component that owns a real cell is preserved and still reaches the hard topology gate. The exact Oxford `80 nm`, `10 nm`-mesh failure is now a regression test. The formerly failing first two steps pass locally.

This is a generic numerical fix in `feature_step_3d.py`, not an Oxford exception.

### Generic species boundary

The same commit added `species_resolved_feature_boundary.py`. It supports:

- arbitrary ion and neutral names;
- scalar or discrete per-species ion-energy measures;
- per-species masses and charge states;
- one common or species-specific two-component IADF;
- deterministic polar/azimuthal quadrature;
- analytic half-Maxwellian neutral flux distributions;
- ion-only or ion-plus-neutral boundaries;
- exact absolute-current reconstruction.

### Species-specific surface response

Commit `2b9bf62` fixes a second aggregation issue. The reduced fluorinated-oxide mechanism previously accepted one bare-oxide, one fluorinated-oxide, and one passivation sputter-yield curve and applied each curve to all energetic species. It can now accept a fail-closed map from ion name to yield law. A positive ion missing from any declared species map makes the surface step refuse execution; it is never silently assigned another ion's coefficient.

Legacy scalar decks retain their previous behavior and provenance format. Rate-normalized bridge mechanisms can also explicitly declare inert neutrals, allowing full reactor inventories to pass through conditional mask laws without pretending an undeclared neutral was consumed.

The focused validation after this change is `105 passed`; the narrower surface-kernel group is `47 passed`.

## TiO2 surface science

The surface topology is substantially clearer than it was at pause:

- Choi's bias/XPS board requires Ti-F formation, energy-dependent bond breaking, and ion-assisted product desorption.
- Choi's nonmonotonic oxygen response requires competing fluorination, oxygen blocking/cleanup, and neutral supply; one monotonic scalar etch coefficient is insufficient.
- Ji's RF morphology board shows positive retained/deposited volume and mask/passivation evolution are required; a removal-only solid cannot reproduce strict interfeature-gap narrowing.
- Ji's spacing board requires pattern-dependent transport and surface-state feedback.
- The Ha/Janissen oxygen profile series now quantitatively binds TiO2 rate, Cr rate, selectivity, and profile-class transitions.
- Depla and Van Bever's checksum-bound TiO2/Ar sputter curve gives a bare physical reference of about `0.192143 TiO2 formula units/ion` at `276 eV`.

Under the stronger conditional Oxford reactor boundary, clearing 700 nm in 1200 s requires roughly `0.981--1.253 TiO2 formula units per incident positive ion` before feature attenuation. Bare physical Ar sputtering is therefore only a minority of the necessary effective removal. Reactive fluorination plus ion-assisted desorption must dominate if the conditional dose is near reality.

The physical reduced TiO2 kernel already implements bounded fluorination, oxygen blocking/cleanup, passivation inventory, energetic removal, passivation growth/recession, conservative state remapping, and moving material geometry. What is missing is not another code branch; it is target-relevant numerical evidence for:

- surface site density;
- ALD-film formula-unit density within the measured material bracket;
- fluorination/complex-formation probability by radical;
- passivation sticking and bulk density;
- oxygen cleanup and blocking probabilities;
- species- and energy-resolved bare/fluorinated TiO2 yields;
- passivation sputter yields;
- Cr surface-response law under the same chemistry;
- chemistry-dependent roughness evolution.

Cross-process plasma recipes identify response signs and topology but generally do not report simultaneous absolute wafer fluxes and state-resolved surface inventories, so they do not uniquely identify those probabilities. The code continues to fail closed rather than promote cross-machine values to Oxford constants.

## Profile state

The earlier frozen conditional atlas remains valid within its declared one-moving-material/rate-normalized scope. It is not the final moving-mask board.

A new production campaign is currently running on Vast instance `48118439` (`RTX 3090`, label `petch-zhu-moving-cr-20260819`). It executes 56 independent trajectories:

- seven widths: `80, 120, 160, 200, 240, 280, 320 nm`;
- four ion energy/angular cases: low/high energy crossed with tail fraction `0/0.65`;
- two TiO2:Cr selectivities: `14.0` and `18.0167`;
- each trajectory emits both independent TiO2 rate endpoints, `34.125` and `43.4667 nm/min`;
- total intended endpoint count: 112;
- `10 nm` production mesh;
- TiO2 and Cr both move; fused silica remains pinned;
- 12 spawned deterministic workers; CUDA is used for the common transport path.

The campaign was restarted only after the exact production topology failure passed locally. At this write-up it is alive, with all 12 workers consuming CPU and the GPU context initialized. It writes a trajectory cache only when a complete trajectory finishes, so the absence of new production cache files during the first minutes is not a failure signal. The instance must be monitored, results copied and checked locally, and then stopped/destroyed so it does not continue billing.

The moving board remains explicitly cross-machine conditional: its absolute TiO2 rates and TiO2:Cr selectivities come from independent Janissen/Ha witnesses, not Freddie's SEM. No target outcome selected a coefficient.

## Krueger depth

Krueger remains unresolved rather than falsely matched.

- Published-boundary profile depth: approximately `346.833 nm`.
- Reported paper depth: `825 nm`.
- The earlier apparent `790--811 nm` agreement was retracted because two implementation errors canceled.
- The proposed `13.75 nm/s` blanket anchor is invalid: that number is a feature-average rate, not a published blanket measurement.
- The surface model agrees with direct beam evidence at feature-relevant energy, including the independently graded `4.7%` beam-measurement result.
- Reaching `825 nm` from the published aggregate flux would require a removal yield inconsistent with the direct physical ceiling used in the certified comparison.

The remaining honest route is a validated C4F6 reactor boundary or the authors' species-resolved HPEM/PCMCM wafer flux and IEAD. The repository has already advanced the C4F6 chemistry with NIST ion-fragmentation evidence, Benck radical/voltage boards, and a differential-loss closure, but those do not uniquely recover Krueger's absolute wafer boundary.

## Validation status

The last full all-files suite before today's changes was `2101 passed, 1 skipped` at the continuation checkpoint. Today's changes have not yet been followed by another full 20-minute suite. Focused results are:

- generalized boundary/topology group: `122 passed`;
- species-boundary audit group: `6 passed`;
- material/feature/surface group after species-specific yield support: `105 passed`;
- exact local Oxford production reproduction: passed, two accepted steps, requested-duration endpoint.

All today's completed checkpoints are committed and pushed. The only intentionally uncommitted paths are the two unrelated pre-existing artifacts named at the top of this document.

## What “generalizable” means here

To move from Oxford/TiO2 to another chemistry or material, the shared engine does not change. A new case supplies:

1. a gas-phase species/reaction/collision deck;
2. reactor geometry, power coupling, wall state, pressure/flow control, and diagnostics;
3. sheath and wafer-transfer evidence;
4. a material surface deck with species-resolved probabilities/yields and density;
5. mask/substrate material laws and geometry;
6. untouched measurements for grading.

The same deterministic transport, charging, conservative state update, multi-material routing, level-set motion, remapping, convergence checks, and implicit-differentiation contracts remain shared. This is the intended process-EDA architecture.

General does not mean universal constants. Chemistry and material response belong in evidence-bearing decks. Tool serial number and chamber state belong in boundary providers. Mixing those layers would make transfer look easier while destroying predictive meaning.

## Atomic-accuracy claim boundary

No atomic-accuracy claim is currently supportable.

- The best Oxford refinement rung is nanometre-scale, not atom-by-atom.
- Gas and surface states are reduced kinetic populations, not explicit atomic configurations.
- Several target surface coefficients remain unmeasured.
- The final target SEM has not been scored.

The architecture can approach atom-counted mass conservation and sub-nanometre process increments, but profile accuracy must be demonstrated against measurements at the relevant scale. DFT or molecular dynamics may help supply selected reaction barriers/yields for emerging materials; they do not replace reactor diagnostics, surface-beam experiments, or blind profile validation.

## Highest-value next actions

1. Monitor the current 56-trajectory moving-Cr campaign; retrieve, checksum, and locally certify all 112 endpoints; render the common-scale atlas; stop and destroy the Vast instance.
2. Complete a species-resolved physical-surface full-stack sentinel using the new per-ion yield-map contract. It must remain a sensitivity until every TiO2 coefficient has evidence.
3. Run the full repository suite after the current campaign artifacts and sentinel land.
4. Ask Freddie for the exact-run achieved DC self-bias/electrode waveform, blanket TiO2 loss, remaining Cr thickness, actual GDS/sample radius, and then the top-down/cross-section SEMs as the blind answer key.
5. Use those independent same-run observables to condition the boundary/surface deck; do not fit directly to the target profile.
6. For Krueger, obtain the species-resolved HPEM/PCMCM boundary or build and validate a C4F6 reactor against independent flux/diagnostic data before reopening the `825 nm` grade.

The main scientific accomplishment of this resumption is not a cosmetic depth match. It is that the full-stack interfaces now preserve species, charge, energy, angle, material identity, and evidence status without silent aggregation—and the production geometry bug that blocked the moving-mask board is reproduced, fixed, tested, committed, and running again.
