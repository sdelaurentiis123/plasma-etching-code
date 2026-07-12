# Claude Code handoff: petch unified physics engine

Date: 2026-07-12  
Repository: `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code`  
Branch: `codex/unified-engine-root-fixes`  
Handoff base: `4281c8f`  
Remote policy: **do not pull and do not push**. Make local commits only.

## Read this first

The user wants one fast, GPU-oriented, first-principles feature-scale plasma-etch engine that accepts
physical geometry, material, plasma-boundary, and chemistry inputs and predicts profiles and failure
modes across aspect ratio, geometry, and supported chemistries. Unknown surface-reaction data must remain
explicit calibrated/measured/atomistic inputs with provenance and uncertainty. Known conservation laws,
transport, electrostatics, and interface physics must not be replaced by fitted profile formulas.

The near-term product is not “arbitrary recipe name in, exact wafer out.” It is:

> physical feature geometry + dimensional plasma boundary + material mechanisms in → profile, surface
> state, flux/energy/angle, charge/field, convergence, uncertainty, and validity/refusal out.

The commercial proof is calibration of declared uncertain physical inputs on structure N followed by a
held-out prediction on structure N+1. Differentiable calibration and inverse recipe design are later,
after the forward engine is verified.

Do not start by running a large charging or AR campaign. First establish which common-engine capability
and external observable the run would close. Every investigation needs one causal hypothesis, one
external or analytic judge, a numerical error budget, and a compute ceiling.

## Non-negotiable operating rules

1. Work only in the repository and branch above.
2. Before editing, run `git status --short --branch`, `git log --oneline -12`, and `git reflog -5`.
3. Confirm no other Codex/Claude session is writing this repository. Two concurrent sessions previously
   edited and committed into the same linear branch. During this reconciliation, a resumed Codex CLI and
   a second Codex app node session were found and terminated. Do not spawn another writer or subagent.
4. Do not reset, rewrite, squash, pull, push, or discard local commits. If a bounded subsystem is wrong,
   replace or revert it with a new local commit.
5. Use `apply_patch` for edits. Preserve unrelated user work.
6. Run focused tests while iterating and the full suite before each promoted checkpoint.
7. Never describe agreement with ViennaPS, Hwang–Giapis, MCFPM, MD, or a manufactured solution as
   experimental validation.
8. No benchmark name, aspect-ratio band, expected answer, or surface region may select governing physics.
9. Adaptive numerical effort may depend on measured error, not on “floor,” “corner,” “AR40,” or a target
   answer.
10. Do not silently provide a Maxwellian, cosine distribution, reaction probability, product branch, or
    material conductivity when the physical source did not specify it. Refuse or expose an explicit
    closure.

## Current Git and verification state

The worktree was clean at handoff creation except for this report. The last verified full suite on the
`4281c8f` tree was:

```text
296 passed, 1 skipped in 73.53s
```

Run:

```bash
pytest -q
```

before trusting the tree after any further edits.

The repository has one worktree and a linear local history. `main` is far ahead of `origin/main`; the
current branch begins at local `main` commit `109badb` and contains the reconciliation/root-fix commits.
There was no hidden competing branch to merge. `git fsck --unreachable --no-reflogs` found only an older
Jeon README blob, not a lost code commit.

### Commits added in the reconciliation arc

| Commit | Meaning |
|---|---|
| `72f425d` | Conserved periodic transport geometry and interface motion; added CR-2/redistancing gates |
| `72a005e` | Repeated subcell interface-motion regression gate |
| `a92c713` | Reproduced the pinned-mask subcell-motion defect |
| `b16dbec` | Qualified finite-extrusion evidence rather than overstating it |
| `fda21a3` | Evolve exposed material level sets independently; stop mask/substrate union pinning |
| `7842f23` | Validate layered geometry and null transport atoms |
| `9f419ee` | Reset the Jeon ladder after the geometry fix; did not preserve obsolete favorable numbers |
| `bb65f34` | Reconciled legacy and common engine paths; public legacy results now identify their engine |
| `a60b77f` | Recorded local history/lineage audit and withdrew misleading 2-D charging claims |
| `ae0c066` | Added explicit `PhysicalProcess` public entry into `feature-3d-v1` with validity/provenance |
| `95ab945` | Added exact surface material-exchange ledger and unresolved-product refusal |
| `eb4cc32` | Added globally conservative surface-emission form-factor transport |
| `4281c8f` | Added named product populations, generic physical sputtering, pinned DeepMD Si sputtering, and common-engine product readiness diagnostics |

Read [LOCAL_HISTORY_RECONCILIATION.md](LOCAL_HISTORY_RECONCILIATION.md) for the complete history audit and
[PHYSICS_PRODUCT_CAPABILITY_MAP.md](PHYSICS_PRODUCT_CAPABILITY_MAP.md) for the current physics/product map.

## The central reconciliation: there are three lineages

Do not call these one validated engine yet.

### 1. Legacy compatibility/demo engine

Entry point:

```text
src/petch/api.py::Process.run
    -> src/petch/threed.py::run_etch_3d
```

Result engine id: `legacy-threed-v1`.

This path preserves historical SF6/O2, ViennaPS parity, ARDE, charging/notching, Bosch, cryogenic, ALE,
and speed demonstrations. It is a roughly 2,000-line flag-selected monolith. It contains useful physics,
but also calibrated closures and at least one unacceptable region exception:

```text
src/petch/threed.py around the redeposition block
```

When charging and an etch-stop are enabled, it suppresses redeposition in a geometry-defined notch-foot
band. Do not port this. Replace it with species/material product generation, transport, deposition, and
resputtering.

The legacy high-level API was not redirected because doing so would silently lose mechanisms. Its
`Result.engine` now reports `legacy-threed-v1`.

### 2. Experimental 2-D charging/adjoint lineage

Principal files:

- `src/petch/charging_backward.py`
- `src/petch/boundary_transport.py`
- `src/petch/charging_nodal.py`
- `src/petch/charging_nodal_fixed_point.py`
- `scripts/charging_nodal_campaign.py`

Purpose: cheaply derive, test, and falsify adjoint/current-balance estimators on 2-D cell/nodal geometry.
This is not the product solver. Its analytic gates are useful; its reference-emulation conventions must
not leak into production.

Important correction: the old module header claimed the backward solver converged a validated dipole
across AR4/8/15. That statement was withdrawn. Current-balance convergence, nonuniform-field reciprocity,
and the HG reference do not all close in one demonstrated configuration. A `Vs^-0.35` RF phase weight is
HG-reference-specific, not a first-principles production source.

### 3. Common dimensional 3-D engine

Entry point:

```text
src/petch/physical_api.py::PhysicalProcess
    -> src/petch/feature_step_3d.py::solve_feature_3d
```

Result engine id: `feature-3d-v1`.

This is the only target for new governing feature physics.

Core files:

- `src/petch/boundary_state.py`: immutable SI plasma-to-feature boundary distribution contract.
- `src/petch/boundary_transport_3d.py`: arbitrary-triangle first-hit, field transport, deterministic face
  gather, and form-factor geometry.
- `src/petch/charging_poisson_3d.py`: Q1 variable-permittivity Poisson operator and compatible sheet-charge
  projection.
- `src/petch/charging_coupled_3d.py`: physical-time and safeguarded quasi-static dielectric charging.
- `src/petch/surface_kinetics.py`: reduced stateful SiO2/fluorocarbon mechanism with conserved sites and
  polymer inventory.
- `src/petch/tabulated_chemistry.py`: narrow sourced Si-Cl2-Ar+ and Si physical-sputter mechanisms.
- `src/petch/surface_interaction_table.py`: versioned interaction tables, uncertainty, interpolation audit,
  and extrapolation refusal.
- `src/petch/feature_step_3d.py`: material-local surface state, conservative remap, additive material level
  sets, interface motion, validity, and optional charging/radiosity.
- `src/petch/physical_sputtering.py`: general configurable physical-sputter mechanism.
- `src/petch/surface_exchange.py`: exact removed/outgoing/unresolved/deposited material ledger and named
  product-population contract.
- `src/petch/neutral_radiosity_3d.py`: conservative diffuse neutral and surface-emitted product transport.

The common engine may reuse independently gated numerical primitives from `threed.py`—marching cubes,
velocity extension, level-set advection, and redistancing. That reuse does not make the legacy monolith the
common engine.

## What is now physically true in the common engine

### Boundary state

`PlasmaBoundaryState` carries dimensional species flux, signed charge, mass, joint velocity/energy/angle,
phase/position where available, a physical reference plane, and provenance. Analytic Maxwellian,
finite-transit RF sheath, histogram, mixture, and sampled representations share this contract.

The source law is upstream input. Transport must not know whether it came from a reactor model, PIC,
diagnostic reconstruction, or learned surrogate.

### Transport

The common 3-D adapters preserve species-resolved flux and energetic hit events on arbitrary triangles.
Charged particles can traverse a nodal electrostatic field. The deterministic face gather is useful for
reproducibility but does not by itself make the evolving geometry differentiable; hard visibility/hit
boundaries remain.

Diffuse neutral radiosity uses the required area reciprocity:

```text
B[i,j] = A[j] F[j->i] / A[i]
```

and closes source = reacted + escaped. Surface-emitted populations now use the same operator but account
separately for direct escape before first impact, first impact, repeated reaction/re-emission, and final
escape.

### Electrostatics/charging

The common 3-D path has variable permittivity, compatible nodal charge deposition, electrostatic work,
physical-time current deposition, and a safeguarded steady dielectric solve that refuses nonconvergence.
It is quasi-static across evolving geometries: each new geometry rebuilds the material operator and solves
a new root. Transient moving-surface charge remap, mixed conductors/dielectrics, finite conductivity,
secondary-electron emission, and experimentally closed HARC charging remain open.

### Surface state and chemistry

The reduced SiO2 mechanism conserves:

- accessible complex sites;
- finite fluorocarbon film inventory;
- complex formation/removal;
- bare/complex energetic oxide removal;
- fluorocarbon deposition;
- oxygen/energetic film removal.

Its current test parameters are manufactured/nonpredictive. Crosslinking, resolved fluorination sequence,
complex stoichiometry, product branching, conductivity, and ACL chemistry are open.

The pinned CC-BY DeepMD data in `data/surface_interactions/kounis_melas_2024/` provides narrow sourced
Si-Cl2-Ar+ evidence. It is MD evidence, not experiment and not SiO2 chemistry.

### Material exchange and products

Every reduced SiO2 and tabulated Si step now returns `SurfaceMaterialExchange`:

```text
removed material = routed outgoing material + chemically unresolved material
```

face-by-face to floating-point tolerance. Unresolved reactive products cannot be transported or
redeposited.

`SurfaceProductPopulation` maps named particles back to the material-origin ledger. A product is
transport-ready only when it has both:

- an angular law supported by the chosen transport backend; and
- a complete parameterized energy law.

The generic physical-sputter mechanism can emit a fully declared product. The pinned Si DeepMD
energy→yield table routes every removed Si atom to named neutral Si, but correctly reports that it is not
transport-ready because the source table lacks a differential emitted energy-angle distribution.

The common feature step reports:

- `product_routing_complete`;
- `product_population_count`;
- `product_transport_ready`;
- material/product limitations in its validity result.

This is not full redeposition yet. Deposited-product state, growth/interface coupling, and resputtering are
still open.

### Geometry/interface

The common engine uses additive per-material level sets. Only exposed etchable material fields move. This
fixed a real defect where a pinned mask level set won the union beneath a narrow opening and stopped a
physically bombarded substrate after roughly one cell.

CR-2 constrained redistancing and material-local conservative state remap are gated. Topology-changing
events are refused rather than silently remapped. The current interface scheme remains first-order and
needs full grid/time/profile convergence.

## Evidence classes and current claim discipline

Use the definitions in [EXPERIMENTAL_VALIDATION_MATRIX.md](EXPERIMENTAL_VALIDATION_MATRIX.md):

- `A`: analytic/manufactured/numerical invariant;
- `S`: published simulation or independent-code comparison;
- `E2`: qualitative/normalized/partial experimental observable;
- `E3`: calibrated experiment;
- `E4`: held-out experiment.

“Validated” without qualification is reserved for E3/E4 and must identify fitted parameters and held-out
observables.

Current honest state:

| Capability | Best historical evidence | `feature-3d-v1` state |
|---|---|---|
| Boundary/ballistic transport | A | substantial analytic/numerical gates pass |
| Interface/remap/material geometry | A | substantial gates pass; full convergence open |
| Neutral radiosity | A | conservation and refinement gates pass |
| Physical surface-emission transport | A | conservation gates pass; deposited-state feedback open |
| Reduced SiO2 chemistry | A + literature topology | executable but current parameter set is nonpredictive |
| Si-Cl2-Ar+ RIE | S/MD | narrow sourced table passes through the common engine |
| Si physical sputtering | S/MD | yield/removal/routing pass; emission distribution missing |
| de Boer Si SF6/O2 ARDE | historical E3/E4 split | **not re-earned** |
| HG/Hwang–Giapis charging | S, low-AR simulation | **not re-earned** |
| Fujiwara/Nozawa notching | historical E2 | **not re-earned** |
| Jeon/Jeong SiO2 transfer | extracted experiment + preregistered score | no passing common-engine held-out prediction |
| Bosch DRIE | historical E2/E3 reduced gates | not migrated |
| cryogenic etch | historical narrow E3 reduced gate | not migrated |
| Si-Cl2-Ar+ ALE | reduced 0-D S/E2/E3 | transient spatial ALE state not migrated |
| bowing/necking | open/partial historical | product/mask feedback open |
| stochastic twist/systematic tilt | open | common evolving ensemble gate open |
| differentiable feature calibration | A only in reduced ALE | open in common feature solver |
| reactor-feature two-way coupling | open | inbound boundary exists; feedback schema/open validation remain |

## The historical benchmarks that must be re-earned

The user specifically asked for de Boer/ARDE, HG, Bosch, cryogenic, ALE, and all prior useful behavior to
become one coherent engine. Green legacy tests are regression protection, not migration proof.

### de Boer / Blauw ARDE

Historical claim: a wall-loss parameter calibrated at AR10/20 predicted held-out AR40 normalized rate.
Reported full-curve RMSE was roughly 0.031–0.043 by seed.

Problems:

- the result is not zero-calibration first principles;
- absolute rate/profile contours and independent chemistry are missing;
- several old scripts contain AR-shaped closure parameters;
- no canonical de Boer campaign currently runs through `feature-3d-v1` in `pytest`.

Migration requirement:

1. Use one physical boundary distribution and one stateful Si/SF6/O2 mechanism.
2. Use common ballistic/radiosity transport without AR-dependent formulas.
3. Calibrate only declared wall/material probabilities on the preregistered calibration points.
4. Predict held-out AR40 and preferably an independent profile/time dataset.
5. Report grid, face quadrature/form-factor, time-step, calibration, digitization, and model-form error.

Relevant files:

- `scripts/deboer_*`
- `scripts/validate_experiment_arde.py`
- `src/petch/knudsen.py`
- `src/petch/params.py`
- `EXPERIMENTAL_VALIDATION_MATRIX.md`
- `FINDINGS.md`

Do not port `ar_pass`, fresh-band, or expected-floor formulas as governing physics.

### Hwang–Giapis / Fujiwara / Nozawa charging and notching

The local HG PDFs are:

- `refs/HG_jap97.pdf`
- `refs/HG_jvstb97.pdf`
- `refs/HG_deep_read.md`

These are low-AR poly-Si-over-oxide overetch/notch studies, not modern HARC SiO2 validation. HG solves a
specific cell/Laplace charging model and ignores several effects, including surface currents, SEE, and
tunneling. Use it as a causal/unit/reference gate.

Migration requirement:

1. Reconstruct the exact material geometry and conductor/dielectric topology.
2. Use the common physical boundary and 3-D Poisson/current deposition.
3. Require converged current balance, grid/sample/init invariance, and transport uncertainty below the
   profile error budget.
4. Demonstrate charging-off/on notch causality.
5. Compare absolute matched profile/energy observables where the paper supplies them.
6. Keep HG simulation evidence (`S`) separate from Fujiwara/Nozawa experimental trend evidence (`E2`).

Do not import `ion_ied_phase_exponent=0.35` into production. It is reference emulation.

### Bosch, cryogenic, and ALE

These currently live as separate reduced modules:

- `src/petch/bosch.py`
- `src/petch/cryo.py`
- `src/petch/ale.py`
- `src/petch/ale_diff.py`

They must become material mechanisms/time schedules on the unchanged common boundary, transport,
surface-state, remap, geometry, validity, and result contracts. The core must not acquire a `bosch=True`,
`cryo=True`, or `ale=True` branch.

For ALE, transient dose/history, mixed/damaged layer, pulse subdivision invariance, and product branching
must be explicit state. A 0-D cycle curve and gradient do not prove spatial profile evolution.

## OSTI and Oehrlein objectives

The relevant sources and their operational interpretation are:

1. **DOE/FES semiconductor plasma report (2023)**  
   `https://science.osti.gov/-/media/fes/pdf/2023/DOE_FES_PlasmaScience_Semiconductors_Final.pdf`
   - Feature-scale reaction probabilities are often selected empirically and calibrated to SEMs, limiting
     quantitative range.
   - Reactor and feature models should exchange incident flux distributions and effective surface
     reaction/sticking/product feedback.
   - No single mesh spans sub-nanometres to a half-metre chamber; use physically coupled scale-specific
     solvers.

2. **DOE workshop roadmap / OSTI 2349231**  
   `https://www.osti.gov/servlets/purl/2349231`
   - Recipe-development speed and process-window extension are key industry constraints.
   - The forward solver is necessary but not the final inverse/design loop.

3. **Oehrlein et al., JVST B 2024 roadmap**  
   DOI `10.1116/6.0003579`  
   `https://cpseg.eecs.umich.edu/pub/articles/JVSTB_42_041501_2024.pdf`
   - Calls for predictive multiscale plasma-surface models, experiments, UQ, hybrid physics/ML, and digital
     twins under increasingly three-dimensional, atomic-scale, selective, low-damage processing.
   - “Digital twin” is a destination, not a claim this repository has earned.

4. **Huang/Shim/Nam/Kushner pattern-dependent distortion**  
   DOI `10.1116/1.5132800`
   - At extreme AR, stochastic arrivals/charge produce random twist with zero ensemble mean in a symmetric
     array.
   - Dense/sparse nonlocal electrostatics produce systematic tilt.
   - A single smooth or isolated trench cannot validate either.

5. **Krueger et al., JVST A 2024 autonomous hybrid optimization**  
   DOI `10.1116/6.0003554`  
   `https://cpseg.eecs.umich.edu/pub/articles/JVSTA_42_043008_2024.pdf`
   - MCFPM parameters are optimized against a SEM, then tested on process trends.
   - Roughly 20 CPU-hours for one feature and partial transfer define a real competitive baseline.
   - petch must beat it on held-out transfer, data efficiency, and compute—not merely training error.

6. **K-SPEED physical sputtering, Physics of Plasmas 2026**  
   DOI `10.1063/5.0304157`
   - Fits time-resolved trench profiles from 0–60 min.
   - Reproduces mask facets, entrance necking, taper, and transfers to hole geometry.
   - Trench/hole differences arise from geometry-dependent redeposition.
   - This is a stronger validation form than matching only the final depth.

7. **Transport-mediated necking, Vacuum 2026**  
   DOI `10.1016/j.vacuum.2026.115515`
   - Couples 3-D transport and surface state under matched C4F8/Ar/O2 conditions.
   - Necking comes from local polymer deposition versus ion-assisted removal near the hardmask entrance.
   - The effective opening then feeds back on deep-feature transport and process-window width.

The unified acceptance map is therefore:

| Roadmap objective | Required engine evidence |
|---|---|
| Predict plasma-surface interactions | dimensional boundary → transport → state → interface, with sourced/explicit unknown inputs |
| HAR transport and ARDE | same operator across wafer/trench/hole and AR ladder; grid/ray/domain convergence; held-out experiment |
| Differential charging | compatible Poisson/current deposition; converged root; charging-off/on absolute profile |
| Dynamic surface chemistry | conserved state/history; sourced mechanism; time-resolved and held-out transfer |
| UQ/validity | separated numerical/measurement/parameter/model error and explicit refusal |
| Reactor-feature coupling | inbound flux distribution plus outbound effective surface/product feedback |
| Recipe acceleration | verified forward transfer first, then calibration/inverse loop and comparison to DOE/derivative-free baselines |
| AI-native product | learn unknown closures/reduced reactor mappings only; preserve known physics/conservation and provenance |

## Recommended next execution sequence

### Task 0: create the backend-qualified capability manifest

This was about to be implemented when the user requested this handoff. No incomplete manifest code is
left in the tree.

Create `verification/capability_manifest_v1.json` and a strict loader/verifier. Each capability record
should contain:

- gate id and name;
- roadmap objective ids;
- required engine (`feature-3d-v1`);
- current common-engine status;
- backend-qualified evidence records (`legacy-threed-v1`, `research-2d`, `reduced-0d`, or
  `feature-3d-v1`);
- evidence class (`A`, `S`, `E2`, `E3`, `E4`);
- repository-owned artifacts;
- promotion evidence required;
- next blocker.

Required entries: common boundary transport, interface/remap, SiO2 state, de Boer ARDE, HG charging,
notching, Jeon transfer, Bosch, cryo, ALE, Si-Cl2-Ar+ RIE, physical sputtering, product redeposition,
bow/neck/twist/tilt, reactor feedback, differentiability, performance, and public product behavior.

The verifier must reject an external capability marked passed unless a `feature-3d-v1` E3/E4 record and
existing artifacts support it.

### Task 1: re-earn de Boer/ARDE through the common engine

This is the first historical migration because it stresses geometry/neutral transport before charging.
Do not start with a long evolving campaign. First build static open-wafer/trench/hole flux gates using the
same boundary and material mechanism. Establish grid/form-factor/face-quadrature convergence and inspect
which physical channel creates the measured AR trend.

Only then run the preregistered calibration/held-out curve. Keep the legacy score beside the common score.

### Task 2: implement product deposition as material state

The outgoing ledger and conservative transport exist. Add:

- sourced product branching;
- sourced/explicit launch distribution;
- material/state-dependent sticking or reaction;
- deposited inventory on the receiving material;
- volume/interface growth or passivation effect;
- resputtering and exact material balance across steps.

Use one common surface interaction API. Do not subtract a scalar redeposition flux from velocity. Gate with
manufactured two-face conservation, then time-resolved necking/taper and trench→hole transfer.

### Task 3: energetic ion reflection

Replace the legacy coned-cosine closure with an interaction-table kernel over material, projectile,
energy, incidence, and roughness. Cagomoc et al. 2023 is a primary starting point for inert-gas scattering
on smooth/rough Si and SiO2. Expose missing dimensions/uncertainty and refuse extrapolation.

### Task 4: mask as material chemistry

Use additive material level sets and a mechanism registry. ACL oxidation/sputter/deposition must evolve the
mask as a material, not via a pinned region. This is required for entrance widening, facets, bow, and
necking.

### Task 5: Jeon/Jeong held-out transfer

The repository contains 54 checksum-verified trench depths, plasma controls, electron/bias controls, and
a frozen scorecard. The paper does not provide complete IEDF, absolute species fluxes, or etch duration;
keep those as explicit closures. Only the 20% C4F8 continuous-wave curve is calibration. Predict held-out
width ratios and pulse-response reversals without refit.

Relevant files:

- `data/experimental/jeon_2022/`
- `src/petch/experimental_data.py`
- `src/petch/experimental_boundary.py`
- `src/petch/validation_demo.py`
- `scripts/jeon_unified_baseline.py`

### Task 6: common 3-D charging/HG causality

After mean transport/material state is stable, construct the HG reference in the common engine. Do not use
HG as the modern HARC product demo. Then select a matched dielectric experiment where charging changes the
profile above the measurement/numerical error budget.

### Task 7: migrate Bosch, cryo, and ALE schedules

Represent process sequencing as time-dependent boundary/mechanism inputs. Reuse the same solver. Validate
history and pulse-subdivision invariance before profile claims.

### Task 8: differentiation and calibration

Only after the forward transfer gates. The current deterministic architecture is necessary but not
sufficient. Hard visibility/hit boundaries and adaptive estimator selection invalidate naive tracer
autodiff. Require finite-difference agreement through transport, converged charging, surface state, remap,
and profile loss. Then compare calibration cost and held-out transfer against Krueger/derivative-free
baselines.

## Tests and useful commands

Baseline:

```bash
git status --short --branch
git log --oneline --decorate -12
pytest -q
```

Focused common-engine tests:

```bash
pytest -q tests/test_boundary_state.py
pytest -q tests/test_boundary_transport_3d.py
pytest -q tests/test_charging_poisson_3d.py tests/test_charging_coupled_3d.py
pytest -q tests/test_surface_kinetics.py tests/test_tabulated_chemistry.py
pytest -q tests/test_surface_exchange.py tests/test_physical_sputtering.py
pytest -q tests/test_neutral_radiosity_3d.py
pytest -q tests/test_feature_step_3d.py tests/test_reinitialization_3d.py
pytest -q tests/test_experimental_data.py tests/test_experimental_boundary.py
pytest -q tests/test_validation_demo.py
```

Legacy regression tests must remain green, but their pass status does not close common-engine migration.

## Stop rules

Stop and reassess rather than adding compute when:

- the gating observable worsens non-monotonically under the same solver method;
- numerical uncertainty is already below measurement/model error;
- a missing physical input dominates the answer;
- a proposed fix needs an AR/benchmark/region conditional;
- a profile match improves only by degrading held-out transfer;
- a transport correction violates reciprocity, probability, energy, or material balance;
- charging residuals improve but the target profile is insensitive within its error budget;
- a source table lacks the dimension being extrapolated.

## Final claim boundary at handoff

It is accurate to say:

> petch now has an explicit dimensional common 3-D feature engine with physical plasma boundary states,
> species-resolved transport, stateful surface mechanisms, conservative material/interface evolution,
> optional self-consistent dielectric charging, exact surface material ledgers, and conservative emitted-
> product transport primitives. A narrow sourced second chemistry and physical-sputter table run through
> the same contracts. Unsupported product and chemistry inputs are exposed or refused.

It is not yet accurate to say:

- all historical capabilities run through one engine;
- de Boer ARDE is revalidated in the common engine;
- HG/notching is common-engine experimentally validated;
- SiO2 reactive product branching/redeposition is complete;
- ACL mask evolution, bowing, necking, twisting, and tilt are predictive;
- the feature solver is differentiable end-to-end;
- reactor-to-feature coupling is two-way;
- petch is an experimentally validated digital twin or profile-to-recipe product.

The next Claude Code session should preserve this honesty boundary while methodically moving each row from
legacy/research evidence to `feature-3d-v1` evidence.
