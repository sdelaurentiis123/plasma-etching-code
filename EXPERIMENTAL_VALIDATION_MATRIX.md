# Experimental validation matrix

Audit date: 2026-07-11. This is the claim-control document for petch. A numerical invariant, agreement
with another simulator, agreement with a published simulation, and agreement with a physical experiment
are different evidence classes and must not be described interchangeably.

## Evidence classes

- **E4 — held-out experiment:** physical data not used to choose the mechanism or its parameters.
- **E3 — calibrated experiment:** physical data matched after fitting declared physical parameters.
- **E2 — experimental trend/partial observable:** qualitative, normalized, cross-scale, or incomplete
  comparison to physical data.
- **S — simulation reference:** comparison to ViennaPS, MCFPM, Hwang–Giapis simulation, MD, or another
  computational result. Useful, but not experimental validation.
- **A — analytic/numerical invariant:** conservation, limiting law, convergence, manufactured solution,
  or finite-difference consistency.
- **Open:** no adequate gate.

“Validated” without a qualifier should be reserved for E3/E4 and must name the calibrated quantities.

## Current scorecard

| Physics/module | Best evidence | Current status | What is actually established | What remains open |
|---|---:|---|---|---|
| Backward electron gather | A | passing on current gates | Open-wafer normalization, Langmuir `exp(V/Te)` retardation, positive-potential saturation, fixed-seed reproducibility, and independent forward/backward scoring of a refined trench at 0 and +9.2 V floor potential are automated. | Reciprocity in a nonuniform self-consistent field; explicit spatial-quadrature convergence; smooth/verified parameter sensitivities. A coarse version of the trench gate differed by 16% at zero field and about 6% at +9.2 V, showing face-center discretization materially biases flux. |
| Backward ion gather | A | passing on current gates | Open uncharged wafer returns unit normalized flux; independent forward/backward scoring passes on a refined trench at zero field and a +5 V retarding floor; a manufactured uniform-field orbit bounds production-step impact-energy error below 0.7%. The core now defaults to analytic uniform RF phase. | IEDF moment replay and retarding-potential curve; nonuniform self-consistent-field reciprocity; reflection/re-emission; replace the reduced instantaneous sheath with a derived time-dependent sheath when its regime requires it. The `p=0.35` phase weight is explicitly HG-benchmark-only. |
| Self-consistent backward charging | A partial; S failing when converged | **conservation bug fixed; reference mismatch exposed** | Corner-face current pooling restores residual decay; deterministic sub-face integration removes point collocation; electron/ion sample ladders are converged; the 3.7 µm physical source-plane correction is inert; W32 does not move toward the reference; and 250/500/1000 Laplace sweeps give the same iteration-10 floor (34.10–34.11 V) with field RMS residual ~6–9×10⁻⁵. | The published 34 V match is exactly an iteration-10 transient while charging residual RMS is ~0.61. Continuing toward current balance drives AR4 toward floor ~51 V and foot energy ~11.5 eV, not 34.0 V/21.6 eV. Remaining suspects are source/current normalization, conductor/insulator physics, boundary geometry, or missing charge transport—not particles, source height, fixed Laplace sweeps, or simple W16→W32 refinement. |
| Backward notch-foot ion energy | S | failing | The observable is computed. | Replayed 2026-07-11: 10.5 eV RMSE and correlation -0.323 versus `_PETCH_FOOT_E`. The “face convention” note does not close this gate. Do not call it validated. |
| HG floor-flux charging curve | S | passing in legacy/closure lineage; mixed in newer solvers | The code can reproduce the Hwang–Giapis published computational curve under documented configurations. | HG is a simulation benchmark, not wafer measurement. The general/backward estimator must pass without geometry-specific overrides and with convergence evidence. |
| Charging-driven notch mechanism | E2 | partial | Charging-on creates a localized foot notch while charging-off does not; normalized notch trend is compared with published notch measurements. | Absolute geometry/process matched experiment; width dependence; uncertainty in digitization; foot-energy mechanism inconsistency above. |
| Neutral/ion transport and de Boer ARDE | E3/E4 split | strongest feature-scale result | A declared wall-loss parameter is calibrated on AR10/20; the evolving model predicts the held-out AR40 normalized rate around the measured value. Full-curve RMSE is reported as 0.031--0.043 by seed. | Independent chemistry/wafer dataset; absolute rate and profile contours; uncertainty in digitized ARDE; disentangle reaction model from transport. This is not zero-calibration first principles. |
| ViennaPS transport parity | S | passing for selected gates | Ballistic/radiosity and profile comparisons characterize implementation differences and speed. | No experimental authority. Agreement with ViennaPS cannot validate shared approximations. |
| Continuous SF6/O2 surface chemistry | E3 + literature parameters | partial | Belen/Gomez-derived mechanism and de Boer calibration reproduce selected rates/ARDE. | Energy-angle-flux-ratio yields across independent experiments; surface state memory; product branching; material generality. |
| Si-Cl2-Ar+ ALE reduced model | S + E2/E3 | passing at 0-D process level | Published MD-derived reduced model reproduces ALE window, dose behavior, synergy, and selected experimental etch-per-cycle anchors; autodiff matches finite differences. | Spatial surface-state field, mixed-layer/damage-depth evolution, moving-interface conservation, independent held-out cycles. |
| Cryogenic etch rate model | E3 | narrow calibrated pass | Langmuir form with measured adsorption energy and two declared calibrated parameters reproduces a published 1.6x rate anchor. | Correct chemistry match, independent temperatures/pressures, condensed-layer transport/conductivity, profile evolution. |
| Bosch cyclic etch | E3/E2 | partial | Existing Ayon/Tillocher gates constrain depth/scallop/process regimes. | Independent profile dataset, reactor drift, spatial wafer uniformity, mask selectivity, and held-out cycle recipes. The newly imported Zenodo data begins the wafer/reactor layer. |
| 3-D feature evolution | S + A partial | partial | Level-set/transport operation and selected parity/speed tests exist. | Experimental 3-D contours, grid/time/ray convergence as a suite, charging in 3-D, stochastic twist distribution. |
| Pattern-dependent charging/tilt | Open | open | Literature-defined sign expectations are documented. | Symmetric zero-mean lateral field, dense-to-sparse field sign, neighbor/domain convergence, stochastic ensemble, evolving 3-D tilt. |
| Reactor-to-feature coupling | Open + public data acquired | open | Physical source-plane boundary is conceptually clear; a CC-BY Bosch wafer dataset is now locally available with verified provenance. | Common source object, machine/OES ingestion, reactor model or surrogate, flux-distribution inference, surface feedback, held-out wafer prediction. |
| End-to-end differentiable calibration | A at ALE only | open for feature solver | A small 0-D ALE inverse example works. | Gradients through transport, charging fixed point, chemistry, and level set; finite-difference gates; contour loss; comparison against derivative-free calibration. |

## Immediate physics blockers

### 1. Charging convergence is not yet demonstrated

`self_consistent_backward()` historically performed a fixed number of damped voltage updates. It now
reports per-surface and pooled-conductor `log(Gi/Ge)` residuals and supports opt-in tolerance stopping.
The first measurement falsified the prior convergence claim. It also exposed a concrete bug: a corner
cell was assigned one independent floating balance per exposed face, so opposing face currents applied
multiple contradictory updates to the same potential. Pooling faces by physical cell restores residual
decay, but the converged direction now misses the forward reference badly. Estimator/source/field parity
must be diagnosed before early stopping or any charging-validation claim.

Subsequent kill tests narrowed the cause. Electron and ion fluxes are stable from `n_log2=9` through 12;
using the physical 3.7 µm source plane instead of the reduced 2.5 µm plane is inert; W32 does not move the
solution toward 34 V; and tripling/quadrupling Laplace sweeps leaves the voltage unchanged with small PDE
residual. The next work must audit the actual current/source normalization and material charge physics.

The actual-field reciprocity audit found an additional geometry bug: backward rays were launched 1.5
cells away from the interface with an asymmetric origin formula. Launching just outside the true cell face
reduces W16 electron reciprocity error from +43.5% to +15.3%. Residual electron/ion errors have opposite
signs and are step-size invariant from 0.15 to 0.04, so field/interface discretization remains open.
At W32 the electron error falls to +5.3%, confirming spatial convergence, while the ion error remains
−21.3%. The ion launch-energy map and nearest-cell electric-field integration are therefore not yet a
discrete conservative pair on the actual nonuniform field; fix the discretization or derive an exit-state
weight rather than introducing a fitted correction.

The audit is now reproducible as `scripts/backward_actual_field_reciprocity.py`. Its W16 production-
statistics check with analytic uniform RF phase reports backward/forward floor-flux errors of +8.6% for
electrons and −16.9% for ions at iteration 10 (field residual RMS 1.0e-4; charging residual RMS 0.75).
The fixed iteration is deliberately not called converged; the script freezes one field solely to test the
transport reciprocity invariant. Run its W16/W32 ladder before accepting any particle-mover replacement.

A Liouville-consistent opt-in ion estimator now evaluates the declared RF-arcsine/Gaussian source density
at the actual traced plasma-exit velocity. This supplies the phase-space-density ratio missing from the old
1-D launch-energy map, without a fitted coefficient. At W16 it reduces the ion error from −16.9% to −9.0%;
at W32 it reduces −20.7% to −6.0%, while electron error falls to +1.1%. A 4x smaller timestep leaves the
W16 residual unchanged. Thus the source-score derivation fixes the dominant ion error, and the remaining
error is spatial/interface convergence. The new estimator remains opt-in until that interface gate closes.

The remaining gap was then traced to proposal support, not a fitted physical term. The 1-D proposal only
launches surface vertical energies in the shifted RF interval, but a 2-D field can exchange lateral and
vertical kinetic energy. An opt-in multiple-importance proposal now mixes that efficient analytic stratum
with a broad, known-density surface-energy stratum; the exact mixture density appears in the Liouville
weight. At production statistics it gives W16 ion error +1.4% and W32 +2.7% (electron +0.6%), closing the
4% reciprocity gate. W32 one Sobol level lower overshoots by 10%, exposing slow convergence near the
integrable RF-arcsine horns. Therefore this is a correct-but-statistics-sensitive experimental path, not
yet the self-consistent default; sample convergence must be automated before promotion.

The corrected estimator is also wired through `self_consistent_backward()` behind explicit
`ion_exit_state_weight` / `ion_exit_energy_mixture` options. A W16, iteration-10, `n_log2=10` diagnostic
moves the floor from 38.8 V to 36.4 V, reduces charging-residual RMS from about 0.75 to 0.61, and gives
−0.5% ion reciprocity on the newly generated frozen field. This is a mechanism consequence, not validation
against HG: iteration 10 is not current-balanced, W16 electron reciprocity is still +13%, and the solver
default remains unchanged until spatial and sample convergence are enforced together.

A standalone boundary-fitted nodal Laplace/Q1 tracer now passes analytic field, exact-face absorption,
no-tunnelling, and manufactured electron/ion reciprocity gates. On the high-stat W32 frozen AR4 audit it
closes corrected-ion reciprocity to +1.6% but leaves electron reciprocity at -7.6%; therefore it is not a
production replacement and is not wired into the fixed point. This separates the remaining frontier:
forward electron source-plane quadrature/domain convergence. A support-complete natural/barrier-shifted
Liouville proposal resolves the +39 V rare event and agrees with the legacy electron gather, ruling out
proposal support as the dominant residual. Exact lateral reflection fixes a smaller orbit error. Four
forward scrambles reduce the W16 electron discrepancy to -1.4% (shifted Liouville -1.9%) and quantify
forward standard error as 0.00144; the former single-scramble comparison overstated electron bias.

The corrected-ion variance was traced to the uniform-energy broad proposal at RF-arcsine horns. Replacing
it with an exact Chebyshev/arcsine broad density reduces W16 discrepancy to +0.45% in a four-scramble run,
but replicate tails remain and a universal per-element adaptive controller correctly refuses stringent
convergence even after refining all elements. This is not a reason for named corner sampling. A new
collisionless finite-transit-time RF sheath model now passes static-energy, Child-thickness, and high-
frequency phase-mixing gates; it is the upstream path toward a nonsingular physical boundary state.

`PlasmaBoundaryState` now provides the single immutable plasma-to-feature contract: normalized weighted
joint velocity-energy samples, phase/position, signed charge, mass, absolute flux, reference plane, and
provenance. Instantaneous and finite-transit sheath constructors pass the same contract and current-density
invariant. Production transport does not consume it yet; embedded source laws remain an explicit migration
blocker.

The first boundary-state transport adapter now preserves probability, absolute flux, and full 3-D kinetic
energy and passes the same open-wafer transport gate for charged ions and neutral CF2. This demonstrates a
species-independent interface, not universal chemistry validation. Production source migration and an
adjoint density representation remain open.

A normalized joint velocity-histogram density now supplies the common adjoint scoring contract for
reactor/PIC/diagnostic distributions. The unified boundary transport passes the same vertical-particle gate
at AR 1/4/16 with geometry as the only change. This is an architecture/no-tunnelling invariant, not evidence
that charging, collisions, chemistry, or profile evolution are validated at AR16.

A generic boundary-state adjoint floor gather now passes the same open-surface normalization for charged
Ar+ and neutral CF2. It contains no species source equation and applies the common density ratio and
Liouville normal-velocity Jacobian without weight clipping. Nonuniform-field and arbitrary-face gates are
still required before production migration.

The ion source audit also separated reference emulation from first principles. The formerly hard-coded
`Vs^-0.35` phase weight came from the Hwang-Giapis simulated IEDF horn ratio. The backward core now uses
uniform RF phase by default; `ion_ied_phase_exponent=0.35` is passed explicitly only by the HG benchmark
script. Experimental/reference curves must not silently define the production transport law.
With analytic uniform phase, the W32 AR4 iteration-10 floor is 38.6 V rather than the benchmark-shaped
value near 34 V. That movement is honest: restoring the high-energy IEDF horn increases penetrating ion
current. A time-dependent sheath model must determine the real phase distribution from physical inputs.

### 2. Adjoint reciprocity is asserted more broadly than it is gated

An independent forward-launch/backward-gather gate now passes on a refined trench at zero field and with
a physically consistent +9.2 V floor potential. Coarse point-collocated versions differed by 16% and 6%,
respectively, so the gather now integrates deterministic Sobol samples across each finite face. Grid
convergence is still required and cannot be hidden by sample count. The
ion reciprocity also passes at zero and +5 V retarding floor potential. The next reciprocity gate needs
a nonuniform self-consistent field, with confidence/error estimates.

### 3. Pattern electrostatics is the first direct attack on academic state of the art

Huang et al. identify two distinct observables: symmetric arrays have zero ensemble-mean systematic tilt,
while asymmetric dense/sparse patterns create a mean horizontal field from dense to sparse. A frozen-pattern
field/sign gate is cheaper and more diagnostic than immediately evolving a twisting hole.

### 4. Surface chemistry needs state, not more scalar knobs

The Vella/Graves work shows that top-layer coverage, subsurface mixing, damage depth, and fluence history
matter. Existing ALE state is zero-dimensional. The correct lift is a conserved surface-state field whose
uniform-flux limit reproduces the current solver.

## New experimental/data targets from the 2026 search

### Acquired: Bosch reactor/wafer dataset (2025)

Sayyed et al., Zenodo 17122442, CC BY 4.0:

- synchronized OES (3648 wavelengths at 25 Hz), 31 machine variables at 5 Hz, and wafer measurements;
- 76 wafers in the nine-point table, with silicon depth and oxide-mask loss;
- designed around chamber conditioning and within-lot drift;
- useful for reactor-state inference, selectivity, uniformity, and held-out wafer prediction;
- not a feature-profile or charging validation dataset.

The small nine-point measurement table is vendored under `data/experimental/zenodo_17122442/` with its
source checksum and a strict loader. Larger NetCDF files remain remote until a reactor experiment is
defined.

### High-value public feature-profile experiment: VLSet-AE DRIE (2025/2026)

The published experiment reports a 16-run orthogonal Bosch design and 1000 cross-sectional SEM images,
with etch/passivation times, profile angle, scallop depth/width, and trench depth. The paper exposes the
condition table and aggregate measurements; availability and reuse terms for the image corpus must be
confirmed before ingestion. This is a strong held-out profile-transfer target because petch already has a
Bosch module and the recipes vary physical cycle times.

### High-value experimental-profile paper: autonomous hybrid optimization (2024)

Krüger et al., JVST A 42, 043008 (2024), couples HPEM/MCFPM optimization to experimental SEMs. The paper
reports an unusually valuable failure map:

- O2/C4F6 variation reproduces clogging/necking trends but misses mask erosion quantitatively;
- maximum experimental depth near O2/C4F6 = 1.5 is only marginally reproduced;
- doubling low-frequency power from 4 to 8 kW changes experimental depth little, indicating neutral
  transport limitation above the ion-energy threshold;
- simulation depth still increases (635, 715, 720 nm), suggesting ion-energy influence is overestimated;
- the authors explicitly condition their conclusions on HPEM wafer fluxes being accurate.

This is nearly the exact commercial wedge: fit physical parameters on one condition, predict trends under
new gas ratio and power, and expose missing mechanisms instead of silently refitting every structure.

### Public 3-D reactor/etch profiles: magnetic sheath tailoring (2024)

Jüngling et al. release CC-BY-4.0 profilometer traces for silicon and a-C:H with and without transverse
magnetic field, paired with PIC/MCC work. This is a reactor/sheath-to-profile coupling target rather than
a conventional vertical-feature etch gate. It becomes useful after `PlasmaBoundaryState` exists.

### Public materials breadth: binary III-V etch-rate compilation (2025)

Clarke releases a CC-BY-4.0 literature-derived dataset of etch rates, materials, techniques, chemistries,
and sources. It is appropriate for mechanism/data coverage analysis and experimental design, but heterogeneous
literature conditions make it a poor direct solver gate unless subsets have comparable reactor conditions.

## State-of-the-art academic reading

The strongest current academic pattern is not “replace physics with AI.” It is:

1. high-fidelity reactor simulation or diagnostics produce wafer flux distributions;
2. feature Monte Carlo evolves chemistry and geometry;
3. physical parameters are optimized against multiple feature observables;
4. transfer to unseen process conditions tests whether parameters are mechanisms or mere fit coefficients;
5. discrepancies identify missing reactions, transport limits, or source errors.

The 2026 Kushner SEE paper sharpens the charging frontier: secondary electrons redistribute charge and
lower in-feature positive potential, but can reduce the incremental benefit of anisotropic high-energy
electrons because charge redistribution becomes dominant. This reinforces the need to treat SEE as a
conservative particle/charge channel, not a voltage correction.

For petch, the opportunity is to execute this loop deterministically, with verified gradients, explicit
uncertainty, open mechanisms, and far lower compute—then demonstrate transfer on data that were not used
to fit the model.

## Rules for adding experimental data

1. Record DOI/URL, authors, license, acquisition date, source checksum, units, and any digitization method.
2. Preserve raw source data unchanged; transformations live in scripts and generate derived artifacts.
3. Never silently impute missing provenance or measurements.
4. Separate calibration and held-out splits before parameter selection.
5. Track experimental uncertainty, digitization uncertainty, numerical uncertainty, and model discrepancy
   separately wherever the source permits.
6. Do not redistribute copyrighted figures or paper text; store structured measurements only when license
   or applicable data reuse permits, otherwise store citations and reproducible digitization instructions.
