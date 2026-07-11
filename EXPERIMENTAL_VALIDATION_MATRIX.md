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
