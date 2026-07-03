# Robust physics model plan: Kushner + Graves pass

Research pass: 2026-07-03. Scope: use the Kushner HPEM/PCMCM/MCFPM lineage and
Graves/Princeton/PPPL plasma-surface work to decide what petch should build next, without tuning
the Hwang-Giapis curves by hand.

Working rule: no new scalar knobs to match HG. Each added model must have a published source,
an isolated reproducer, and a numerical gate before it is allowed into the coupled charging or
chemistry path.

## 1. What the literature says

### Kushner layer: source + feature transport + charging

Kushner's stack is not a single-trench analytic-source model. It is a coupled, staged model:

1. HPEM computes reactor-scale plasma properties, charged-particle fluxes, and wafer energy-angle
   distributions. In the 2026 SEE charging paper, HPEM feeds wafer fluxes and EADs to MCFPM, which
   simulates the feature charging process.
2. PCMCM is the Monte-Carlo transport layer above the wafer. The official CPSEG page says PCMCM uses
   HPEM fields, sheath potential/thickness, chemistry, and species distributions to produce energy and
   angular distributions at surfaces. These are exactly the boundary distributions petch currently
   approximates analytically.
3. MCFPM is a 2-D/3-D voxel kinetic feature model. The official CPSEG page lists the needed inputs:
   material mesh, species/chemistry file with energy and angular dependence, and gas species flux
   distributions. The 2024 VWT paper states that MCFPM launches pseudoparticles with energies and
   angles sampled from HPEM, advances charged trajectories in feature-charging fields, and then uses
   the incident particle energy/angle to select a surface reaction.
4. Kushner's 2020 pattern-distortion paper says pattern charging is nonlocal: asymmetric dense/open
   patterns create horizontal fields from dense, more positive regions toward sparse/open regions.
   These fields act on ions in adjacent features. That is directly relevant to the HG edge-line/open-area
   problem: a single periodic trench cannot represent the published mechanism.
5. Kushner's 2026 SEE charging paper says SEE redistributes charge inside HAR features, lowers positive
   potential, and can dominate the effect of anisotropic electron remediation. This supports keeping SEE
   as a first-class transport/charge channel, not a post-hoc voltage correction.
6. The MCFPM feature model does not stop at a scalar charging table: charged pseudo-particles deposit
   charge on material-tagged feature surfaces, the feature electric field is updated, and later charged
   trajectories are bent by that field. The current petch 3-D HG hook is a table closure; the end state
   is one shared electrostatic interface that the Monte Carlo, radiosity, Knudsen, and DDA paths can all
   query without duplicating charging logic.

### Graves layer: surface state, not just a scalar yield

Graves' MD work points to a richer surface-response state:

1. The 2009 Graves-Brault review frames MD as the way to resolve bond breaking/forming, ion impact
   chemistry, products, damage, roughness, and energy transfer. The review also makes clear why MD is
   not an inline feature solver: it must be amortized into yield/state laws.
2. The 2025 DeepMD Si-Cl-Ar paper validates a learned interatomic potential against the exact objects
   petch needs: Ar+ physical sputter yields, amorphous-crystalline interface depth, Cl2/Ar+ etch yield
   vs flux ratio and ion energy, ALE etch products, and spontaneous Cl etch probability. It uses DP-GEN
   active learning plus a ZBL short-range correction for energetic collisions.
3. The 2025 Si-Cl2-Ar+ ALE-window paper uses MD plus a reduced-order transient site-balance model.
   The ROM splits the near surface into top, mixed, and crystalline layers; evolves top-layer Cl coverage
   during chlorination; evolves top/mixed-layer Cl during ion bombardment; and predicts product fluxes.
   It finds a narrow normal-incidence ALE window around 15-20 eV, sub-monolayer EPC around the
   self-limiting regime, and strong fluence dependence outside the window.

## 2. Implications for petch

The current edge-array WIP is the right direction but not the final architecture:

- Keep: explicit edge/open/nonperiodic geometry, separate edge and neighboring conductor potentials,
  source launched at the feature mouth, right-hand neighboring trench, and the boolean-index clip fix.
  These are physical corrections, not fit knobs.
- Treat as diagnostic only for now: the line-of-sight edge boundary. Kushner pattern papers point to
  explicit pattern electrostatics, not a scalar auxiliary current. The top-up formulation is safer than
  additive current, but the final model should prefer explicit multi-feature pattern domains.
- Do not claim the HG gates are closed from the edge-array WIP. Current reduced runs still miss at least
  the floor-flux RMSE and/or neighbor-poly potential gates, depending on source model and statistics.
- Do not replace the chemistry with a black-box network first. Graves' work says the first useful step is
  a physics-state interface that can accept MD/DeepMD tables later.
- Do not promote the 2-D diagnostic closure into the 3-D production path. Production needs material
  electrical labels, charge deposition, a Poisson/Laplace residual check, conductor equipotentials, and
  charge-conserving SEE before the table closure can be retired.

## 3. Implementation order

### W0 - Lock current WIP numerics and diagnostics

Files:

- `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/src/petch/charging2d.py`
- `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/scripts/charging_gate.py`
- `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/scripts/notching_gate.py`
- `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/tests/test_charging_edge_open.py`

Actions:

1. Commit the real edge-array fixes locally after tests: correct boolean-index clipping, source at
   feature mouth, right-side neighboring trench, and explicit-vs-boundary current accounting.
2. Add a reduced gate runner that prints the exact env/config into the output and refuses to overwrite
   saved figures unless `PETCH_WRITE_RESULTS=1`.
3. Add diagnostics for edge, neighbor, PR, trench-floor, open-floor, and right-trench current residuals.
4. Add a confidence interval over seeds for each gate metric. A pass must survive seed spread.

Kill criterion: if reduced gates remain residual-limited above 0.08 after the clip/source fixes, do not
run high-stat gates. Fix charge conservation first.

Current local reduced benchmarks, no-write, 2026-07-03:

- `python scripts/charging_reduced_gate.py`
  (`edge_array`, `line_of_sight`, `source_model=analytic`, `W=16`, `mouth=80`, `n=1200`, `it=100`):
  floor RMSE `0.077` fail, survivor max `0.0000` pass, residual max `0.329` fail, foot-energy
  max relative error `0.592` fail, foot-flux ratio `1.45` pass, poly-potential max relative error
  `0.504` fail, Matsui pass, 0-D closure pass.
- `PETCH_SOURCE_MODEL=sheath_electrons python scripts/charging_reduced_gate.py`
  (same geometry/statistics): floor RMSE `0.069` fail, survivor max `0.0000` pass, residual max
  `0.437` fail, foot-energy max relative error `0.251` pass, foot-flux ratio `1.38` pass,
  poly-potential max relative error `0.386` fail, Matsui pass, 0-D closure pass.

Interpretation: particle loss is not the current limiter at this reduced setting; current residual and
neighbor-poly/open-pattern electrostatics are. Do not spend GPU/Vast time on high-stat repeats until the
explicit charge-conservation and pattern-domain pieces are implemented.

### W1 - Replace analytic wafer source with a table source

Files:

- `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/src/petch/sheath1d.py`
- `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/src/petch/charging2d.py`

Actions:

1. Define an explicit boundary-table interface: species, charge, weight, energy, theta, RF phase or
   sheath snapshot, and source-plane x distribution.
2. Make current analytic and sheath samplers write/read this table format. No physics change first.
3. Add an HPEM/PCMCM-compatible import path for future tabulated EADs.
4. Gate table replay against current analytic/sheath samplers bit-for-bit statistically.

Gate: table replay changes HG metrics by less than Monte-Carlo seed uncertainty.

### W2 - Explicit pattern electrostatics

Files:

- `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/src/petch/charging2d.py`
- new tests under `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/tests/`

Actions:

1. Generalize edge-array geometry to a small pattern-domain builder: open field, edge line, one or more
   trenches, neighboring lines, optional symmetric/dense/open patterns.
2. Solve all conductors as connected components; solve insulators as local surface charges.
3. Remove scalar edge boundary current from the production path once explicit open field and neighboring
   trenches are in the domain.
4. Gate against Kushner pattern facts: symmetric patterns have no systematic lateral field on average;
   asymmetric/open-side patterns create a dense-to-open horizontal field and systematic edge tilt.

HG gate expectation: neighbor poly potential and foot energy should become geometry-derived rather than
boundary-current-derived.

### W3 - SEE and charge-transfer cascades in the explicit pattern model

Files:

- `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/src/petch/charging2d.py`

Actions:

1. Port the existing PR-sidewall PMMA SEE branch into the edge/pattern tracer.
2. Add material-tagged SEE curves for SiO2 and poly-Si only from published curves; keep PR-only as the
   first isolated gate.
3. Track emitted-electron source and landing material separately. Charge conservation must include
   emission site positive charge and landing site negative charge.
4. Reproduce the literature effect direction: SEE must lower positive potential by redistributing
   electrons inside the feature. No rescaling to hit HG.

Gate: in an AR4 fixed-geometry run, SEE must reduce positive potential and conserve total charge within
statistical error. If it increases potential or only changes results below seed noise, stop and diagnose.

### W4 - Graves-style surface state layer for Cl2/HBr and ALE

Files:

- `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/src/petch/chemistry.py`
- new plugin/state files under `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/src/petch/`

Actions:

1. Add a `SurfaceState` abstraction with at least: material id, top-layer Cl coverage, mixed-layer Cl
   coverage, mixed/amorphous depth, accumulated Cl and ion fluence, product counters, deposited film
   thickness, and optional roughness moments.
2. Implement Cl2/poly-Si as the first non-SF6 chemistry using published Chang-Sawin yields as the
   immediate gate, but design the state shape to accept Graves/DeepMD tables.
3. Implement Si-Cl2-Ar+ ALE as a cycle scheduler and transient site-balance model:
   top-layer chlorination, mixed-layer state, ion-bombardment product yields, and fluence saturation.
4. Gate against Graves/Vella: ALE window around 15-20 eV for normal incidence Ar+, sub-monolayer
   self-limiting EPC, product partition, and fluence saturation.

Quantitative anchors to digitize:

- DeepMD atomic-Cl paper: thermal Cl exposure at 300 K gives steady Cl coverage near 1.25 ML.
- Damage/mixing papers: start/end Si yields under Ar bombardment are about 0.017/0.002 at 25 eV,
  0.115/0.039 at 100 eV, and 0.257/0.148 at 215 eV; 70 eV Ar+ produces about a 1 nm amorphous region.
- ALE products: early Ar fluence is dominated by Cl and SiClx removal; late fluence trends toward
  physical Si sputter, with atomic Si staying near the physical sputter yield.
- DP-GEN provenance: four-model ensemble, force-deviation uncertainty, candidate bounds around
  0.1-0.5 eV/A, final training set 43,682 frames, and ZBL interpolation for short-range ion collisions.

This is independent of HG charging closure, but it is the path to real Resona-relevant chemistry and
future MD-trained operators.

### W5 - 3-D charging/electrostatics production interface

Files:

- `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/src/petch/threed.py`
- new shared charging/electrostatics module under
  `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/src/petch/`

Actions:

1. Add material/electrical labels to the 3-D surface representation: grounded conductor, floating
   connected conductor, dielectric charge storage, resist, oxide, poly-Si, and open boundary.
2. Define one shared charging interface: deposit charged-particle current, solve the electrostatic
   update, expose field interpolation to charged-particle tracers, and report charge-conservation plus
   Poisson/Laplace residual diagnostics.
3. Route all transport modes through the same interface. The existing HG table can remain as a fallback
   validation closure, but `mc`, `radiosity`, `knudsen`, and `dda` must not each grow separate charging
   behavior.
4. Gate simple 3-D electrostatics before chemistry coupling: grounded conductive Si floors bypass
   dielectric charging; symmetric patterns have zero seed-averaged lateral tilt; asymmetric open/dense
   patterns produce the Kushner-sign horizontal field.

Kill criterion: if the 3-D solver cannot pass charge conservation and residual checks on fixed geometry,
do not couple it to level-set evolution. Fix the electrostatics first.

## 4. Gate matrix

| Layer | Gate | Pass condition |
|---|---|---|
| W0 current conservation | residual max | < 0.08 reduced, < 0.05 high-stat target |
| W0 survivors | ion/electron survivor frac | < 0.001 |
| HG floor flux | RMSE vs 8 HG points | <= 0.05 |
| HG AR4 potential | Vc | 33 V +/- 40%, with no floor-flux regression |
| HG notching | foot energy | rising 15 -> 27.5 eV, max rel err <= 30% |
| HG notching | foot flux | max/min <= 2 for AR >= 1.6 |
| HG poly line | neighbor/line potential | rising 6 -> 39 V, max rel err <= 30% |
| Kushner pattern sanity | symmetric pattern | no systematic lateral field after seed average |
| Kushner pattern sanity | open/dense asymmetry | dense-to-open horizontal field, edge-line tilt direction correct |
| SEE sanity | charge transfer | total charge conserved; positive potential lowered in comparable HAR geometry |
| Graves Cl chemistry | continuous Cl2 etch | yield vs energy/flux ratio/angle within digitization tolerance |
| Graves ALE | cycle model | 15-20 eV window, sub-monolayer EPC, fluence saturation |
| 3-D electrostatics | fixed-geometry charge solve | charge conserved; Poisson/Laplace residual reported and below tolerance |
| 3-D pattern charging | symmetric/asymmetric arrays | no systematic symmetric tilt; dense-to-open field sign correct |

## 5. Sources read in this pass

- Mark J. Kushner group, "Monte Carlo Feature Profile Model (MCFPM)," CPSEG:
  https://cpseg.eecs.umich.edu/Projects/MCFPM/MCFPM.htm
- Mark J. Kushner group, "Plasma Chemistry Monte Carlo Model (PCMCM)," CPSEG:
  https://cpseg.eecs.umich.edu/Projects/PCMCM/PCMCM.htm
- Chenyao Huang and Mark J. Kushner, "Consequences of secondary electron emission on charging of SiO2
  features in capacitively coupled plasmas having sinusoidal and tailored bias waveforms," JVST A 44,
  023013 (2026): https://cpseg.eecs.umich.edu/pub/articles/JVSTA_44_023013_2026.pdf
- Florian Krueger, Hyunjae Lee, Sang Ki Nam, and Mark J. Kushner, "Voltage waveform tailoring for high
  aspect ratio plasma etching of SiO2 using Ar/CF4/O2 mixtures," Phys. Plasmas 31, 033508 (2024):
  https://cpseg.eecs.umich.edu/pub/articles/PhysPlasmas_31_033508_2024.pdf
- Huang, Shim, Nam, and Kushner, "Pattern dependent profile distortion during plasma etching of high
  aspect ratio features in SiO2," JVST A 38, 023001 (2020):
  https://www.osti.gov/biblio/1802573
- Kounis-Melas, Vella, Panagiotopoulos, and Graves, "Deep potential molecular dynamics simulations of
  low-temperature plasma-surface interactions," JVST A 43, 012603 (2025):
  https://www.osti.gov/biblio/2514378
- Kounis-Melas, Panagiotopoulos, and Graves, "Deep-potential molecular-dynamics simulations of
  ion-enhanced etching of silicon by atomic chlorine," JVST A 43, 063204 (2025):
  https://pubs.aip.org/avs/jva/article/43/6/063204/3368098/Deep-potential-molecular-dynamics-simulations-of
- Vella and Graves, "Si-Cl2-Ar+ Atomic Layer Etching Window: A Fundamental Study Using Molecular
  Dynamics Simulations and a Reduced Order Model," J. Phys. Chem. B (2025):
  https://www.osti.gov/servlets/purl/2586627
- Vella, Hao, Donnelly, and Graves, "Dynamics of plasma atomic layer etching: Molecular dynamics
  simulations and optical emission spectroscopy," JVST A 41 (2023):
  https://www.osti.gov/biblio/2248044
- Vella, Hao, Elgarhy, Donnelly, and Graves, "A transient site balance model for atomic layer etching,"
  Plasma Sources Sci. Technol. 33 (2024):
  https://www.osti.gov/pages/biblio/2406001
- Vella and Graves, "Near-surface damage and mixing in Si-Cl2-Ar atomic layer etching processes:
  Insights from molecular dynamics simulations," JVST A 41 (2023):
  https://www.osti.gov/biblio/1999799
- Graves and Brault, "Molecular dynamics for low-temperature plasma-surface interaction studies,"
  J. Phys. D (2009): https://arxiv.org/abs/0902.2695
