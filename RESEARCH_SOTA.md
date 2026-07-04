# SOTA of feature-scale plasma-etch charging simulation (research pass 2026-07-03)

Motivating question: are we right to build a GENERAL material-grid particle-in-field charging
engine, or is the hardcoded single-feature ("edge array") cell justified? Three parallel
deep-research passes over ViennaPS, Kushner's MCFPM lineage, and the commercial/academic/open
landscape. Bottom line up front: **the engine should be general; the single-feature *domain* is a
legitimate, near-universal simplification. Those are different things, and we conflated them.**

## The one distinction that matters

- **General ENGINE** = the solver takes geometry as *data* (a material-tagged grid or a level-set)
  and runs the same physics on any arrangement. This is universal SOTA. Nobody hardcodes feature
  shapes into the solver.
- **Single-feature DOMAIN** = the *problem* is reduced to one representative feature (a half-trench,
  a 2D infinite-trench cross-section, or a periodic unit cell with reflecting/periodic BCs). This is
  near-universal and *well-motivated*: real device arrays (gratings, DRAM, 3D-NAND) genuinely are
  periodic; full-die 3D at fine voxel resolution is tens of hours per feature; long trenches are
  translationally invariant; and the canonical benchmarks (Coburn-Winters, de Boer, Hwang-Giapis)
  are themselves single or periodic structures. (Huard thesis periodic sub-array; half-trench Cl2
  studies; DRIE replicated-domain 20h/20cycles.)

Our mistake was **not** using a single periodic edge-array cell — that is standard and fine. It was
**hardcoding that geometry into the engine** (column indices `edge0/edge1`, tracer labels like
"hit == 3 means edge outer wall"). The fix is a general engine that takes the edge-array as input
(`charging_general.py`), then runs the same periodic single-cell problem. Instinct validated;
precise target clarified.

## What each SOTA code actually does

### Kushner MCFPM (Michigan) — the reference charging engine
Voxel / cubic-cell grid, **one material tag per cell**. Charged pseudoparticles deposit charge on
tagged cells; the field is obtained by **implicitly solving Poisson's equation ∇·ε∇Φ = ρ with
finite volume + red-black SOR**; the resulting field bends subsequent charged trajectories via a
Lorentz force in the equations of motion.
- **Conductors/insulators/dielectrics are NOT equipotential hacks** — every cell carries a
  dielectric constant ε and charge mobilities μ±. Conductor = high mobility (charge relaxes/spreads);
  insulator = low mobility (charge stays local); ε varies as material ids change. Behaviour emerges
  from per-cell material properties. (Huang, JVSTA 37, 031304, 2019.)
- Poisson updated only every ~400 charged-particle strikes (cost control); millions launched.
- Captures **inter-feature field penetration** (charge in one feature bends ions in the next) and
  pattern-dependent systematic tilt (Huang, JVSTA 38, 023001, 2020) — the general form of the exact
  edge/neighbour split we've been chasing.
- Initialized from "arbitrary" multi-material geometry; same engine runs finFET gate etch, HAR
  contacts, DTI, focus-ring charging, SEE (2024 GEC). This is our roadmap.

### ViennaPS (TU Wien) — general geometry + GPU flux, but NO charging
Sparse-field **level-set** surface (ViennaLS/HRLE) advected by fluxes from **top-down Monte-Carlo
ballistic ray tracing** (ViennaRay, Embree on CPU / **NVIDIA OptiX on GPU**). General geometry, no
hardcoded features.
- **No charging module exists**: its 27 process-model headers contain zero charging/potential/
  Laplace/electron code; the SF6O2 paper states the model is kinetic chemistry "without
  electrostatic charging effects."
- **Why (structural):** a level-set is a scalar signed-distance field that "does not hold
  directional information" and carries no volumetric state (charge, potential). A charging field
  solve needs a *separate* volumetric representation — hence ViennaCS (a voxel-like cell set) exists
  alongside it. Level-sets are great at moving surfaces (topology change, pinch-off), bad at holding
  volumetric charge.
- Charging on the ViennaPS/level-set lineage is a *research add-on* by an independent group
  (Zhai, …, Filipovic, JAP 137, 063302, 2025), not the core release.

### Kokkoris/Gogolides (NTUA), Radjenović (Belgrade)
Level-set surface + Monte-Carlo ion/electron flux + **Laplace/Poisson solved by finite elements**
for the in-feature potential. General surface, single-trench setups. Confirms the level-set +
separate-FEM-field-solve pattern.

### Commercial
- **Coventor/Lam SEMULATOR3D**: physics-driven **voxel** engine + hybrid level-set; MC sheath flux
  with shadowing. No documented charging.
- **Synopsys Sentaurus Topography**: **level-set** surface + MC particle tracking; a Synopsys patent
  also describes a discrete-VOF **voxel** engine. No documented charging.
- **Silvaco Victory Process**: **level-set** on Cartesian meshes + MC etch; ions as "charged,
  accelerated, non-reacting particles" — no local charging field.
- No commercial tool documents a feature-charging / electrostatic ion-deflection model (medium
  confidence; some datasheets were image-only).

### Hwang-Giapis 1997 (our benchmark) — the obsolete predecessor
Static single 2D trench, **Laplace** (charge-free interior), hardcoded **equipotential** conductors,
sidewall reactions frozen. This is exactly the bespoke single-feature model the field moved *away*
from. Reproducing its Fig. 3/4/6 is our validation target, not our architecture.

## The frontier — and petch's defensible whitespace

Three documented cutting edges: (1) GPU/OptiX ray tracing on level-set surfaces for flux (ViennaPS,
Riedel 2026); (2) self-consistent in-feature charging MC (Kushner MCFPM, Krüger thesis 2024); (3)
ML surrogates replacing the profile evolver (CAE+RNN, cascade-RNN).

Documented **whitespace** (things that do NOT exist in the literature/repos):
- **No open-source feature-scale CHARGING code exists.** ViennaPS — the only production open
  framework — has no charging. The capability lives only in closed tools (MCFPM, SEMULATOR3D,
  Silvaco). Genuine open-source gap.
- **No fully differentiable / adjoint feature-scale etch-profile solver exists.** The clearest gap;
  adjoint methods are mature only in adjacent fields.
- **No mature GPU in-feature charged-particle transport** (GPU flux exists; GPU *charging* MC does
  not).

petch hits all three at once: a **general MCFPM-style voxel charging engine (Poisson + per-cell
ε/μ), GPU-resident (Warp), and differentiable (autodiff)**, coupled to petch's existing level-set
for profile evolution. That is a defensible, literature-confirmed position — not another hand-built
trench.

## Concrete architecture recommendation for petch

1. **Keep the general engine** (`charging_general.py` direction): material-tagged grid, geometry is
   input data. The edge array / any device is just an input.
2. **Upgrade the field solve to the MCFPM formulation**: replace Laplace + connected-component
   equipotential conductors with **Poisson ∇·ε∇Φ = ρ, per-cell ε and mobility μ±**. Conductor/
   insulator/dielectric become material properties, not special cases. This also removes the
   equipotential awkwardness and is what lets charge redistribute physically.
3. **Trace both species through the field** (ions directional, electrons isotropic) — the
   electrostatic focusing HG demands is then automatic (the piece the view-factor shortcut lost).
4. **Single representative / periodic feature domain is fine** for validation (it's standard); use
   periodic/reflecting BCs as MCFPM does.
5. **GPU-port the tracer to Warp** (61× precedent from the 3D `_trace3d` kernel) once the physics
   converges; keep it **differentiable**.
6. **Couple to the existing level-set** for profile evolution (level-set surface ⊕ voxel charging is
   how the field splits the two representations; petch already has the level-set half).

## Primary sources
- Huang et al., "Plasma etching of HAR features in SiO2 (Ar/C4F8/O2)," JVSTA 37, 031304 (2019) —
  MCFPM voxel + Poisson(FV, red-black SOR) + per-cell ε/μ charging.
- Huang et al., "Pattern dependent profile distortion ... HAR SiO2," JVSTA 38, 023001 (2020) —
  inter-feature field penetration, pattern tilt.
- Krüger et al., VWT series, JVSTA 41, 013006 (2023); Phys. Plasmas 31, 033508 (2024).
- Litch et al., JVSTA 43, 033001 (2025) — focus-ring/DTI charging. Kushner GEC 2024 — SEE in MCFPM.
- Reiter & Filipovic, "ViennaPS," SoftwareX 32, 102453 (2025); ViennaPS/ViennaLS/ViennaRay GitHub;
  Manstetten PhD 2018 (level-set limitation: scalar field, no directional/volumetric state).
- Zhai, …, Filipovic, JAP 137, 063302 (2025) — charging on a level-set etch model (independent).
- Kokkoris & Gogolides (NTUA); Radjenović & Radmilović-Radjenović (Belgrade) — level-set + FEM Poisson.
- Giapis/Hwang (Caltech), Economou (Houston), Fujiwara/Nozawa (Panasonic), Matsui (TEL),
  Hamaguchi (Osaka) — single-trench MC + charge-site/cell arrays.
- Hwang & Giapis, "Aspect-ratio-dependent charging," JAP 82, 566 (1997) — our benchmark; static
  single trench, Laplace, equipotential conductors.
- Reviews: Oehrlein/Kanarik et al., JVST B 42, 041501 (2024); JJAP 63 methods review (2024,
  "voxel + MC dominant; level-set for etch-front; smart voxels"); Coventor SEMULATOR3D; Synopsys
  Sentaurus Topography; Silvaco Victory Process. Open repos: ViennaTools/ViennaPS (no charging),
  cococastano/PythonDryEtchModel, sparta/sparta (DSMC, no charging).
