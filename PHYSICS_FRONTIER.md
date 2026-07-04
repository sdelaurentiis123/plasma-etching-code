# The physics frontier: MD-derived surface-state chemistry (Graves) as petch's moat

Research pass 2026-07-03 (David Graves program + the atomistic->feature-scale gap). The engineering
frontier (GPU + differentiable + open-source charging) is necessary but copyable. The *durable* moat
is owning a PHYSICS frontier no incumbent has folded into a full solver. This is that frontier.

## The gap, stated sharply

**Every feature-scale engine uses INSTANTANEOUS EMPIRICAL surface chemistry.** ViennaPS and Kushner
MCFPM both run a steady-state Langmuir-Hinshelwood site balance with empirically-calibrated
sticking/yield coefficients (Belen SF6O2, Chang-Sawin). The surface is assumed to be in *instant
chemical steady state*: no memory of fluence history, no evolving mixed-layer/damage depth, no ALE
transient. That is structurally unable to represent:
- ALE self-limiting saturation (the plateau that defines atomic-layer etching),
- the etch yield *decaying* from a chemically-enhanced peak toward the physical sputter yield as an
  amorphous mixed layer builds (Graves MD proves this is real),
- damage/amorphization depth accumulation,
- energy-AND-angle-resolved yields + product partitioning.

**Graves' group generates exactly this physics — but stops at 0-D.** DeepMD/ML-potential MD now
gives energy-angle-resolved Si etch yields, SiClx product partitioning, mixed-layer thickness, and
1.25-ML steady Cl coverage (Kounis-Melas 2025). It's amortized into reduced-order models: a
transient two-ODE Cl site-balance (Vella 2024), a depth-resolved diffusion-reaction Cl field (2026),
an ALE-window ROM (2025), and a physics-constrained **Neural Master Equation** whose transition
rates are NNs trained on MD (Nath-Vella-Graves-Mesbah 2025). But every one of these is **0-D or
1-D-in-depth reactor-observable physics — never coupled to a spatial feature-scale profile or to
charging.**

**Nobody bridges them, and it's institutionally hard.** The MD/ML-potential community (Graves) stops
at 0-D observables; the feature-scale/charging community (ViennaPS, Kushner) stays with empirical
chemistry; and neither builds autodiff solvers. Differentiable simulation is thriving in adjacent
fields (JAX-BTE phonons, differentiable Maxwell, ICF) but has never been applied to plasma etch.
Owning BOTH sides of that scale gap, differentiably, is the moat.

## petch's unique position

A **differentiable, GPU, multiscale solver** spanning: MD/ML-potential surface chemistry ->
per-voxel surface-state memory -> feature-scale transport + charging -> profile evolution. Because
it's differentiable end-to-end, you can INVERT/CALIBRATE the whole chain against experiment (fit
MD-informed chemistry + charging jointly to measured profiles via autodiff) -- something no
incumbent can do. petch already has the charging engine (this repo), the level-set profile, the GPU
Warp kernels, and autodiff. The missing piece is the surface-state chemistry layer.

## Modules to own (each with its Graves anchor)

- **(a) MD/ML-potential energy+angle-resolved yield & product tables** — Y(E, theta, flux-ratio,
  coverage) + product-branching vector, replacing scalar Belen/Chang-Sawin. Anchor: Kounis-Melas,
  JVST A 43(6) 2025 (OSTI 3001885); seed numbers: Cl coverage ~1.25 ML, yield decays to physical
  sputter yield.
- **(b) Per-voxel SURFACE-STATE layer** {top-Cl coverage, subsurface-Cl, amorphous/damage depth,
  cumulative fluence} advanced by an MD-parameterized ODE/PDE each timestep. THE biggest departure
  from ViennaPS/MCFPM. Anchors: Vella 2024 transient site-balance (OSTI 2406001); Vella 2026
  diffusion-reaction depth-resolved Cl (JVST A 44(2) 022602); Vella-Graves 2023 damage/mixing
  (OSTI 1999799).
- **(c) ALE cycle scheduling with self-limiting saturation** — native mod/removal half-cycles
  reproducing the saturation plateau + synergy %. Anchors: ALE-window ROM (JPCB 2025); Kanarik
  synergy framing (surface binding E0 ~4.7 eV Si).
- **(d) SEE / electron-driven surface chemistry** — most relevant to CHARGING coupling, but the
  weakest-anchored (Graves names SEE as important but has no recent MD SEE-yield dataset; needs
  external SEE data). Higher risk.
- **(e) Roughness/damage evolution** per-voxel, coupling back into yield. Anchor: bias-pulsed SiC
  ALE 99% synergy / sub-angstrom roughness (JVST A 41(3) 032607 2023).

## First buildable piece (de-risked)

Replace the instantaneous Langmuir-Hinshelwood site balance with the **Vella 2024 two-ODE transient
Cl site-balance** (top-monolayer Cl coverage + perfectly-mixed subsurface Cl, advanced vs fluence),
seeded with the **Kounis-Melas 2025 energy-angle yield table and 1.25-ML coverage**. It's a
0-D-per-voxel surface state that the existing charging framework can already carry, and it
immediately gives ALE saturation + damage memory that no competitor's chemistry model can express.
Because petch is differentiable, the ODE rate parameters become fit/invertible against experiment.

## The full-solver vision (what no one else has)

```
 MD / ML-potential (DeepMD, fine-tuned MACE)   <- physics source (Graves)
        |  energy-angle yields, product branching, rate params
        v
 per-voxel SURFACE STATE (coverage, mixed-layer, damage depth, fluence)   <- (b), the moat
        |  local etch rate + charge deposition
        v
 feature-scale TRANSPORT + CHARGING (Poisson eps/mu, ion+electron MC)   <- this repo (general engine)
        |  surface velocity + fields
        v
 LEVEL-SET profile evolution   <- petch already has it
        |
        v  (all GPU/Warp, all autodiff)
 END-TO-END DIFFERENTIABLE: invert the whole chain against measured profiles
```

Incumbents have pieces: ViennaPS (level-set + GPU flux, empirical chemistry, no charging, not
differentiable); MCFPM (voxel + charging, empirical chemistry, not differentiable); Graves NME
(differentiable MD-chemistry, but 0-D, no transport/charging/profile). No single solver spans the
whole chain differentiably. That is the thing to build.

## Primary sources
- Kounis-Melas, Vella, Panagiotopoulos, Graves, JVST A 43(1) 2025 (OSTI 2514378) — DeepMD Si-Cl-Ar.
- Kounis-Melas, Panagiotopoulos, Graves, JVST A 43(6) 2025 (OSTI 3001885) — energy-angle yields, 1.25 ML.
- Vella, Hao, Elgarhy, Donnelly, Graves, PSST 33(7) 2024 (OSTI 2406001) — transient two-ODE site balance.
- Vella et al., JVST A 44(2) 022602 2026 — diffusion-reaction depth-resolved Cl.
- Vella & Graves, JVST A 41(4) 042601 2023 (OSTI 1999799) — near-surface damage/mixing.
- Nath, Vella, Graves, Mesbah, npj Comput. Mater. 11(1) 2025 (10.1038/s41524-025-01677-4) — Neural Master Equation.
- 2026 arxiv 2606.21632 — fine-tuning MACE-MP-0 (universal MLIP) for plasma-surface chemistry.
- Kanarik et al., JPC Lett. 2018 (10.1021/acs.jpclett.8b00997) — ALE synergy, surface binding energy.
- (Adjacent differentiable-sim proof) JAX-BTE arxiv 2503.23657; ICF arxiv 2606.08827.
