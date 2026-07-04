# petch: deep breakdown — how it all works, where we are, where it goes

Written 2026-07-04, after closing the Hwang-Giapis charging benchmark and three frontier research
passes. This is the map of the whole system: the physics, the engine, the honest state, the
landscape, the frontier, and the build plan. Companion docs: RESEARCH_SOTA.md (landscape),
NEW_PHYSICS_FRONTIER.md (physics frontier), METHODS_FRONTIER.md (methods frontier),
NEIGHBOR_MECHANISM.md (the benchmark-closing mechanism), PHYSICS_FRONTIER.md (Graves/chemistry moat).

---

## PART 1 — The physics: how feature charging actually works

**The setup.** Etching a chip feature (trench, hole, gate) means firing plasma at a wafer through a
mask. Two species arrive:
- **Ions** — heavy, fast, *directional* (the sheath accelerates them straight down). Energy set by
  the RF bias: V_s = 37 + 30·sin(ωt), so ~7-67 eV. Angular spread ~4°.
- **Electrons** — light, thermal (~4 eV), *isotropic* (arrive from all angles). They only reach the
  wafer during the brief RF-phase window when the sheath collapses (RF bursts).

**Why anything charges.** On an insulator (SiO2 floor, photoresist), charge can't flow away. The
directional ions reach the trench bottom easily; the isotropic electrons are shadowed by the walls.
More ions than electrons at the bottom → the floor charges **positive** until it repels enough ions
that ion current = electron current (an insulator sits at zero net current at steady state). This is
the master rule: **every insulator cell floats to local current balance; every connected conductor
floats as one equipotential to whole-line current balance.**

**The notch mechanism (the hard part).** In a line-and-space pattern, the **edge line** (bordering
the open area) keeps its outer wall fed by electrons from the open reservoir → stays LOW (~7 V). The
**neighbor line** (walled in by trenches on both sides) is electron-starved → charges HIGH (~39 V).
That 32 V difference makes a sideways field across the trench that **deflects ions into the foot of
the edge line's inner wall** — the deflected-ion bombardment carves the notch. Higher aspect ratio →
bigger split → deeper notch.

**The subtlety that cost us weeks (now understood).** The neighbor rises to 39 V through a *slow*
self-consistent runaway: as the trench charges, ions deflect onto the neighbor, charging it more,
deflecting more. Critically, **the potential keeps evolving long after the currents balance** — HG's
own paper: currents balance at ~1500 charging steps, but the potential needs ~7000 steps to reach
steady state (a 5× shortfall). We were stopping at current balance and reading 23 V; running to
potential steady state, the neighbor climbs straight to 39 V. Not a missing mechanism — premature
convergence. (SEE is the wrong sign — it *lowers* potential; hot neutrals deposit no charge; narrower
source EAD is second-order. All ruled out by primary sources.)

---

## PART 2 — The petch engine: architecture

Five layers, of which we own four and the fifth is the moat:

```
 REACTOR (HPEM/Kushner, or a sheath model)  -> ion+electron energy-angle distributions at the mouth
        |  [our SOURCE — currently analytic, should ingest sheath-derived EADs]
        v
 ATOMISTIC (Graves MD/DeepMD)  -> surface reaction physics (yields, coverage, damage)
        |  [our SURFACE CHEMISTRY — the moat, not yet built]
        v
 SURFACE STATE (per-voxel: coverage, mixed-layer, damage, fluence)      <- Graves bridge, unbuilt
        |
        v
 FEATURE TRANSPORT + CHARGING  <- THIS REPO. ions+electrons in self-consistent field, charge deposit,
        |                          Poisson/Laplace solve, iterate to steady state
        v
 LEVEL-SET PROFILE EVOLUTION   <- petch already has it (EO advect + FSM reinit)
        |
        v  (all GPU/Warp, all autodiff)
 END-TO-END DIFFERENTIABLE  -> invert the whole chain against measured profiles
```

**What's actually in the repo now:**
- `charging_general.py` — the general geometry-agnostic charging engine (material grid → launch
  ions+electrons → trace in field → deposit → Laplace/Poisson → iterate). Conductors auto-detected as
  connected components. This is the MCFPM-family voxel particle-in-field engine.
- `charging_gpu.py` — the Warp GPU particle tracer (24.6M particles/s on CUDA, parity 0.9995 vs CPU).
  The kinetic core.
- `charging2d.py` — the earlier hardcoded edge-array solver + sky view factors (diagnostic-rich,
  validated reference).
- `threed.py` — 3D transport (mc/radiosity/knudsen/dda) + level-set, GPU Warp kernels (61× on the 3D
  ray tracer).
- The level-set profile evolution, the SF6O2 Belen chemistry, marching cubes, etc.

---

## PART 3 — Honest current state (what's matched, what's a fudge)

**Hwang-Giapis charging benchmark (JAP 82,566) — SUBSTANTIALLY CLOSED:**
| metric (AR 1→4) | petch | HG | status |
|---|---|---|---|
| floor ion flux | 0.60→0.16 | 0.59→0.22 | RMSE ~0.05 ✓ |
| foot ion energy | 17→26 | 10→28 | right trend ✓ |
| floor center V | 16→48 | 8→33 | right trend, calibratable |
| edge line V | ~5 | 2→7 | ✓ (open-wall physics) |
| **neighbor line V** | **38.5 @ it1000** | **39** | **✓ once run to potential steady state** |

The split — the last holdout — closes by running to potential steady state (no fudge). The `vf_focus`
factor was a red herring; the real fix was convergence. Remaining polish: the neighbor slightly
over-shoots past 39 (44.9, still climbing at it1500) — the clean plateau needs either the equipotential
foot-charge-redistribution BC or a d(V)/d(iter) convergence criterion instead of a fixed iteration count.

**Other validated milestones (from the reconciliation arc):** wafer-gate PASS vs de Boer (evolving
Knudsen held-out AR40 prediction, RMSE 0.031-0.043); radiosity+GMRES vs measured ViennaPS static
(RMSE 0.043); ~14× speed vs ViennaPS-GPU on a specific 3D benchmark (honest, with caveats); the notch
shape gate (corr 0.92).

**Known honest limits:** absolute depth uncalibrated on sub-micron grids; the Poisson-through-dielectric
mode needs a grounded substrate (currently gas-only Laplace is the clean default); the source is
analytic (should be sheath-derived EADs); no surface-state chemistry yet.

---

## PART 4 — The landscape (who does what)

- **ViennaPS** (TU Wien, open): general level-set + GPU MC ray-trace flux. **No charging at all** —
  structural (a level-set is a scalar field with no volumetric state). The open baseline we leapfrog
  by simply having charging.
- **Kushner MCFPM** (Michigan, closed): the reference charging engine — voxel material grid + Poisson
  (∇·ε∇Φ=ρ, per-cell ε/mobility) + Lorentz-bent trajectories. CPU, non-differentiable. The physics
  benchmark to reproduce and beat on openness/differentiability/speed. Fed by HPEM (reactor) → PCMCM
  (sheath EADs).
- **Kokkoris/Gogolides** (NTUA): level-set + MC flux + FEM-Laplace charging. General but CPU.
- **Commercial** (SEMULATOR3D voxel, Sentaurus/Victory level-set): general engines, MC flux, **no
  documented charging**.
- **Graves** (Princeton/PPPL): atomistic MD/DeepMD surface chemistry + reduced-order models + Neural
  Master Equation. Stops at 0-D; never coupled to a feature-scale charging engine.

**The gap that is our entire opening:** open + GPU + differentiable + charging exists NOWHERE, and
nobody bridges the reactor (Kushner) and the atoms (Graves) at the feature scale.

---

## PART 5 — The physics frontier (new modules to own)

Ranked (full detail + citations in NEW_PHYSICS_FRONTIER.md):
1. **Stochastic differential-charging → twisting kernel.** The accepted #1 frontier defect at extreme
   AR. Deposit charge as a discrete Poisson *process*, deflect ions, backprop through the level-set to
   *invert* for the twist-suppressing process window. We're ~80% there.
2. **Cryo surface-conductivity charging — hits OUR deep-AR floor collapse.** At −60 °C a condensed
   HF/H2O layer raises insulator surface conductivity 3-6 orders of magnitude and *shorts* feature
   charge. A T-gated lateral-conduction term modulating the exact field we already compute. Kushner+
   Samsung just did it closed (2026); every open tool has neither temperature nor charging.
3. **Curvature field-enhancement at corners (FEE).** Closed-form function of local curvature, trivially
   differentiable, never in any solver — sharpens our notching regime.
Build-fourth: phase-resolved pulsed charging (glow/afterglow relaxation ODE). Materials: MRAM redep
shunt, Ru/Mo volatility, GST stoichiometry drift, ferroelectric-HfO2 vacancy charge.

**The Graves chemistry moat** (PHYSICS_FRONTIER.md): every incumbent uses *instantaneous empirical*
surface chemistry. A per-voxel MD-derived **surface-state layer** (coverage + mixed-layer/damage depth
+ fluence history + energy-angle yields) enables ALE saturation + damage memory nobody can express.
First piece: Vella 2024 two-ODE transient Cl site-balance + Kounis-Melas 2025 energy-angle yield table.

---

## PART 6 — The methods frontier (how to build it fast + differentiable)

Top 3 (full detail + citations in METHODS_FRONTIER.md):
1. **Path-Replay Backprop + boundary-term gradients.** Our MC-charging-over-a-moving-level-set *is*
   differentiable SDF rendering (a solved, production problem in graphics — Mitsuba3/Dr.Jit).
   Constant-memory autodiff; recovers the d(charge)/d(surface) gradient naive AD silently drops.
2. **Implicit / deep-equilibrium differentiation of the steady-state charging.** It's a fixed point →
   implicit function theorem → exact gradients in O(1) memory from one linear solve (vs unrolling
   hundreds of relaxation steps).
3. **Neural-preconditioned CG for the Poisson inner loop** (~30× iteration reduction, HYBRID not
   surrogate — outer CG keeps true residual to tolerance).
GPU patterns (WarpX): cell-sort particles (7.5× deposit), direct atomics, fused gather+push, device-
resident, fp16 SoA, Warp tile + cuFFTDx FFT-Poisson. Collision gradients: adjoint DSMC. Chemistry
bridge: Neural Master Equation (Graves lineage). Precedent that GPU+diff beats CPU incumbents: FDTDX
(415× vs Ceviche).

---

## PART 7 — The moat, in one paragraph

We are not out-racing ViennaPS's ray tracer, rebuilding Kushner's reactor, or redoing Graves's MD. We
are building **the one thing all three point toward and none own: the fast, open, differentiable,
GPU-native feature-scale kinetic charging engine that couples the reactor's EADs and the atoms'
chemistry.** MCFPM has the physics but is closed/CPU/non-differentiable; ViennaPS is open but has no
charging; Graves has the chemistry but stops at 0-D. Differentiability through the full charge → field
→ trajectory → profile loop (proven next door, 140× inverse fits, unclaimed for etch) plus GPU speed
is the defensible position — it turns the simulator from a forward model into an *invertible* one you
calibrate and design against.

---

## PART 8 — The build plan (ordered)

**Now (benchmark polish):** add a d(V)/d(iter) potential-steady-state convergence criterion (kills the
over-shoot, closes HG cleanly with zero fudge); remove the vf_focus stand-in once electrons trace
correctly to steady state.

**Next (the GPU kinetic engine — the load-bearing build):**
1. Warp GPU particle tracer — DONE (`charging_gpu.py`, 24.6M particles/s).
2. GPU Poisson solve (Warp tile + cuFFTDx or neural-preconditioned CG).
3. Keep the whole loop device-resident (cell-sort, direct atomics, fused kernels) — the real speedup.
4. Trace BOTH species self-consistently from a proper sheath-derived source (kill the analytic source).
5. Crank grid resolution (GPU affords it) so floor vs sidewall-foot separate physically.

**Then (differentiable):** wrap the loop with Path-Replay Backprop + implicit-fixed-point gradients →
end-to-end autodiff → calibrate cal_F / IADF / sheath against de Boer + CD-SEM by gradient descent.

**Then (frontier physics, pick by leverage):** cryo surface-conductivity module (our floor collapse) →
stochastic twisting kernel (the #1 frontier defect) → FEE corner charging → Graves surface-state
chemistry (the moat).

The through-line: **close HG cleanly (almost done) → GPU-native kinetic engine → differentiable →
frontier physics nobody else has.** Each step is buildable, and the first three are de-risked.
