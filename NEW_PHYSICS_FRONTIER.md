# New-physics frontier for a GPU differentiable feature-charging engine (research 2026-07-04)

Frontier hunt across HAR/3D-NAND charging, cryo etch, pulsed plasma, stochastic effects, new
materials, and novel mechanisms (2022-2026). Full citations below. The through-line:

**The accepted #1 frontier defect is stochastic differential charging -> twisting/tilting/bowing at
extreme AR. Exactly one group models it well (Kushner MCFPM, Michigan+Samsung): closed-source, CPU,
non-differentiable. The open baseline (ViennaPS) has NO charging at all. The gap between "the
accepted frontier mechanism" and "any open / GPU / differentiable implementation of it" is the entire
opening.** Field roadmap that legitimizes this: Oehrlein, Kushner, Donnelly, Economou, "Future of
plasma etching for microelectronics," JVST B 42, 041501 (2024),
https://cpseg.eecs.umich.edu/pub/articles/JVSTB_42_041501_2024.pdf

## Hard calibration targets the frontier now has (use these as gates)
- Vertical charging gradient 200-300 V top-to-bottom in a HAR hole (Lam twisting patent WO2022271526A1).
- Only ~2.5% of neutral flux reaches the bottom of a 50:1 hole; selectivity drops ~50% over AR 40->140;
  IADF must be <0.5 deg 1-sigma to keep the etch front flat (Shen/Lill, Lam, JJAP 62, SI0801, 2023).
- Cryo: HF/H2O condensed layer raises SiO2 SURFACE CONDUCTIVITY 3-6 orders of magnitude at -60C,
  dissipating feature charge and un-bending ions (APL 123, 212106, 2023).
- Ultra-low-Te cuts charging 5 V -> ~0 V, eliminates notching (Chung/Hanyang, PSST 34, 045009, 2025).
- Pulsed: 50 Hz DC-pulse bias cut TiN loss 0.9 -> 0.4 nm (Hitachi, JJAP 2026); phase-controlled CCP
  +32% bottom CD (Samsung, PSST 2025).

## Top 3 new-physics modules to build (ranked)

### #1 Stochastic differential-charging -> twisting kernel (discrete-charge, differentiable, GPU)
The accepted #1 frontier defect. Deposit charge as a discrete/Poisson PROCESS (not a smooth density),
solve sidewall potential, deflect ions, and BACKPROP through the level-set to invert for the process
window that suppresses twisting. Nobody has this: Kushner's is closed/CPU/coarse-discreteness;
Kokkoris has charging and stochastic roughness as SEPARATE uncoupled modules; ViennaPS has no
charging. Our charging + notching code is ~80% of the way there. Anchors: Kushner/Wang twisting model
IEEE 11017954 (2025); pattern-dependent distortion OSTI 1802573.

### #2 Cryo charging module: temperature-dependent condensed-layer SURFACE CONDUCTIVITY  <-- hits OUR pain point
This is the highest-upside module for OUR specific unsolved problem. The known petch failure is the
deep-AR (AR>20) FLOOR COLLAPSE / over-charging. The cryo mechanism is a direct physical fix: below
~-60C an HF/H2O layer condenses and raises the insulator SURFACE CONDUCTIVITY by 3-6 orders of
magnitude, which SHORTS the accumulated feature charge and un-bends the ion trajectories. As a module:
a temperature-gated surface-conduction term (charge relaxes laterally along surfaces with a T-dependent
mobility) that MODULATES THE EXACT ion-bending field our charging solver already computes -- a small,
differentiable term. Kushner+Samsung JUST coupled cryo+charging+HAR (JVST A 44, 033006, May 2026) --
closed and non-differentiable; the frontier is moving fast. Every open tool has neither temperature nor
charging. Cryo is also the real industry direction (Lam Cryo 3.0, TEL 400+ layer). Anchors:
APL 123, 212106 (2023, surface conductivity); Small Methods 2024 (pseudo-wet H2O autocatalysis);
JVST A 44, 033006 (2026, Kushner cryo+charging).

### #3 Curvature field-enhancement (FEE) + curvature-dependent sputter yield (corner-charging module)
Sheath/charge fields are amplified at sharp geometry (lightning-rod effect); ions steer into corner
hotspots with size set by a field-enhancement factor that is a CLOSED-FORM function of local curvature
-> trivially differentiable, cheap, directly sharpens notch/corner fidelity (our notching regime).
Demonstrated experimentally 2025, NEVER embedded in any feature-scale solver -- unoccupied. Lowest
risk, high novelty. Anchors: Chang/DTU "Field enhancement effect in reactive ion etching," Materials &
Design 254, 114144 (2025), open PDF
https://backend.orbit.dtu.dk/ws/files/402109473/1-s2.0-S0264127525005647-main.pdf ;
Bradley-Hobler curvature-dependent yield, J. Appl. Phys. 133, 065303 (2023).

### Build-fourth: phase-resolved pulsed charging (glow vs afterglow)
Pure whitespace: charge deposits in the glow phase, drains in the afterglow -- a pulse-phase-indexed
flux BC + per-cycle surface-charge relaxation ODE. Differentiability lets you invert for the pulse
frequency/duty/phase that minimizes bottom charging (5V->0V, +32% bottom CD targets to gate against).

## New materials (each a per-material surface-state add-on, all unmodeled at feature scale)
- MRAM/MTJ: conductive redeposition electrically SHORTS the tunnel junction (most redep-dominated device
  in the industry); no simulator couples angle-dependent multi-material sputter + conductive-shunt.
- Ru/Mo interconnect: volatile RuO4 vs non-volatile RuO2 redep branch (Mo word-lines, 300+ layer NAND).
- GST/chalcogenides: per-element preferential-sputter stoichiometry drift (no code tracks a per-element
  surface stoichiometry vector).
- Ferroelectric HfO2: etch damage governed by charged oxygen-vacancy state (arXiv:2412.06416).

## Two cross-cutting moats nobody has combined
1. Differentiability through the charge -> field -> trajectory -> profile loop. Method PROVEN in
   adjacent plasma physics (differentiable programming through kinetic plasma sims, 140x speedups,
   >1000-param inverse fits, arXiv:2603.11231, 2026) but UNCLAIMED for feature-scale etch charging.
2. GPU speed at feature scale (our Warp tracer, 24.6M particles/s; the 3D kernel hit 61x).

Benchmark to reproduce+beat: Kushner MCFPM (physics) on openness+differentiability+speed. Leapfrog
ViennaPS by simply HAVING charging.

## Skip / scope-first (don't waste time)
- SKIP: charge-driven Joule/thermal runaway (feature currents ~pA, negligible); "plasma catalysis" for
  etch (no feature-scale support -- the real "5x etch" is cryo autocatalytic H2O, not catalysis).
- SCOPE FIRST (back-of-envelope the field): Fowler-Nordheim tunneling corner-discharge -- only if the
  corner field reaches ~V/nm.
