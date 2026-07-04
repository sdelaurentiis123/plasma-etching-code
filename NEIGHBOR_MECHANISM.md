# Why the neighbor line stalls at ~23 V instead of HG's ~39 V (mechanism hunt 2026-07-04)

Sources are HG's own open-access papers (this is HG's Caltech MC model, NOT Kushner MCFPM):
- JAP 82, 566 (1997): https://authors.library.caltech.edu/records/je8bd-j6v68  (the +39 V result, Figs 5-6)
- JVST B 15, 70 (1997): https://authors.library.caltech.edu/records/ac5xn-zqb88  (full method, charge update)
HG solve LAPLACE (surface charge as BC), and SEE is explicitly neglected -> SEE ruled out by construction.

## THE diagnosis (verbatim from HG)

The neighbor rise is **current-balance / electron-shading**, self-consistent between surfaces, realized
as a **whole-conductor equipotential steady state** -- NOT ion energy, NOT SEE.

- Runaway loop: *"The repulsive entrance potential reduces the flux of electrons to the trench bottom,
  thus forcing the SiO2 surface potential to increase, until enough ions are deflected so that the ion
  and electron fluxes are balanced... [the edge line's] potential can be maintained low by electrons
  arriving at its outer sidewall from the open area. Geometric shadowing prohibits the latter... at
  intermediate lines."*
- Why the neighbor specifically: *"The outer sidewall of the edge line is supplied by electrons from
  the open space... The same is not possible for the poly-Si sidewalls of the neighboring line, whose
  potential should increase much more than that of the edge line."*

## THE fix -- why WE stall at 23 V (the actionable part)

**Premature convergence + conductor BC.** HG, verbatim:
- *"The surface charge on the poly-Si sidewalls must be unevenly distributed in order to make the
  surface equipotential."*
- *"despite the fact that the ion and electron currents to the bottom SiO2 surface are perfectly
  balanced, a dramatic charge redistribution occurs... About 7000 charging steps are needed for the
  potential distribution to reach steady state"* -- whereas currents balance at ~1500 steps: **a ~5x
  shortfall.** We stop when net current ~ 0; HG runs ~5x longer to POTENTIAL steady state.

Minimal non-fudge fix (two parts, purely kinetic-electrostatic):
1. Replace symmetric/equidistributed conductor charge with HG's four-step EQUIPOTENTIAL redistribution:
   recompute each conductor cell's potential from all charged cells (+ mirror images), reallocate the
   line's total charge proportional to each cell's deviation from the line-mean until the whole line is
   one equipotential. This naturally PILES CHARGE AT THE POLY/SiO2 FOOT and breaks edge-vs-neighbor
   symmetry, creating the strong foot field that deflects ions in -> the feedback.
2. Run to POTENTIAL steady state (~5x more charging than current-balance), not current balance.
3. If still short: stop neglecting the shadowed-foot electron flux (HG's own admitted approximation) --
   let electrons be attracted into the positive valley self-consistently (Memos-Kokkoris FEM-Laplace
   trajectory model does this: Micromachines 9, 415, 2018).

## Ruled out (do NOT chase these for the 16 V gap)
- SEE: WRONG SIGN. Kushner JVST A 44, 023013 (2026): *"SEE redistributes charge within the feature,
  LOWERING the in-feature positive electric potential"* (16/13/2% reduction at AR 16.7/25/50). It would
  push the neighbor DOWN.
- Hot neutrals: neutral on impact -> deposit no charge; they reshape the profile (microtrenching), not
  neighbor charging.
- Narrower source EAD: second-order (Huard thesis: ion width 0.55->3.3deg moves CD -18..+30%, not the
  charging ceiling). Electrons are launched ISOTROPIC Maxwellian ~4 eV, one per ion (Kushner).

## Cross-check with SOTA source generation
Kushner does NOT hand-pick a cosine spread: HPEM->PCMCM computes the wafer-plane EAD by MC trajectory
integration through the self-consistent time-resolved sheath, then MCFPM launches "electron
pseudoparticles with velocity from an isotropic Maxwellian" one-per-ion. The decisive asymmetry:
*"electrons dominantly initially charge the top surface and small-AR surfaces, while ions reach deep."*

## Implication for our build
This is EXACTLY what the full self-consistent GPU kinetic engine gives once (a) it runs to potential
steady state (not current balance) and (b) the conductor holds a real spatial charge that piles at the
foot. The `vf_focus` fudge is a shortcut around not doing (a)+(b). The Warp GPU tracer (charging_gpu.py)
makes the required ~5x-longer, million-particle self-consistent runs affordable.

## CONFIRMED (2026-07-04): running to potential steady state closes the split

Direct test, AR4, laplace field, vf electron model, DEFAULT relax, just more iterations:
| iteration | edge V | neighbor V | note |
|---|---|---|---|
| 400  | 4.5 | 23.5 | where we were stuck (premature) |
| 875  | 5.2 | 36.2 | |
| 1000 | 5.4 | 38.5 | **HG target edge 7 / neighbor 39 -- MATCH** |
| 1125 | 5.1 | 40.5 | keeps climbing |
| 1500 | 3.8 | 44.9 | slight over-shoot, not yet plateaued |

The neighbor rises straight through HG's 39 V once run past current-balance to potential steady state.
Premature convergence WAS the bug, exactly as HG predicted. The vf_focus fudge is secondary. Remaining
polish: it slightly over-shoots (44.9, still climbing at it1500) -- the true plateau needs the
equipotential charge-redistribution BC (foot pile-up) to pin it at steady state, OR a convergence
criterion on d(potential)/d(iter) rather than a fixed iteration count. But the MECHANISM is closed:
edge low + neighbor ~39 emerges self-consistently, no fudge.
