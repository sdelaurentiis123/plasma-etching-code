# Floor over-charge: mechanistic diagnosis + ruled-out fixes (2026-07-06)

The last open HG-charging residual is the **floor over-charge**: in the geometry-agnostic
engine (`charging_general.solve_charging`) the AR-4 trench floor charges to ~46 V and repels
ions down to a floor flux ~0.16, versus Hwang–Giapis floor V 33 / floor flux 0.22. This note
records what the over-charge actually is, the three fixes that were tested and **refuted**, and
the one path that remains.

## What it is: the floor is electron-GEOMETRIC-PINNED (sub-geometric, even)

Direct instrumentation at AR 4 (`scripts`-style diag, `electron_model="trace"`, nit 800):

| floor electron flux | value | note |
|---|---|---|
| USED (fed to charge balance) | 0.161 | = max(traced, geometric) |
| TRACED (raw kinetic)         | 0.115 | what the trajectory integrator actually delivers |
| GEOMETRIC (sky-view × boost) | 0.124 | analytic solid-angle shadowing |

**traced (0.115) < geometric (0.124).** The reduced 2-D trace delivers *fewer* electrons to
the floor than pure geometric shadowing would — a **negative** anti-shadowing (an entrance-barrier
deficit), the exact opposite of HG's *positive* electrostatic anti-shadowing ("electrostatics
decreases the geometric shadowing", JAP 82,566). Because collection is floored by
`np.maximum(traced, geometric)`, the floor sits at ~geometric shadowing, which is not enough
electron flux to hold the floor at 33 V, so it over-charges to 46 V. The floor-flux deficit and
the floor-V excess are the same fact seen twice.

## Three fixes tested and REFUTED

1. **Electron energy** (hypothesis: our electrons too hot). Sweep `electron_Te` 4→1 eV at fixed
   angle: floor flux barely moves and moves the *wrong* way (0.154→0.171). Refuted — and HG Fig 10
   agrees (their 5 eV electrons reach the floor; 3 eV ones are repelled at the entrance).

2. **Electron angular breadth** (the research agents' leading fix: launch HG's broader `cos^0.6`
   EADF). Added `e_flux_power` (flux ∝ cos^p θ) and swept p 2.0→0.05 (isotropic→grazing): floor
   flux is flat-to-declining (0.182→0.152). **Refuted** — the launch distribution never even sets
   the floor number, because the geometric `max()` overwrites the traced value at the floor. The
   knob is still correct physics (HG's measured EADF) so it is kept, opt-in, default off.

3. **Scalar V-dependent focusing closure** (`electron_model="vf"`, collection ∝ vf·(1+k·V)). Sweep
   `vf_focus_pot` 0→0.10:

   | vf_focus_pot | floorFlux | floorV | edge | neigh |
   |---|---|---|---|---|
   | 0.00 | 0.310 | 67.0 | 32.4 | 45.1 |
   | 0.03 | 0.297 | 41.5 |  9.0 | 27.0 |
   | 0.06 | 0.342 | 29.0 |  5.6 | 20.4 |
   | 0.10 | 0.394 | 21.0 |  2.3 | 14.7 |
   | HG   | 0.22  | 33   |  7   | 39   |

   As focusing rises, floorV sweeps down through 33 and edge through 7 (good) — but floorFlux goes
   the *wrong way* (always ≥0.30) and the **neighbor collapses** (45→15, away from 39). **No single
   knob hits floor + edge + neighbor together.**

## Why the closures fail (the real physics)

HG's anti-shadowing is **spatially selective**: the floor (a broad positive insulator at the
bottom) collects *more* than geometric, while the neighbor line's poly sidewalls stay *shadowed*
and rise to 39 V. A scalar "focus ∝ local V" term cannot distinguish those — it focuses electrons
onto every positive surface, so it feeds the neighbor and destroys the edge/neighbor split that
drives the notch. Only genuine trajectory integration in the self-consistent field separates
"focus into the floor" from "focus onto the neighbor wall."

## The one remaining path

The floor over-charge is a **model-form limit of the reduced 2-D trace**, not a mis-set knob. At
W=16 with a 64-cell-deep trench, the +46 V floor field is too weak a lever over the drop for the
coarse trace to bend oblique electrons in — so it under-focuses to sub-geometric. Closing it needs
the **fine-grid, device-resident GPU kinetic engine** (`charging_gpu.py` is the tracer core):
resolve the in-trench field finely, launch the phase-resolved Boltzmann-gated source, and let real
trajectories produce the spatially-selective focusing. That is a CUDA-box build (cell-sort +
CUDA-graph per the architecture research), not a CPU knob-sweep.

## Status of the committed benchmark page

`docs/charging.html`'s floor-flux curve (0.62→0.22, RMSE 0.016) comes from the **older**
`charging2d.solve_edge_array_charging` at a tuned closure config — a different code path from the
general engine diagnosed here. The page is not over-stated by this finding (it already lists the
high-AR floor over-charge as the open frontier), but the two solvers should be reconciled before
the page is edited: the general engine is the geometry-agnostic future and currently under-performs
the old solver on floor flux for exactly the geometric-pinning reason above.

## Code left in place (all opt-in, default behavior unchanged; 4/4 charging tests pass)
- `sample_electrons(..., flux_power=p)` — cos^p θ flux launch (HG cos^0.6 EADF). Default None.
- `solve_charging(..., e_flux_power=None, electron_Te=None)` — thread the above + electron energy.
- return dict gains `electron_traced` / `electron_geom` (trace model) for floor-flux instrumentation.
