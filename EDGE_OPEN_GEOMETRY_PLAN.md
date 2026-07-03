# Explicit HG edge/open geometry plan

Status: planned 2026-07-03 after `edge_open_model="line_of_sight"` was measured and failed the full
HG gates. The auxiliary boundary matched HG Fig. 3 gross outer electron flux, but still failed floor
flux, foot energy, and neighbor-poly potential. Therefore the next implementation target is the
actual nonperiodic geometry/electrostatics, not another scalar current correction.

## Target geometry

Build a 2-D edge-array charging solver with:

- open area on the left;
- edge poly-Si line with outer and inner sidewalls tied to one conductor potential;
- edge trench between edge and neighboring line;
- neighboring poly-Si line as a second conductor;
- optional right-side continuation/buffer so the neighbor is not artificially grounded by the
  numerical boundary;
- photoresist sidewalls above 0.3 um poly-Si, oxide floor below.

HG observables must be side-resolved:

- oxide floor ion flux in the edge trench;
- ion flux and mean incident energy on the edge-line **inner** poly sidewall;
- edge-line and neighboring-line poly equipotentials separately;
- gross outer edge-line electron flux for Fig. 3 comparison;
- current residual per insulator segment and per conductor.

## Code strategy

Do this as a new solver path first, not as a rewrite of `solve_trench_charging`.

Primary file:

- `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/src/petch/charging2d.py`

Add:

- `solve_edge_array_charging(...)` with the same source/SEE/integrator options as
  `solve_trench_charging`;
- a nonperiodic Laplace helper: x boundaries should be Neumann/open-buffer, not `np.roll`;
- an occupancy/material grid for solids instead of the implicit `x < pad or x >= pad+W` wall test;
- a segment map that labels each exposed surface as PR insulator, oxide floor, edge poly outer,
  edge poly inner, neighbor poly, top mask, or escape boundary;
- a tracer variant that terminates on the first solid-cell crossing and returns the segment id.

Keep the existing periodic solver as the regression baseline. Do not replace
`charging_floor_profile(AR)` until the new edge-array solver passes gates.

## Minimal implementation slices

1. Geometry rasterizer
   - Inputs: `W`, `open_width_um`, `poly_um`, `feature_w_um`, `AR`, `n_lines=2`.
   - Output: `solid`, `material_id`, `segment_id`, conductor masks, and a small geometry dict.
   - Unit test: segment counts are nonzero and conductor ids connect outer+inner edge sidewalls.

2. Nonperiodic Laplace
   - Replace x-periodic `np.roll` with neighbor averaging using mirrored boundary cells or a
     sufficiently wide grounded/open buffer.
   - Unit test: with symmetric two-trench geometry and no open side, left/right potentials remain
     symmetric.

3. Segment tracer
   - First version can be scalar/Numba but must preserve the adaptive step logic and survivor gate.
   - Unit test: vertical ions hit floor, shallow-angle particles hit expected sidewalls.

4. Current relaxation
   - Insulators update per segment.
   - Conductors update per conductor id from total ion-electron current over all tied segments.
   - Gate: current residual max <= 0.08 and survivor max < 0.001 in reduced AR4.

5. Official gates and figures
   - Extend `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/scripts/charging_gate.py`
     and `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/scripts/notching_gate.py` with
     `PETCH_CHARGING_GEOMETRY=edge_array`.
   - Add/extend figures:
     - `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/viz/edge_open_current.png`
     - `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/viz/charging_hg.png`
     - `/Users/stanislavdelaurentiis/chip-etch/plasma-etching-code/viz/notching.png`

## Pass/fail gates

- HG floor flux RMSE <= 0.05.
- AR4 bottom-center potential within 33 V +/- 40%.
- Edge-line inner poly foot energy rises 15 -> 27.5 eV with max error <= 30%.
- Edge-line foot flux max/min <= 2 for AR >= 1.6.
- Neighbor poly potential rises 6 -> 39 V with max error <= 30%.
- Survivor fraction < 0.001.
- Current residual max < 0.08.
- Matsui 300 eV floor remains open at AR4.

## Kill criteria

- If the edge-array solver still matches HG Fig. 3 outer gross electron flux but fails neighbor-poly
  potential by >30%, the residual is the x-boundary/electrostatic domain size. Sweep open/buffer
  widths once, document convergence, then stop.
- If it passes potentials but fails foot energy, the remaining error is the sheath-source IEAD/EEAD
  approximation; move to the sheath-MC source gate.
- If it passes charging but not notch depth, move to material etch/overetch coupling rather than
  changing charging constants.
