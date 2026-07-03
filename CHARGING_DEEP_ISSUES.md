# Charging deep issues after W1/W2

Status: researched and measured 2026-07-03. This memo is the handoff after the adaptive
trajectory integrator and PR-sidewall SEE work. It narrows the remaining Hwang-Giapis miss.

## Measured facts

- W1 fixed the numerical drop: survivor fractions are below 0.1% in the gates.
- PR-sidewall PMMA SEE is real and sign-correct, but not the missing lever.
  - First full SEE gate (`charging_gate_see_result.npz`, `notching_gate_see_result.npz`):
    floor-flux RMSE 0.091 fail; AR4 `V_c = 40.6 V`; AR4 foot energy 16.45 eV vs HG 27.5 eV.
  - Corrected cascade AR4 probe (8000 x 110, seed 10): default foot energy 15.77 eV;
    `see_model="pmma_pr", see_generations=3` gives 16.70 eV, with 2187 emitted electrons and
    zero survivor leakage.
- Poly conductor charge sharing is probably not the dominant residual: `V_poly` passes the HG
  curve while deep foot energy still fails.
- Scalar field knobs are not enough. Relaxing the insulator minimum moves `V_c` and floor flux but
  leaves AR4 foot energy near 22 eV in reduced probes.
- The foot-hit distribution has an energetic tail but the mean is dominated by low-energy hits:
  full AR4 default p50 11.56 eV, p90 29.12 eV; corrected PR SEE p50 12.79 eV, p90 32.61 eV.
- The primary HG paper defines the notching energy on the **inner poly-Si sidewall of the edge
  line**, and explains the driver as the potential difference between the edge line and neighboring
  line. The previous petch mechanism cell tied both poly sidewalls to one equipotential, so that
  lateral line-to-line tilt was absent by construction.

## Deep issue

The remaining miss is source/trajectory selection into the poly foot, not missing termination,
not the scalar poly potential, and not PR-sidewall SEE alone. HG's method used sheath-MC-derived
joint arrival distributions. petch still uses analytic source shortcuts, especially for the joint
ion/electron energy-angle-phase structure.

The sharper issue after reading HG is geometry/current-balance: the notching gate is an edge-line
quantity, while petch's mechanism cell is symmetric/periodic. A split-conductor diagnostic with
left/right poly lines floating independently stayed symmetric at AR4 (left/right potentials differed
by only ~0.24 V and foot energies were equal), proving that independent conductors alone are not
enough; the open-area electron supply to the edge line must be represented. An imposed edge/neighbor
poly bias diagnostic moved the AR4 relevant-side foot energy strongly:

| imposed line-to-line bias | left poly Emean | right poly Emean | avg foot E |
|---:|---:|---:|---:|
| 0 V | 19.1 eV | 19.0 eV | 19.1 eV |
| 5 V | 20.5 eV | 18.8 eV | 19.9 eV |
| 10 V | 22.5 eV | 19.2 eV | 21.9 eV |
| 15 V | 24.3 eV | 22.2 eV | 24.1 eV |

This is the first diagnostic that moves deep-AR foot energy by multiple eV without relying on a
yield knob. It points to an edge-line geometry/source-current implementation as the next real fix.

## Falsified or deprioritized

- **PR-sidewall SEE alone:** falsified by corrected cascade AR4 probe.
- **Insulator floating clip as main cause:** falsified by reduced `insul_vmin_Te` sweep; foot energy
  barely moves while `V_c` and flux move strongly.
- **Simple high-energy/high-angle ion coupling:** reduced `ion_angle_energy_corr="positive"` made
  AR4 foot energy worse, not better.
- **Single conductor charge sharing:** falsified as insufficient. Split left/right conductors in the
  symmetric periodic cell remain nearly equal; HG needs edge-line/open-area asymmetry.

## Next diagnostics, in order

1. **Recollection audit for SEE.** Count emitted electrons by source surface and final absorbing
   surface. If PR-sidewall emitted electrons mostly re-hit upper PR or escape, SEE cannot fix the
   deep gate without source changes.
2. **Edge-line geometry/current balance.** Replace the periodic one-trench/tied-poly mechanism cell
   with the HG edge-line cell: edge poly line plus neighboring poly line, separate equipotentials,
   and an open-side electron supply on the outer edge-line sidewall. This is now higher priority
   than source tuning.
3. **Material SEE sign test.** Reuse the current SEE machinery as a binary diagnostic on PR-only,
   oxide-floor-only, poly-only, and all-wall surfaces. Do not treat PMMA yields on oxide/poly as a
   calibrated model; this is only a sign/magnitude test.
4. **Sheath-source implementation.** Keep the new `source_model="sheath_mc"` interface for A/B, but
   do not expect it to close the gate until the edge-line geometry exists.
5. **Ion foot-hit phase audit.** For ions that hit the poly foot, record launch phase, initial
   energy, initial angle, impact `z`, and impact energy. The question is whether HG's rising foot
   energy requires a population our analytic source under-samples.
6. **Mesh/geometry control.** Repeat the AR4 foot diagnostics at finer `W`/`D` only after the source
   audit; geometry is lower priority unless source diagnostics are flat.

## Primary references

- Hwang and Giapis, JAP 82, 566 (1997), notching mechanism and charging curves:
  https://authors.library.caltech.edu/records/je8bd-j6v68
- Memos, Lidorikis, and Kokkoris, Micromachines 9, 415 (2018), PMMA SEEE model and adaptive
  trajectory treatment: https://www.mdpi.com/2072-666X/9/8/415
- Donko et al., RF sheath ion distribution context: https://arxiv.org/abs/1809.06779
