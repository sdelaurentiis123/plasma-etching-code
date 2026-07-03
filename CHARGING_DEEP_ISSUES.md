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

## Deep issue

The remaining miss is source/trajectory selection into the poly foot, not missing termination,
not the scalar poly potential, and not PR-sidewall SEE alone. HG's method used sheath-MC-derived
joint arrival distributions. petch still uses analytic source shortcuts, especially for the joint
ion/electron energy-angle-phase structure.

## Falsified or deprioritized

- **PR-sidewall SEE alone:** falsified by corrected cascade AR4 probe.
- **Insulator floating clip as main cause:** falsified by reduced `insul_vmin_Te` sweep; foot energy
  barely moves while `V_c` and flux move strongly.
- **Simple high-energy/high-angle ion coupling:** reduced `ion_angle_energy_corr="positive"` made
  AR4 foot energy worse, not better.
- **Conductor charge sharing:** deprioritized because `V_poly` passes while foot energy fails.

## Next diagnostics, in order

1. **Recollection audit for SEE.** Count emitted electrons by source surface and final absorbing
   surface. If PR-sidewall emitted electrons mostly re-hit upper PR or escape, SEE cannot fix the
   deep gate without source changes.
2. **Material SEE sign test.** Reuse the current SEE machinery as a binary diagnostic on PR-only,
   oxide-floor-only, poly-only, and all-wall surfaces. Do not treat PMMA yields on oxide/poly as a
   calibrated model; this is only a sign/magnitude test.
3. **Sheath-source implementation.** Build the W3 source interface before another tuning pass:
   sampled joint `f_e(E, theta, phase)` and `f_i(E, theta, phase)` from a 1-D RF sheath MC, with the
   current analytic source retained for A/B.
4. **Ion foot-hit phase audit.** For ions that hit the poly foot, record launch phase, initial
   energy, initial angle, impact `z`, and impact energy. The question is whether HG's rising foot
   energy requires a population our analytic source under-samples.
5. **Mesh/geometry control.** Repeat the AR4 foot diagnostics at finer `W`/`D` only after the source
   audit; geometry is lower priority unless source diagnostics are flat.

## Primary references

- Hwang and Giapis, JAP 82, 566 (1997), notching mechanism and charging curves:
  https://authors.library.caltech.edu/records/je8bd-j6v68
- Memos, Lidorikis, and Kokkoris, Micromachines 9, 415 (2018), PMMA SEEE model and adaptive
  trajectory treatment: https://www.mdpi.com/2072-666X/9/8/415
- Donko et al., RF sheath ion distribution context: https://arxiv.org/abs/1809.06779
