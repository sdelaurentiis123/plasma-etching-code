"""Forecast both base gates under the crosslink bond multiplicity.

Krueger thesis sec. 2.2.3 states the physical basis verbatim: "Each material
has a maximum number of crosslink partners associated with it, which is based
on the number of available bonds (3 in the example depicted in Figure 2.2)."

For a fluorocarbon cell CF_x the available (non-fluorine) bonds on the carbon
are 4 - x, so the partner count follows the film's OWN composition ledger --
an F-rich film is chain-terminated and crosslinks little, a C-rich film forms
a network.  Nothing is fitted: 4 is carbon's valence and x = F/C is measured
from the layer inventory.

Measures the two observables the base gates ride on, at the audited
deliveries, so a change can be forecast before any box spend:
  LIP   -> x_xl and film growth (the mouth closure driver)
  FLOOR -> oxide recession (the depth driver)
"""
import sys
sys.path.insert(0, "src")

import numpy as np

from petch.mixed_layer_mechanism import build_krueger_2024_mixed_layer_mechanisms
from petch.surface_kinetics import EnergeticFlux, SurfaceFluxes

# Krueger Table-I source fluxes (m^-2 s^-1).
SRC = {"CF": 4.4e20, "CF2": 9.4e20, "C2F3": 6.8e20, "CF3": 8.4e19,
       "O": 7.7e20, "C3F4": 9.5e20}
ION_SRC, E_ION = 9.6e19, 1500.0

# Audited top-band (lip) delivery and wall tilt -- RESULTS_LIP_CROSSLINK doc.
LIP_DELIV = 0.37189910507849583
LIP_TILT_DEG = 0.4722228422825369
LIP_VIS = 0.7417599287854645

# ml19 measured feature rates (results/curated/.../ml19-depxl-60s/audit.json).
ML19_RATE_AR1 = 21.91      # nm/s at t=3.6 s, AR ~ 1 (open field)
ML19_RATE_AR15 = 16.71     # nm/s at t=42.9 s, AR ~ 15.4
KRUEGER_AVG = 825.0 / 60.0  # 13.75 nm/s run average


def relax(mech, neutral_scale, cos_inc, ion_scale, dt=2.0, steps=600):
    """Relax a surface to steady state at a prescribed delivery."""
    neutral = {k: v * neutral_scale for k, v in SRC.items()}
    ion = EnergeticFlux(name="Ar+", flux_m2_s=ION_SRC * ion_scale,
                        energy_eV=np.array([E_ION]),
                        cosine_incidence=np.array([max(cos_inc, 1e-4)]),
                        weight=np.array([1.0]))
    fluxes = SurfaceFluxes(neutral_flux_m2_s=neutral, energetic_fluxes=(ion,))
    state = mech.initial_state(())
    result = None
    for _ in range(steps):
        result = mech.advance(state, fluxes, dt)
        state = result.state
    return state, result


def report(label, state, result):
    n_c = float(np.asarray(state.n_c_film))
    n_f = float(np.asarray(state.n_f_film))
    total = n_c + n_f
    x_xl = float(np.asarray(state.n_xl_film)) / max(total, 1e-300)
    fc = n_f / max(n_c, 1e-300)
    etch = float(np.asarray(result.etch_velocity_m_s)) * 1e9
    growth = float(np.asarray(result.normal_growth_velocity_m_s)) * 1e9
    print(f"  {label:32s} x_xl={x_xl:6.3f}  F/C={fc:6.3f}  "
          f"film={total:9.3e}  etch={etch:8.3f}  growth={growth:7.3f} nm/s")
    return dict(x_xl=x_xl, fc=fc, film=total, etch=etch, growth=growth)


def main():
    oxide, mask = build_krueger_2024_mixed_layer_mechanisms()

    print("FLOOR (oxide, normal incidence) -- the depth gate driver")
    blanket = report("blanket (AR~0, full delivery)",
                     *relax(oxide, 1.0, 1.0, 1.0))
    # AR ~ 15 delivery: ions survive the collimated-beam wall loss (~0.82),
    # thermal neutrals are Clausing-throttled by the aperture.
    deep = report("floor at AR~15 (ion 0.82)",
                  *relax(oxide, 0.10, 1.0, 0.82))

    print("\nLIP (mask film, near-vertical wall) -- the mouth gate driver")
    cos_lip = float(np.sin(np.deg2rad(LIP_TILT_DEG)))
    lip = report("top band (tilt 0.47deg)",
                 *relax(mask, LIP_DELIV, cos_lip, cos_lip * LIP_VIS))

    print("\nConsistency vs the ml19 feature run")
    print(f"  0-D blanket etch      {blanket['etch']:8.3f} nm/s")
    print(f"  ml19 measured AR~1    {ML19_RATE_AR1:8.3f} nm/s")
    print(f"  ml19 measured AR~15   {ML19_RATE_AR15:8.3f} nm/s")
    print(f"  0-D floor AR~15       {deep['etch']:8.3f} nm/s")
    print(f"  Krueger run average   {KRUEGER_AVG:8.3f} nm/s")
    if deep['etch'] > 0:
        print(f"  depth overshoot vs Krueger avg: "
              f"{ML19_RATE_AR15 / KRUEGER_AVG:.3f}x (measured feature)")
    print(f"\n  lip growth {lip['growth']:.3f} nm/s "
          f"(Krueger per-side closure 0.427 nm/s -> "
          f"{lip['growth'] / 0.427:.2f}x)")


if __name__ == "__main__":
    main()
