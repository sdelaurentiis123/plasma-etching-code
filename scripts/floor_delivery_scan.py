"""Is the ml19 floor etch rate reachable by any physical delivery?

ml19 etches 21.9 nm/s at its shallowest recorded point and 16.7 nm/s at its
deepest -- but the etch front sits below an 850 nm mask, so even the first
record is AR ~ 11.7 and the last is AR ~ 32.  A shadowed floor cannot receive
more NEUTRAL flux than the open field, so scanning the delivery plane bounds
what the chemistry can produce and isolates whether the feature's rate implies
an energetic over-delivery (the funnelled hot-neutral cascade).
"""
import sys
sys.path.insert(0, "src")

import numpy as np

from petch.mixed_layer_mechanism import build_krueger_2024_mixed_layer_mechanisms
from petch.surface_kinetics import EnergeticFlux, SurfaceFluxes

SRC = {"CF": 4.4e20, "CF2": 9.4e20, "C2F3": 6.8e20, "CF3": 8.4e19,
       "O": 7.7e20, "C3F4": 9.5e20}
ION_SRC, E_ION = 9.6e19, 1500.0


def etch_rate(oxide, neutral_scale, ion_scale, dt=2.0, steps=400):
    neutral = {k: v * neutral_scale for k, v in SRC.items()}
    ion = EnergeticFlux(name="Ar+", flux_m2_s=ION_SRC * ion_scale,
                        energy_eV=np.array([E_ION]),
                        cosine_incidence=np.array([1.0]),
                        weight=np.array([1.0]))
    fluxes = SurfaceFluxes(neutral_flux_m2_s=neutral, energetic_fluxes=(ion,))
    state = oxide.initial_state(())
    result = None
    for _ in range(steps):
        result = oxide.advance(state, fluxes, dt)
        state = result.state
    return float(np.asarray(result.etch_velocity_m_s)) * 1e9


def main():
    oxide, _ = build_krueger_2024_mixed_layer_mechanisms()
    ion_scales = [1.0, 2.0, 3.0]
    neutral_scales = [0.02, 0.05, 0.1, 0.3, 0.6, 1.0]
    print("oxide recession (nm/s); rows = ion flux / source, cols = neutral "
          "delivery / source")
    header = "  ion\\neu " + "".join(f"{n:>8.2f}" for n in neutral_scales)
    print(header)
    for i in ion_scales:
        row = [etch_rate(oxide, n, i) for n in neutral_scales]
        print(f"  {i:6.2f}  " + "".join(f"{v:8.2f}" for v in row))
    print("\n  ml19 measured: 21.91 nm/s at AR~11.7, 16.71 nm/s at AR~32")
    print("  Krueger run average: 13.75 nm/s")


if __name__ == "__main__":
    main()
