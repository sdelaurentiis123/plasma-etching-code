"""Forecast: what declared ion-flux normalization reproduces the 825 nm depth?

Frozen-geometry floor rates at HAR delivery, scanning the ion-only and uniform
scalings.  The distinction matters and is the reason ion-only is the mode
implemented: uniform scaling raises deposition in lockstep with removal, so the
closure/etch ratio (measured 2.33x Krueger's) does not improve; ion-only raises
removal against a fixed depositor flux, so it acts on the depth gate and the
closure ratio together.

Yields are never touched here -- only the boundary flux, which is the quantity
with no measurement behind it.  See RESULTS_BLANKET_ANCHOR_2026-08-06.md.
"""
import sys

import numpy as np

sys.path.insert(0, "src")

from petch.chemistry_deck import build_mixed_layer_mechanisms_from_deck
from petch.surface_kinetics import EnergeticFlux, SurfaceFluxes

_ION_PUBLISHED = 1.2e20
_NEUTRALS = {"CF": 4.4e20, "CF2": 9.4e20, "C2F3": 6.8e20, "CF3": 8.4e19,
             "O": 7.7e20, "C3F4": 9.5e20}
# Measured board (RESULTS_SCORECARD_ENDPOINT_2026-08-06.md).
_MEASURED_DEPTH_NM = 346.8
_TARGET_DEPTH_NM = 825.0


def floor_rate_nm_s(ion_scale=1.0, neutral_scale=1.0, *, neutral_delivery=0.10,
                    ion_delivery=0.70, energy_eV=3406.0, steps=120, dt=2.0):
    oxide, _ = build_mixed_layer_mechanisms_from_deck()
    neutrals = {k: v * neutral_delivery * neutral_scale
                for k, v in _NEUTRALS.items()}
    ion = EnergeticFlux(
        name="Ar+", flux_m2_s=_ION_PUBLISHED * ion_delivery * ion_scale,
        energy_eV=np.array([energy_eV]), cosine_incidence=np.array([1.0]),
        weight=np.array([1.0]))
    fluxes = SurfaceFluxes(neutral_flux_m2_s=neutrals, energetic_fluxes=(ion,))
    state = oxide.initial_state(())
    result = None
    for _ in range(steps):
        result = oxide.advance(state, fluxes, dt)
        state = result.state
    return float(np.asarray(result.etch_velocity_m_s)) * 1e9


def main():
    need = _TARGET_DEPTH_NM / _MEASURED_DEPTH_NM
    base = floor_rate_nm_s()
    print(f"measured 60 s depth {_MEASURED_DEPTH_NM} nm -> need {need:.3f}x "
          f"to reach {_TARGET_DEPTH_NM} nm")
    print(f"floor rate at published flux: {base:.3f} nm/s\n")
    print(f"{'scale':>6} {'ion-only':>10} {'ratio':>7} {'uniform':>10} {'ratio':>7}")
    rows = []
    for s in (1.0, 1.5, 2.0, 2.4, 2.8, 3.2, 4.0):
        io = floor_rate_nm_s(ion_scale=s)
        un = floor_rate_nm_s(ion_scale=s, neutral_scale=s)
        rows.append((s, io / base))
        print(f"{s:6.2f} {io:10.3f} {io/base:7.3f} {un:10.3f} {un/base:7.3f}")
    xs = [r[0] for r in rows]
    ys = [r[1] for r in rows]
    print(f"\nion-only scale for {need:.3f}x floor rate: "
          f"{np.interp(need, ys, xs):.2f}x")
    print("NOTE: frozen-geometry estimate. The coupled run has mouth feedback "
          "(a wider mouth raises delivery), so the realized depth response is "
          "expected to exceed this; the run measures it.")


if __name__ == "__main__":
    main()
