import sys; sys.path.insert(0, "src")
import numpy as np
from petch.mixed_layer_mechanism import build_krueger_2024_mixed_layer_mechanisms
from petch.surface_kinetics import EnergeticFlux, SurfaceFluxes

# Krueger Table-I source fluxes (m^-2 s^-1) and the audited top-band delivery.
SRC = {"CF": 4.4e20, "CF2": 9.4e20, "C2F3": 6.8e20, "CF3": 8.4e19,
       "O": 7.7e20, "C3F4": 9.5e20}
ION_SRC, E_ION = 9.6e19, 1500.0
DELIV = 0.37189910507849583          # audited top-band depositor delivery
TILT_DEG = 0.4722228422825369        # audited top-band wall tilt

def run(delivery, cos_inc, ion_scale, label, dt=2.0, steps=4000):
    oxide, mask = build_krueger_2024_mixed_layer_mechanisms()
    neutral = {k: v * delivery for k, v in SRC.items()}
    ion = EnergeticFlux(name="Ar+", flux_m2_s=ION_SRC * ion_scale,
                        energy_eV=np.array([E_ION]),
                        cosine_incidence=np.array([max(cos_inc, 1e-4)]),
                        weight=np.array([1.0]))
    fx = SurfaceFluxes(neutral_flux_m2_s=neutral, energetic_fluxes=(ion,))
    st = mask.initial_state(())
    for _ in range(steps):
        r = mask.advance(st, fx, dt); st = r.state
    tot = float(np.asarray(st.n_c_film)) + float(np.asarray(st.n_f_film))
    x = float(np.asarray(st.n_xl_film)) / max(tot, 1e-300)
    growth = float(np.asarray(r.normal_growth_velocity_m_s)) * 1e9
    print(f"{label:34s} x_xl={x:6.3f}  film={tot:8.3e}  growth={growth:7.3f} nm/s")
    return x, growth

# lip: near-vertical wall, ion areal flux collapses as sin(tilt) x visibility
cos_lip = np.sin(np.deg2rad(TILT_DEG))
run(DELIV, cos_lip, cos_lip * 0.7417599287854645, "LIP (top band, tilt 0.47deg)")
run(1.0, 1.0, 1.0, "BLANKET (normal incidence)")
