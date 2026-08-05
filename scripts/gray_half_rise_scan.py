import numpy as np, sys
sys.path.insert(0,"src")
from petch import mixed_layer as ml
# Rung-0 closed form: theta_F = J/(J + 4*capacity) -> half-rise at J = 4*capacity.
# Sticking s scales J_effective, so half-rise F/Ar+ = 4*capacity/(s*ion_flux).
for s in (1.0, 0.3, 0.1, 0.06, 0.03):
    ml._THERMAL_F_STICKING = s
    ion = 1.0e19
    def Y(r):
        fx = ml.SurfaceFluxes(precursor_flux=0.0, fluorine_flux=r*ion, oxygen_flux=0.0,
                              ion_flux=ion, ion_energy_eV=350.0, cosine_incidence=1.0)
        return float(np.asarray(ml.steady_state(fx).substrate_removal_rate))/ion
    lo, hi = Y(0.0), Y(500.0)
    rs = np.logspace(-1, 2.7, 24); ys = np.array([Y(r) for r in rs])
    half = float(rs[np.argmax(ys >= 0.5*(lo+hi))])
    print(f"s={s:5.2f}  dyn_range={lo/hi:.3f}  half_rise={half:6.2f}")
ml._THERMAL_F_STICKING = 1.0
print("Gray 1993 target: dyn_range 0.20-0.30, half-rise 27 +/- 8")
