import sys, json; sys.path.insert(0,'src'); sys.path.insert(0,'scripts')
import numpy as np
from petch.mixed_layer import steady_state, SurfaceFluxes
recs = json.load(open('results/curated/e8_thermalized_return/floor_composition.json'))
KR_NEUTRALS = {"C3F4":9.5e20,"C2F3":6.8e20,"CF":4.4e20,"CF2":9.4e20,
               "CF3":8.4e19,"O":7.7e20}   # Krueger Table 6.1 (cm-2s-1 -> m-2s-1)
print("Krueger target floor rate for 825nm/60s = 13.75 nm/s; ml23 gives ~7 nm/s\n")
print("  f_FC   floor C-flux(m-2s-1)   rate(nm/s)")
for rec in recs:
    if rec["ar"] < 5: continue
    e8 = rec["e8_source"]           # per-area thermalized source at floor
    ion = rec["direct_ion"] + rec["hot_neutral"]
    for f in (0.0, 0.15, 0.30, 0.50, 1.0):
        # thermalized FC ions return as precursor carbon+fluorine
        prec = rec["plasma_neutral"] + f*e8
        fx = SurfaceFluxes(precursor_flux=prec, fluorine_flux=0.0,
                           oxygen_flux=0.0, ion_flux=ion, ion_energy_eV=3406.0)
        r = steady_state(fx)
        v = float(np.asarray(r.recession_velocity_m_s))*1e9
        print(f"  {f:4.2f}   {prec:.4e}          {v:8.3f}")
