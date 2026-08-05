import sys, numpy as np
sys.path.insert(0,"src")
from petch import mixed_layer as ml
ION=1.0e19
def Y(r, E=350.0):
    fx = ml.SurfaceFluxes(precursor_flux=0.0, fluorine_flux=r*ION, oxygen_flux=0.0,
                          ion_flux=ION, ion_energy_eV=E, cosine_incidence=1.0)
    return float(np.asarray(ml.steady_state(fx).substrate_removal_rate))/ION
print("Gray 1993 (350 eV Ar+ on SiO2, via Kwon Fig 3.4): floor 0.28, plateau 1.10")
print(f"petch:  floor {Y(0.0):.3f}   plateau {Y(500.0):.3f}")
print(f"ratio petch/Gray:  floor {Y(0.0)/0.28:.2f}x   plateau {Y(500.0)/1.10:.2f}x")
