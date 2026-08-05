import importlib.util, numpy as np
spec = importlib.util.spec_from_file_location("solve", "scripts/ion_channel_model_solve.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
ml = m.ml
o_tpy, o_cef = ml._threshold_power_yield, ml._complex_energy_factor
# Gray thesis (MIT, 1993): Table 5-1  Y_sputter = 0.0139 (sqrt(E) - sqrt(18))
#                          Eq. 5-35   beta_e    = 0.053  (sqrt(E) - sqrt(4))
def tpy(E,p0,eth,er,q):
    if abs(float(p0)-0.0852)<1e-12:
        E=np.asarray(E,dtype=float)
        return 0.0139*np.maximum(np.sqrt(np.maximum(E,0.0))-np.sqrt(18.0),0.0)
    return o_tpy(E,p0,eth,er,q)
def cef(e,eps,ref):
    E=np.asarray(e,dtype=float)
    # kernel_complex multiplies by 0.1471, so return beta_e/0.1471
    return 0.053*np.maximum(np.sqrt(np.maximum(E,0.0))-np.sqrt(4.0),0.0)/0.1471
base,_ = m.coupled_rate(1.0, m.ION_SRC, m.E_FRONT)
ml._threshold_power_yield, ml._complex_energy_factor = tpy, cef
print("Gray's own laws (Table 5-1 + Eq 5-35), absolute, no free parameters\n")
for E in (350.0, 2000.0, 3406.0):
    b = 0.0139*(np.sqrt(E)-np.sqrt(18.0)); c = 0.053*(np.sqrt(E)-np.sqrt(4.0))
    print(f"  {E:6.0f} eV : bare {b:6.3f}   complex {c:6.3f}   ratio {b/c:5.3f}")
dyn,_ = m.gray_curve()
front,_ = m.coupled_rate(1.0, m.ION_SRC, m.E_FRONT)
lo,_ = m.coupled_rate(0.05, m.ION_SRC, m.E_FRONT)
mid,_ = m.coupled_rate(0.3, m.ION_SRC, m.E_FRONT)
arde,_ = m.arde_ratio()
print(f"\n  Gray dyn range {dyn:.3f} (measured 0.255)")
print(f"  coupled front {front:.2f} nm/s   depth factor {front/base:.3f}   (need 0.735-0.812)")
print(f"  neutral response 0.05x -> {lo/front:.3f} ; 0.3x -> {mid/front:.3f}   (was 1.017 flat)")
print(f"  ARDE (AR16/AR0) {arde:.3f}   (was 0.972)")
ml._threshold_power_yield, ml._complex_energy_factor = o_tpy, o_cef
