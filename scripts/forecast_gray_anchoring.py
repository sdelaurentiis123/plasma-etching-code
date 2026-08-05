import importlib.util, numpy as np
spec = importlib.util.spec_from_file_location("solve", "scripts/ion_channel_model_solve.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
ml = m.ml
orig_tpy, orig_cef = ml._threshold_power_yield, ml._complex_energy_factor
BARE = 0.28/0.341      # -> Gray floor
CPLX = 1.10/0.390      # -> Gray plateau (dyn range then follows as a CHECK)
print(f"bare x{BARE:.3f}  complex x{CPLX:.2f}   (dynamic range is then a check, not a fit)\n")
def tpy(E,p0,eth,er,q):
    v=orig_tpy(E,p0,eth,er,q)
    return v*BARE if abs(float(p0)-0.0852)<1e-12 else v
def cef(e,eps,ref): return orig_cef(e,eps,ref)*CPLX
base_front,_ = m.coupled_rate(1.0, m.ION_SRC, m.E_FRONT)
ml._threshold_power_yield, ml._complex_energy_factor = tpy, cef
dyn,_ = m.gray_curve()
front,_ = m.coupled_rate(1.0, m.ION_SRC, m.E_FRONT)
lo,_ = m.coupled_rate(0.05, m.ION_SRC, m.E_FRONT)
mid,_ = m.coupled_rate(0.3, m.ION_SRC, m.E_FRONT)
arde,_ = m.arde_ratio()
print(f"  Gray floor {m.beam_yield(0.0):.3f} (0.28)  plateau {m.beam_yield(500.0):.3f} (1.10)")
print(f"  Gray dynamic range {dyn:.3f}  vs measured 0.255   <-- CHECK, follows from the two anchors")
print(f"  coupled front rate {front:.2f} nm/s   depth factor {front/base_front:.3f}  (need 0.735-0.812)")
print(f"  neutral response: 0.05x -> {lo/front:.3f} ; 0.3x -> {mid/front:.3f}  (was 1.017 = flat)")
print(f"  ARDE (AR16/AR0): {arde:.3f}  (was 0.972)")
ml._threshold_power_yield, ml._complex_energy_factor = orig_tpy, orig_cef
