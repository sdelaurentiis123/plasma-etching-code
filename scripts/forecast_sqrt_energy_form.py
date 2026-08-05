import importlib.util, numpy as np, sys
spec = importlib.util.spec_from_file_location("solve", "scripts/ion_channel_model_solve.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
ml = m.ml
orig_tpy = ml._threshold_power_yield
orig_cef = ml._complex_energy_factor

def sqrt_tpy(E, p0, eth, er, q):
    # n=0.5 (Steinbruchel) ONLY on the bare-oxide row (p0=0.0852); all other
    # rows untouched so polymer/mask results are unaffected.
    if abs(float(p0) - 0.0852) < 1e-12:
        E = np.asarray(E, dtype=float)
        return p0 * np.maximum(np.sqrt(np.maximum(E, 0.0)) - np.sqrt(eth), 0.0) \
               / (np.sqrt(er) - np.sqrt(eth))
    return orig_tpy(E, p0, eth, er, q)

def sqrt_cef(e_iface, eps_dep, ref):
    E = np.asarray(e_iface, dtype=float)
    return np.maximum(np.sqrt(np.maximum(E, 0.0)) - np.sqrt(35.0), 0.0) / (np.sqrt(140.0) - np.sqrt(35.0))

for label, patch in (("current (linear/ZBL)", False), ("Steinbruchel sqrt(E) both rows", True)):
    if patch:
        ml._threshold_power_yield = sqrt_tpy; ml._complex_energy_factor = sqrt_cef
    dyn, _ = m.gray_curve()
    plateau = m.beam_yield(500.0); floor = m.beam_yield(0.0)
    blanket, _ = m.coupled_rate(1.0, m.ION_SRC, m.E_BLANKET)
    front, film = m.coupled_rate(1.0, m.ION_SRC, m.E_FRONT)
    lo, _ = m.coupled_rate(0.05, m.ION_SRC, m.E_FRONT)   # radical-starved
    arde, _ = m.arde_ratio()
    print(f"\n=== {label}")
    print(f"  Gray 350eV: floor {floor:.3f} (meas 0.28)  plateau {plateau:.3f} (meas 1.10)  dyn {dyn:.3f} (meas 0.255)")
    print(f"  coupled: blanket {blanket:.2f} nm/s   front {front:.2f} nm/s   film {film:.4f} nm")
    print(f"  neutral response: rate at 0.05x radicals {lo:.2f} vs 1.0x {front:.2f}  ->  ratio {lo/front:.3f}")
    print(f"  ARDE ratio (AR16/AR0): {arde:.3f}")
ml._threshold_power_yield = orig_tpy; ml._complex_energy_factor = orig_cef
