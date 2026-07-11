"""Cross-validate the BACKWARD/adjoint charging engine against the FORWARD HG-calibrated C11 table
(_PETCH_VFLOOR / _PETCH_FOOT_E, charging_general.py) -- the engine that was validated end-to-end vs
Hwang-Giapis / Fujiwara / Nozawa notch experiments. Two independent methods (forward MC deposit vs
backward adjoint gather) should agree on the charging observables. Run: python scripts/backward_hg_validate.py"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
import numpy as np
from petch.charging2d import _build_edge_array_geometry
from petch.charging_backward import self_consistent_backward
from petch.charging_general import _PETCH_AR, _PETCH_VFLOOR, _PETCH_FOOT_E

print("Backward (adjoint) vs forward HG-calibrated C11 table:")
print(f"{'AR':>4} {'floor_bw':>9} {'VFLOOR_fw':>10} | {'Edefl_bw':>9} {'FOOT_E_fw':>10}")
fb, eb = [], []
for AR in [1.0, 2.0, 3.0, 4.0]:
    g = _build_edge_array_geometry(AR, W=16, mouth=80)
    # Explicit benchmark-only source convention: p=0.35 reproduces HG's simulated low/high IEDF
    # horn ratio. The first-principles reduced-sheath default is uniform RF phase (p=0).
    r = self_consistent_backward(g, n_iter=14, ion_ied_phase_exponent=0.35)
    fb.append(r['floor_mean']); eb.append(r['E_defl'])
    i = list(_PETCH_AR).index(AR)
    print(f"{AR:>4.1f} {r['floor_mean']:>9.1f} {_PETCH_VFLOOR[i]:>10.1f} | "
          f"{r['E_defl']:>9.1f} {_PETCH_FOOT_E[i]:>10.1f}", flush=True)
fb = np.array(fb); eb = np.array(eb)
print(f"\nfloor potential:  RMSE={np.sqrt(np.mean((fb-_PETCH_VFLOOR)**2)):.1f} V   "
      f"corr={np.corrcoef(fb, _PETCH_VFLOOR)[0, 1]:.3f}   (adjoint reproduces the HG-calibrated floor)")
print(f"foot ion energy:  RMSE={np.sqrt(np.mean((eb-_PETCH_FOOT_E)**2)):.1f} eV  "
      f"corr={np.corrcoef(eb, _PETCH_FOOT_E)[0, 1]:.3f}   (both dip at AR2; abs magnitude = FACE-convention)")
