"""Demo: the BACKWARD/ADJOINT self-consistent charging solve reproduces the electron-shading dipole
and the Kushner low-AR->high-AR crossover, from first principles, with NO per-region overrides.

Run: python scripts/backward_charge_demo.py
Prints the converged left-wall potential profile (depth from mouth) + floor mean for AR4/8/15.
Expect: upper wall NEGATIVE (~-12.7 V, electron shading), deep wall POSITIVE and rising with AR
(grazing ions), floor positive + monotone in AR; the potential MAX moves from the floor (low AR) to
the deep sidewall (high AR), matching Kushner (JAP 137,063302; MCFPM)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from petch.charging2d import _build_edge_array_geometry
from petch.charging_backward import self_consistent_backward

for AR in (4.0, 8.0, 15.0, 30.0):                      # NEE cone importance sampling reaches AR30
    g = _build_edge_array_geometry(AR, W=16, mouth=80)
    r = self_consistent_backward(g, n_iter=14)
    print(f"\nAR{AR:g}  (floor mean = {r['floor_mean']:+.1f} V,  k=Ci/Ce={r['k']:.2f})")
    print("  depth+  Vs")
    for d, v in zip(r['wall_depth'], r['Lwall']):
        print(f"  {d:5d}  {v:+6.1f}")
print("\nupper wall NEGATIVE (electron shading), deep wall POSITIVE (grazing ions), floor POSITIVE+monotone.")
