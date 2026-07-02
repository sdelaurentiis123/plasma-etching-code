#!/usr/bin/env python3
"""Notching mechanism gate: the DEFLECTED-ion population at the sidewall foot vs Hwang-Giapis.

The notch is dug by ions deflected off the charged floor into the sidewall foot (HG JVST B 15,70;
JAP 82,566). Full notch-DEPTH evolution vs Nozawa/Fujiwara needs multi-material etch-stop (poly on
oxide) which petch does not have yet -- so the gate here is the published INTERMEDIATE observable,
HG's own tabulated deflected-ion numbers (JAP 82,566 Fig. 4, right axis):

  GATE A: mean impact energy of ions on the inner sidewall (foot region) RISES 15 -> 27.5 eV over
          AR 1 -> 4 (8 tabulated points; tolerance +-30% -- first wiring, digitized reference).
  GATE B: the deflected-ion flux to the foot is ~AR-INDEPENDENT (HG: "the deflected-ion flux to
          the notch sidewall stays ~constant with AR"); tolerance: max/min <= 2 over AR 1.6-4
          (the shallow AR<1.6 features barely charge, so little deflection -- excluded per HG).

Everything from the same nothing-tuned solver that passed the floor-flux gate (RMSE 0.039).
Writes notching_gate_result.npz for viz/notching.png (scripts/plot_notching.py).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from petch.charging2d import solve_trench_charging

HG_AR = np.array([1.0, 1.2, 1.6, 2.0, 2.6, 3.0, 3.6, 4.0])
HG_FOOT_E = np.array([15.0, 16.5, 17.5, 20.0, 23.0, 25.0, 26.5, 27.5])   # eV, inner-sidewall ions

Em, Fl = [], []
for ar in HG_AR:
    r = solve_trench_charging(ar, n_per_iter=8000, n_iter=120, seed=3)
    Em.append(r["foot_ion_Emean"]); Fl.append(r["foot_ion_flux"])
    print(f"  AR {ar:3.1f}:  foot E_mean = {Em[-1]:5.1f} eV (HG {HG_FOOT_E[list(HG_AR).index(ar)]:.1f})"
          f"   foot flux = {Fl[-1]:.3f}", flush=True)
Em = np.array(Em); Fl = np.array(Fl)

relerr = np.abs(Em - HG_FOOT_E) / HG_FOOT_E
okA = bool((relerr <= 0.30).all()) and Em[-1] > Em[0]
print(f"GATE A (foot ion energy 15->27.5 eV, +-30%): max rel err = {relerr.max()*100:.0f}%  "
      f"rising = {Em[-1] > Em[0]}   [{'PASS' if okA else 'fail'}]")
sel = HG_AR >= 1.6
ratio = Fl[sel].max() / max(Fl[sel].min(), 1e-9)
okB = bool(ratio <= 2.0)
print(f"GATE B (foot flux ~AR-independent, AR>=1.6): max/min = {ratio:.2f}  [{'PASS' if okB else 'fail'}]")
np.savez(os.path.join(os.path.dirname(__file__), "..", "notching_gate_result.npz"),
         ar=HG_AR, foot_E=Em, hg_foot_E=HG_FOOT_E, foot_flux=Fl, okA=okA, okB=okB)
print("DONE")
