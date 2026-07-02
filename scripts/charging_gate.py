#!/usr/bin/env python3
"""Charging gate: the Hwang-Giapis current-balance closure vs the published curves.

GATE 1 (primary): floor ion flux vs AR — HG JAP 82, 566 Fig. 4 (digitized +-0.02),
Cl2 HDP, V_s = 37 + 30 sin(wt), T_e = 4 eV, EADF ~ cos^0.6. All constants from the paper;
NOTHING tuned. Gate: RMSE <= 0.05 over the 8 points.
GATE 2: floor potential V_f(AR) in ~8 -> 33 V over AR 1 -> 4 (HG text).
GATE 3 (asymptote): with 300 eV ions (Matsui APL 78, 883), the floor must NOT cut off at AR 4-6
(full cutoff only above AR ~ 7).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from petch.charging import floor_balance, ied_survival

HG_AR = np.array([1.0, 1.2, 1.6, 2.0, 2.6, 3.0, 3.6, 4.0])
HG_FLUX = np.array([0.59, 0.55, 0.47, 0.40, 0.34, 0.30, 0.26, 0.22])

print("=== GATE 1: floor ion flux vs AR (HG JAP 82,566 Fig.4; Te=4eV, Vs=37+30sin) ===")
pred, vfs = [], []
for ar in HG_AR:
    vf, flux = floor_balance(ar, Te=4.0, V_dc=37.0, V_rf=30.0, cos_power=0.6, n=400000)
    pred.append(flux); vfs.append(vf)
    print(f"  AR {ar:3.1f}:  model flux={flux:.3f}  HG={HG_FLUX[list(HG_AR).index(ar)]:.2f}   V_f={vf:5.1f} V")
pred = np.array(pred)
rmse = float(np.sqrt(np.mean((pred - HG_FLUX) ** 2)))
print(f"  GATE 1 RMSE = {rmse:.3f}  [{'PASS' if rmse <= 0.05 else 'fail'}] (gate 0.05)")

print("=== GATE 2: floor potential range (HG: ~8 V @ AR1 -> ~33 V @ AR4) ===")
ok2 = 4.0 <= vfs[0] <= 14.0 and 26.0 <= vfs[-1] <= 40.0
print(f"  V_f(AR1) = {vfs[0]:.1f} V (HG ~8)   V_f(AR4) = {vfs[-1]:.1f} V (HG ~33)   [{'PASS' if ok2 else 'fail'}]")

print("=== GATE 3: Matsui asymptote (300 eV ions: no cutoff at AR 4-6; full cutoff ~AR7+) ===")
f4 = floor_balance(4.0, Te=4.0, V_dc=300.0, V_rf=30.0, cos_power=0.6, n=400000)[1]
f6 = floor_balance(6.0, Te=4.0, V_dc=300.0, V_rf=30.0, cos_power=0.6, n=400000)[1]
ok3 = f4 > 0.05 and f6 > 0.02
print(f"  300 eV ions: floor flux AR4 = {f4:.3f}, AR6 = {f6:.3f}   [{'PASS' if ok3 else 'fail'}] (must stay > 0)")
print("DONE")
