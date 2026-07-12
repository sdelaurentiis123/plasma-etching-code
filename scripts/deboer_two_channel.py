#!/usr/bin/env python3
"""de Boer cryo SF6/O2 ARDE: two-channel (radical + directional ion) reduced model.

The de Boer/Blauw normalized-rate experiment is [1.0, 0.43, 0.29, 0.20] at AR [0, 10, 20, 40]. The
VALIDATED radical-transport channel (scripts/arde_mc_reference.mc_reactive_transmission, agreeing
with the common engine's radiosity and an analytic view factor) reproduces the low-AR knee but
collapses too fast at high AR (best radical-only ~0.035 at AR40, max_bounce-converged -> a genuine
physics gap, not MC under-sampling). de Boer SF6/O2 is ion-ASSISTED: directional ions sustain the
high-AR floor that broad radicals cannot. This script adds a reduced directional-ion channel and
fits the experiment, calibrating on the AR10/20 knee and predicting the AR40 floor held-out.

CAVEATS (honest scope): the ion channel here is a REDUCED analytic model (Gaussian cross-slot angle
+ absorbing walls), not the full common-engine ion transport; the rate law is an additive
radical+ion assumption. The radical channel is the validated one. Best fit: s=0.06, ion IAD sigma
~1 deg, beta~0.4 -> NR [1.0, 0.43, 0.279, 0.169] vs [1.0, 0.43, 0.29, 0.20], knee RMSE 0.008,
held-out AR40 error 0.031. The AR40 residual (0.169 vs 0.20) and the near-sub-degree ion IAD are the
remaining frontier (sub-degree IADF and/or charging). Run: python scripts/deboer_two_channel.py
"""
from __future__ import annotations

import sys
import os

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from arde_mc_reference import mc_reactive_transmission


def ion_transmission(ar, sigma_deg, *, opening=0.10, mask=0.05, n=1 << 20, seed=5):
    """Directional-ion floor transmission: ions enter the opening with Gaussian cross-slot angle
    sigma_deg and are absorbed on any sidewall hit (s=1). Straight-line drift over the mask+trench
    tube depth; reaches the floor iff it stays within the opening."""
    etched = ar * opening
    tube = mask + etched
    Wc = 2.0 * opening
    ox0 = (Wc - opening) / 2.0
    ox1 = ox0 + opening
    rng = np.random.default_rng(seed)
    xe = rng.uniform(ox0, ox1, n)
    th = rng.normal(0.0, np.deg2rad(sigma_deg), n)
    xf = xe + tube * np.tan(th)
    return float(np.mean((xf >= ox0) & (xf <= ox1)))


def main():
    exp_ar = [0.0, 10.0, 20.0, 40.0]
    exp = np.array([1.0, 0.43, 0.29, 0.20])
    best = None
    for s in (0.03, 0.04, 0.05, 0.06, 0.08):
        rad = {a: mc_reactive_transmission(max(a, 0.0), s, log2n=18, max_bounce=2000) for a in exp_ar}
        for sig in (0.75, 1.0, 1.5, 2.0):
            ion = {a: ion_transmission(a, sig) for a in exp_ar}
            for beta in np.linspace(0.2, 8.0, 40):
                nr = np.array([(rad[a] + beta * ion[a]) / (rad[0.0] + beta * ion[0.0]) for a in exp_ar])
                cal = np.sqrt(np.mean((nr[1:3] - exp[1:3]) ** 2))   # calibrate on AR10, AR20
                if best is None or cal < best[0]:
                    best = (cal, s, sig, beta, nr)
    cal, s, sig, beta, nr = best
    print("de Boer two-channel reduced model (calibrate AR10,20; predict AR40 held-out):")
    print(f"  radical sticking s={s}, ion IAD sigma={sig} deg, ion/radical strength beta={beta:.2f}")
    print(f"  model NR   = {[round(float(v), 3) for v in nr]}")
    print(f"  experiment = {list(exp)}")
    print(f"  calibration RMSE(AR10,20) = {cal:.3f}")
    print(f"  held-out AR40: predicted {nr[3]:.3f} vs experiment 0.20  (error {abs(nr[3]-0.20):.3f})")
    print("  NOTE reduced ion model + additive rate law; radical channel is the validated one.")


if __name__ == "__main__":
    main()
