"""Smoke + parity tests for the petch package.

- test_small_run: fast sanity (small grid) — the simulator runs and produces a sane profile.
- test_baseline_parity: full width-8 reproduces the PoC baseline depth (~9.18 um) with the
  langmuir default and rate_scale=0.29. This is the parity anchor for the package refactor.
"""
import numpy as np
import petch


POC = petch.Flags(chemistry="langmuir", yield_angular="cosine")   # original PoC config


def test_small_run():
    par = dict(petch.PAR); par['cal_F'] = 1.0   # PoC (uncalibrated) for the parity-style check
    r = petch.run_etch(W=12, H=12, dx=0.5, trench_width=5, mask_thickness=2.0,
                       sub_top=9, t_end=0.5, n_steps=3, par=par, flags=POC, verbose=False)
    assert len(r['segs']) > 0
    d = petch.center_depth(r)
    assert 0.0 < d < 9.0          # some etching happened, not absurd
    assert r['timings']['total'] > 0.0


def test_baseline_parity():
    OURS = dict(petch.PAR); OURS['rate_scale'] = 0.29; OURS['cal_F'] = 1.0   # uncalibrated PoC
    r = petch.run_etch(W=20.0, H=24.0, dx=0.25, trench_width=8.0, mask_thickness=2.0,
                       sub_top=18.0, t_end=3.0, n_steps=60, par=OURS, flags=POC, verbose=False)
    d = petch.center_depth(r)
    # PoC baseline on this machine = 9.18 um; allow margin for MC variance.
    assert abs(d - 9.18) < 0.4, f"baseline depth {d:.3f} drifted from 9.18"


def test_calibrated_match():
    """Regression guard: the calibrated default config lands the width-8 depth near ViennaPS
    (10.05) -- catches accidental breakage of the cal_F / chemistry / angular-yield calibration."""
    par = dict(petch.PAR); par['rate_scale'] = 0.034   # cal_F=12 default
    r = petch.run_etch(W=20.0, H=24.0, dx=0.25, trench_width=8.0, mask_thickness=2.0,
                       sub_top=18.0, t_end=3.0, n_steps=60, par=par, verbose=False)  # default flags
    d = petch.center_depth(r)
    assert 8.5 < d < 12.0, f"calibrated width-8 depth {d:.2f} drifted from ~10 um"


def test_3d_loop():
    """3D etch loop runs end-to-end (Warp flux kernel) and deepens a feature."""
    from petch import threed as t3
    par = dict(petch.PAR); par['rate_scale'] = 0.05   # cal_F=12 default -> keep tiny grid sane
    geo = t3.run_etch_3d(Lx=8, Ly=4, Lz=12, dx=0.5, trench_width=3, mask_th=2, sub_top=8,
                         t_end=1.5, n_steps=6, par=par, flags=petch.Flags(chemistry="langmuir"),
                         n_ion=4000, n_neu=4000, verbose=False)
    d = t3.center_depth_3d(geo)
    assert d > 0.1, f"3D etch produced no depth ({d})"


if __name__ == "__main__":
    test_small_run(); print("test_small_run OK")
    test_baseline_parity(); print("test_baseline_parity OK")
    test_3d_loop(); print("test_3d_loop OK")
