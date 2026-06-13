"""Smoke + parity tests for the petch package.

- test_small_run: fast sanity (small grid) — the simulator runs and produces a sane profile.
- test_baseline_parity: full width-8 reproduces the PoC baseline depth (~9.18 um) with the
  langmuir default and rate_scale=0.29. This is the parity anchor for the package refactor.
"""
import numpy as np
import petch


def test_small_run():
    r = petch.run_etch(W=12, H=12, dx=0.5, trench_width=5, mask_thickness=2.0,
                       sub_top=9, t_end=0.5, n_steps=3, verbose=False)
    assert len(r['segs']) > 0
    d = petch.center_depth(r)
    assert 0.0 < d < 9.0          # some etching happened, not absurd
    assert r['timings']['total'] > 0.0


def test_baseline_parity():
    OURS = dict(petch.PAR); OURS['rate_scale'] = 0.29
    r = petch.run_etch(W=20.0, H=24.0, dx=0.25, trench_width=8.0, mask_thickness=2.0,
                       sub_top=18.0, t_end=3.0, n_steps=60, par=OURS, verbose=False)
    d = petch.center_depth(r)
    # PoC baseline on this machine = 9.18 um; allow margin for MC variance.
    assert abs(d - 9.18) < 0.4, f"baseline depth {d:.3f} drifted from 9.18"


if __name__ == "__main__":
    test_small_run(); print("test_small_run OK")
    test_baseline_parity(); print("test_baseline_parity OK")
