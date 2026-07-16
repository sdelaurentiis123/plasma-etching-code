"""Smoke + regression tests for the petch package. Run: `pytest tests/`.

Covers: the 2D reference path (parity anchor), the 3D production engine (the real simulator),
and the high-level ViennaPS-shaped API (Domain/Process/Result). All run CPU-only (slower) or GPU.
"""
from pathlib import Path
import re

import numpy as np
import petch


POC = petch.Flags(chemistry="langmuir", yield_angular="cosine")   # original 2D proof-of-concept config


def test_runtime_version_matches_project_metadata():
    """The import-visible version must agree with the version stamped into built wheels."""
    project_text = (Path(__file__).parents[1] / "pyproject.toml").read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', project_text, flags=re.MULTILINE)
    assert match is not None, "pyproject.toml has no declared project version"
    assert petch.__version__ == match.group(1)


# ----------------------------- 2D reference path (legacy parity anchor) -----------------------------

def test_small_run():
    """2D sanity: the reference simulator runs and produces a sane profile."""
    r = petch.run_etch(W=12, H=12, dx=0.5, trench_width=5, mask_thickness=2.0,
                       sub_top=9, t_end=0.5, n_steps=3, par=dict(petch.PAR), flags=POC, verbose=False)
    assert len(r['segs']) > 0
    d = petch.center_depth(r)
    assert 0.0 < d < 9.0          # some etching happened, not absurd
    assert r['timings']['total'] > 0.0


def test_baseline_parity():
    """2D PoC parity: width-8 reproduces the proof-of-concept baseline depth (~9.18 um) at rate_scale=0.29.
    The 2D path is the historical reference implementation; production use is the 3D API below."""
    OURS = dict(petch.PAR); OURS['rate_scale'] = 0.29
    r = petch.run_etch(W=20.0, H=24.0, dx=0.25, trench_width=8.0, mask_thickness=2.0,
                       sub_top=18.0, t_end=3.0, n_steps=60, par=OURS, flags=POC, verbose=False)
    d = petch.center_depth(r)
    assert abs(d - 9.18) < 0.5, f"baseline depth {d:.3f} drifted from 9.18"   # margin for MC variance


# ----------------------------- 3D legacy compatibility engine -----------------------------

def test_3d_loop():
    """3D etch runs end-to-end with the faithful ViennaPS config (belen coverages + ion reflection)
    and deepens the feature. Asserts on max_depth (footprint metric) + no NaN blowup."""
    from petch import threed as t3
    par = dict(petch.PAR); par['rate_scale'] = 0.1; par['periodic_y'] = 1
    fl = petch.Flags(chemistry="belen", yield_angular="viennaps", coverage_sticking=True,
                     warm_start_coverage=True, sampling="sobol", ion_reflection=True)
    geo = t3.run_etch_3d(Lx=3.0, Ly=0.9, Lz=6.5, dx=0.15, trench_width=0.9, mask_th=0.3, sub_top=6.0,
                         t_end=1.5, n_steps=20, par=par, flags=fl, n_ion=8000, n_neu=8000,
                         reinit_method="fsm", verbose=False)
    md = t3.max_depth_3d(geo); cd = t3.center_depth_3d(geo)
    assert np.isfinite(geo['phi']).all(), "level set blew up (non-finite phi)"
    assert md > 0.3, f"3D etch produced no depth (max_depth={md})"
    assert cd > 0.1, f"center did not etch (center_depth={cd}, max_depth={md})"


# ----------------------------- high-level public API -----------------------------

def test_api_trench():
    """The ViennaPS-shaped public API (Domain/SF6O2/Process/Result) runs and reports a sane etch."""
    dom = petch.Domain.trench(extent=3.0, dx=0.15, width=0.9, mask=0.3, depth=6.0)
    res = petch.Process(dom, petch.SF6O2(rate_scale=0.1), duration=1.5).run(steps=20)
    assert res.engine == "legacy-threed-v1"
    assert res.max_depth > 0.3, f"API etch shallow (max_depth={res.max_depth})"
    assert res.aspect_ratio > 0.3
    assert res.wall_time > 0.0
    v, f = res.mesh
    assert len(v) > 0 and len(f) > 0          # surface mesh extracted


if __name__ == "__main__":
    test_small_run(); print("test_small_run OK")
    test_baseline_parity(); print("test_baseline_parity OK")
    test_3d_loop(); print("test_3d_loop OK")
    test_api_trench(); print("test_api_trench OK")
