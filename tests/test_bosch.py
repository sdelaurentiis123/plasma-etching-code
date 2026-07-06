"""Gate the Bosch DRIE cycle driver against published SEM measurements (BOSCH_BENCHMARK_SPEC.md).
Config R = Ayon et al. 1999 (STS-HRM, via the McVittie NNIN deck): 2 um trench, 65 x 3.5 s cycles.
Two per-cycle rates are calibrated from published endpoints (r_iso from the undercut/arc, d_dir from
the pitch); depth, pitch, scallop and undercut are then all EMERGENT and scored with no further
tuning. The sequential punch-then-iso semantics is load-bearing: the simultaneous model provably
gives s ~ 32 nm (4x under Ayon's measured 140 nm)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from petch.bosch import run_bosch


def _config_r():
    return run_bosch(width_um=2.0, n_cycles=65, r_iso_um=0.238, d_dir_um=0.196, dx_um=0.02)


def test_config_r_all_four_gates():
    """Ayon 1999: D=28.2+-10%, p=434+-10%, s=140+-35 nm, U=250+-50 nm — all four simultaneously."""
    r = _config_r()
    assert abs(r["depth"] - 28.2) < 2.9, r["depth"]
    assert abs(r["pitch"] * 1000 - 434) < 44, r["pitch"]
    assert abs(r["scallop"] * 1000 - 140) < 35, r["scallop"]
    assert abs(r["undercut"] * 1000 - 250) < 50, r["undercut"]


def test_scallop_scales_with_etch_time():
    """Cross-config physics: scallop depth grows with the per-cycle etch (roughly linearly).
    Halving the per-cycle rates must shrink the scallop by >=1.5x (kills accidentally-right models)."""
    r_full = _config_r()
    r_half = run_bosch(width_um=2.0, n_cycles=40, r_iso_um=0.119, d_dir_um=0.098, dx_um=0.02)
    assert r_half["scallop"] < r_full["scallop"] / 1.5


def test_sequential_semantics_analytic():
    """The sequential-mechanics arc formula: s = r - sqrt(r^2 - (p/2)^2) with r=238, p=434 -> 140 nm.
    The grid must reproduce its own analytic geometry (mechanics-to-math consistency)."""
    r = _config_r()
    p = r["pitch"] * 1000
    s_pred = 238.0 - np.sqrt(238.0 ** 2 - (p / 2) ** 2)
    assert abs(r["scallop"] * 1000 - s_pred) < 35
