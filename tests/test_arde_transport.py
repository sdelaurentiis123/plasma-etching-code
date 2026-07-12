"""Aspect-ratio-dependent-etching (ARDE) neutral-transport gate for the common feature-3d engine.

This is the first ARDE physics gate in the suite. The pre-existing "aspect ratio ladder" test
(test_boundary_transport.py) fires a vertical monodirectional beam into a straight trench and asserts
floor flux == 1 at every AR -- true by geometry, and by construction unable to exhibit ARDE.

Here the common engine's ballistic transport (forward first-hit tracer + a half-Maxwellian *flux*
cosine source, no fitted angular closure) is driven at s=1 (pure line-of-sight shadowing) and its
floor transmission vs aspect ratio is cross-checked against an INDEPENDENT pure-numpy analytic
ray-trace of the same mask+trench aperture (`reference_floor_transmission`, no engine transport code).

Physics under test:
  1. ARDE: floor transmission decreases monotonically with aspect ratio (free-molecular shadowing).
  2. First-principles agreement: the engine reproduces the independent geometric reference (the
     residual is grid error that shrinks with refinement; at dx=0.02 it is ~3-12%, worse at low AR).
  3. AMR necessity: a fixed, coarse angular quadrature OVER-predicts the deep-feature transmission
     because it aliases the ~arctan(1/A) acceptance cone; the QMC-refined source does not.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from deboer_arde_static import floor_transmission, reference_floor_transmission  # noqa: E402
from arde_mc_reference import mc_reactive_transmission  # noqa: E402


def _engine_reactive(ar, s, dx=0.0125, log2=15):
    return floor_transmission(
        ar, s, opening_um=0.10, dx_um=dx, rays_per_face=64,
        n_transverse=0, n_normal=0, transport_method="forward", n_position=32,
        source_method="qmc", log2_samples=log2)["floor_transmission"]


def _engine(ar, dx=0.02, log2=15):
    return floor_transmission(
        ar, 1.0, opening_um=0.10, dx_um=dx, rays_per_face=64,
        n_transverse=0, n_normal=0, transport_method="forward", n_position=16,
        source_method="qmc", log2_samples=log2)["floor_transmission"]


def test_arde_floor_transmission_is_monotone_in_aspect_ratio():
    t4 = _engine(4.0)
    t8 = _engine(8.0)
    assert 0.0 < t8 < t4 < 1.0, (t4, t8)


@pytest.mark.parametrize("ar", [4.0, 8.0])
def test_engine_reproduces_independent_geometric_reference(ar):
    ref = reference_floor_transmission(ar, dx_um=0.02, log2_samples=18)
    eng = _engine(ar)
    # dx=0.02 grid error is worse at low AR; 15% brackets it and shrinks under refinement.
    assert eng == pytest.approx(ref, rel=0.15), (ar, ref, eng)


@pytest.mark.parametrize("n_transverse", [4, 5])
def test_fixed_coarse_quadrature_fails_deep_feature_without_amr(n_transverse):
    # A fixed coarse angular quadrature cannot resolve the ~arctan(1/A) acceptance cone at high AR.
    # It fails in one of two ways: an odd-node grid has a perfectly vertical atom and OVER-predicts
    # (flatlines high); an even-node grid has none and UNDER-predicts (the cone falls between nodes,
    # artificial ~zero). Either way it is badly wrong, while the QMC-refined (AMR) source matches the
    # independent reference. This encodes the adaptive-phase-space-refinement requirement.
    ref8 = reference_floor_transmission(8.0, dx_um=0.02, log2_samples=18)
    fixed = floor_transmission(
        8.0, 1.0, opening_um=0.10, dx_um=0.02, rays_per_face=64,
        n_transverse=n_transverse, n_normal=8)["floor_transmission"]
    assert (fixed > 2.0 * ref8) or (fixed < 0.5 * ref8), (n_transverse, fixed, ref8)
    assert _engine(8.0) == pytest.approx(ref8, rel=0.15), (ref8,)


@pytest.mark.parametrize("ar,s", [(1.0, 0.5), (2.0, 0.3)])
def test_reactive_radiosity_matches_independent_particle_mc(ar, s):
    # The engine's deterministic diffuse radiosity (re-emission with sticking s) must agree with a
    # completely independent stochastic particle Monte Carlo of the same box. Two different methods
    # agreeing validates the reactive transport. Tolerance covers grid (continuum MC vs dx mesh),
    # QMC, and MC noise.
    mc = mc_reactive_transmission(ar, s, dx=0.0125, log2n=16)
    eng = _engine_reactive(ar, s)
    assert eng == pytest.approx(mc, rel=0.08), (ar, s, mc, eng)


def test_reactive_family_monotone_in_ar_and_ordered_in_sticking():
    # Coburn-Winters: floor transmission decreases with aspect ratio at fixed sticking, and decreases
    # with sticking at fixed aspect ratio (more wall consumption -> less flux reaches the floor).
    t = {(ar, s): _engine_reactive(ar, s, dx=0.02, log2=14) for ar in (1.0, 4.0) for s in (0.1, 0.5)}
    for s in (0.1, 0.5):
        assert t[(4.0, s)] < t[(1.0, s)], (s, t)
    for ar in (1.0, 4.0):
        assert t[(ar, 0.5)] < t[(ar, 0.1)], (ar, t)
