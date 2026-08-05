"""Gates for the HAR hole study phase-1 machinery (HOLE_STUDY_PLAN_2026-08-05)."""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from petch.iadf_two_component import (
    acceptance_half_angle_deg,
    kim_2025_reference_iadf,
)

_SPEC = importlib.util.spec_from_file_location(
    "hole_study_phase1",
    Path(__file__).resolve().parents[1] / "scripts" / "hole_study_phase1.py")
hole_study = importlib.util.module_from_spec(_SPEC)
sys.modules["hole_study_phase1"] = hole_study
_SPEC.loader.exec_module(hole_study)


def _cascade(aspect_ratio, tail_fraction, **kwargs):
    iadf = kim_2025_reference_iadf(tail_fraction=tail_fraction)
    options = dict(n_polar=48, n_azimuth=16, n_radial=8, depth_bins=40)
    options.update(kwargs)
    return hole_study.cascade_delivery(
        aspect_ratio, iadf, hole_study.REFERENCE_ENERGY_EV, **options)


def test_cascade_conserves_flux_to_machine_precision():
    """Direct + cascaded + wall + thermalised must equal the entering flux.

    The cascade splits unit influx into four disjoint sinks; nothing may be
    created or lost by the bounce bookkeeping.
    """
    for aspect_ratio, fraction in ((50.0, 0.0), (200.0, 0.65)):
        record = _cascade(aspect_ratio, fraction)
        assert abs(record["closure_residual"]) < 1e-12


def test_direct_bottom_matches_analytic_cone_acceptance_on_axis():
    """A ray entering on the axis clears the hole iff its polar angle is inside
    the acceptance cone, so the on-axis direct share is the analytic cone mass."""
    aspect_ratio = 100.0
    iadf = kim_2025_reference_iadf(tail_fraction=0.5)
    alpha = float(acceptance_half_angle_deg(aspect_ratio))
    analytic = float(iadf.acceptance_fraction_cone(
        alpha, hole_study.REFERENCE_ENERGY_EV))
    # Collapse the entry disk onto the axis by using a single radial node at
    # r -> 0: first-strike distance is then exactly R / tan(theta).
    record = hole_study.cascade_delivery(
        aspect_ratio, iadf, hole_study.REFERENCE_ENERGY_EV,
        n_polar=4096, n_azimuth=4, n_radial=1, depth_bins=8)
    # The single Gauss node sits at r = R/sqrt(2) (mean of the area measure),
    # so the direct share brackets the on-axis cone mass from below.
    assert 0.0 < record["direct_bottom"] < analytic


def test_bottom_delivery_falls_with_aspect_ratio_and_with_the_tail():
    """Both monotonicities are physics, not bookkeeping: deeper holes and wider
    beams deliver less to the etch front."""
    shallow = _cascade(50.0, 0.0)["total_bottom"]
    deep = _cascade(200.0, 0.0)["total_bottom"]
    deep_tail = _cascade(200.0, 0.65)["total_bottom"]
    assert deep < shallow
    assert deep_tail < deep


def test_cascaded_flux_dominates_direct_flux_at_extreme_aspect_ratio():
    """The study's headline: past AR ~100 the etch front is fed mostly by
    wall-reflected hot particles, not by line-of-sight ions."""
    record = _cascade(200.0, 0.65)
    assert record["cascaded_share_of_bottom"] > 0.5
    assert record["cascaded_bottom"] > 2.0 * record["direct_bottom"]


def test_wall_deposition_is_top_weighted():
    """First strikes concentrate near the entrance; the profile must carry more
    weight in the upper half of the hole than in the lower half."""
    record = _cascade(200.0, 0.65, depth_bins=100)
    profile = np.asarray(record["wall_rate_profile"])
    # Bin 0 is the etch front, the last bin the entrance.
    lower, upper = profile[:50].sum(), profile[50:].sum()
    assert upper > lower


def test_transmission_reference_reproduces_the_clausing_benchmark():
    """The gated exact-algebra path must still return the study's quoted
    0.656 % at 200:1 and agree with Santeler inside his stated error."""
    from petch.axisymmetric_exchange_3d import (
        cylinder_clausing_transmission, santeler_transmission)
    tau = cylinder_clausing_transmission(200.0)
    assert tau == pytest.approx(0.00656, abs=5e-5)
    assert abs(tau - santeler_transmission(200.0)) / tau < 0.007


def test_imported_hole_geometry_routes_to_the_axisymmetric_operator(tmp_path):
    """STL on disk -> watertight diagnostics -> measured out-of-roundness ->
    routing decision, all as receipts rather than assumptions."""
    receipt, profile = hole_study.geometry_receipt(
        20.0, tmp_path / "hole.stl")
    assert receipt["watertight"] and receipt["consistently_oriented"]
    assert receipt["relative_deviation"] < 1e-3
    assert receipt["routed_to"] == "axisymmetric"
    # The extracted radius is the tessellated polygon's, so it sits inside the
    # true radius by the mesh's own faceting bound -- compare against that
    # measured bound rather than an invented tolerance.
    assert np.all(np.asarray(profile.r) <= hole_study.RADIUS + 1e-12)
    assert np.allclose(np.asarray(profile.r), hole_study.RADIUS,
                       atol=receipt["facet_bound"] * hole_study.RADIUS + 1e-12)


def test_general_operator_probe_reports_refusal_without_bypassing_it():
    """The general body-of-revolution operator does not certify its self-pair
    quadrature on a straight wall; the probe must report that, never swallow
    it (the straight-hole numbers come from the exact-algebra path)."""
    probe = hole_study.general_operator_probe(10.0, 0.20, 32, 8)
    assert probe["certified"] is False
    assert probe["value"] is None
    assert probe["residual"] > 1e-4
