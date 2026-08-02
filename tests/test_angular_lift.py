"""Gates for the axisymmetric angular lift (P1a).

The published Krüger IEAD supplies the ion angle *projected into one plane*.
Lifting it to 3-D must invert that marginalization; the pre-2026-08-02 closure
identified the planar angle with the polar angle and so discarded exactly
sqrt(2) of the angular width (RESULTS_ANGULAR_CONVERGENCE_P0_2026-08-02, EXP C).
"""

import os

import numpy as np
import pytest

from petch.angular_lift import (
    axisymmetric_polar_weights,
    invert_planar_marginal,
    shell_projection_mass,
)
from petch.reactor_boundary import (
    build_krueger_2024_development_boundary,
    load_krueger_2024_digitized_iead,
)

KRUEGER_DATA = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "experimental", "krueger_2024")

# The repo's own digitization band for the published planar width; the P0
# harness measured 0.5893 deg from the broken lift against 0.8334 published.
DIGITIZATION_BAND_DEG = (0.822, 0.860)


def _gaussian_planar(sigma_deg, *, span=6.0, samples=200001):
    angle = np.linspace(-span * sigma_deg, span * sigma_deg, samples)
    weight = np.exp(-angle ** 2 / (2.0 * sigma_deg ** 2))
    return angle, weight / weight.sum()


def _tan_rms_deg(angle_deg, weight):
    """Width in the variable the projection law is exact in (the slope)."""
    tangent = np.tan(np.deg2rad(np.asarray(angle_deg, dtype=float)))
    return np.degrees(np.arctan(np.sqrt(
        np.average(tangent ** 2, weights=np.asarray(weight, dtype=float)))))


def _node_widths(species):
    velocity = np.asarray(species.velocity_sqrt_eV, dtype=float)
    weight = np.asarray(species.weight, dtype=float)
    weight = weight / weight.sum()
    planar = velocity[:, 0] / velocity[:, 2]
    polar = np.hypot(velocity[:, 0], velocity[:, 1]) / velocity[:, 2]
    return (float(np.average(planar ** 2, weights=weight)),
            float(np.average(polar ** 2, weights=weight)))


def test_projection_law_is_the_arcsine_measure():
    """A shell projects onto |theta_x| <= theta with the circle-projection law."""
    # Whole shell lands inside its own radius, nothing beyond it.
    assert shell_projection_mass(1.0, 0.0, 1.0) == pytest.approx(1.0, abs=1e-12)
    assert shell_projection_mass(1.0, 1.0, 5.0) == pytest.approx(0.0, abs=1e-12)
    # Median of the folded arcsine sits at sin(pi/4) = 1/sqrt(2) of the radius.
    half = np.degrees(np.arctan(np.tan(np.deg2rad(2.0)) / np.sqrt(2.0)))
    assert shell_projection_mass(2.0, 0.0, half) == pytest.approx(0.5, abs=1e-12)
    # Bins tile the shell exactly.
    edges = np.linspace(0.0, 3.0, 25)
    assert float(shell_projection_mass(3.0, edges[:-1], edges[1:]).sum()) == (
        pytest.approx(1.0, abs=1e-12))


def test_onion_peel_round_trip_is_exact():
    """The peel is exact back-substitution: projecting the shells reproduces
    the planar histogram to machine precision, with no free parameter."""
    angle, weight = _gaussian_planar(0.8334)
    edges = np.arange(0.0, 5.5, 0.25)
    index = np.clip(np.searchsorted(edges, np.abs(angle), side="right") - 1,
                    0, edges.size - 2)
    planar = np.bincount(index, weights=weight, minlength=edges.size - 1)

    centres, shells, diagnostics = invert_planar_marginal(edges, planar)
    kernel = shell_projection_mass(
        centres[None, :], edges[:-1, None], edges[1:, None])
    assert np.abs(kernel @ shells - planar).max() < 1e-12
    assert diagnostics["clamped_negative_fraction"] <= 0.0
    assert shells.sum() == pytest.approx(planar.sum(), rel=1e-12)
    # Inverting widens: the polar measure carries more weight at large angle.
    assert _tan_rms_deg(centres, shells) > _tan_rms_deg(centres, planar)


def test_lift_preserves_the_exact_factor_two_second_moment():
    """E[tan^2 polar] = 2 E[tan^2 planar] holds for ANY axisymmetric measure —
    it is the analytic content of the sqrt(2), not a Gaussian coincidence.
    The discrete lift honours it exactly because a uniform azimuthal ring
    averages cos^2 to 1/2 bitwise for order >= 3."""
    for kwargs in ({}, {"ion_energy_bin_eV": 500.0, "ion_angle_bin_deg": 0.5}):
        boundary = build_krueger_2024_development_boundary(
            KRUEGER_DATA, reference_plane_m=1.0e-6,
            ion_azimuthal_closure="axisymmetric_uniform", **kwargs)
        planar_moment, polar_moment = _node_widths(boundary.get("ions"))
        assert polar_moment / planar_moment == pytest.approx(2.0, abs=1e-12)


def test_gaussian_planar_marginal_inverts_toward_the_rayleigh_width():
    """A Gaussian planar marginal must invert to the Rayleigh polar measure
    whose rms is sqrt(2) larger, converging as the peel grid refines.

    The residual is the Abel discretization (shells sit at bin centres), and
    it floors near 1e-4 because inverting the projection is ill-posed —
    refining past the source resolution amplifies noise instead of removing
    bias, which is why the production path peels on the 0.25 deg source grid.
    """
    sigma = 0.8334
    angle, weight = _gaussian_planar(sigma)
    expected = _tan_rms_deg(angle, weight) * np.sqrt(2.0)

    coarse = _tan_rms_deg(*axisymmetric_polar_weights(
        angle, weight, bin_deg=0.1)[:2])
    fine = _tan_rms_deg(*axisymmetric_polar_weights(
        angle, weight, bin_deg=0.025)[:2])
    assert coarse == pytest.approx(expected, rel=5e-3)
    assert fine == pytest.approx(expected, rel=1e-3)
    assert abs(fine - expected) < abs(coarse - expected)


def test_real_iead_lift_lands_in_the_published_digitization_band():
    """The gate that matters: the lifted 3-D beam's own planar marginal must
    reproduce the published planar width. Pre-fix this read 0.5893 deg."""
    iead = load_krueger_2024_digitized_iead(KRUEGER_DATA)
    published = _tan_rms_deg(iead.signed_angle_deg, iead.probability_weight)
    assert DIGITIZATION_BAND_DEG[0] <= published <= DIGITIZATION_BAND_DEG[1]

    boundary = build_krueger_2024_development_boundary(
        KRUEGER_DATA, reference_plane_m=1.0e-6,
        ion_azimuthal_closure="axisymmetric_uniform")
    planar_moment, _ = _node_widths(boundary.get("ions"))
    lifted = np.degrees(np.arctan(np.sqrt(planar_moment)))
    assert DIGITIZATION_BAND_DEG[0] <= lifted <= DIGITIZATION_BAND_DEG[1]
    # The P0 deficit (published / lifted = 1.4141) is retired.
    assert published / lifted == pytest.approx(1.0, rel=0.02)


def test_lift_conserves_probability_and_support():
    iead = load_krueger_2024_digitized_iead(KRUEGER_DATA)
    angle = iead.signed_angle_deg
    weight = iead.probability_weight
    polar, lifted, diagnostics = axisymmetric_polar_weights(angle, weight)

    assert float(lifted.sum()) == pytest.approx(float(weight.sum()), rel=1e-12)
    assert np.all(lifted >= 0.0)
    # A shell projects onto its own radius, so the lift reweights the published
    # angular support without extending it.
    assert float(polar.max()) == pytest.approx(float(np.abs(angle).max()),
                                               rel=1e-12)
    assert diagnostics["bin_deg"] == 0.25
    assert diagnostics["clamped_negative_fraction"] <= 0.0

    boundary = build_krueger_2024_development_boundary(
        KRUEGER_DATA, reference_plane_m=1.0e-6,
        ion_azimuthal_closure="axisymmetric_uniform")
    ion = boundary.get("ions")
    assert float(np.sum(ion.weight)) == pytest.approx(1.0, rel=1e-12)
    assert "onion-peel" in ion.provenance["three_dimensional_polar_inversion"]
    assert ion.provenance[
        "three_dimensional_polar_inversion_diagnostics"]["bin_deg"] == 0.25


def test_over_fine_peel_is_refused_not_silently_clamped():
    """Abel inversion below the source resolution diverges into negative
    shells; the guard refuses rather than clamping a manufactured answer."""
    iead = load_krueger_2024_digitized_iead(KRUEGER_DATA)
    with pytest.raises(ValueError, match="diverged"):
        axisymmetric_polar_weights(
            iead.signed_angle_deg, iead.probability_weight, bin_deg=0.10)
    # The same call is admissible once the guard is lifted deliberately.
    _, _, diagnostics = axisymmetric_polar_weights(
        iead.signed_angle_deg, iead.probability_weight, bin_deg=0.10,
        max_clamped_fraction=None)
    assert diagnostics["clamped_negative_fraction"] > 1.0


def test_degenerate_and_invalid_inputs():
    zero = np.zeros(4)
    polar, lifted, diagnostics = axisymmetric_polar_weights(
        zero, np.full(4, 0.25))
    assert diagnostics["degenerate"] is True
    assert np.all(polar == 0.0) and lifted.sum() == pytest.approx(1.0)

    with pytest.raises(ValueError, match="bin_deg"):
        axisymmetric_polar_weights(np.array([1.0]), np.array([1.0]), bin_deg=0.0)
    with pytest.raises(ValueError, match="90 degrees"):
        axisymmetric_polar_weights(np.array([91.0]), np.array([1.0]))
    with pytest.raises(ValueError, match="share a shape"):
        axisymmetric_polar_weights(np.array([1.0, 2.0]), np.array([1.0]))
    with pytest.raises(ValueError, match="ascend"):
        invert_planar_marginal(np.array([0.0, 1.0, 0.5]), np.array([0.5, 0.5]))
    with pytest.raises(ValueError, match="one entry per bin"):
        invert_planar_marginal(np.array([0.0, 1.0, 2.0]), np.array([1.0]))
