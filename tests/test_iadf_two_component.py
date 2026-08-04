"""Gates for the two-component (core + collisional tail) ion angular distribution.

The physics targets are measurements, not simulation outputs:
Kim et al., Jpn. J. Appl. Phys. 64, 05SP15 (2025) for the widths, and the repo's
own digitized Krueger Figure-4 marginal (with its 0.25 degree digitization band)
for representability of a published bi-Gaussian.
"""

import numpy as np
import pytest
from scipy.integrate import quad

from petch.iadf_two_component import (
    AngularComponent,
    KIM_2025_CORE_TEMPERATURE_EV,
    KIM_2025_TAIL_TEMPERATURE_EV,
    KRUEGER_2024_FIGURE4_BIGAUSSIAN,
    TwoComponentIADF,
    acceptance_half_angle_deg,
    build_two_component_boundary,
    kim_2025_reference_iadf,
    krueger_2024_figure4_iadf,
)

# Krueger base-case mean ion energy; the energy at which the research document's
# section A.7 acceptance table is tabulated.
_KRUEGER_MEAN_ENERGY_EV = 3465.0


def _numeric_cone_acceptance(iadf, half_angle_deg, energy_eV):
    """P(theta <= alpha) by quadrature over the Rayleigh density in tan-space."""
    upper = np.tan(np.deg2rad(half_angle_deg))
    total = 0.0
    for component in iadf.components:
        sigma = float(np.asarray(component.sigma_planar_rad(energy_eV)).reshape(()))
        value, _ = quad(lambda r: (r / sigma ** 2) * np.exp(-0.5 * (r / sigma) ** 2),
                        0.0, upper, epsabs=1e-13, epsrel=1e-13, limit=400)
        total += component.fraction * value
    return total


def _numeric_planar_acceptance(iadf, half_angle_deg, energy_eV):
    """P(|theta_x| <= alpha) by quadrature over the Gaussian planar marginal."""
    upper = np.tan(np.deg2rad(half_angle_deg))
    total = 0.0
    for component in iadf.components:
        sigma = float(np.asarray(component.sigma_planar_rad(energy_eV)).reshape(()))
        norm = 1.0 / (sigma * np.sqrt(2.0 * np.pi))
        half, _ = quad(lambda t: norm * np.exp(-0.5 * (t / sigma) ** 2),
                       0.0, upper, epsabs=1e-13, epsrel=1e-13, limit=400)
        total += component.fraction * 2.0 * half
    return total


# --- gate (a): analytic acceptance == numerical quadrature -------------------

def test_analytic_acceptance_matches_quadrature():
    """Both closed forms must equal direct integration of their own density."""
    iadf = kim_2025_reference_iadf()
    energy = _KRUEGER_MEAN_ENERGY_EV
    for half_angle in (0.1, 0.2, 0.286, 0.5, 1.0, 2.0, 5.0, 10.0):
        cone = float(iadf.acceptance_fraction_cone(half_angle, energy))
        planar = float(iadf.acceptance_fraction_planar(half_angle, energy))
        assert cone == pytest.approx(
            _numeric_cone_acceptance(iadf, half_angle, energy), abs=1e-10)
        assert planar == pytest.approx(
            _numeric_planar_acceptance(iadf, half_angle, energy), abs=1e-10)
        # A round hole rejects more than the planar projection -- both
        # transverse components must be small simultaneously -- strictly so
        # wherever the acceptance is not saturated at 1 in double precision.
        assert cone <= planar
        if planar < 1.0 - 1e-12:
            assert cone < planar


def test_cone_and_planar_conventions_are_distinct_at_ar200():
    """The conventions must not be silently interchangeable at the AR-200 cone."""
    iadf = kim_2025_reference_iadf()
    alpha = acceptance_half_angle_deg(200.0)
    cone = float(iadf.acceptance_fraction_cone(alpha, _KRUEGER_MEAN_ENERGY_EV))
    planar = float(iadf.acceptance_fraction_planar(alpha, _KRUEGER_MEAN_ENERGY_EV))
    assert planar - cone > 0.1


def test_section_a7_acceptance_table_reproduced():
    """Regression gate: the published section A.7 table, to its printed digits."""
    energy = _KRUEGER_MEAN_ENERGY_EV
    expected = {
        0.50: (1.000, 0.986, 0.865, 0.686),
        0.65: (1.000, 0.982, 0.824, 0.606),
    }
    for tail_fraction, row in expected.items():
        iadf = kim_2025_reference_iadf(tail_fraction)
        for aspect, target in zip((30.0, 50.0, 100.0, 200.0), row):
            alpha = acceptance_half_angle_deg(aspect)
            assert float(iadf.acceptance_fraction_planar(alpha, energy)) == pytest.approx(
                target, abs=5e-4)
    krueger = krueger_2024_figure4_iadf()
    for aspect, target in zip((30.0, 50.0, 100.0, 200.0),
                              (0.971, 0.834, 0.527, 0.283)):
        alpha = acceptance_half_angle_deg(aspect)
        assert float(krueger.acceptance_fraction_planar(alpha, energy)) == pytest.approx(
            target, abs=5e-4)


# --- gate (b): reproduce the Kim 2025 measured widths ------------------------

def test_kim_2025_measured_widths_reproduced():
    """Declared temperatures must reproduce the published widths and ratios.

    Sheath-speed conversion: the sheath field is normal and does no work on the
    transverse velocity, so the angular width is set by the ratio of the
    unchanged transverse thermal energy to the *impact* (normal) energy the
    sheath voltage delivers -- ``tan(theta) = v_perp / v_z`` with
    ``v_z = sqrt(2 e V_sheath / M)``, giving ``s = sqrt(T_perp / (2 E))`` and the
    2-D radial spread ``theta_th = sqrt(T_perp / E)``.  Both are mass
    independent, which is why the measurement reports the same main-component
    width for Ar+ and Kr+.
    """
    iadf = kim_2025_reference_iadf()
    # Khrabrov & Kaganovich, arXiv:2604.04214v2 section 6.2, verbatim: for
    # T_perp = 0.044 eV and E_b = 1 keV, "theta_th = sqrt(T_perp/E_b) rad or
    # approximately 0.4 degrees".
    core_theta_th = float(iadf.theta_thermal_deg(1000.0)[0])
    assert core_theta_th == pytest.approx(0.380, abs=5e-3)
    assert 0.35 < core_theta_th < 0.45

    # "tail component 0.57 eV" is 13x hotter than the 0.044 eV core, hence
    # sqrt(13) = 3.6x wider -- the published ratio.
    assert (KIM_2025_TAIL_TEMPERATURE_EV / KIM_2025_CORE_TEMPERATURE_EV) == pytest.approx(
        12.95, abs=0.05)
    sigma = iadf.sigma_planar_deg(_KRUEGER_MEAN_ENERGY_EV)
    assert float(sigma[1] / sigma[0]) == pytest.approx(3.6, abs=0.02)

    # Width falls as E**-0.5, monotone across the measurement's 1.4-2.0 keV band.
    energies = np.linspace(1400.0, 2000.0, 13)
    widths = np.array([float(iadf.sigma_planar_deg(value)[0]) for value in energies])
    assert np.all(np.diff(widths) < 0.0)
    ratio = widths / (energies ** -0.5)
    assert np.allclose(ratio, ratio[0], rtol=1e-12)

    # The core temperature is the neutral gas temperature, 0.044 eV = 511 K.
    assert KIM_2025_CORE_TEMPERATURE_EV / 8.6173e-5 == pytest.approx(511.0, abs=2.0)


# --- gate (c): single-component limit == existing Gaussian machinery ---------

def test_single_component_limit_matches_gauss_hermite_machinery():
    """f_tail = 0 must reproduce the transverse Gauss-Hermite second moment.

    ``collisionless_sheath_boundary_state`` builds its transverse tensor as
    ``sqrt(T) * hermgauss nodes`` with weights ``w / sqrt(pi)``; the resulting
    ``E[tan^2 theta_x]`` is ``T / (2 E)``, which is exactly ``s^2`` of the
    two-component model in its single-component limit.
    """
    temperature = KIM_2025_CORE_TEMPERATURE_EV
    energy = _KRUEGER_MEAN_ENERGY_EV
    nodes, weights = np.polynomial.hermite.hermgauss(3)
    transverse = np.sqrt(temperature) * nodes
    weights = weights / np.sqrt(np.pi)
    gh_second_moment = float(np.dot(weights, transverse ** 2) / energy)

    single = TwoComponentIADF(
        components=(AngularComponent(fraction=1.0, temperature_eV=temperature,
                                     label="thermal_core"),),
        provenance={"source": "single-component limit gate"})
    sigma = float(single.sigma_planar_deg(energy)[0])
    assert np.deg2rad(sigma) ** 2 == pytest.approx(gh_second_moment, rel=1e-12)

    # Same limit reached through the two-component constructor with zero tail.
    zero_tail = kim_2025_reference_iadf(0.0)
    assert float(zero_tail.marginal_sigma_planar_deg(energy)) == pytest.approx(
        sigma, rel=1e-12)
    alpha = acceptance_half_angle_deg(200.0)
    assert float(zero_tail.acceptance_fraction_planar(alpha, energy)) == pytest.approx(
        0.953, abs=5e-4)  # section A.7 "thermal only, measured core" row


def test_polar_rms_exceeds_planar_by_sqrt_two():
    """The P1a identity, rederived analytically rather than by inversion."""
    iadf = kim_2025_reference_iadf()
    energy = _KRUEGER_MEAN_ENERGY_EV
    planar = float(iadf.marginal_sigma_planar_deg(energy))
    assert float(iadf.polar_rms_deg(energy)) == pytest.approx(
        np.sqrt(2.0) * planar, rel=1e-14)


# --- gate (d): Krueger digitized bi-Gaussian round-trip ---------------------

def test_krueger_figure4_bigaussian_round_trips_within_digitization_band():
    """The fitted bi-Gaussian must be representable and land in the repo band."""
    iadf = krueger_2024_figure4_iadf()
    energy = _KRUEGER_MEAN_ENERGY_EV
    low, high = KRUEGER_2024_FIGURE4_BIGAUSSIAN["marginal_sigma_band_deg"]
    marginal = float(iadf.marginal_sigma_planar_deg(energy))
    assert low <= marginal <= high
    assert marginal == pytest.approx(0.8412, abs=1e-3)
    fractions = [component.fraction for component in iadf.components]
    assert fractions[1] == pytest.approx(0.65, abs=1e-12)

    # Round trip through the discrete export: the quadrature must recover the
    # analytic polar second moment of the same measure.
    polar_deg, mass = iadf.polar_quadrature(energy, n_polar=4096, max_sigma=12.0)
    recovered = np.sqrt(np.dot(mass, np.deg2rad(polar_deg) ** 2) / mass.sum())
    assert np.rad2deg(recovered) == pytest.approx(
        float(iadf.polar_rms_deg(energy)), rel=2e-3)


def test_digitized_widths_are_energy_independent():
    """A digitized width is a fixed angle, not a temperature: no E dependence."""
    iadf = krueger_2024_figure4_iadf()
    assert float(iadf.marginal_sigma_planar_deg(1000.0)) == pytest.approx(
        float(iadf.marginal_sigma_planar_deg(4000.0)), rel=1e-14)


# --- gate (e): flux conservation --------------------------------------------

def test_quadrature_conserves_flux():
    iadf = kim_2025_reference_iadf()
    for n_polar in (16, 64, 256):
        _, mass = iadf.polar_quadrature(_KRUEGER_MEAN_ENERGY_EV, n_polar=n_polar)
        assert mass.sum() == pytest.approx(1.0, abs=1e-12)
        assert np.all(mass >= 0.0)
    energy, polar, weight = iadf.discrete_nodes(
        np.array([2000.0, 3465.0, 5000.0]), np.array([0.2, 0.5, 0.3]), n_polar=32)
    assert weight.sum() == pytest.approx(1.0, abs=1e-12)
    assert energy.shape == polar.shape == weight.shape


def test_boundary_state_export_is_valid_and_conserves_weight():
    iadf = kim_2025_reference_iadf()
    boundary = build_two_component_boundary(
        iadf, 1.0e19, np.array([3000.0, 3465.0]), energy_weight=np.array([0.4, 0.6]),
        n_polar=24, azimuthal_order=8)
    species = boundary.species[0]
    assert species.weight.sum() == pytest.approx(1.0, abs=1e-12)
    assert species.velocity_sqrt_eV.shape[0] == species.weight.size == 24 * 2 * 8
    assert np.all(species.velocity_sqrt_eV[:, 2] > 0.0)
    energies = species.kinetic_energy_eV
    assert float(energies.min()) == pytest.approx(3000.0, rel=1e-9)
    assert species.mean_energy_eV == pytest.approx(0.4 * 3000.0 + 0.6 * 3465.0, rel=1e-6)
    assert species.provenance["model"] == "bi_gaussian_core_plus_collisional_tail"
    assert "[VERIFY]" in species.provenance["iadf_provenance"]["tail_fraction"]
    # Azimuth is closed uniformly: the transverse mean must vanish.
    transverse_mean = np.dot(species.weight, species.velocity_sqrt_eV[:, :2])
    assert np.allclose(transverse_mean, 0.0, atol=1e-12)


# --- declaration discipline --------------------------------------------------

def test_fractions_must_be_declared_not_normalized():
    with pytest.raises(ValueError, match="sum to 1"):
        TwoComponentIADF(
            components=(
                AngularComponent(fraction=0.4, temperature_eV=0.044),
                AngularComponent(fraction=0.4, temperature_eV=0.57),
            ),
            provenance={"source": "gate"})


def test_provenance_is_mandatory():
    with pytest.raises(ValueError, match="provenance"):
        TwoComponentIADF(
            components=(AngularComponent(fraction=1.0, temperature_eV=0.044),),
            provenance={})


def test_component_requires_exactly_one_width_source():
    with pytest.raises(ValueError, match="exactly one"):
        AngularComponent(fraction=1.0, temperature_eV=0.044, sigma_planar_deg=0.6)
    with pytest.raises(ValueError, match="exactly one"):
        AngularComponent(fraction=1.0)
