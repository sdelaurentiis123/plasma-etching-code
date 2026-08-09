import numpy as np
import pytest

from petch.chlorine_species_resolved_si import (
    SpeciesResolvedChlorineSiMechanism,
)
from petch.surface_kinetics import EnergeticFlux, SurfaceFluxes


def _fluxes(*, sicl2=0.0, cl2_angle_deg=0.0):
    return SurfaceFluxes(
        {"Cl": 2.0e22, "SiCl2": sicl2},
        (
            EnergeticFlux("Cl+", 1.0e20, [35.0], [1.0], [1.0]),
            EnergeticFlux(
                "Cl2+",
                2.0e20,
                [100.0],
                [np.cos(np.deg2rad(cl2_angle_deg))],
                [1.0],
            ),
        ),
    )


def test_species_resolved_join_sums_measured_clplus_and_cl2plus_rates():
    mechanism = SpeciesResolvedChlorineSiMechanism()
    result = mechanism.advance(
        mechanism.initial_state(), _fluxes(), 2.0, strict=False
    )
    theta = (0.18 * 200.0 + 1.0) / (0.18 * 200.0 + 1.0 + 4.0 * 1.14)
    clplus_rate = 1.0e20 * 1.14 * theta
    cl2_yield = 0.22 * (np.sqrt(100.0) - np.sqrt(25.998846756576185))
    cl2plus_rate = 2.0e20 * cl2_yield

    assert result.chlorination_fraction == pytest.approx(theta)
    assert result.clplus_removal_rate_si_m2_s == pytest.approx(clplus_rate)
    assert result.cl2plus_removal_rate_si_m2_s == pytest.approx(cl2plus_rate)
    assert result.etch_velocity_m_s == pytest.approx(
        (clplus_rate + cl2plus_rate) / 5.0e28
    )
    assert result.removed_si_atoms_m2 == pytest.approx(
        2.0 * (clplus_rate + cl2plus_rate)
    )
    assert result.material_exchange.product_routing_complete


def test_species_resolved_join_applies_eq5_6_before_clplus_removal():
    mechanism = SpeciesResolvedChlorineSiMechanism()
    clean = mechanism.advance(
        mechanism.initial_state(), _fluxes(), 1.0, strict=False
    )
    returned = mechanism.advance(
        mechanism.initial_state(), _fluxes(sicl2=1.0e21), 1.0, strict=False
    )

    expected_theta = (36.0 + 1.0 + 3.0) / (
        36.0 + 1.0 + 3.0 + 4.0 * 1.14 + 90.0
    )
    assert returned.chlorination_fraction == pytest.approx(expected_theta)
    assert returned.sicl2_to_clplus_flux_ratio == pytest.approx(10.0)
    assert returned.clplus_removal_rate_si_m2_s < (
        clean.clplus_removal_rate_si_m2_s
    )
    # The Balooch branch is independently measured on a saturated surface;
    # Eq. 5.6 does not supply a coverage-rescaling law for it.
    assert returned.cl2plus_removal_rate_si_m2_s == pytest.approx(
        clean.cl2plus_removal_rate_si_m2_s
    )


def test_source_bounded_cl2plus_coverage_sensitivity_uses_printed_endpoints():
    mechanism = SpeciesResolvedChlorineSiMechanism(
        cl2plus_coverage_mode="source_bounded_linear",
    )
    clean = mechanism.advance(
        mechanism.initial_state(), _fluxes(), 1.0, strict=False
    )
    returned = mechanism.advance(
        mechanism.initial_state(), _fluxes(sicl2=1.0e21), 1.0, strict=False
    )
    saturated_yield = 0.22 * (
        np.sqrt(100.0) - np.sqrt(25.998846756576185)
    )
    bare_fraction = np.sqrt(2.0) * 0.06 / 0.22
    expected_scale = bare_fraction + (
        1.0 - bare_fraction
    ) * returned.chlorination_fraction

    assert returned.cl2plus_removal_rate_si_m2_s == pytest.approx(
        2.0e20 * saturated_yield * expected_scale
    )
    assert returned.cl2plus_removal_rate_si_m2_s < (
        clean.cl2plus_removal_rate_si_m2_s
    )
    assert "cl2plus_coverage_interpolation" in (
        returned.validity.nonpredictive_parameters
    )


def test_species_resolved_join_fails_closed_on_unknown_cl2plus_angle():
    mechanism = SpeciesResolvedChlorineSiMechanism()
    with pytest.raises(ValueError, match="normal incidence"):
        mechanism.advance(
            mechanism.initial_state(), _fluxes(cl2_angle_deg=20.0), 1.0
        )
    sensitivity = mechanism.advance(
        mechanism.initial_state(),
        _fluxes(cl2_angle_deg=20.0),
        1.0,
        strict=False,
    )
    assert not sensitivity.validity.within_declared_scope
    assert "cl2plus_angular_response" in (
        sensitivity.validity.nonpredictive_parameters
    )
