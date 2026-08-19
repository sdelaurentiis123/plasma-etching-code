import numpy as np
import pytest

from petch.rate_normalized_removal import (
    RateNormalizedRemovalMechanism,
    RateNormalizedRemovalParameters,
)
from petch.surface_kinetics import (
    FaceResolvedEnergeticFlux,
    ParameterEvidence,
    SurfaceFluxes,
)


def _mechanism(*, predictive_rate=False):
    return RateNormalizedRemovalMechanism(RateNormalizedRemovalParameters(
        material_name="TiO2",
        material_inventory_name="tio2_formula_units",
        projectile_species=("ion",),
        reference_projectile_flux_m2_s=2.0e19,
        blanket_removal_velocity_m_s=40.0e-9 / 60.0,
        bulk_material_unit_density_m3=3.2e28,
        evidence={
            "reference_projectile_flux_m2_s": ParameterEvidence(
                "conserved reactor boundary", "calculated",
                supports_prediction_within_declared_domain=True,
            ),
            "blanket_removal_velocity_m_s": ParameterEvidence(
                "cross-machine process analog", "cross_machine_analog",
                supports_prediction_within_declared_domain=predictive_rate,
            ),
            "bulk_material_unit_density_m3": ParameterEvidence(
                "measured ALD density range", "measured_cross_process",
                supports_prediction_within_declared_domain=True,
            ),
        },
    ))


def _flux(values, name="ion"):
    values = np.asarray(values, dtype=float)
    return SurfaceFluxes({}, (FaceResolvedEnergeticFlux(
        name,
        len(values),
        np.arange(len(values)),
        values,
        np.full(len(values), 200.0),
        np.ones(len(values)),
    ),))


def test_rate_normalized_law_recovers_blanket_speed_and_local_dose_ratio():
    mechanism = _mechanism()
    state = mechanism.initial_state((3,))
    result = mechanism.advance(state, _flux([2.0e19, 1.0e19, 0.0]), 3.0)

    expected = 40.0e-9 / 60.0
    assert np.allclose(result.etch_velocity_m_s, [expected, 0.5 * expected, 0.0])
    assert np.allclose(
        result.removed_material_units_m2,
        result.etch_velocity_m_s * 3.2e28 * 3.0,
    )
    assert not result.material_exchange.product_routing_complete
    assert not result.validity.parameter_evidence_supports_prediction
    assert result.validity.nonpredictive_parameters == (
        "blanket_removal_velocity_m_s",
    )


def test_rate_normalized_law_refuses_undeclared_positive_population():
    mechanism = _mechanism()
    with pytest.raises(ValueError, match="no rate-normalized removal channel"):
        mechanism.advance(mechanism.initial_state((1,)), _flux([1.0], "other"), 1.0)


def test_rate_normalized_parameters_require_complete_evidence():
    with pytest.raises(ValueError, match="cover every physical input"):
        RateNormalizedRemovalParameters(
            material_name="TiO2",
            material_inventory_name="tio2_formula_units",
            projectile_species=("ion",),
            reference_projectile_flux_m2_s=1.0,
            blanket_removal_velocity_m_s=1.0,
            bulk_material_unit_density_m3=1.0,
            evidence={},
        )
