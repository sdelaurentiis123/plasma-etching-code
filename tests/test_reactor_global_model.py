import numpy as np
import pytest

from petch.reactor_global import (
    ArgonGlobalCondition,
    ArgonTransportState,
    CylindricalReactor,
    FixedArgonTransportProvider,
    LeeLiebermanArgonGlobalModel,
    PASCAL_PER_MTORR,
)


def _condition(*, wall_evidence="assumed", power_evidence="assumed"):
    return ArgonGlobalCondition(
        condition_id="manufactured-argon",
        absorbed_power_W=1000.0,
        pressure_Pa=10.0 * PASCAL_PER_MTORR,
        gas_temperature_K=600.0,
        geometry=CylindricalReactor(radius_m=0.1525, length_m=0.075),
        ion_wall_energy_factor_Te=5.0,
        ion_wall_energy_source="lee-lieberman-1994-global stated range 5--8 Te",
        ion_wall_energy_evidence=wall_evidence,
        absorbed_power_source="manufactured absorbed power",
        absorbed_power_evidence=power_evidence,
        absorbed_power_boundary_kind="manufactured_test",
    )


def _provider(*, evidence="assumed"):
    return FixedArgonTransportProvider(ArgonTransportState(
        ion_mean_free_path_m=0.004,
        ambipolar_diffusion_m2_s=50.0,
        metastable_effective_diffusion_m2_s=1.0,
        source="manufactured transport closure",
        evidence_kind=evidence,
    ))


def test_argon_global_model_closes_particle_and_power_balances():
    solution = LeeLiebermanArgonGlobalModel(_provider()).solve(_condition())
    assert 0.1 < solution.electron_temperature_eV < 100.0
    assert 0.0 < solution.electron_density_m3 < solution.ground_density_m3
    assert 0.0 < solution.metastable_density_m3 < solution.ground_density_m3
    assert solution.ion_density_m3 == solution.electron_density_m3
    assert solution.axial_ion_flux_m2_s > 0.0
    assert solution.radial_ion_flux_m2_s > 0.0
    assert solution.maximum_normalized_residual <= 1.0e-8
    assert np.isclose(
        solution.modeled_power_loss_W,
        solution.absorbed_power_W,
        rtol=1.0e-8,
        atol=0.0,
    )
    assert not solution.supports_prediction


def test_predictive_flag_requires_both_transport_and_wall_energy_evidence():
    modeled_transport = LeeLiebermanArgonGlobalModel(
        _provider(evidence="validated_model"))
    assert not modeled_transport.solve(_condition()).supports_prediction
    predictive = modeled_transport.solve(
        _condition(
            wall_evidence="validated_model",
            power_evidence="validated_model",
        ))
    assert predictive.supports_prediction


def test_residual_api_rejects_nonpositive_candidate_state():
    model = LeeLiebermanArgonGlobalModel(_provider())
    with pytest.raises(ValueError, match="must be positive"):
        model.residuals(
            electron_temperature_eV=3.0,
            electron_density_m3=-1.0,
            metastable_density_m3=1.0e15,
            condition=_condition(),
        )


def test_solver_fails_closed_when_residual_gate_is_impossible():
    model = LeeLiebermanArgonGlobalModel(_provider())
    with pytest.raises(RuntimeError, match="failed conservation gate"):
        model.solve(
            _condition(),
            residual_tolerance=1.0e-30,
            maximum_evaluations=100,
        )
