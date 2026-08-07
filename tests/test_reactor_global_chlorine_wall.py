import numpy as np
import pytest

from petch.reactor_global import (
    CHLORINE_ATOM_MASS_AMU,
    ChlorineIncidentVelocityState,
    ChlorineWallRecombinationBoundary,
    chlorine_atom_mean_thermal_speed_m_s,
    thermalized_chlorine_incident_velocity_state,
)
from petch.reactor_global.chlorine_wall import BOLTZMANN_J_K
from petch.reactor_global.transport import ATOMIC_MASS_UNIT_KG


def _boundary(**changes):
    arguments = {
        "recombination_probability": 0.02,
        "surface_state": "conditioned silica-rich stainless steel",
        "source": "stafford-2010-cl-wall Figure 8 measured envelope member",
        "evidence_kind": "published_range_member",
        "valid_cl_to_cl2_ratio": (0.1, 0.8),
        "valid_pressure_Pa": (1.25 * 0.1333223684, 20 * 0.1333223684),
        "valid_icp_power_W": (100.0, 600.0),
        "valid_gas_temperature_K": (280.0, 350.0),
        "relative_measurement_uncertainty": None,
    }
    arguments.update(changes)
    return ChlorineWallRecombinationBoundary(**arguments)


def test_chlorine_mean_speed_is_exact_maxwellian_kinetic_theory():
    temperature = 300.0
    expected = np.sqrt(
        8.0 * BOLTZMANN_J_K * temperature
        / (
            np.pi
            * CHLORINE_ATOM_MASS_AMU
            * ATOMIC_MASS_UNIT_KG
        )
    )
    assert chlorine_atom_mean_thermal_speed_m_s(
        temperature) == pytest.approx(expected, rel=1.0e-15)


def test_local_wall_flux_uses_quarter_n_vbar_and_conserves_chlorine():
    boundary = _boundary()
    density = 2.0e20
    flux = boundary.evaluate(
        chlorine_atom_density_m3=density,
        incident_velocity_state=(
            thermalized_chlorine_incident_velocity_state(
                350.0,
                source="assumed thermalized chlorine population",
                evidence_kind="assumed",
                relative_uncertainty=None,
            )
        ),
        gas_temperature_K=350.0,
        cl_to_cl2_ratio=0.35,
        pressure_Pa=5.0 * 0.1333223684,
        icp_power_W=300.0,
    )
    assert flux.incident_cl_atom_flux_m2_s == pytest.approx(
        0.25 * density * flux.mean_cl_speed_m_s)
    assert flux.recombined_cl_atom_flux_m2_s == pytest.approx(
        0.02 * flux.incident_cl_atom_flux_m2_s)
    assert flux.returned_cl2_molecule_flux_m2_s == pytest.approx(
        0.5 * flux.recombined_cl_atom_flux_m2_s)
    assert flux.chlorine_atom_inventory_residual_m2_s == 0.0
    assert not flux.reactor_volume_closure_ready
    assert not flux.supports_local_prediction


def test_measured_nonthermal_isotropic_speed_is_used_exactly():
    boundary = _boundary(
        evidence_kind="measured",
        relative_measurement_uncertainty=0.2,
    )
    velocity = ChlorineIncidentVelocityState(
        mean_speed_m_s=1500.0,
        distribution_kind="measured_isotropic",
        source="hypothetical velocity-distribution measurement",
        evidence_kind="measured",
        relative_uncertainty=0.1,
    )
    flux = boundary.evaluate(
        chlorine_atom_density_m3=2.0e20,
        incident_velocity_state=velocity,
        gas_temperature_K=300.0,
        cl_to_cl2_ratio=0.35,
        pressure_Pa=5.0 * 0.1333223684,
        icp_power_W=300.0,
    )
    assert flux.incident_cl_atom_flux_m2_s == pytest.approx(7.5e22)
    assert flux.mean_cl_speed_m_s == 1500.0
    assert flux.supports_local_prediction


def test_thermalized_state_refuses_a_different_reactor_temperature():
    velocity = thermalized_chlorine_incident_velocity_state(
        300.0,
        source="thermalization sensitivity",
        evidence_kind="sensitivity",
        relative_uncertainty=None,
    )
    with pytest.raises(ValueError, match="reference temperature"):
        _boundary().evaluate(
            chlorine_atom_density_m3=1.0e20,
            incident_velocity_state=velocity,
            gas_temperature_K=301.0,
            cl_to_cl2_ratio=0.35,
            pressure_Pa=5.0 * 0.1333223684,
            icp_power_W=300.0,
        )


@pytest.mark.parametrize(
    ("field", "value"),
    (
        ("recombination_probability", -0.01),
        ("recombination_probability", 1.01),
        ("surface_state", ""),
        ("evidence_kind", "fitted_to_depth"),
        ("valid_cl_to_cl2_ratio", (0.8, 0.1)),
        ("valid_pressure_Pa", (0.0, 1.0)),
        ("valid_icp_power_W", (600.0, 100.0)),
        ("valid_gas_temperature_K", (350.0, 280.0)),
    ),
)
def test_wall_boundary_rejects_unphysical_or_unscoped_inputs(field, value):
    with pytest.raises(ValueError):
        _boundary(**{field: value})


@pytest.mark.parametrize(
    "condition",
    (
        {
            "cl_to_cl2_ratio": 0.05,
            "pressure_Pa": 5.0 * 0.1333223684,
            "icp_power_W": 300.0,
            "gas_temperature_K": 300.0,
        },
        {
            "cl_to_cl2_ratio": 0.35,
            "pressure_Pa": 30.0 * 0.1333223684,
            "icp_power_W": 300.0,
            "gas_temperature_K": 300.0,
        },
        {
            "cl_to_cl2_ratio": 0.35,
            "pressure_Pa": 5.0 * 0.1333223684,
            "icp_power_W": 700.0,
            "gas_temperature_K": 300.0,
        },
        {
            "cl_to_cl2_ratio": 0.35,
            "pressure_Pa": 5.0 * 0.1333223684,
            "icp_power_W": 300.0,
            "gas_temperature_K": 400.0,
        },
    ),
)
def test_wall_boundary_fails_outside_declared_evidence_domain(condition):
    with pytest.raises(ValueError, match="outside"):
        _boundary().evaluate(
            chlorine_atom_density_m3=1.0e20,
            incident_velocity_state=(
                thermalized_chlorine_incident_velocity_state(
                    condition["gas_temperature_K"],
                    source="assumed thermalized chlorine population",
                    evidence_kind="assumed",
                    relative_uncertainty=None,
                )
            ),
            **condition,
        )


def test_local_prediction_requires_predictive_evidence_and_uncertainty():
    assert not _boundary().supports_local_prediction
    assert not _boundary(
        evidence_kind="measured",
        relative_measurement_uncertainty=None,
    ).supports_local_prediction
    assert _boundary(
        evidence_kind="measured",
        relative_measurement_uncertainty=0.2,
    ).supports_local_prediction
