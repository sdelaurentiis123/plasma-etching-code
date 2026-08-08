import csv
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import least_squares

from petch.reactor_global import (
    CHLORINE_ATOM_MASS_AMU,
    ChlorineIncidentVelocityState,
    ChlorineWallRecombinationBoundary,
    chlorine_atom_mean_thermal_speed_m_s,
    stafford_2010_bounded_hill_wall_recombination_provider,
    stafford_2010_conditioned_wall_recombination_provider,
    thermalized_chlorine_incident_velocity_state,
)
from petch.reactor_global.chlorine_wall import BOLTZMANN_J_K
from petch.reactor_global.transport import ATOMIC_MASS_UNIT_KG

ROOT = Path(__file__).resolve().parents[1]
STAFFORD_FIGURE8 = (
    ROOT / "data" / "experimental" / "stafford_2010"
    / "figure8_chlorine_wall_recombination.csv"
)


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


@pytest.mark.parametrize(
    "material", ("anodized_aluminum", "stainless_steel"))
def test_stafford_wall_fit_is_exactly_reproducible_from_digitized_markers(
        material):
    with STAFFORD_FIGURE8.open(newline="", encoding="utf-8") as stream:
        rows = [
            row for row in csv.DictReader(stream)
            if row["material"] == material
        ]
    ratio = np.asarray([
        float(row["cl_to_cl2_density_ratio"]) for row in rows])
    log_probability = np.log10([
        float(row["cl_recombination_probability"]) for row in rows])
    slope, intercept = np.polyfit(ratio, log_probability, 1)
    residual = log_probability - (intercept + slope * ratio)
    leave_one_out_residual = []
    for index in range(len(rows)):
        loo_slope, loo_intercept = np.polyfit(
            np.delete(ratio, index),
            np.delete(log_probability, index),
            1,
        )
        leave_one_out_residual.append(
            log_probability[index]
            - (loo_intercept + loo_slope * ratio[index])
        )
    leave_one_out_residual = np.asarray(leave_one_out_residual)
    provider = stafford_2010_conditioned_wall_recombination_provider(
        material)

    assert provider.marker_count == len(rows)
    assert provider.slope_per_ratio == pytest.approx(slope, rel=1.0e-15)
    assert provider.intercept_log10 == pytest.approx(
        intercept, rel=1.0e-15)
    assert provider.fit_rmse_log10 == pytest.approx(
        np.sqrt(np.mean(residual ** 2)), rel=1.0e-15)
    assert provider.fit_maximum_absolute_residual_log10 == pytest.approx(
        np.max(np.abs(residual)), rel=1.0e-15)
    assert provider.leave_one_out_rmse_log10 == pytest.approx(
        np.sqrt(np.mean(leave_one_out_residual ** 2)), rel=1.0e-15)
    assert (
        provider.leave_one_out_maximum_absolute_residual_log10
        == pytest.approx(
            np.max(np.abs(leave_one_out_residual)), rel=1.0e-15)
    )


def test_stafford_wall_provider_is_ratio_dependent_and_fail_closed():
    provider = stafford_2010_conditioned_wall_recombination_provider(
        "anodized_aluminum")
    low = provider.predict(
        cl_to_cl2_ratio=0.2,
        pressure_Pa=5.0 * 0.1333223684,
        icp_power_W=300.0,
        gas_temperature_K=300.0,
    )
    high = provider.predict(
        cl_to_cl2_ratio=0.7,
        pressure_Pa=5.0 * 0.1333223684,
        icp_power_W=300.0,
        gas_temperature_K=300.0,
    )
    assert high.recombination_probability > low.recombination_probability
    assert high.evidence_kind == "regressed"
    assert high.relative_measurement_uncertainty is None
    assert not high.supports_local_prediction
    assert not provider.supports_prediction
    assert high.provenance["coefficient_selection_target"].endswith(
        "no reactor, feature, or depth observable")

    with pytest.raises(ValueError, match="Cl/Cl2 ratio"):
        provider.predict(
            cl_to_cl2_ratio=0.05,
            pressure_Pa=5.0 * 0.1333223684,
            icp_power_W=300.0,
            gas_temperature_K=300.0,
        )


@pytest.mark.parametrize(
    "material", ("anodized_aluminum", "stainless_steel"))
def test_bounded_hill_wall_fit_replays_direct_markers(material):
    with STAFFORD_FIGURE8.open(newline="", encoding="utf-8") as stream:
        rows = [
            row for row in csv.DictReader(stream)
            if row["material"] == material
        ]
    ratio = np.asarray([
        float(row["cl_to_cl2_density_ratio"]) for row in rows])
    probability = np.asarray([
        float(row["cl_recombination_probability"]) for row in rows])
    asymptote = float(np.max(probability))

    def residual(log_parameters):
        half_saturation, hill_exponent = np.exp(log_parameters)
        modeled = asymptote / (
            1.0 + (half_saturation / ratio) ** hill_exponent)
        return np.log10(modeled) - np.log10(probability)

    fit = least_squares(
        residual,
        np.log((0.7, 1.3)),
        xtol=1.0e-14,
        ftol=1.0e-14,
        gtol=1.0e-14,
        max_nfev=100_000,
    )
    provider = stafford_2010_bounded_hill_wall_recombination_provider(
        material)
    assert provider.asymptotic_probability == asymptote
    assert provider.half_saturation_ratio == pytest.approx(
        np.exp(fit.x[0]), rel=1.0e-8)
    assert provider.hill_exponent == pytest.approx(
        np.exp(fit.x[1]), rel=1.0e-8)
    assert provider.fit_rmse_log10 == pytest.approx(
        np.sqrt(np.mean(residual(fit.x) ** 2)), rel=1.0e-12)
    assert provider.recombination_probability(1.0e6) < asymptote
    assert provider.recombination_probability(1.0e6) == pytest.approx(
        asymptote, rel=2.0e-7)
    assert not provider.extends_direct_evidence


def test_bounded_hill_wall_transfer_is_bounded_and_sensitivity_only():
    provider = stafford_2010_bounded_hill_wall_recombination_provider(
        "anodized_aluminum",
        valid_cl_to_cl2_ratio=(1.0e-5, 30.0),
        valid_gas_temperature_K=(300.0, 333.0),
        transfer_source=(
            "declared Malyshev anodized-wall ratio/temperature sensitivity"),
    )
    boundary = provider.predict(
        cl_to_cl2_ratio=3.0,
        pressure_Pa=2.0 * 0.1333223684,
        icp_power_W=500.0,
        gas_temperature_K=333.0,
    )
    assert 0.0 < boundary.recombination_probability < (
        provider.asymptotic_probability)
    assert boundary.evidence_kind == "sensitivity"
    assert boundary.provenance["runtime_extension_is_sensitivity"] is True
    assert boundary.provenance["coefficient_selection_target"].endswith(
        "no reactor, feature, or depth observable")
    assert not provider.supports_prediction
    with pytest.raises(ValueError, match="transfer source"):
        stafford_2010_bounded_hill_wall_recombination_provider(
            "anodized_aluminum",
            valid_cl_to_cl2_ratio=(1.0e-5, 30.0),
        )
