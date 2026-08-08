import pytest

from petch.reactor_global import (
    ChlorineWallRecombinationBoundary,
    CylindricalReactor,
    thermalized_chlorine_incident_velocity_state,
)
from petch.reactor_global.model import BOLTZMANN_J_K
from petch.reactor_global.chlorine_transport import (
    ECONOMOU_CHLORINE_REDUCED_DIFFUSIVITY_M_INV_S,
    lymberopoulos_economou_1995_chlorine_diffusivity,
    malyshev_1998_chlorine_in_chlorine_diffusivity,
    ramamurthi_economou_2002_chlorine_diffusivity,
    solve_chlorine_neutral_wall_transport,
)


def _wall_boundary():
    return ChlorineWallRecombinationBoundary(
        recombination_probability=0.02,
        surface_state="conditioned chlorine reactor wall",
        source="hypothetical scoped wall measurement",
        evidence_kind="measured",
        valid_cl_to_cl2_ratio=(0.1, 0.8),
        valid_pressure_Pa=(1.0, 2.0),
        valid_icp_power_W=(100.0, 600.0),
        valid_gas_temperature_K=(500.0, 500.0),
        relative_measurement_uncertainty=0.2,
    )


def test_published_chlorine_diffusivities_retain_temperature_conflict():
    original = lymberopoulos_economou_1995_chlorine_diffusivity()
    reused = ramamurthi_economou_2002_chlorine_diffusivity()
    assert (
        original.reduced_diffusivity_m_inv_s
        == reused.reduced_diffusivity_m_inv_s
        == ECONOMOU_CHLORINE_REDUCED_DIFFUSIVITY_M_INV_S
    )
    assert original.reference_temperature_K == 500.0
    assert reused.reference_temperature_K == 300.0
    assert "temperature_conflict" in original.provenance
    assert "temperature_conflict" in reused.provenance
    assert not original.supports_prediction
    assert not reused.supports_prediction


def test_malyshev_chapman_enskog_diffusivity_replays_source_correction():
    model = malyshev_1998_chlorine_in_chlorine_diffusivity()
    assert model.mixed_sigma_angstrom == pytest.approx(3.8315)
    assert model.mixed_epsilon_over_k_K == pytest.approx(
        163.63068171953572)
    assert model.collision_integral(333.0) == pytest.approx(
        1.0689196299076194)
    assert model.diffusivity_cm2_s_at_pressure(
        gas_temperature_K=298.15,
        pressure_Pa=101325.0,
    ) == pytest.approx(0.15167523928085008)
    assert model.diffusivity_cm2_s_at_pressure(
        gas_temperature_K=333.0,
        pressure_Pa=101325.0,
    ) == pytest.approx(0.1862154115776554)
    assert model.source_correction_factor == 1.25
    assert model.provenance["source_bibkey"] == "malyshev-1998-lam-cl2"
    assert model.provenance["collision_integral_bibkey"] == (
        "neufeld-1972-collision-integrals")
    assert model.provenance["coefficient_selection_target"] is None
    assert model.provenance["feature_depth_target"] is None
    assert not model.supports_prediction


def test_malyshev_diffusivity_is_temperature_consistent_and_composable():
    pressure = 10.0 * 0.1333223684
    temperature = 333.0
    density = pressure / (BOLTZMANN_J_K * temperature)
    model = malyshev_1998_chlorine_in_chlorine_diffusivity()
    state = model.evaluate(
        total_neutral_density_m3=density,
        gas_temperature_K=temperature,
    )
    assert state.diffusivity_m2_s * density == pytest.approx(
        4.103975103414651e20)
    assert state.provenance["omega_11"] == pytest.approx(
        1.0689196299076194)
    assert state.diffusivity_m2_s * density < (
        ECONOMOU_CHLORINE_REDUCED_DIFFUSIVITY_M_INV_S)
    assert not state.supports_prediction

    boundary = ChlorineWallRecombinationBoundary(
        recombination_probability=0.035,
        surface_state="source-model anodized aluminum",
        source="Malyshev 1998 fitted reactor value",
        evidence_kind="regressed",
        valid_cl_to_cl2_ratio=(0.1, 1.0),
        valid_pressure_Pa=(pressure, pressure),
        valid_icp_power_W=(500.0, 500.0),
        valid_gas_temperature_K=(temperature, temperature),
        relative_measurement_uncertainty=None,
    )
    transport = solve_chlorine_neutral_wall_transport(
        geometry=CylindricalReactor(radius_m=0.215, length_m=0.11),
        wall_boundary=boundary,
        incident_velocity_state=(
            thermalized_chlorine_incident_velocity_state(
                temperature,
                source="source-model initial gas temperature",
                evidence_kind="assumed",
                relative_uncertainty=None,
            )
        ),
        diffusivity_model=model,
        total_neutral_density_m3=density,
        gas_temperature_K=temperature,
        cl_to_cl2_ratio=0.5,
        pressure_Pa=pressure,
        icp_power_W=500.0,
    )
    assert transport.wall_loss.numerical_closure_passes
    assert transport.diffusivity == state
    assert not transport.supports_prediction


def test_composed_chlorine_wall_transport_is_conserved_but_not_predictive():
    pressure = 1.333223684
    temperature = 500.0
    total_density = pressure / (BOLTZMANN_J_K * temperature)
    transport = solve_chlorine_neutral_wall_transport(
        geometry=CylindricalReactor(radius_m=0.14, length_m=0.10),
        wall_boundary=_wall_boundary(),
        incident_velocity_state=(
            thermalized_chlorine_incident_velocity_state(
                temperature,
                source="published-model thermal population",
                evidence_kind="assumed",
                relative_uncertainty=None,
            )
        ),
        diffusivity_model=(
            lymberopoulos_economou_1995_chlorine_diffusivity()),
        total_neutral_density_m3=total_density,
        gas_temperature_K=temperature,
        cl_to_cl2_ratio=0.3,
        pressure_Pa=pressure,
        icp_power_W=300.0,
    )
    assert transport.diffusivity.diffusivity_m2_s == pytest.approx(
        ECONOMOU_CHLORINE_REDUCED_DIFFUSIVITY_M_INV_S
        / total_density
    )
    assert transport.geometry == CylindricalReactor(
        radius_m=0.14, length_m=0.10)
    assert transport.wall_loss.numerical_closure_passes
    assert not transport.supports_prediction
    assert not transport.incident_velocity_state.supports_prediction

    rates = transport.evaluate_volume_rates(2.0e20)
    assert rates.chlorine_atom_loss_m3_s > 0.0
    assert rates.chlorine_molecule_return_m3_s == pytest.approx(
        0.5 * rates.chlorine_atom_loss_m3_s)
    assert rates.chlorine_atom_inventory_residual_m3_s == 0.0


def test_composed_transport_refuses_cross_temperature_coefficient_use():
    with pytest.raises(ValueError, match="temperature"):
        solve_chlorine_neutral_wall_transport(
            geometry=CylindricalReactor(radius_m=0.14, length_m=0.10),
            wall_boundary=_wall_boundary(),
            incident_velocity_state=(
                thermalized_chlorine_incident_velocity_state(
                    500.0,
                    source="published-model thermal population",
                    evidence_kind="assumed",
                    relative_uncertainty=None,
                )
            ),
            diffusivity_model=(
                ramamurthi_economou_2002_chlorine_diffusivity()),
            total_neutral_density_m3=1.0e21,
            gas_temperature_K=500.0,
            cl_to_cl2_ratio=0.3,
            pressure_Pa=1.333223684,
            icp_power_W=300.0,
        )


def test_transport_refuses_mismatched_maxwellian_temperature():
    with pytest.raises(ValueError, match="reference temperature"):
        solve_chlorine_neutral_wall_transport(
            geometry=CylindricalReactor(radius_m=0.14, length_m=0.10),
            wall_boundary=_wall_boundary(),
            incident_velocity_state=(
                thermalized_chlorine_incident_velocity_state(
                    499.0,
                    source="temperature mismatch test",
                    evidence_kind="assumed",
                    relative_uncertainty=None,
                )
            ),
            diffusivity_model=(
                lymberopoulos_economou_1995_chlorine_diffusivity()),
            total_neutral_density_m3=1.0e21,
            gas_temperature_K=500.0,
            cl_to_cl2_ratio=0.3,
            pressure_Pa=1.333223684,
            icp_power_W=300.0,
        )
