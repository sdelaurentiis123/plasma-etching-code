import numpy as np
import pytest

from petch.reactor_global import (
    ChlorineChargedTransportState,
    ChlorineFixedPressureCondition,
    ChlorineWallRecombinationBoundary,
    CylindricalReactor,
    FixedElectronTemperatureChlorineParticleModel,
    FixedChlorineChargedTransportProvider,
    FixedChlorineNeutralWallTransportProvider,
    LeeEconomouChlorineChargedTransportProvider,
    PositiveIonWallTransport,
    ReactorScalarInput,
    build_hamilton_dissociation_chlorine_particle_network,
    build_lee_lieberman_chlorine_particle_network,
    lymberopoulos_economou_1995_chlorine_diffusivity,
    lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities,
    solve_chlorine_neutral_wall_transport,
    standard_volume_flow_molecules_s,
    thermalized_chlorine_incident_velocity_state,
)
from petch.reactor_global.chlorine_particle_model import BOLTZMANN_J_K


def _scalar(value, unit, *, evidence_kind="assumed", uncertainty=None):
    return ReactorScalarInput(
        value=value,
        unit=unit,
        source="test input independent of solver output",
        evidence_kind=evidence_kind,
        relative_uncertainty=uncertainty,
    )


def _condition(geometry):
    flow = standard_volume_flow_molecules_s(
        35.0,
        standard_temperature_K=273.15,
        standard_pressure_Pa=101325.0,
    )
    return ChlorineFixedPressureCondition(
        condition_id="lee-chlorine-particle-reproduction",
        geometry=geometry,
        pressure=_scalar(1.333223684, "Pa"),
        gas_temperature=_scalar(500.0, "K"),
        electron_temperature=_scalar(3.0, "eV"),
        chlorine_molecule_feed=_scalar(flow, "molecule s^-1"),
        source_power=_scalar(1000.0, "W"),
    )


def _wall_boundary(*, ratio_domain=(0.001, 100.0)):
    return ChlorineWallRecombinationBoundary(
        recombination_probability=0.02,
        surface_state="conditioned chlorine reactor wall",
        source="hypothetical direct wall measurement",
        evidence_kind="measured",
        valid_cl_to_cl2_ratio=ratio_domain,
        valid_pressure_Pa=(1.0, 2.0),
        valid_icp_power_W=(100.0, 1200.0),
        valid_gas_temperature_K=(500.0, 500.0),
        relative_measurement_uncertainty=0.2,
    )


def _neutral_transport(condition, *, ratio_domain=(0.001, 100.0)):
    return solve_chlorine_neutral_wall_transport(
        geometry=condition.geometry,
        wall_boundary=_wall_boundary(ratio_domain=ratio_domain),
        incident_velocity_state=(
            thermalized_chlorine_incident_velocity_state(
                condition.gas_temperature.value,
                source="assumed thermalized source-model chlorine",
                evidence_kind="assumed",
                relative_uncertainty=None,
            )
        ),
        diffusivity_model=(
            lymberopoulos_economou_1995_chlorine_diffusivity()),
        total_neutral_density_m3=condition.target_neutral_density_m3,
        gas_temperature_K=condition.gas_temperature.value,
        cl_to_cl2_ratio=1.0,
        pressure_Pa=condition.pressure.value,
        icp_power_W=condition.source_power.value,
    )


def _neutral_provider(condition, *, ratio_domain=(0.001, 100.0)):
    return FixedChlorineNeutralWallTransportProvider(
        _neutral_transport(condition, ratio_domain=ratio_domain))


def _charged_transport(geometry):
    return ChlorineChargedTransportState(
        geometry=geometry,
        positive_ion_transport={
            "Cl2+": PositiveIonWallTransport(
                axial_flux_velocity_m_s=1000.0,
                radial_flux_velocity_m_s=500.0,
                source="declared fixed source-model transport",
                evidence_kind="assumed",
            ),
            "Cl+": PositiveIonWallTransport(
                axial_flux_velocity_m_s=1400.0,
                radial_flux_velocity_m_s=700.0,
                source="declared fixed source-model transport",
                evidence_kind="assumed",
            ),
        },
        negative_ion_confinement_source=(
            "lee-lieberman-1994 parabolic negative-ion profile with "
            "zero sheath-edge density"
        ),
        negative_ion_confinement_evidence="published_model",
    )


def test_sccm_conversion_requires_and_uses_explicit_standard_state():
    flow = standard_volume_flow_molecules_s(
        1.0,
        standard_temperature_K=273.15,
        standard_pressure_Pa=101325.0,
    )
    expected = (
        101325.0 / (BOLTZMANN_J_K * 273.15) * 1.0e-6 / 60.0)
    assert flow == pytest.approx(expected, rel=1.0e-15)
    assert flow == pytest.approx(4.47796e17, rel=2.0e-6)
    with pytest.raises(ValueError):
        standard_volume_flow_molecules_s(
            1.0,
            standard_temperature_K=0.0,
            standard_pressure_Pa=101325.0,
        )


def test_fixed_pressure_chlorine_particle_solver_closes_every_ledger():
    geometry = CylindricalReactor(radius_m=0.1525, length_m=0.075)
    condition = _condition(geometry)
    solution = FixedElectronTemperatureChlorineParticleModel(
        build_lee_lieberman_chlorine_particle_network()
    ).solve(
        condition,
        charged_transport_provider=FixedChlorineChargedTransportProvider(
            _charged_transport(geometry)),
        neutral_wall_transport_provider=_neutral_provider(condition),
    )

    assert solution.maximum_normalized_residual < 1.0e-12
    assert solution.supports_particle_reproduction
    assert not solution.supports_prediction
    assert "electron_power_balance" in solution.missing_prediction_closures
    assert solution.chlorine_atom_dissociation_fraction == pytest.approx(
        0.6125981926, rel=2.0e-8)
    assert (
        solution.densities_m3["Cl2"] + solution.densities_m3["Cl"]
        == pytest.approx(condition.target_neutral_density_m3, rel=2.0e-13)
    )
    feed_atom_rate_density = (
        2.0 * condition.chlorine_molecule_feed.value / geometry.volume_m3)
    assert abs(solution.chlorine_atom_inventory_residual_m3_s) <= (
        2.0e-12 * feed_atom_rate_density)
    electron_wall_loss = sum(
        solution.positive_ion_wall_loss_m3_s.values())
    assert abs(solution.electron_current_balance_residual_m3_s) <= (
        2.0e-12 * electron_wall_loss)
    assert solution.axial_positive_ion_flux_m2_s["Cl2+"] == pytest.approx(
        1000.0 * solution.densities_m3["Cl2+"])
    assert solution.axial_positive_ion_flux_m2_s["Cl+"] == pytest.approx(
        1400.0 * solution.densities_m3["Cl+"])


def test_pressure_controller_exhaust_is_solved_after_dissociation():
    geometry = CylindricalReactor(radius_m=0.1525, length_m=0.075)
    condition = _condition(geometry)
    solution = FixedElectronTemperatureChlorineParticleModel(
        build_lee_lieberman_chlorine_particle_network()
    ).solve(
        condition,
        charged_transport_provider=FixedChlorineChargedTransportProvider(
            _charged_transport(geometry)),
        neutral_wall_transport_provider=_neutral_provider(condition),
    )
    naive_feed_over_neutral_inventory = (
        condition.chlorine_molecule_feed.value
        / (condition.target_neutral_density_m3 * geometry.volume_m3)
    )
    assert solution.exhaust_loss_frequency_s_inv > (
        1.5 * naive_feed_over_neutral_inventory)


def test_solution_refuses_wall_probability_outside_final_ratio_domain():
    geometry = CylindricalReactor(radius_m=0.1525, length_m=0.075)
    condition = _condition(geometry)
    with pytest.raises(ValueError, match="Cl/Cl2 ratio"):
        FixedElectronTemperatureChlorineParticleModel(
            build_lee_lieberman_chlorine_particle_network()
        ).solve(
            condition,
            charged_transport_provider=FixedChlorineChargedTransportProvider(
                _charged_transport(geometry)),
            neutral_wall_transport_provider=_neutral_provider(
                condition, ratio_domain=(0.9, 1.1)),
        )


def test_solver_rejects_transport_from_a_different_geometry():
    geometry = CylindricalReactor(radius_m=0.1525, length_m=0.075)
    condition = _condition(geometry)
    other_geometry = CylindricalReactor(radius_m=0.14, length_m=0.075)
    with pytest.raises(ValueError, match="charged transport geometry"):
        FixedElectronTemperatureChlorineParticleModel(
            build_lee_lieberman_chlorine_particle_network()
        ).solve(
            condition,
            charged_transport_provider=FixedChlorineChargedTransportProvider(
                _charged_transport(other_geometry)),
            neutral_wall_transport_provider=_neutral_provider(condition),
        )


def test_charged_transport_requires_both_positive_ion_species():
    with pytest.raises(ValueError, match="charged-transport"):
        ChlorineChargedTransportState(
            geometry=CylindricalReactor(radius_m=0.15, length_m=0.075),
            positive_ion_transport={
                "Cl+": PositiveIonWallTransport(
                    axial_flux_velocity_m_s=1000.0,
                    radial_flux_velocity_m_s=500.0,
                    source="incomplete test state",
                    evidence_kind="assumed",
                )
            },
            negative_ion_confinement_source="test",
            negative_ion_confinement_evidence="assumed",
        )


def test_no_untracked_units_enter_fixed_pressure_condition():
    geometry = CylindricalReactor(radius_m=0.15, length_m=0.075)
    condition = _condition(geometry)
    with pytest.raises(ValueError, match="explicit unit Pa"):
        ChlorineFixedPressureCondition(
            condition_id="wrong-pressure-unit",
            geometry=geometry,
            pressure=_scalar(condition.pressure.value, "mTorr"),
            gas_temperature=condition.gas_temperature,
            electron_temperature=condition.electron_temperature,
            chlorine_molecule_feed=condition.chlorine_molecule_feed,
            source_power=condition.source_power,
        )


def test_state_dependent_charged_transport_closes_inside_particle_solve():
    geometry = CylindricalReactor(radius_m=0.1525, length_m=0.075)
    condition = _condition(geometry)
    mobilities = (
        lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities())
    provider = LeeEconomouChlorineChargedTransportProvider(
        reduced_mobilities={
            species: mobilities[species] for species in ("Cl2+", "Cl+")
        },
        ion_temperature=_scalar(0.12, "eV"),
    )
    solution = FixedElectronTemperatureChlorineParticleModel(
        build_lee_lieberman_chlorine_particle_network()
    ).solve(
        condition,
        charged_transport_provider=provider,
        neutral_wall_transport_provider=_neutral_provider(condition),
    )

    assert solution.maximum_normalized_residual < 1.0e-11
    assert solution.supports_particle_reproduction
    assert not solution.charged_transport_supports_prediction
    assert solution.total_axial_positive_ion_flux_m2_s > 0.0
    assert abs(solution.electron_current_balance_residual_m3_s) < (
        2.0e-11 * sum(solution.positive_ion_wall_loss_m3_s.values())
    )


def test_hamilton_dissociation_deck_closes_particle_solve_without_depth_fit():
    geometry = CylindricalReactor(radius_m=0.1525, length_m=0.075)
    condition = _condition(geometry)
    mobilities = (
        lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities())
    provider = LeeEconomouChlorineChargedTransportProvider(
        reduced_mobilities={
            species: mobilities[species] for species in ("Cl2+", "Cl+")
        },
        ion_temperature=_scalar(0.12, "eV"),
    )
    neutral_transport_provider = _neutral_provider(condition)
    upgraded = FixedElectronTemperatureChlorineParticleModel(
        build_hamilton_dissociation_chlorine_particle_network()
    ).solve(
        condition,
        charged_transport_provider=provider,
        neutral_wall_transport_provider=neutral_transport_provider,
    )
    legacy = FixedElectronTemperatureChlorineParticleModel(
        build_lee_lieberman_chlorine_particle_network()
    ).solve(
        condition,
        charged_transport_provider=provider,
        neutral_wall_transport_provider=neutral_transport_provider,
    )

    assert upgraded.maximum_normalized_residual < 1.0e-11
    assert upgraded.supports_particle_reproduction
    assert not upgraded.supports_prediction
    assert upgraded.chlorine_atom_dissociation_fraction < (
        legacy.chlorine_atom_dissociation_fraction)
