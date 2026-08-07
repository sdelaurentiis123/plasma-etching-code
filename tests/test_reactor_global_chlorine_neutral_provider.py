import pytest

from petch.reactor_global import (
    ChlorineFixedPressureCondition,
    CylindricalReactor,
    FixedElectronTemperatureChlorineParticleModel,
    LeeEconomouChlorineChargedTransportProvider,
    LogLinearChlorineWallRecombinationProvider,
    ReactorScalarInput,
    StateDependentChlorineNeutralTransportProvider,
    build_hamilton_dissociation_chlorine_particle_network,
    lymberopoulos_economou_1995_chlorine_diffusivity,
    lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities,
    ramamurthi_economou_2002_chlorine_diffusivity,
    stafford_2010_conditioned_wall_recombination_provider,
    standard_volume_flow_molecules_s,
    thermalized_chlorine_incident_velocity_state,
)


def _scalar(value, unit, *, evidence_kind="published_model"):
    return ReactorScalarInput(
        value=value,
        unit=unit,
        source="independent neutral-transport regression input",
        evidence_kind=evidence_kind,
        relative_uncertainty=None,
    )


def _condition(*, temperature_K, pressure_Pa, source_power_W):
    geometry = CylindricalReactor(radius_m=0.1525, length_m=0.075)
    return ChlorineFixedPressureCondition(
        condition_id="state-dependent-neutral-source-reproduction",
        geometry=geometry,
        neutral_control_volume=_scalar(geometry.volume_m3, "m3"),
        pressure=_scalar(pressure_Pa, "Pa"),
        gas_temperature=_scalar(temperature_K, "K"),
        electron_temperature=_scalar(3.0, "eV"),
        chlorine_molecule_feed=_scalar(
            standard_volume_flow_molecules_s(
                35.0,
                standard_temperature_K=273.15,
                standard_pressure_Pa=101325.0,
            ),
            "molecule s^-1",
        ),
        source_power=_scalar(source_power_W, "W"),
    )


def _densities(*, ratio):
    chlorine_molecule_density = 1.0e20
    return {
        "e": 1.0e16,
        "Cl2": chlorine_molecule_density,
        "Cl": ratio * chlorine_molecule_density,
        "Cl2+": 1.2e16,
        "Cl+": 0.8e16,
        "Cl-": 1.0e16,
    }


def test_stafford_wall_law_recomputes_exact_transport_from_density_ratio():
    condition = _condition(
        temperature_K=300.0,
        pressure_Pa=5.0 * 0.1333223684,
        source_power_W=300.0,
    )
    provider = StateDependentChlorineNeutralTransportProvider(
        wall_recombination_provider=(
            stafford_2010_conditioned_wall_recombination_provider(
                "anodized_aluminum")),
        incident_velocity_state=(
            thermalized_chlorine_incident_velocity_state(
                300.0,
                source="Stafford wall-flux 300 K assumption",
                evidence_kind="assumed",
                relative_uncertainty=None,
            )
        ),
        diffusivity_model=ramamurthi_economou_2002_chlorine_diffusivity(),
    )
    low = provider.predict(condition, _densities(ratio=0.2))
    high = provider.predict(condition, _densities(ratio=0.7))

    assert high.wall_boundary.recombination_probability > (
        low.wall_boundary.recombination_probability)
    assert high.wall_loss.exact_loss_frequency_s_inv > (
        low.wall_loss.exact_loss_frequency_s_inv)
    assert high.wall_loss.numerical_closure_passes
    assert not high.supports_prediction


def test_stafford_provider_refuses_lam_wall_temperature_transplant():
    condition = _condition(
        temperature_K=333.0,
        pressure_Pa=5.0 * 0.1333223684,
        source_power_W=300.0,
    )
    provider = StateDependentChlorineNeutralTransportProvider(
        wall_recombination_provider=(
            stafford_2010_conditioned_wall_recombination_provider(
                "anodized_aluminum")),
        incident_velocity_state=(
            thermalized_chlorine_incident_velocity_state(
                333.0,
                source="Lam wall-temperature sensitivity",
                evidence_kind="sensitivity",
                relative_uncertainty=None,
            )
        ),
        diffusivity_model=ramamurthi_economou_2002_chlorine_diffusivity(),
    )
    with pytest.raises(ValueError, match="gas temperature"):
        provider.predict(condition, _densities(ratio=0.4))


def test_dynamic_neutral_and_charged_transport_close_together_in_solver():
    condition = _condition(
        temperature_K=500.0,
        pressure_Pa=1.333223684,
        source_power_W=1000.0,
    )
    wall_provider = LogLinearChlorineWallRecombinationProvider(
        slope_per_ratio=0.02,
        intercept_log10=-1.7,
        surface_state="manufactured state-dependent wall regression",
        source="manufactured transport-coupling regression",
        valid_cl_to_cl2_ratio=(0.001, 10.0),
        valid_pressure_Pa=(1.0, 2.0),
        valid_icp_power_W=(100.0, 1200.0),
        valid_gas_temperature_K=(500.0, 500.0),
        marker_count=10,
        fit_rmse_log10=0.0,
        fit_maximum_absolute_residual_log10=0.0,
        leave_one_out_rmse_log10=0.0,
        leave_one_out_maximum_absolute_residual_log10=0.0,
    )
    neutral_provider = StateDependentChlorineNeutralTransportProvider(
        wall_recombination_provider=wall_provider,
        incident_velocity_state=(
            thermalized_chlorine_incident_velocity_state(
                500.0,
                source="manufactured thermalized velocity state",
                evidence_kind="assumed",
                relative_uncertainty=None,
            )
        ),
        diffusivity_model=lymberopoulos_economou_1995_chlorine_diffusivity(),
    )
    mobilities = (
        lymberopoulos_economou_1995_chlorine_reduced_ion_mobilities())
    charged_provider = LeeEconomouChlorineChargedTransportProvider(
        reduced_mobilities={
            species: mobilities[species] for species in ("Cl2+", "Cl+")
        },
        ion_temperature=_scalar(0.12, "eV"),
    )
    solution = FixedElectronTemperatureChlorineParticleModel(
        build_hamilton_dissociation_chlorine_particle_network()
    ).solve(
        condition,
        charged_transport_provider=charged_provider,
        neutral_wall_transport_provider=neutral_provider,
    )
    final_transport = neutral_provider.predict(condition, solution.densities_m3)
    expected_atom_loss = (
        final_transport.wall_loss.exact_loss_frequency_s_inv
        * solution.densities_m3["Cl"]
    )

    assert solution.maximum_normalized_residual < 1.0e-11
    assert solution.chlorine_wall_atom_loss_m3_s == pytest.approx(
        expected_atom_loss, rel=1.0e-14)
    assert not solution.neutral_transport_supports_prediction
    assert not solution.supports_prediction
