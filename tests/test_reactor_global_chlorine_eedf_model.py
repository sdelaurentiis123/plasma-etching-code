from dataclasses import replace

import pytest

from petch.reactor_global.chlorine import (
    build_lee_lieberman_chlorine_particle_network,
)
from petch.reactor_global.chlorine_eedf_model import (
    EEDFChlorineAbsorbedPowerModel,
    EEDFChlorineCondition,
    FixedPositiveIonWallEnergyProvider,
    PositiveIonWallEnergyState,
)
from petch.reactor_global.chlorine_particle_model import (
    ChlorineChargedTransportState,
    FixedChlorineChargedTransportProvider,
    FixedChlorineNeutralWallTransportProvider,
    PositiveIonWallTransport,
    ReactorScalarInput,
    standard_volume_flow_molecules_s,
)
from petch.reactor_global.chlorine_transport import (
    lymberopoulos_economou_1995_chlorine_diffusivity,
    solve_chlorine_neutral_wall_transport,
)
from petch.reactor_global.chlorine_wall import (
    ChlorineWallRecombinationBoundary,
    thermalized_chlorine_incident_velocity_state,
)
from petch.reactor_global.electron_collision_chemistry import (
    ElectronCollisionChemistry,
    ElectronCollisionHeavyMapping,
)
from petch.reactor_global.electron_collision_deck import (
    ElectronCollisionDeck,
    ElectronCollisionProcess,
)
from petch.reactor_global.electron_kinetics import (
    DeterministicTwoTermBoltzmannSolver,
    ElectronEnergyGrid,
)
from petch.reactor_global.geometry import CylindricalReactor
from petch.reactor_global.network import ReactionNetwork
from petch.reactor_global.power import AbsorbedPowerEstimate


def _scalar(value, unit):
    return ReactorScalarInput(
        value=value,
        unit=unit,
        source="manufactured global-model gate",
        evidence_kind="assumed",
    )


def _deck():
    common = {
        "target": "Cl2",
        "electron_energy_eV": (0.0, 100.0),
    }
    return ElectronCollisionDeck(
        processes=(
            ElectronCollisionProcess(
                kind="ELASTIC",
                product=None,
                cross_section_m2=(2.0e-20, 2.0e-20),
                mass_ratio=7.68e-6,
                **common,
            ),
            ElectronCollisionProcess(
                kind="EXCITATION",
                product="Cl2(v=1)",
                electron_energy_eV=(0.0, 2.0, 5.0, 100.0),
                cross_section_m2=(0.0, 0.0, 1.0e-20, 1.0e-20),
                energy_loss_eV=2.0,
                target="Cl2",
            ),
            ElectronCollisionProcess(
                kind="ATTACHMENT",
                product="Cl- + Cl",
                electron_energy_eV=(0.0, 0.2, 1.0, 100.0),
                cross_section_m2=(2.0e-20, 2.0e-20, 0.0, 0.0),
                energy_loss_eV=0.0,
                target="Cl2",
            ),
            ElectronCollisionProcess(
                kind="IONIZATION",
                product="Cl2+",
                electron_energy_eV=(0.0, 10.0, 15.0, 100.0),
                cross_section_m2=(0.0, 0.0, 1.4e-20, 1.4e-20),
                energy_loss_eV=10.0,
                target="Cl2",
            ),
            ElectronCollisionProcess(
                kind="IONIZATION",
                product="Cl+ + Cl",
                electron_energy_eV=(0.0, 12.0, 18.0, 100.0),
                cross_section_m2=(0.0, 0.0, 8.0e-21, 8.0e-21),
                energy_loss_eV=12.0,
                target="Cl2",
            ),
        ),
        payload_sha256="e" * 64,
        source_database="manufactured global model deck",
        retrieved_at="2026-08-08",
        source_reference="test only",
    )


def _chemistry(deck, species):
    common = {"source": "manufactured", "evidence_kind": "manufactured"}
    return ElectronCollisionChemistry(
        deck,
        species,
        (
            ElectronCollisionHeavyMapping(
                1, "excitation", {"Cl2": 1}, {"Cl2": 1}, **common),
            ElectronCollisionHeavyMapping(
                2, "attachment", {"Cl2": 1}, {"Cl-": 1, "Cl": 1},
                **common,
            ),
            ElectronCollisionHeavyMapping(
                3, "molecular_ionization", {"Cl2": 1}, {"Cl2+": 1},
                **common,
            ),
            ElectronCollisionHeavyMapping(
                4, "dissociative_ionization", {"Cl2": 1},
                {"Cl+": 1, "Cl": 1}, **common,
            ),
        ),
    )


def _condition(geometry):
    flow = standard_volume_flow_molecules_s(
        20.0,
        standard_temperature_K=273.15,
        standard_pressure_Pa=101325.0,
    )
    return EEDFChlorineCondition(
        condition_id="manufactured-eepf-power",
        geometry=geometry,
        neutral_control_volume=_scalar(geometry.volume_m3, "m3"),
        pressure=_scalar(1.333223684, "Pa"),
        gas_temperature=_scalar(500.0, "K"),
        chlorine_molecule_feed=_scalar(flow, "molecule s^-1"),
        source_power=_scalar(300.0, "W"),
        absorbed_power=AbsorbedPowerEstimate(
            lower_W=20.0,
            upper_W=20.0,
            point_W=20.0,
            boundary_kind="manufactured direct absorbed power",
            measurement_source="manufactured test",
            loss_source="manufactured test",
            measurement_evidence="assumed",
            loss_evidence="assumed",
        ),
        reduced_field_bounds_Td=(20.0, 250.0),
    )


def _providers(condition):
    geometry = condition.geometry
    wall = ChlorineWallRecombinationBoundary(
        recombination_probability=0.02,
        surface_state="manufactured conditioned wall",
        source="manufactured test",
        evidence_kind="measured",
        valid_cl_to_cl2_ratio=(1.0e-5, 1.0e5),
        valid_pressure_Pa=(1.0, 2.0),
        valid_icp_power_W=(100.0, 500.0),
        valid_gas_temperature_K=(500.0, 500.0),
        relative_measurement_uncertainty=0.2,
    )
    neutral = solve_chlorine_neutral_wall_transport(
        geometry=geometry,
        wall_boundary=wall,
        incident_velocity_state=thermalized_chlorine_incident_velocity_state(
            500.0,
            source="manufactured thermal wall state",
            evidence_kind="assumed",
            relative_uncertainty=None,
        ),
        diffusivity_model=lymberopoulos_economou_1995_chlorine_diffusivity(),
        total_neutral_density_m3=condition.target_neutral_density_m3,
        gas_temperature_K=500.0,
        cl_to_cl2_ratio=1.0,
        pressure_Pa=condition.pressure.value,
        icp_power_W=condition.source_power.value,
    )
    charged = ChlorineChargedTransportState(
        geometry=geometry,
        positive_ion_transport={
            "Cl2+": PositiveIonWallTransport(
                800.0, 400.0, "manufactured", "assumed"),
            "Cl+": PositiveIonWallTransport(
                1100.0, 550.0, "manufactured", "assumed"),
        },
        negative_ion_confinement_source="manufactured confinement",
        negative_ion_confinement_evidence="assumed",
    )
    return (
        FixedChlorineChargedTransportProvider(charged),
        FixedChlorineNeutralWallTransportProvider(neutral),
        FixedPositiveIonWallEnergyProvider(PositiveIonWallEnergyState(
            energy_eV_per_lost_ion={"Cl2+": 12.0, "Cl+": 14.0},
            source="manufactured wall energy",
            evidence_kind="assumed",
        )),
    )


def test_absorbed_power_eepf_model_closes_knobs_to_species_fluxes():
    geometry = CylindricalReactor(radius_m=0.10, length_m=0.08)
    condition = _condition(geometry)
    lee = build_lee_lieberman_chlorine_particle_network()
    heavy = ReactionNetwork(species=lee.species, reactions=lee.reactions[6:8])
    deck = _deck()
    model = EEDFChlorineAbsorbedPowerModel(
        DeterministicTwoTermBoltzmannSolver(
            ElectronEnergyGrid.linear(80.0, 320), deck),
        _chemistry(deck, heavy.species),
        heavy,
    )
    charged, neutral, wall_energy = _providers(condition)
    solution = model.solve(
        condition,
        charged_transport_provider=charged,
        neutral_wall_transport_provider=neutral,
        wall_energy_provider=wall_energy,
        initial_reduced_electric_field_Td=100.0,
        maximum_tail_population_fraction=1.0e-5,
        residual_tolerance=2.0e-7,
        maximum_evaluations=800,
    )
    assert solution.maximum_normalized_residual < 2.0e-7
    assert solution.modeled_power_density_W_m3 == pytest.approx(
        solution.absorbed_power_density_W_m3, rel=2.0e-7)
    assert 20.0 < solution.reduced_electric_field_Td < 250.0
    assert solution.mean_electron_energy_eV > 0.0
    assert solution.axial_positive_ion_flux_m2_s["Cl2+"] > 0.0
    assert solution.axial_positive_ion_flux_m2_s["Cl+"] > 0.0
    assert (
        solution.collision_chemistry_state.relative_electron_growth_closure
        < 2.0e-15
    )
    assert not solution.supports_implicit_differentiation
    assert not solution.supports_reactor_state_prediction
    assert not solution.supports_wafer_flux
    assert not solution.supports_feature_depth


def test_eq11_condition_constrains_chlorine_nuclei_not_particle_count():
    condition = replace(
        _condition(CylindricalReactor(radius_m=0.10, length_m=0.08)),
        neutral_density_constraint="chlorine_nuclei_equivalent_molecules",
    )
    target = condition.target_neutral_density_m3
    densities = {"Cl2": 0.6 * target, "Cl": 0.8 * target}
    assert condition.constrained_neutral_density_m3(densities) == (
        pytest.approx(target))
    assert densities["Cl2"] + densities["Cl"] == pytest.approx(1.4 * target)
    assert condition.angular_field_frequency_over_density(
        1.4 * target) == pytest.approx(0.0)
