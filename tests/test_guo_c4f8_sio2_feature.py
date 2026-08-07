import numpy as np
import pytest

from petch.guo_c4f8_sio2 import (
    GuoC4F8ArSiO2Mechanism,
    GuoIncidentComposition,
    GuoIonQuadrature,
    GuoSourceLawUnderspecified,
    GuoTmlState,
    nearest_neighbor_probability,
)
from petch.guo_c4f8_sio2_feature import (
    GuoC4F8ArSiO2FeatureMechanism,
    GuoC4F8ArSiO2FeatureState,
)
from petch.surface_kinetics import (
    EnergeticFlux,
    FaceResolvedEnergeticFlux,
    SurfaceFluxes,
)


def test_neutral_collision_probabilities_are_source_products_not_raw_coefficients():
    state = GuoC4F8ArSiO2FeatureState(
        0.39, 0.36, 0.12, 0.13, 0.05,
    )
    mechanism = GuoC4F8ArSiO2FeatureMechanism(
        neutral_species=("F", "CF2"),
        ion_species_mapping={"Ar+": "Ar"},
    )
    probability = mechanism.neutral_reaction_probability(state)
    scalar = GuoTmlState(0.39, 0.36, 0.12, 0.13, 0.05)

    assert probability["F"] == pytest.approx(
        20.0 * nearest_neighbor_probability(scalar, "Si", "V"))
    assert probability["CF2"] == pytest.approx(
        3.5 * nearest_neighbor_probability(scalar, "C", "V")
        + 1.8 * nearest_neighbor_probability(scalar, "O", "V"))
    assert 0.0 < probability["F"] < 1.0
    assert 0.0 < probability["CF2"] < 1.0


def test_collision_probability_refuses_unpublished_saturation_closure():
    mechanism = GuoC4F8ArSiO2FeatureMechanism(
        neutral_species=("F",),
        ion_species_mapping={"Ar+": "Ar"},
    )
    with pytest.raises(GuoSourceLawUnderspecified, match=r"leave \[0,1\]"):
        mechanism.neutral_reaction_probability(
            GuoC4F8ArSiO2FeatureState(
                0.95, 0.02, 0.01, 0.02, 10.0,
            )
        )


def test_planar_feature_adapter_replays_the_same_guo_steady_state_and_yield():
    neutral_ratio = {"CF": 3.7, "CF2": 7.8, "CF3": 0.7, "O": 6.4}
    ion_flux = 1.2e20
    fluxes = SurfaceFluxes(
        {
            name: ratio * ion_flux
            for name, ratio in neutral_ratio.items()
        },
        (EnergeticFlux(
            "Ar+", ion_flux, np.asarray([350.0]),
            np.asarray([1.0]), np.asarray([1.0]),
        ),),
    )
    feature = GuoC4F8ArSiO2FeatureMechanism(
        neutral_species=tuple(neutral_ratio),
        ion_species_mapping={"Ar+": "Ar"},
    )
    observed = feature.advance(feature.initial_state(), fluxes, 2.0)
    reference = GuoC4F8ArSiO2Mechanism(
        GuoIncidentComposition(neutral_ratio, {}),
        GuoIonQuadrature.monoenergetic(350.0),
    ).solve_steady_state_algebraic()

    assert observed.validity.within_declared_scope
    assert observed.sio2_yield_per_ion == pytest.approx(
        reference.sio2_yield_per_ion)
    assert observed.state.si == pytest.approx(reference.state.si)
    assert observed.state.vacancy == pytest.approx(reference.state.vacancy)
    assert observed.etch_velocity_m_s == pytest.approx(
        reference.sio2_yield_per_ion * ion_flux / 2.2e28)
    assert observed.removed_sio2_formula_units_m2 == pytest.approx(
        reference.sio2_yield_per_ion * ion_flux * 2.0)
    assert abs(observed.atom_ledger_residual_atoms_per_ion) < 1.0e-12


def test_exact_face_event_measure_drives_independent_local_steady_states():
    mechanism = GuoC4F8ArSiO2FeatureMechanism(
        neutral_species=("CF", "CF2", "CF3", "O"),
        ion_species_mapping={"Ar+": "Ar"},
        allow_out_of_board_transfer_audit=True,
    )
    events = FaceResolvedEnergeticFlux(
        "Ar+",
        2,
        event_face=np.asarray([0, 0, 1]),
        event_flux_m2_s=np.asarray([0.4e20, 0.6e20, 2.0e20]),
        event_energy_eV=np.asarray([500.0, 1500.0, 3500.0]),
        event_cosine_incidence=np.asarray([1.0, 0.95, 0.8]),
    )
    result = mechanism.advance(
        mechanism.initial_state((2,)),
        SurfaceFluxes(
            {
                "CF": np.asarray([3.7e20, 3.7e20]),
                "CF2": np.asarray([7.8e20, 7.8e20]),
                "CF3": np.asarray([0.7e20, 0.7e20]),
                "O": np.asarray([6.4e20, 6.4e20]),
            },
            (events,),
        ),
        0.25,
    )

    assert result.sio2_yield_per_ion.shape == (2,)
    assert result.sio2_yield_per_ion[0] != result.sio2_yield_per_ion[1]
    assert np.all(result.etch_velocity_m_s > 0.0)
    assert np.max(np.abs(
        result.atom_ledger_residual_atoms_per_ion)) < 1.0e-12
    assert not result.validity.within_declared_scope
    assert any(
        "<=370 eV" in reason for reason in result.validity.reasons)


def test_krueger_transfer_mode_executes_but_does_not_upgrade_evidence():
    mechanism = GuoC4F8ArSiO2FeatureMechanism.krueger_2024_transfer_audit()
    ion_flux = 1.2e20
    result = mechanism.advance(
        mechanism.initial_state(),
        SurfaceFluxes(
            {
                "C3F4": 9.5e20,
                "C2F3": 6.8e20,
                "CF": 4.4e20,
                "CF2": 9.4e20,
                "CF3": 8.4e19,
                "O": 7.7e20,
            },
            (EnergeticFlux(
                "ions", ion_flux, [3500.0], [1.0], [1.0]),),
        ),
        1.0,
    )

    assert result.sio2_yield_per_ion > 0.0
    assert not result.validity.within_declared_scope
    assert not result.validity.parameter_evidence_supports_prediction
    assert "feature_depth_used" in mechanism.provenance["calibration"]
    assert mechanism.provenance["calibration"]["feature_depth_used"] is False


def test_feature_adapter_resolves_source_deposition_complementarity_face():
    mechanism = GuoC4F8ArSiO2FeatureMechanism(
        neutral_species=("C3F4", "C2F3", "CF", "CF2", "CF3", "O"),
        ion_species_mapping={"Ar+": "Ar"},
        allow_out_of_board_transfer_audit=True,
    )
    ion_flux = 1.0e20
    ratio = {
        "C3F4": 45.291916210116725,
        "C2F3": 17.647329508552286,
        "CF": 7.59935150762662,
        "CF2": 16.23497822083869,
        "CF3": 1.4507852878196275,
        "O": 16.4880332386623,
    }
    result = mechanism.advance(
        mechanism.initial_state(),
        SurfaceFluxes(
            {name: value * ion_flux for name, value in ratio.items()},
            (EnergeticFlux(
                "Ar+", ion_flux,
                np.asarray([500.0, 1500.0, 3500.0]),
                np.cos(np.deg2rad([4.0, 12.0, 25.0])),
                np.asarray([0.2, 0.3, 0.5]),
            ),),
        ),
        1.0,
    )

    assert result.bdf_fallback_face_count == 0
    assert result.net_movement_atoms_per_ion < 0.0
    assert result.normal_growth_velocity_m_s > 0.0
    assert result.etch_velocity_m_s == 0.0
    assert result.state.c + result.state.f > 0.9
    assert result.steady_state_residual < 2.0e-8
