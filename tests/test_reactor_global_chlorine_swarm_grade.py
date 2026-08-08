from dataclasses import replace

import pytest

from petch.reactor_global.chlorine_swarm import (
    GonzalezMaganaPureChlorineSwarmBoard,
)
from petch.reactor_global.chlorine_swarm_grade import (
    ElectronSwarmPrediction,
    grade_gonzalez_magana_pure_cl2_swarm,
)


DECK_HASH = "a" * 64


def _exact_predictions(board):
    return tuple(
        ElectronSwarmPrediction(
            observation_id=item.observation_id,
            transport_definition=item.transport_definition,
            value_si=item.value_si,
            si_unit=item.si_unit,
            solver_id="manufactured-exact-solver-v1",
            collision_deck_sha256=DECK_HASH,
        )
        for item in board.measurements
    )


def test_measurement_board_exposes_noninterchangeable_transport_definitions():
    board = GonzalezMaganaPureChlorineSwarmBoard.from_package_data()
    definitions = {
        item.observable: item.transport_definition
        for item in board.measurements
    }
    assert definitions == {
        "electron_drift_velocity": (
            "pulsed_townsend_mean_arrival_time_drift_velocity"
        ),
        "effective_ionization_coefficient": (
            "steady_state_townsend_effective_ionization_coefficient"
        ),
        "density_normalized_longitudinal_diffusion": (
            "pulsed_townsend_spatiotemporal_longitudinal_diffusion"
        ),
    }


def test_exact_manufactured_predictions_close_board_without_depth_authority():
    board = GonzalezMaganaPureChlorineSwarmBoard.from_package_data()
    grade = grade_gonzalez_magana_pure_cl2_swarm(
        board, _exact_predictions(board))
    assert grade.all_inside_source_wide_typical_intervals
    assert grade.inside_interval_counts == {
        "density_normalized_longitudinal_diffusion": (8, 8),
        "effective_ionization_coefficient": (21, 21),
        "electron_drift_velocity": (23, 23),
    }
    assert grade.supports_collision_set_validation
    assert not grade.supports_reactor_state_prediction
    assert not grade.supports_wafer_flux
    assert not grade.supports_feature_depth


def test_flux_mobility_cannot_masquerade_as_mean_arrival_time_drift():
    board = GonzalezMaganaPureChlorineSwarmBoard.from_package_data()
    predictions = list(_exact_predictions(board))
    drift_index = next(
        index for index, item in enumerate(board.measurements)
        if item.observable == "electron_drift_velocity"
    )
    predictions[drift_index] = replace(
        predictions[drift_index],
        transport_definition="flux_drift_velocity",
    )
    with pytest.raises(ValueError, match="transport definition mismatch"):
        grade_gonzalez_magana_pure_cl2_swarm(
            board, tuple(predictions))


def test_grade_refuses_missing_markers_mixed_decks_and_outside_interval():
    board = GonzalezMaganaPureChlorineSwarmBoard.from_package_data()
    predictions = _exact_predictions(board)
    with pytest.raises(ValueError, match="exactly one prediction"):
        grade_gonzalez_magana_pure_cl2_swarm(board, predictions[:-1])

    mixed = list(predictions)
    mixed[-1] = replace(mixed[-1], collision_deck_sha256="b" * 64)
    with pytest.raises(ValueError, match="one solver and collision deck"):
        grade_gonzalez_magana_pure_cl2_swarm(board, tuple(mixed))

    displaced = list(predictions)
    displaced[0] = replace(
        displaced[0], value_si=1.2 * displaced[0].value_si)
    grade = grade_gonzalez_magana_pure_cl2_swarm(
        board, tuple(displaced))
    assert not grade.all_inside_source_wide_typical_intervals
    assert grade.inside_interval_counts["electron_drift_velocity"] == (22, 23)
