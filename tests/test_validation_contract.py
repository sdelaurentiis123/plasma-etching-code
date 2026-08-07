import pytest

from petch.validation_contract import (
    CalibrationReveal,
    HeldOutPrediction,
    ObservationValueReveal,
    PreRegisteredValidationProtocol,
    ValidationObservation,
    ValidationParameter,
    ValidationProtocol,
    ValidationTargetCommitment,
    reveal_preregistered_observations,
    score_held_out_predictions,
)


SOURCE_SHA = "a" * 64
RUN_SHA = "b" * 64
OPERATOR_EPOCH_SHA = "c" * 64
OTHER_OPERATOR_EPOCH_SHA = "d" * 64


def _observation(observation_id, split, value, measurement_uncertainty=0.1):
    return ValidationObservation(
        observation_id=observation_id,
        condition_id=observation_id,
        observable="etch_depth",
        value=value,
        unit="um",
        split=split,
        source="manufactured experimental gate",
        source_sha256=SOURCE_SHA,
        digitization_uncertainty=0.02,
        measurement_uncertainty=measurement_uncertainty,
    )


def _protocol(*, held_measurement_uncertainty=0.1, predictive_parameter=True):
    return ValidationProtocol(
        protocol_id="manufactured-transfer",
        revision="R1",
        intended_use="predict etch depth at a held-out process condition",
        observations=(
            _observation("cal", "calibration", 1.0),
            _observation(
                "held", "held_out_transfer", 2.0,
                measurement_uncertainty=held_measurement_uncertainty,
            ),
        ),
        parameters=(
            ValidationParameter(
                name="reaction_probability",
                role="physical_mechanism",
                unit="1",
                lower=0.0,
                upper=1.0,
                source="bounded literature prior",
                calibration_allowed=True,
                supports_prediction=predictive_parameter,
            ),
            ValidationParameter(
                name="sample_count",
                role="numerical_control",
                unit="count",
                lower=128,
                upper=65536,
                source="numerical refinement",
                calibration_allowed=False,
                supports_prediction=True,
            ),
        ),
        maximum_calibrated_parameters=1,
        mechanism_id="mechanism-v1",
        boundary_provider_id="reactor-v1",
    )


def _prediction(
    reveal,
    *,
    value=2.05,
    boundary=True,
    mechanism=True,
    operator_epoch_sha256=OPERATOR_EPOCH_SHA,
):
    return HeldOutPrediction(
        observation_id="held",
        predicted_value=value,
        numerical_uncertainty=0.02,
        parameter_uncertainty=0.03,
        model_discrepancy_uncertainty=0.01,
        within_declared_scope=True,
        boundary_supports_prediction=boundary,
        mechanism_supports_prediction=mechanism,
        run_manifest_sha256=RUN_SHA,
        calibration_reveal_sha256=reveal.reveal_sha256,
        operator_epoch_sha256=operator_epoch_sha256,
    )


def _target(observation_id, split, chemistry_family, condition_id=None):
    return ValidationTargetCommitment(
        observation_id=observation_id,
        chemistry_family=chemistry_family,
        material="Si",
        condition_id=condition_id or observation_id,
        observable="etch_depth",
        unit="um",
        split=split,
        boundary_evidence_tier="A_species_energy_angle_measured",
        source="manufactured source panel",
        source_locator=f"figure-1:{observation_id}",
        source_sha256=SOURCE_SHA,
    )


def _preregistration():
    return PreRegisteredValidationProtocol(
        protocol_id="value-blind-cross-chemistry",
        revision="R1",
        intended_use="predict held-out depth in two chemistry families",
        targets=(
            _target("cal", "calibration", "chlorine"),
            _target("held_cl", "held_out_transfer", "chlorine"),
            _target("held_f", "held_out_transfer", "fluorine"),
        ),
        parameters=_protocol().parameters,
        maximum_calibrated_parameters=1,
        mechanism_id="species-resolved-mechanisms-v1",
        boundary_provider_id="measured-beam-boundaries-v1",
        minimum_held_out_chemistry_families=2,
    )


def test_numerical_controls_cannot_be_calibration_parameters():
    with pytest.raises(ValueError, match="validation parameter"):
        ValidationParameter(
            name="grid_spacing",
            role="numerical_control",
            unit="m",
            lower=1.0e-9,
            upper=1.0e-6,
            source="solver control",
            calibration_allowed=True,
            supports_prediction=True,
        )


def test_value_blind_preregistration_freezes_sources_splits_and_chemistries():
    preregistration = _preregistration()

    assert len(preregistration.commit_sha256) == 64
    assert {
        item.chemistry_family for item in preregistration.targets
        if item.split == "held_out_transfer"
    } == {"chlorine", "fluorine"}
    assert all(not hasattr(item, "value") for item in preregistration.targets)


def test_observation_reveal_requires_exactly_the_frozen_target_set():
    preregistration = _preregistration()
    reveals = (
        ObservationValueReveal("cal", 1.0, 0.02, 0.1),
        ObservationValueReveal("held_cl", 2.0, 0.02, 0.1),
    )

    with pytest.raises(ValueError, match="every and only preregistered"):
        reveal_preregistered_observations(preregistration, reveals)


def test_observation_reveal_preserves_committed_metadata_and_binds_protocol():
    preregistration = _preregistration()
    protocol = reveal_preregistered_observations(
        preregistration,
        (
            ObservationValueReveal("cal", 1.0, 0.02, 0.1),
            ObservationValueReveal("held_cl", 2.0, 0.02, 0.1),
            ObservationValueReveal("held_f", 3.0, 0.03, 0.2),
        ),
    )

    assert protocol.preregistration_sha256 == preregistration.commit_sha256
    by_id = {item.observation_id: item for item in protocol.observations}
    assert by_id["held_f"].chemistry_family == "fluorine"
    assert by_id["held_f"].source_locator == "figure-1:held_f"
    assert by_id["held_f"].source_sha256 == SOURCE_SHA
    assert by_id["held_f"].split == "held_out_transfer"
    assert by_id["held_f"].value == 3.0


def test_preregistration_refuses_condition_leakage_between_splits():
    with pytest.raises(ValueError, match="preregistered"):
        PreRegisteredValidationProtocol(
            protocol_id="leaky",
            revision="R1",
            intended_use="reject calibration and held-out reuse",
            targets=(
                _target("cal", "calibration", "chlorine", "same-condition"),
                _target(
                    "held", "held_out_transfer", "chlorine", "same-condition"),
            ),
            parameters=_protocol().parameters,
            maximum_calibrated_parameters=1,
            mechanism_id="mechanism-v1",
            boundary_provider_id="boundary-v1",
            minimum_held_out_chemistry_families=1,
        )


def test_preregistration_enforces_held_out_cross_chemistry_count():
    with pytest.raises(ValueError, match="preregistered"):
        PreRegisteredValidationProtocol(
            protocol_id="single-chemistry",
            revision="R1",
            intended_use="reject an undersized held-out chemistry set",
            targets=(
                _target("cal", "calibration", "chlorine"),
                _target("held", "held_out_transfer", "chlorine"),
            ),
            parameters=_protocol().parameters,
            maximum_calibrated_parameters=1,
            mechanism_id="mechanism-v1",
            boundary_provider_id="boundary-v1",
            minimum_held_out_chemistry_families=2,
        )


def test_protocol_and_reveal_freeze_exact_parameter_and_data_splits():
    protocol = _protocol()
    reveal = CalibrationReveal.from_protocol(
        protocol,
        {"reaction_probability": 0.3},
        operator_epoch_sha256=OPERATOR_EPOCH_SHA,
    )

    assert len(protocol.calibration_observations) == 1
    assert len(protocol.held_out_observations) == 1
    assert reveal.protocol_commit_sha256 == protocol.commit_sha256
    with pytest.raises(ValueError, match="exactly"):
        CalibrationReveal.from_protocol(
            protocol,
            {"reaction_probability": 0.3, "sample_count": 4096},
            operator_epoch_sha256=OPERATOR_EPOCH_SHA,
        )


def test_score_accepts_every_and_only_held_out_rows():
    protocol = _protocol()
    reveal = CalibrationReveal.from_protocol(
        protocol,
        {"reaction_probability": 0.3},
        operator_epoch_sha256=OPERATOR_EPOCH_SHA,
    )
    with pytest.raises(ValueError, match="every and only held-out"):
        score_held_out_predictions(protocol, reveal, ())

    score = score_held_out_predictions(
        protocol, reveal, (_prediction(reveal),))
    assert score.formal_validation_passed
    assert score.uncertainty_complete
    assert score.absolute_error["held"] == pytest.approx(0.05)
    assert score.operator_epoch_sha256 == OPERATOR_EPOCH_SHA


def test_score_refuses_prediction_from_a_different_operator_epoch():
    protocol = _protocol()
    reveal = CalibrationReveal.from_protocol(
        protocol,
        {"reaction_probability": 0.3},
        operator_epoch_sha256=OPERATOR_EPOCH_SHA,
    )

    with pytest.raises(ValueError, match="operator epoch"):
        score_held_out_predictions(
            protocol,
            reveal,
            (_prediction(
                reveal,
                operator_epoch_sha256=OTHER_OPERATOR_EPOCH_SHA,
            ),),
        )


def test_operator_epoch_is_bound_into_calibration_reveal_checksum():
    protocol = _protocol()
    values = {"reaction_probability": 0.3}
    first = CalibrationReveal.from_protocol(
        protocol,
        values,
        operator_epoch_sha256=OPERATOR_EPOCH_SHA,
    )
    second = CalibrationReveal.from_protocol(
        protocol,
        values,
        operator_epoch_sha256=OTHER_OPERATOR_EPOCH_SHA,
    )

    assert first.reveal_sha256 != second.reveal_sha256


def test_formal_score_refuses_missing_uncertainty_or_nonpredictive_inputs():
    missing_uncertainty = _protocol(held_measurement_uncertainty=None)
    reveal = CalibrationReveal.from_protocol(
        missing_uncertainty,
        {"reaction_probability": 0.3},
        operator_epoch_sha256=OPERATOR_EPOCH_SHA,
    )
    score = score_held_out_predictions(
        missing_uncertainty, reveal, (_prediction(reveal),))
    assert not score.formal_validation_passed
    assert not score.uncertainty_complete
    assert any("measurement uncertainty" in reason for reason in score.reasons)

    nonpredictive = _protocol(predictive_parameter=False)
    reveal = CalibrationReveal.from_protocol(
        nonpredictive,
        {"reaction_probability": 0.3},
        operator_epoch_sha256=OPERATOR_EPOCH_SHA,
    )
    score = score_held_out_predictions(
        nonpredictive,
        reveal,
        (_prediction(reveal, boundary=False),),
    )
    assert not score.formal_validation_passed
    assert any("reactor boundary" in reason for reason in score.reasons)
    assert any("lack predictive evidence" in reason for reason in score.reasons)
