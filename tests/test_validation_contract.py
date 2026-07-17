import pytest

from petch.validation_contract import (
    CalibrationReveal,
    HeldOutPrediction,
    ValidationObservation,
    ValidationParameter,
    ValidationProtocol,
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
