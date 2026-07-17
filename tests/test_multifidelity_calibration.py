import numpy as np
import pytest

from petch.multifidelity_calibration import (
    assess_high_fidelity_trial,
    build_corrected_multifidelity_model,
    fit_affine_response,
    predict_affine,
    propose_trust_region_step,
    scaled_l2_merit,
)


def _responses(parameters):
    parameters = np.asarray(parameters, dtype=float)
    low = np.column_stack([
        10.0 + 2.0 * parameters[:, 0] + parameters[:, 1],
        20.0 - parameters[:, 0] + 3.0 * parameters[:, 1],
    ])
    discrepancy = np.column_stack([
        1.0 + 0.5 * parameters[:, 0],
        -2.0 + 0.25 * parameters[:, 1],
    ])
    return low, low + discrepancy


def test_affine_multifidelity_model_recovers_value_and_gradient_and_proposes_target():
    center = np.asarray([0.5, 0.5])
    low_parameters = np.asarray([
        [0.5, 0.5], [0.3, 0.5], [0.7, 0.5], [0.5, 0.3], [0.5, 0.7]])
    low, _ = _responses(low_parameters)
    paired_parameters = np.asarray([[0.5, 0.5], [0.6, 0.5], [0.5, 0.6]])
    paired_low, paired_high = _responses(paired_parameters)
    model = build_corrected_multifidelity_model(
        low_parameters, low, paired_parameters, paired_low, paired_high,
        center=center, parameter_scale=[0.1, 0.1], response_scale=[1.0, 1.0])

    expected_center = _responses(center[None, :])[1][0]
    assert model["first_order_discrepancy_identified"] is True
    assert model["value"] == pytest.approx(expected_center)
    assert predict_affine(model, [0.55, 0.45]) == pytest.approx(
        _responses(np.asarray([[0.55, 0.45]]))[1][0])

    target_parameters = np.asarray([0.55, 0.45])
    target = _responses(target_parameters[None, :])[1][0]
    proposal = propose_trust_region_step(
        model, target, [[0.0, 1.0], [0.0, 1.0]], radius=1.0)
    assert proposal["candidate"] == pytest.approx(target_parameters)
    assert proposal["predicted_merit"] == pytest.approx(0.0, abs=1e-12)


def test_center_only_discrepancy_is_a_receipt_not_proposal_authority():
    center = np.asarray([0.5, 0.5])
    parameters = np.asarray([[0.5, 0.5], [0.6, 0.5], [0.5, 0.6]])
    low, high = _responses(parameters)
    model = build_corrected_multifidelity_model(
        parameters, low, parameters[:1], low[:1], high[:1],
        center=center, parameter_scale=[0.1, 0.1], response_scale=[1.0, 1.0])

    assert model["first_order_discrepancy_identified"] is False
    assert model["discrepancy_direction_rank"] == 0
    with pytest.raises(ValueError, match="first-order"):
        propose_trust_region_step(
            model, high[1], [[0.0, 1.0], [0.0, 1.0]], radius=1.0)


def test_affine_fit_refuses_rank_deficient_or_ill_conditioned_design():
    with pytest.raises(ValueError, match="span"):
        fit_affine_response(
            [[0.0, 0.0], [0.1, 0.0], [0.2, 0.0]],
            [[0.0], [1.0], [2.0]], center=[0.0, 0.0], parameter_scale=[1.0, 1.0])


def test_trial_ratio_rejects_bad_high_fidelity_direction_and_shrinks():
    proposal = {
        "current_merit": 4.0,
        "predicted_merit_reduction": 3.0,
        "radius": 1.0,
        "trust_boundary_active": True,
    }
    result = assess_high_fidelity_trial(
        proposal, actual_response=[5.0, 0.0], target=[0.0, 0.0],
        response_scale=[1.0, 1.0])

    assert result["accepted"] is False
    assert result["decision"] == "reject_shrink"
    assert result["trust_ratio_rho"] < 0.0
    assert result["next_radius"] == pytest.approx(0.5)


def test_trial_ratio_accepts_verified_boundary_step_and_grows():
    proposal = {
        "current_merit": 4.0,
        "predicted_merit_reduction": 3.0,
        "radius": 0.5,
        "trust_boundary_active": True,
    }
    result = assess_high_fidelity_trial(
        proposal, actual_response=[1.0, 0.0], target=[0.0, 0.0],
        response_scale=[1.0, 1.0], maximum_radius=2.0)

    assert result["accepted"] is True
    assert result["decision"] == "accept_grow"
    assert result["trust_ratio_rho"] == pytest.approx(1.0)
    assert result["next_radius"] == pytest.approx(1.0)


def test_scaled_merit_requires_declared_positive_scales():
    assert scaled_l2_merit([3.0, 4.0], [0.0, 0.0], [1.0, 1.0]) == 5.0
    with pytest.raises(ValueError, match="positive"):
        scaled_l2_merit([1.0], [0.0], [0.0])
