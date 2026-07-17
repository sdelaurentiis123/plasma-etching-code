"""Evidence-aware affine multi-fidelity trust-region primitives.

The high-fidelity model remains the authority.  A low-fidelity response may
propose a step only after the low/high discrepancy is anchored at the current
center.  By default a proposal additionally requires enough paired directions
to identify a first-order discrepancy correction.  This prevents a cheap model
from steering calibration merely because it is fast.

The routines are deliberately independent of any plasma chemistry or benchmark
data.  They operate on small arrays of parameters and declared observables so
the same controller can serve feature, reactor-boundary, or material-closure
calibration without becoming a second physics engine.
"""
from __future__ import annotations

import numpy as np


def _array(name, value, *, ndim=None):
    output = np.asarray(value, dtype=float)
    if (ndim is not None and output.ndim != int(ndim)) or np.any(~np.isfinite(output)):
        raise ValueError(f"{name} must be a finite {ndim}-D array")
    return output


def scaled_l2_merit(response, target, response_scale):
    """Return the declared-scale L2 calibration residual."""
    response = _array("response", response, ndim=1)
    target = _array("target", target, ndim=1)
    scale = _array("response_scale", response_scale, ndim=1)
    if response.shape != target.shape or scale.shape != response.shape or np.any(scale <= 0.0):
        raise ValueError("response, target, and positive response_scale must match")
    return float(np.linalg.norm((response - target) / scale))


def fit_affine_response(parameters, responses, *, center, parameter_scale,
                        response_standard_uncertainty=None,
                        maximum_design_condition=1e6):
    """Fit a local affine response with rank and conditioning receipts.

    ``response_standard_uncertainty`` supplies positive per-sample/per-response
    standard uncertainties.  Each response column is then fit by weighted least
    squares.  It is an estimator weight, not permission to hide model residuals.
    """
    parameters = _array("parameters", parameters, ndim=2)
    responses = _array("responses", responses, ndim=2)
    center = _array("center", center, ndim=1)
    parameter_scale = _array("parameter_scale", parameter_scale, ndim=1)
    if (parameters.shape[0] != responses.shape[0]
            or parameters.shape[1] != center.size
            or parameter_scale.shape != center.shape
            or np.any(parameter_scale <= 0.0)):
        raise ValueError("affine-response array shapes or parameter scale are invalid")
    sample_count, parameter_count = parameters.shape
    if sample_count < parameter_count + 1:
        raise ValueError("an affine response needs at least n_parameter + 1 samples")
    design = np.column_stack([
        np.ones(sample_count), (parameters - center) / parameter_scale])
    if np.linalg.matrix_rank(design) < parameter_count + 1:
        raise ValueError("affine-response samples do not span parameter space")

    if response_standard_uncertainty is None:
        uncertainty = np.ones_like(responses)
        weighting = "uniform"
    else:
        uncertainty = _array(
            "response_standard_uncertainty", response_standard_uncertainty, ndim=2)
        if uncertainty.shape != responses.shape or np.any(uncertainty <= 0.0):
            raise ValueError("response uncertainty must be positive and match responses")
        weighting = "inverse-variance"

    coefficients = np.empty((parameter_count + 1, responses.shape[1]), dtype=float)
    conditions = []
    for response_index in range(responses.shape[1]):
        weighted_design = design / uncertainty[:, response_index, None]
        weighted_response = responses[:, response_index] / uncertainty[:, response_index]
        condition = float(np.linalg.cond(weighted_design))
        if not np.isfinite(condition) or condition > float(maximum_design_condition):
            raise ValueError("affine-response design is ill-conditioned")
        coefficients[:, response_index] = np.linalg.lstsq(
            weighted_design, weighted_response, rcond=None)[0]
        conditions.append(condition)

    fitted = design @ coefficients
    residual = responses - fitted
    degrees_of_freedom = sample_count - parameter_count - 1
    residual_rms = np.sqrt(np.mean(residual * residual, axis=0))
    return {
        "center": center,
        "value": coefficients[0],
        "jacobian": coefficients[1:].T / parameter_scale[None, :],
        "parameter_scale": parameter_scale,
        "sample_count": int(sample_count),
        "parameter_count": int(parameter_count),
        "response_count": int(responses.shape[1]),
        "design_rank": int(np.linalg.matrix_rank(design)),
        "maximum_weighted_design_condition": float(max(conditions)),
        "weighting": weighting,
        "degrees_of_freedom": int(degrees_of_freedom),
        "residual_rms": residual_rms,
        "maximum_absolute_residual": np.max(np.abs(residual), axis=0),
    }


def build_corrected_multifidelity_model(
        low_parameters, low_responses, paired_parameters,
        paired_low_responses, paired_high_responses, *, center,
        parameter_scale, response_scale, low_standard_uncertainty=None,
        maximum_design_condition=1e6, center_tolerance=1e-12):
    """Build a value-anchored low-plus-discrepancy affine model.

    A paired high/low observation at ``center`` is mandatory.  Additional
    linearly independent paired directions identify the discrepancy gradient.
    With only the center pair the returned model is useful as a readiness
    receipt but is not first-order identified and cannot propose by default.
    """
    center = _array("center", center, ndim=1)
    parameter_scale = _array("parameter_scale", parameter_scale, ndim=1)
    response_scale = _array("response_scale", response_scale, ndim=1)
    if (parameter_scale.shape != center.shape or np.any(parameter_scale <= 0.0)
            or np.any(response_scale <= 0.0)):
        raise ValueError("positive parameter and response scales are required")
    low = fit_affine_response(
        low_parameters, low_responses, center=center,
        parameter_scale=parameter_scale,
        response_standard_uncertainty=low_standard_uncertainty,
        maximum_design_condition=maximum_design_condition)

    paired_parameters = _array("paired_parameters", paired_parameters, ndim=2)
    paired_low = _array("paired_low_responses", paired_low_responses, ndim=2)
    paired_high = _array("paired_high_responses", paired_high_responses, ndim=2)
    if (paired_parameters.shape[1] != center.size
            or paired_low.shape != paired_high.shape
            or paired_low.shape[0] != paired_parameters.shape[0]
            or paired_low.shape[1] != response_scale.size):
        raise ValueError("paired multi-fidelity arrays are incompatible")
    normalized_offset = (paired_parameters - center) / parameter_scale
    center_index = np.flatnonzero(
        np.max(np.abs(normalized_offset), axis=1) <= float(center_tolerance))
    if center_index.size != 1:
        raise ValueError("exactly one paired low/high observation must anchor the center")
    center_index = int(center_index[0])
    discrepancy = paired_high - paired_low
    center_discrepancy = discrepancy[center_index]

    other = np.arange(paired_parameters.shape[0]) != center_index
    direction_design = normalized_offset[other]
    direction_delta = discrepancy[other] - center_discrepancy
    direction_rank = int(np.linalg.matrix_rank(direction_design)) if np.any(other) else 0
    first_order_identified = bool(direction_rank == center.size)
    discrepancy_jacobian = np.zeros((response_scale.size, center.size), dtype=float)
    discrepancy_condition = None
    discrepancy_residual_rms = np.zeros(response_scale.size, dtype=float)
    if first_order_identified:
        discrepancy_condition = float(np.linalg.cond(direction_design))
        if (not np.isfinite(discrepancy_condition)
                or discrepancy_condition > float(maximum_design_condition)):
            raise ValueError("paired discrepancy directions are ill-conditioned")
        normalized_slope = np.linalg.lstsq(
            direction_design, direction_delta, rcond=None)[0]
        discrepancy_jacobian = normalized_slope.T / parameter_scale[None, :]
        discrepancy_residual = direction_delta - direction_design @ normalized_slope
        discrepancy_residual_rms = np.sqrt(np.mean(
            discrepancy_residual * discrepancy_residual, axis=0))

    value = low["value"] + center_discrepancy
    jacobian = low["jacobian"] + discrepancy_jacobian
    scaled_jacobian = (
        jacobian * parameter_scale[None, :] / response_scale[:, None])
    response_condition = float(np.linalg.cond(scaled_jacobian))
    response_rank = int(np.linalg.matrix_rank(scaled_jacobian))
    return {
        "center": center,
        "value": value,
        "jacobian": jacobian,
        "parameter_scale": parameter_scale,
        "response_scale": response_scale,
        "low_model": low,
        "paired_count": int(paired_parameters.shape[0]),
        "center_discrepancy": center_discrepancy,
        "discrepancy_jacobian": discrepancy_jacobian,
        "discrepancy_direction_rank": direction_rank,
        "discrepancy_direction_condition": discrepancy_condition,
        "discrepancy_residual_rms": discrepancy_residual_rms,
        "first_order_discrepancy_identified": first_order_identified,
        "scaled_response_jacobian": scaled_jacobian,
        "scaled_response_jacobian_rank": response_rank,
        "scaled_response_jacobian_condition": response_condition,
    }


def predict_affine(model, parameters):
    parameters = _array("parameters", parameters, ndim=1)
    center = _array("model center", model["center"], ndim=1)
    value = _array("model value", model["value"], ndim=1)
    jacobian = _array("model jacobian", model["jacobian"], ndim=2)
    if parameters.shape != center.shape or jacobian.shape != (value.size, center.size):
        raise ValueError("affine model shapes are invalid")
    return value + jacobian @ (parameters - center)


def propose_trust_region_step(model, target, bounds, *, radius,
                              require_first_order_discrepancy=True,
                              maximum_response_condition=1e6):
    """Propose one bounded response-matching step from a corrected model."""
    if (require_first_order_discrepancy
            and not model.get("first_order_discrepancy_identified", False)):
        raise ValueError("first-order low/high discrepancy is not identified")
    center = _array("model center", model["center"], ndim=1)
    value = _array("model value", model["value"], ndim=1)
    jacobian = _array("model jacobian", model["jacobian"], ndim=2)
    parameter_scale = _array("parameter_scale", model["parameter_scale"], ndim=1)
    response_scale = _array("response_scale", model["response_scale"], ndim=1)
    target = _array("target", target, ndim=1)
    bounds = _array("bounds", bounds, ndim=2)
    if (target.shape != value.shape or bounds.shape != (center.size, 2)
            or np.any(bounds[:, 0] >= bounds[:, 1])
            or np.any(center < bounds[:, 0]) or np.any(center > bounds[:, 1])
            or not np.isfinite(radius) or radius <= 0.0):
        raise ValueError("trust-region target, bounds, center, or radius are invalid")
    scaled_jacobian = jacobian * parameter_scale[None, :] / response_scale[:, None]
    condition = float(np.linalg.cond(scaled_jacobian))
    if (np.linalg.matrix_rank(scaled_jacobian) < center.size
            or not np.isfinite(condition)
            or condition > float(maximum_response_condition)):
        raise ValueError("scaled calibration response is rank deficient or ill-conditioned")
    scaled_residual = (value - target) / response_scale
    full_scaled_step = np.linalg.lstsq(
        scaled_jacobian, -scaled_residual, rcond=None)[0]
    full_norm = float(np.max(np.abs(full_scaled_step)))
    trust_scale = min(1.0, float(radius) / full_norm) if full_norm > 0.0 else 1.0
    scaled_step = trust_scale * full_scaled_step
    raw_step = parameter_scale * scaled_step

    bound_scale = 1.0
    for index, direction in enumerate(raw_step):
        if direction > 0.0:
            bound_scale = min(
                bound_scale, (bounds[index, 1] - center[index]) / direction)
        elif direction < 0.0:
            bound_scale = min(
                bound_scale, (bounds[index, 0] - center[index]) / direction)
    bound_scale = float(np.clip(bound_scale, 0.0, 1.0))
    scaled_step *= bound_scale
    raw_step *= bound_scale
    candidate = center + raw_step
    prediction = value + jacobian @ raw_step
    current_merit = scaled_l2_merit(value, target, response_scale)
    predicted_merit = scaled_l2_merit(prediction, target, response_scale)
    predicted_reduction = current_merit - predicted_merit
    if predicted_reduction <= 0.0:
        raise ValueError("corrected model does not predict a merit reduction")
    return {
        "center": center,
        "candidate": candidate,
        "raw_step": raw_step,
        "scaled_step": scaled_step,
        "full_scaled_step": full_scaled_step,
        "radius": float(radius),
        "trust_direction_scale": float(trust_scale),
        "physical_bound_direction_scale": bound_scale,
        "trust_boundary_active": bool(
            np.isclose(np.max(np.abs(scaled_step)), radius, rtol=1e-12, atol=1e-14)),
        "predicted_response": prediction,
        "current_merit": current_merit,
        "predicted_merit": predicted_merit,
        "predicted_merit_reduction": predicted_reduction,
        "scaled_response_jacobian_condition": condition,
    }


def assess_high_fidelity_trial(proposal, actual_response, target, response_scale, *,
                               acceptance_threshold=0.1,
                               growth_threshold=0.75,
                               shrink_factor=0.5, growth_factor=2.0,
                               maximum_radius=np.inf):
    """Accept/reject one expensive trial and update only the numerical radius."""
    if (not 0.0 < acceptance_threshold < growth_threshold < 1.0
            or not 0.0 < shrink_factor < 1.0 or growth_factor <= 1.0
            or maximum_radius <= 0.0):
        raise ValueError("invalid trust-region update controls")
    actual_merit = scaled_l2_merit(actual_response, target, response_scale)
    current_merit = float(proposal["current_merit"])
    predicted_reduction = float(proposal["predicted_merit_reduction"])
    if not np.isfinite(predicted_reduction) or predicted_reduction <= 0.0:
        raise ValueError("proposal has no positive predicted reduction")
    actual_reduction = current_merit - actual_merit
    rho = actual_reduction / predicted_reduction
    accepted = bool(actual_reduction > 0.0 and rho >= acceptance_threshold)
    radius = float(proposal["radius"])
    if not accepted:
        decision = "reject_shrink"
        next_radius = radius * shrink_factor
    elif rho >= growth_threshold and proposal.get("trust_boundary_active", False):
        decision = "accept_grow"
        next_radius = min(float(maximum_radius), radius * growth_factor)
    else:
        decision = "accept_hold"
        next_radius = radius
    return {
        "decision": decision,
        "accepted": accepted,
        "actual_merit": actual_merit,
        "actual_merit_reduction": actual_reduction,
        "predicted_merit_reduction": predicted_reduction,
        "trust_ratio_rho": float(rho),
        "previous_radius": radius,
        "next_radius": float(next_radius),
    }
