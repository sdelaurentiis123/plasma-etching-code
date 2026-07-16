"""Refinement and isotropy gates for finite-arrival 3-D twist ensembles."""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .profile_observables_3d import FeatureCenterlineEnsemble3D


@dataclass(frozen=True)
class TwistEnsembleRefinementContract3D:
    minimum_realizations: int
    mean_displacement_tolerance_m: float
    standard_deviation_tolerance_m: float
    onset_probability_tolerance: float
    maximum_systematic_z_score: float = 3.0

    def __post_init__(self):
        values = np.asarray([
            self.mean_displacement_tolerance_m,
            self.standard_deviation_tolerance_m,
            self.onset_probability_tolerance,
            self.maximum_systematic_z_score], dtype=float)
        if (int(self.minimum_realizations) != self.minimum_realizations
                or self.minimum_realizations < 4 or np.any(~np.isfinite(values))
                or np.any(values <= 0.0) or self.onset_probability_tolerance > 1.0):
            raise ValueError("invalid twist-ensemble refinement contract")
        object.__setattr__(self, "minimum_realizations", int(self.minimum_realizations))


@dataclass(frozen=True)
class TwistEnsembleRefinement3DResult:
    passed: bool
    n_doubling_mean_change_m: float
    sample_doubling_mean_change_m: float
    n_doubling_standard_deviation_change_m: float
    sample_doubling_standard_deviation_change_m: float
    n_doubling_onset_probability_change: float
    sample_doubling_onset_probability_change: float
    maximum_systematic_z_score: float
    reasons: tuple[str, ...]
    diagnostics: Mapping[str, object]

    def __post_init__(self):
        values = np.asarray([
            self.n_doubling_mean_change_m,
            self.sample_doubling_mean_change_m,
            self.n_doubling_standard_deviation_change_m,
            self.sample_doubling_standard_deviation_change_m,
            self.n_doubling_onset_probability_change,
            self.sample_doubling_onset_probability_change,
            self.maximum_systematic_z_score], dtype=float)
        if np.any(~np.isfinite(values)) or np.any(values < 0.0):
            raise ValueError("invalid twist-ensemble refinement result")
        object.__setattr__(self, "passed", bool(self.passed))
        object.__setattr__(self, "reasons", tuple(self.reasons))
        object.__setattr__(self, "diagnostics", MappingProxyType(dict(self.diagnostics)))


@dataclass(frozen=True)
class TwistConditionCampaign3DResult:
    """One aspect-ratio condition scored at N, 2N, and doubled transport samples."""

    aspect_ratio: float
    base_n: FeatureCenterlineEnsemble3D
    doubled_n: FeatureCenterlineEnsemble3D
    doubled_sample_level: FeatureCenterlineEnsemble3D
    refinement: TwistEnsembleRefinement3DResult
    base_seeds: tuple[int, ...]
    doubled_seeds: tuple[int, ...]
    doubled_sample_seeds: tuple[int, ...]

    def __post_init__(self):
        if (not np.isfinite(self.aspect_ratio) or self.aspect_ratio <= 0.0
                or not isinstance(self.base_n, FeatureCenterlineEnsemble3D)
                or not isinstance(self.doubled_n, FeatureCenterlineEnsemble3D)
                or not isinstance(self.doubled_sample_level, FeatureCenterlineEnsemble3D)
                or not isinstance(self.refinement, TwistEnsembleRefinement3DResult)):
            raise ValueError("invalid twist-condition campaign result")
        base_seeds = tuple(int(value) for value in self.base_seeds)
        doubled_seeds = tuple(int(value) for value in self.doubled_seeds)
        sample_seeds = tuple(int(value) for value in self.doubled_sample_seeds)
        if (len(base_seeds) != len(self.base_n.members)
                or len(doubled_seeds) != len(self.doubled_n.members)
                or len(sample_seeds) != len(self.doubled_sample_level.members)
                or doubled_seeds[:len(base_seeds)] != base_seeds
                or sample_seeds != doubled_seeds
                or len(set(doubled_seeds)) != len(doubled_seeds)):
            raise ValueError("twist-condition seeds must be nested and paired across refinement")
        object.__setattr__(self, "aspect_ratio", float(self.aspect_ratio))
        object.__setattr__(self, "base_seeds", base_seeds)
        object.__setattr__(self, "doubled_seeds", doubled_seeds)
        object.__setattr__(self, "doubled_sample_seeds", sample_seeds)

    @property
    def numerical_refinement_passed(self):
        return self.refinement.passed

    @property
    def twist_probability(self):
        return _onset_probability(self.doubled_sample_level)

    @property
    def mean_maximum_lateral_displacement_m(self):
        return float(np.mean([
            item.maximum_lateral_displacement_m
            for item in self.doubled_sample_level.members]))

    @property
    def statistical_claim_ready(self):
        # One AR condition cannot by itself satisfy the C5 AR-sweep contract.
        return False


@dataclass(frozen=True)
class TwistAspectRatioCampaign3DResult:
    """Numerically refined statistical twist observables over a declared AR sweep."""

    conditions: tuple[TwistConditionCampaign3DResult, ...]
    onset_probability_threshold: float
    minimum_conditions: int
    passed: bool
    reasons: tuple[str, ...]

    def __post_init__(self):
        conditions = tuple(self.conditions)
        if (not conditions
                or any(not isinstance(item, TwistConditionCampaign3DResult)
                       for item in conditions)
                or not np.isfinite(self.onset_probability_threshold)
                or not 0.0 < self.onset_probability_threshold <= 1.0
                or int(self.minimum_conditions) != self.minimum_conditions
                or self.minimum_conditions < 2):
            raise ValueError("invalid twist aspect-ratio campaign result")
        aspect_ratio = np.asarray([item.aspect_ratio for item in conditions])
        if np.any(np.diff(aspect_ratio) <= 0.0):
            raise ValueError("twist campaign conditions must have unique increasing aspect ratios")
        object.__setattr__(self, "conditions", conditions)
        object.__setattr__(self, "minimum_conditions", int(self.minimum_conditions))
        object.__setattr__(self, "passed", bool(self.passed))
        object.__setattr__(self, "reasons", tuple(self.reasons))

    @property
    def aspect_ratio(self):
        values = np.asarray([item.aspect_ratio for item in self.conditions], dtype=float)
        values.setflags(write=False)
        return values

    @property
    def twist_probability(self):
        values = np.asarray([item.twist_probability for item in self.conditions], dtype=float)
        values.setflags(write=False)
        return values

    @property
    def mean_maximum_lateral_displacement_m(self):
        values = np.asarray([
            item.mean_maximum_lateral_displacement_m for item in self.conditions], dtype=float)
        values.setflags(write=False)
        return values

    @property
    def twist_onset_aspect_ratio(self):
        selected = np.flatnonzero(self.twist_probability >= self.onset_probability_threshold)
        return None if not len(selected) else float(self.aspect_ratio[selected[0]])

    @property
    def statistical_claim_ready(self):
        # Numerical C5 readiness does not override the validity/evidence of the physical processes.
        return self.passed


def _common_valid(*ensembles):
    shape = ensembles[0].mean_displacement_xy_m.shape
    if any(item.mean_displacement_xy_m.shape != shape for item in ensembles[1:]):
        raise ValueError("twist ensembles must share one physical depth grid")
    valid = np.ones(shape[0], dtype=bool)
    for item in ensembles:
        valid &= item.valid_realization_count == len(item.members)
    if not np.any(valid):
        raise ValueError("twist ensembles share no fully resolved depth slice")
    return valid


def _maximum_vector_change(first, second, valid):
    return float(np.max(np.linalg.norm(
        second[valid] - first[valid], axis=1)))


def _onset_probability(ensemble):
    return float(np.mean(np.isfinite(ensemble.onset_aspect_ratio)))


def assess_twist_ensemble_refinement_3d(
        base_n: FeatureCenterlineEnsemble3D,
        doubled_n: FeatureCenterlineEnsemble3D,
        doubled_sample_level: FeatureCenterlineEnsemble3D, *,
        base_transport_sample_count, refined_transport_sample_count,
        contract: TwistEnsembleRefinementContract3D):
    """Apply the C5 N/sample/isotropy gates to three preregistered ensembles."""
    if (not isinstance(base_n, FeatureCenterlineEnsemble3D)
            or not isinstance(doubled_n, FeatureCenterlineEnsemble3D)
            or not isinstance(doubled_sample_level, FeatureCenterlineEnsemble3D)
            or not isinstance(contract, TwistEnsembleRefinementContract3D)):
        raise TypeError("twist refinement requires three ensembles and one contract")
    if (len(base_n.members) < contract.minimum_realizations
            or len(doubled_n.members) < 2 * len(base_n.members)
            or len(doubled_sample_level.members) != len(doubled_n.members)
            or int(base_transport_sample_count) != base_transport_sample_count
            or int(refined_transport_sample_count) != refined_transport_sample_count
            or base_transport_sample_count <= 0
            or refined_transport_sample_count < 2 * base_transport_sample_count):
        raise ValueError(
            "C5 requires the declared base N, N doubling, and transport-sample doubling")
    valid = _common_valid(base_n, doubled_n, doubled_sample_level)
    n_mean = _maximum_vector_change(
        base_n.mean_displacement_xy_m, doubled_n.mean_displacement_xy_m, valid)
    sample_mean = _maximum_vector_change(
        doubled_n.mean_displacement_xy_m,
        doubled_sample_level.mean_displacement_xy_m, valid)
    n_std = _maximum_vector_change(
        base_n.standard_deviation_displacement_xy_m,
        doubled_n.standard_deviation_displacement_xy_m, valid)
    sample_std = _maximum_vector_change(
        doubled_n.standard_deviation_displacement_xy_m,
        doubled_sample_level.standard_deviation_displacement_xy_m, valid)
    n_onset = abs(_onset_probability(doubled_n) - _onset_probability(base_n))
    sample_onset = abs(
        _onset_probability(doubled_sample_level) - _onset_probability(doubled_n))
    systematic = max(
        doubled_n.maximum_systematic_z_score,
        doubled_sample_level.maximum_systematic_z_score)
    reasons = []
    if n_mean > contract.mean_displacement_tolerance_m:
        reasons.append("ensemble mean is not stable under N doubling")
    if sample_mean > contract.mean_displacement_tolerance_m:
        reasons.append("ensemble mean is not stable under transport-sample doubling")
    if n_std > contract.standard_deviation_tolerance_m:
        reasons.append("twist variance is not stable under N doubling")
    if sample_std > contract.standard_deviation_tolerance_m:
        reasons.append("twist variance is not stable under transport-sample doubling")
    if n_onset > contract.onset_probability_tolerance:
        reasons.append("twist-onset probability is not stable under N doubling")
    if sample_onset > contract.onset_probability_tolerance:
        reasons.append("twist-onset probability is not stable under transport-sample doubling")
    if systematic > contract.maximum_systematic_z_score:
        reasons.append("symmetric control has a statistically resolved systematic twist direction")
    if (not doubled_n.has_nonzero_twist_variance
            or not doubled_sample_level.has_nonzero_twist_variance):
        reasons.append("symmetric control has no resolved stochastic twist variance")
    diagnostics = dict(
        base_realizations=len(base_n.members),
        doubled_realizations=len(doubled_n.members),
        sample_refined_realizations=len(doubled_sample_level.members),
        base_transport_sample_count=int(base_transport_sample_count),
        refined_transport_sample_count=int(refined_transport_sample_count),
        shared_depth_slice_count=int(np.count_nonzero(valid)),
        base_onset_probability=_onset_probability(base_n),
        doubled_n_onset_probability=_onset_probability(doubled_n),
        doubled_sample_onset_probability=_onset_probability(doubled_sample_level),
        exact_claim="statistical twist prediction only; never a deterministic single profile")
    return TwistEnsembleRefinement3DResult(
        passed=not reasons,
        n_doubling_mean_change_m=n_mean,
        sample_doubling_mean_change_m=sample_mean,
        n_doubling_standard_deviation_change_m=n_std,
        sample_doubling_standard_deviation_change_m=sample_std,
        n_doubling_onset_probability_change=n_onset,
        sample_doubling_onset_probability_change=sample_onset,
        maximum_systematic_z_score=systematic,
        reasons=tuple(reasons), diagnostics=diagnostics)


def score_twist_condition_campaign_3d(
        base_n: FeatureCenterlineEnsemble3D,
        doubled_n: FeatureCenterlineEnsemble3D,
        doubled_sample_level: FeatureCenterlineEnsemble3D, *,
        aspect_ratio, base_transport_sample_count, refined_transport_sample_count,
        base_seeds, doubled_seeds, doubled_sample_seeds,
        contract: TwistEnsembleRefinementContract3D):
    """Score one preregistered AR condition without executing or mutating the engine."""
    refinement = assess_twist_ensemble_refinement_3d(
        base_n, doubled_n, doubled_sample_level,
        base_transport_sample_count=base_transport_sample_count,
        refined_transport_sample_count=refined_transport_sample_count,
        contract=contract)
    return TwistConditionCampaign3DResult(
        aspect_ratio, base_n, doubled_n, doubled_sample_level, refinement,
        tuple(base_seeds), tuple(doubled_seeds), tuple(doubled_sample_seeds))


def assess_twist_aspect_ratio_campaign_3d(
        conditions, *, onset_probability_threshold=0.1, minimum_conditions=3):
    """Aggregate independently refined conditions into the C5 numerical AR-sweep gate."""
    conditions = tuple(sorted(conditions, key=lambda item: item.aspect_ratio))
    if (not conditions
            or any(not isinstance(item, TwistConditionCampaign3DResult)
                   for item in conditions)
            or not np.isfinite(onset_probability_threshold)
            or not 0.0 < onset_probability_threshold <= 1.0
            or int(minimum_conditions) != minimum_conditions or minimum_conditions < 2):
        raise ValueError("invalid twist aspect-ratio campaign inputs")
    reasons = []
    if len(conditions) < int(minimum_conditions):
        reasons.append("aspect-ratio sweep has fewer than the preregistered condition count")
    if any(not item.numerical_refinement_passed for item in conditions):
        reasons.append("one or more aspect-ratio conditions failed N/sample/isotropy refinement")
    return TwistAspectRatioCampaign3DResult(
        conditions, float(onset_probability_threshold), int(minimum_conditions),
        not reasons, tuple(reasons))
