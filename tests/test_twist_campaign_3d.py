import numpy as np

from petch.feature_step_3d import FeatureGeometry3D
from petch.profile_observables_3d import measure_feature_centerline_ensemble_3d
from petch.twist_campaign_3d import (
    TwistEnsembleRefinementContract3D,
    assess_twist_aspect_ratio_campaign_3d,
    assess_twist_ensemble_refinement_3d,
    score_twist_condition_campaign_3d,
)


def _tilted_hole(slope):
    dx = 0.1
    shape = (21, 21, 11)
    x, y, z = (np.arange(size) * dx for size in shape)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    center_x = 1.0 + float(slope) * (1.0 - Z)
    phi = np.sqrt((X - center_x) ** 2 + (Y - 1.0) ** 2) - 0.3
    return FeatureGeometry3D(phi, np.where(phi > 0.0, 1, 0), dx, 1e-6)


def _ensemble(slopes):
    return measure_feature_centerline_ensemble_3d(
        tuple(_tilted_hole(value) for value in slopes),
        lateral_bounds_m=(0.4e-6, 1.6e-6, 0.4e-6, 1.6e-6),
        opening_width_m=0.6e-6, onset_displacement_m=0.08e-6,
        reference_slice_count=1)


def test_twist_refinement_contract_passes_nested_symmetric_manufactured_ensembles():
    base = _ensemble((-0.2, -0.1, 0.1, 0.2))
    doubled = _ensemble((-0.2, -0.1, 0.1, 0.2) * 2)
    refined = _ensemble((-0.2, -0.1, 0.1, 0.2) * 2)
    contract = TwistEnsembleRefinementContract3D(
        minimum_realizations=4,
        mean_displacement_tolerance_m=0.02e-6,
        standard_deviation_tolerance_m=0.02e-6,
        onset_probability_tolerance=0.1,
        maximum_systematic_z_score=3.0)

    result = assess_twist_ensemble_refinement_3d(
        base, doubled, refined, base_transport_sample_count=128,
        refined_transport_sample_count=256, contract=contract)

    assert result.passed, result.reasons
    assert result.diagnostics["doubled_realizations"] == 8
    assert result.diagnostics["exact_claim"].startswith("statistical")


def test_twist_refinement_contract_rejects_systematic_direction():
    base = _ensemble((0.1, 0.2, 0.1, 0.2))
    doubled = _ensemble((0.1, 0.2, 0.1, 0.2) * 2)
    contract = TwistEnsembleRefinementContract3D(
        minimum_realizations=4,
        mean_displacement_tolerance_m=1e-6,
        standard_deviation_tolerance_m=1e-6,
        onset_probability_tolerance=1.0,
        maximum_systematic_z_score=1.0)

    result = assess_twist_ensemble_refinement_3d(
        base, doubled, doubled, base_transport_sample_count=128,
        refined_transport_sample_count=256, contract=contract)

    assert not result.passed
    assert "systematic twist direction" in " ".join(result.reasons)


def test_twist_condition_and_ar_sweep_preserve_nested_paired_seeds():
    base = _ensemble((-0.2, -0.1, 0.1, 0.2))
    doubled = _ensemble((-0.2, -0.1, 0.1, 0.2) * 2)
    contract = TwistEnsembleRefinementContract3D(
        minimum_realizations=4,
        mean_displacement_tolerance_m=0.02e-6,
        standard_deviation_tolerance_m=0.02e-6,
        onset_probability_tolerance=0.1,
        maximum_systematic_z_score=3.0)
    seeds = tuple(range(8))
    conditions = tuple(score_twist_condition_campaign_3d(
        base, doubled, doubled, aspect_ratio=aspect_ratio,
        base_transport_sample_count=128, refined_transport_sample_count=256,
        base_seeds=seeds[:4], doubled_seeds=seeds,
        doubled_sample_seeds=seeds, contract=contract)
        for aspect_ratio in (4.0, 8.0, 12.0))

    campaign = assess_twist_aspect_ratio_campaign_3d(
        conditions, onset_probability_threshold=0.5, minimum_conditions=3)

    assert campaign.passed
    assert campaign.statistical_claim_ready
    assert campaign.twist_onset_aspect_ratio == 4.0
    assert np.array_equal(campaign.aspect_ratio, [4.0, 8.0, 12.0])
