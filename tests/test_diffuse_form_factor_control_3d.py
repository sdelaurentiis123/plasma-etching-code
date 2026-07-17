import numpy as np
import pytest

from petch.diffuse_form_factor_control_3d import (
    ReplicatedDiffuseFormFactors3D,
    estimate_replicated_diffuse_form_factors_3d,
    form_factor_reciprocity_diagnostic_3d,
    solve_replicated_diffuse_neutral_radiosity_3d,
)
from petch.neutral_radiosity_3d import DiffuseFormFactors3D


def _two_face(transfer, *, rays=8):
    transfer = float(transfer)
    return DiffuseFormFactors3D(
        2,
        np.asarray([0, 1]),
        np.asarray([1, 0]),
        np.asarray([transfer, transfer]),
        np.asarray([1.0 - transfer, 1.0 - transfer]),
        rays,
    )


def _ensemble(transfers=(0.2, 0.4, 0.6, 0.8)):
    return ReplicatedDiffuseFormFactors3D(
        tuple(_two_face(value) for value in transfers),
        (13, 7, 19, 3),
        np.ones(2),
    )


def test_replicate_mean_is_deterministic_row_closing_and_immutable():
    ensemble = _ensemble()

    assert ensemble.replicate_seeds == (3, 7, 13, 19)
    assert ensemble.rays_per_replicate == 8
    assert ensemble.total_rays_per_face == 32
    assert np.array_equal(ensemble.mean_form_factors.transfer_fraction, [0.5, 0.5])
    assert np.array_equal(ensemble.mean_form_factors.escape_fraction, [0.5, 0.5])
    assert not ensemble.face_area_m2.flags.writeable
    assert len(ensemble.sha256) == 64
    outgoing = ensemble.mean_form_factors.escape_fraction + np.bincount(
        ensemble.mean_form_factors.source_face,
        weights=ensemble.mean_form_factors.transfer_fraction,
        minlength=2,
    )
    assert np.allclose(outgoing, 1.0, rtol=0.0, atol=5e-13)


def test_reciprocity_diagnostic_measures_without_projecting():
    factors = DiffuseFormFactors3D(
        2,
        np.asarray([0, 1]),
        np.asarray([1, 0]),
        np.asarray([0.5, 0.25]),
        np.asarray([0.5, 0.75]),
        8,
    )
    exact = form_factor_reciprocity_diagnostic_3d(factors, [1.0, 2.0])
    broken = form_factor_reciprocity_diagnostic_3d(factors, [1.0, 1.0])

    assert exact.relative_l1_error == 0.0
    assert exact.one_sided_pair_count == 0
    assert broken.relative_l1_error == pytest.approx(0.5)
    assert np.array_equal(factors.transfer_fraction, [0.5, 0.25])


def test_replicated_radiosity_reports_zero_uncertainty_for_identical_operators():
    ensemble = _ensemble((0.5, 0.5, 0.5, 0.5))
    result = solve_replicated_diffuse_neutral_radiosity_3d(
        [1.0, 0.2], ensemble, [0.3, 0.4],
        relative_tolerance=1.0e-12, maximum_iterations=2000)

    assert np.array_equal(result.incident_standard_error_m2_s, [0.0, 0.0])
    assert np.array_equal(
        result.incident_confidence_half_width_m2_s, [0.0, 0.0])
    assert result.incident_relative_confidence_linf == 0.0
    assert result.authority.relative_balance_error <= 2.0e-12
    assert not result.incident_confidence_half_width_m2_s.flags.writeable


def test_replicated_radiosity_exposes_nonzero_sampling_uncertainty():
    result = solve_replicated_diffuse_neutral_radiosity_3d(
        [1.0, 0.2], _ensemble(), [0.3, 0.4],
        relative_tolerance=1.0e-12, maximum_iterations=2000)

    assert np.any(result.incident_standard_error_m2_s > 0.0)
    assert result.incident_relative_confidence_linf > 0.0
    assert result.incident_area_weighted_relative_confidence_l1 > 0.0
    assert result.authority_to_replicate_mean_incident_relative_linf > 0.0
    assert result.authority.relative_balance_error <= 2.0e-12


def test_open_plane_estimator_replicates_all_escape():
    verts = np.asarray([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ])
    faces = np.asarray([[0, 1, 2]])
    centroids = verts[faces].mean(axis=1)
    ensemble = estimate_replicated_diffuse_form_factors_3d(
        verts, faces, centroids, [[0.0, 0.0, 1.0]], [0.5],
        rays_per_replicate=8, replicate_seeds=(11, 13, 17, 19),
        domain_size=(1.0, 1.0, 1.0), device="cpu")

    assert ensemble.mean_form_factors.source_face.size == 0
    assert np.array_equal(ensemble.mean_form_factors.escape_fraction, [1.0])
    assert ensemble.mean_reciprocity.relative_l1_error == 0.0
    assert ensemble.source_sampling == "triangle_area"
    assert ensemble.construction_identity["sampling_dimension"] == 4
    assert ensemble.construction_identity["visibility_certification"] == (
        "full_event_float64_parity_14664_ray_real_checkpoint")
    assert ensemble.construction_identity["visibility_operator"] == "cellwise_certified"
    assert ensemble.construction_identity["visibility_float64_evaluated_count"] == 0
    assert ensemble.construction_identity["visibility_open_escape_count"] == 32
    assert ensemble.construction_identity["construction_call_count"] == 4


@pytest.mark.parametrize(
    "factors,seeds,match",
    [
        ((_two_face(0.5),) * 3, (1, 2, 3), "at least 4"),
        ((_two_face(0.5),) * 4, (1, 1, 2, 3), "distinct"),
        ((_two_face(0.5),) * 3 + (_two_face(0.5, rays=16),),
         (1, 2, 3, 4), "common nested level"),
    ],
)
def test_replicate_contract_refuses_underpowered_or_incompatible_inputs(
        factors, seeds, match):
    with pytest.raises(ValueError, match=match):
        ReplicatedDiffuseFormFactors3D(factors, seeds, np.ones(2))


def test_estimator_owns_seed_and_ray_controls():
    with pytest.raises(ValueError, match="controlled"):
        estimate_replicated_diffuse_form_factors_3d(
            np.zeros((3, 3)), np.asarray([[0, 1, 2]]), np.zeros((1, 3)),
            [[0.0, 0.0, 1.0]], [1.0], replicate_seeds=(1, 2, 3, 4),
            seed=9)
