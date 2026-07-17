import numpy as np
import pytest

from petch.surface_patch_convergence_3d import (
    aggregate_surface_field_on_physical_patches_3d,
    compare_physical_patch_fields_3d,
    physical_surface_patch_groups_3d,
    physical_surface_patch_keys_3d,
    score_replicated_surface_field_at_physical_scales_3d,
    score_replicated_surface_field_on_physical_patches_3d,
    score_surface_field_refinement_at_physical_scales_3d,
)


def _plane(cell_count):
    x = np.linspace(0.0, 1.0, cell_count + 1)
    vertices = np.asarray([
        [x_value, y_value, 0.0] for x_value in x for y_value in x])
    faces = []
    stride = cell_count + 1
    for ix in range(cell_count):
        for iy in range(cell_count):
            lower = ix * stride + iy
            faces.extend((
                [lower, lower + stride, lower + stride + 1],
                [lower, lower + stride + 1, lower + 1],
            ))
    faces = np.asarray(faces, dtype=int)
    area = np.full(len(faces), 0.5 / cell_count ** 2)
    normal = np.tile([0.0, 0.0, 1.0], (len(faces), 1))
    return vertices, faces, area, normal


def _surface(cell_count, field):
    vertices, faces, area, normal = _plane(cell_count)
    return dict(
        face_field=np.broadcast_to(float(field), len(faces)).copy(),
        face_area_m2=area, verts=vertices, faces=faces,
        face_gas_normals=normal, face_material_id=np.ones(len(faces), dtype=int),
        mesh_length_unit_m=1.0, mesh_origin_m=(0.0, 0.0, 0.0),
        patch_origin_m=(0.0, 0.0, 0.0))


def _disconnected_rectangles(widths):
    """Axis-aligned unit-height rectangles, one in each successive unit patch."""
    vertices = []
    faces = []
    area = []
    for patch_index, width in enumerate(widths):
        lower = float(patch_index)
        first = len(vertices)
        vertices.extend((
            [lower, 0.0, 0.0], [lower + float(width), 0.0, 0.0],
            [lower + float(width), 1.0, 0.0], [lower, 1.0, 0.0],
        ))
        faces.extend((
            [first, first + 1, first + 2],
            [first, first + 2, first + 3],
        ))
        area.extend((0.5 * float(width), 0.5 * float(width)))
    count = len(faces)
    return (
        np.asarray(vertices, dtype=float), np.asarray(faces, dtype=int),
        np.asarray(area, dtype=float),
        np.tile([0.0, 0.0, 1.0], (count, 1)),
    )


def test_triangle_overlap_patch_integral_is_subdivision_invariant():
    coarse = _surface(1, 3.0)
    fine = _surface(4, 3.0)

    result = score_surface_field_refinement_at_physical_scales_3d(
        coarse, fine, (0.5, 1.0), absolute_tolerance=1e-12,
        relative_tolerance=0.05)

    assert len(result) == 2
    assert all(item.maximum_patch_area_relative_error < 3e-15 for item in result)
    assert all(item.integrated_absolute_linf < 5e-16 for item in result)
    assert all(item.mean_absolute_linf < 2e-15 for item in result)
    assert all(item.integrated_mixed_normalized_linf < 1e-2 for item in result)
    assert all(item.mean_mixed_normalized_linf < 1e-2 for item in result)

    coarse_receipt = aggregate_surface_field_on_physical_patches_3d(
        patch_scale_m=0.5, **coarse)
    fine_receipt = aggregate_surface_field_on_physical_patches_3d(
        patch_scale_m=0.5, **fine)
    assert np.allclose(coarse_receipt.patch_projected_support_fraction, 1.0)
    assert np.allclose(
        coarse_receipt.patch_projected_support_area_m2,
        fine_receipt.patch_projected_support_area_m2,
        rtol=0.0, atol=5e-16)


def test_dominant_axis_projection_removes_tilt_area_inflation():
    vertices = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.5],
        [1.0, 1.0, 0.5], [0.0, 1.0, 0.0],
    ])
    faces = np.asarray([[0, 1, 2], [0, 2, 3]], dtype=int)
    face_area = np.full(2, np.sqrt(1.25) / 2.0)
    normal = np.tile([-1.0 / np.sqrt(5.0), 0.0, 2.0 / np.sqrt(5.0)], (2, 1))
    receipt = aggregate_surface_field_on_physical_patches_3d(
        np.ones(2), face_area, vertices, faces, normal, np.ones(2, dtype=int),
        1.0, mesh_length_unit_m=1.0)

    assert receipt.patch_area_m2[0] == pytest.approx(np.sqrt(1.25))
    assert receipt.patch_projected_support_area_m2[0] == pytest.approx(1.0)
    assert receipt.patch_nominal_projected_area_m2[0] == pytest.approx(1.0)
    assert receipt.patch_projected_support_fraction[0] == pytest.approx(1.0)


def test_support_receipt_is_translation_and_patch_origin_stable():
    base = _surface(4, 3.0)
    translated = dict(base)
    translated["mesh_origin_m"] = (7.0, -3.0, 2.0)
    translated["patch_origin_m"] = (7.0, -3.0, 2.0)
    reference = aggregate_surface_field_on_physical_patches_3d(
        patch_scale_m=0.5, **base)
    candidate = aggregate_surface_field_on_physical_patches_3d(
        patch_scale_m=0.5, **translated)

    assert np.array_equal(reference.patch_key, candidate.patch_key)
    assert np.allclose(
        reference.patch_projected_support_fraction,
        candidate.patch_projected_support_fraction)
    assert np.sum(candidate.integrated_field_area) == pytest.approx(3.0)

    shifted_grid = dict(base)
    shifted_grid["patch_origin_m"] = (-0.25, -0.25, 0.0)
    shifted = aggregate_surface_field_on_physical_patches_3d(
        patch_scale_m=0.5, **shifted_grid)
    assert np.sum(shifted.patch_area_m2) == pytest.approx(1.0)
    assert np.sum(shifted.integrated_field_area) == pytest.approx(3.0)
    assert np.min(shifted.patch_projected_support_fraction) == pytest.approx(0.25)


def test_periodic_domain_normalizes_represented_footprint_and_shifted_origin():
    vertices = np.asarray([
        [0.0, 0.0, 0.0], [2.0, 0.0, 0.0],
        [2.0, 1.0, 0.0], [0.0, 1.0, 0.0],
    ])
    faces = np.asarray([[0, 1, 2], [0, 2, 3]], dtype=int)
    area = np.ones(2)
    normal = np.tile([0.0, 0.0, 1.0], (2, 1))
    common = dict(
        face_field=np.ones(2), face_area_m2=area, verts=vertices, faces=faces,
        face_gas_normals=normal, face_material_id=np.ones(2, dtype=int),
        patch_scale_m=2.0, mesh_length_unit_m=1.0)
    uncorrected = aggregate_surface_field_on_physical_patches_3d(**common)
    periodic = aggregate_surface_field_on_physical_patches_3d(
        **common, periodic_domain_lengths_m=(0.0, 1.0, 0.0))
    split_periodic = aggregate_surface_field_on_physical_patches_3d(
        **common, patch_origin_m=(0.0, 0.5, 0.0),
        periodic_domain_lengths_m=(0.0, 1.0, 0.0),
        periodic_domain_origin_m=(0.0, 0.0, 0.0))

    assert uncorrected.patch_projected_support_fraction.tolist() == pytest.approx([0.5])
    assert periodic.patch_projected_support_fraction.tolist() == pytest.approx([1.0])
    assert split_periodic.patch_projected_support_fraction.tolist() == pytest.approx([1.0, 1.0])
    assert np.sum(split_periodic.integrated_field_area) == pytest.approx(2.0)


def test_partial_patch_mean_is_excluded_but_integrated_inventory_still_gates():
    vertices, faces, area, normal = _disconnected_rectangles((1.0, 0.05))
    authority = np.zeros(len(faces))
    replicates = np.zeros((4, len(faces)))
    replicates[:, 2:] = 100.0
    score = score_replicated_surface_field_on_physical_patches_3d(
        authority, replicates, area, vertices, faces, normal,
        np.ones(len(faces), dtype=int), 1.0,
        absolute_tolerance=1.0, relative_tolerance=0.0,
        mesh_length_unit_m=1.0, minimum_mean_support_fraction=0.1)

    assert score.authority.patch_projected_support_fraction.tolist() == pytest.approx([1.0, 0.05])
    assert len(score.authority.patch_key) == 2
    assert np.sum(score.replicate_mean_integrated_field_area) == pytest.approx(5.0)
    assert score.eligible_mean_patch_count == 1
    assert score.excluded_mean_patch_count == 1
    assert score.excluded_mean_surface_area_fraction == pytest.approx(0.05 / 1.05)
    assert score.excluded_mean_projected_support_fraction == pytest.approx(0.05 / 1.05)
    assert score.maximum_mean_combined_mixed_normalized == 0.0
    assert score.maximum_mean_combined_all_patches_mixed_normalized == pytest.approx(100.0)
    # The excluded sliver is never removed from conservation or the fixed-area integrated gate.
    assert score.maximum_integrated_combined_mixed_normalized == pytest.approx(5.0)
    assert not score.all_mixed_tolerances_pass


def test_mean_support_sensitivity_receipts_are_monotone_before_default_choice():
    vertices, faces, area, normal = _disconnected_rectangles((0.2, 0.4, 0.6, 0.8))
    common = dict(
        authority_face_field=np.zeros(len(faces)),
        replicate_face_fields=np.zeros((4, len(faces))),
        face_area_m2=area, verts=vertices, faces=faces,
        face_gas_normals=normal, face_material_id=np.ones(len(faces), dtype=int),
        patch_scale_m=1.0, absolute_tolerance=1.0, relative_tolerance=0.0,
        mesh_length_unit_m=1.0)
    scores = [
        score_replicated_surface_field_on_physical_patches_3d(
            **common, minimum_mean_support_fraction=threshold)
        for threshold in (0.25, 0.5, 0.75)
    ]

    assert [item.eligible_mean_patch_count for item in scores] == [3, 2, 1]
    assert [item.excluded_mean_patch_count for item in scores] == [1, 2, 3]
    assert [item.excluded_mean_surface_area_fraction for item in scores] == pytest.approx(
        [0.1, 0.3, 0.6])
    assert [item.excluded_mean_projected_support_fraction for item in scores] == pytest.approx(
        [0.1, 0.3, 0.6])


def test_mean_gate_refuses_when_every_patch_is_an_unresolved_sliver():
    vertices, faces, area, normal = _disconnected_rectangles((0.01, 0.02))
    with pytest.raises(ValueError, match="no physical patch"):
        score_replicated_surface_field_on_physical_patches_3d(
            np.zeros(len(faces)), np.zeros((4, len(faces))), area,
            vertices, faces, normal, np.ones(len(faces), dtype=int), 1.0,
            absolute_tolerance=1.0, relative_tolerance=0.0,
            mesh_length_unit_m=1.0, minimum_mean_support_fraction=0.1)


def test_patch_overlap_closes_exact_face_area_and_integrated_inventory():
    vertices, faces, area, normal = _plane(2)
    field = np.arange(1, len(faces) + 1, dtype=float)
    receipt = aggregate_surface_field_on_physical_patches_3d(
        field, area, vertices, faces, normal, np.ones(len(faces), dtype=int),
        0.3, mesh_length_unit_m=1.0)

    per_face = np.bincount(
        receipt.contribution_face_index, weights=receipt.contribution_area_m2,
        minlength=len(faces))
    assert np.allclose(per_face, area, rtol=0.0, atol=5e-17)
    assert np.sum(receipt.patch_area_m2) == pytest.approx(np.sum(area), abs=2e-16)
    assert np.sum(receipt.integrated_field_area) == pytest.approx(
        np.sum(area * field), abs=2e-15)
    assert len(receipt.scheme_sha256) == 64
    assert not receipt.contribution_area_m2.flags.writeable


def test_exact_patch_boundary_has_no_roundoff_sliver_and_is_ulp_stable():
    boundary = 0.5
    offsets = (
        np.nextafter(boundary, -np.inf) - boundary,
        0.0,
        np.nextafter(boundary, np.inf) - boundary,
    )
    receipts = []
    for offset in offsets:
        vertices = np.asarray([
            [boundary + offset, 0.1, 0.1],
            [boundary + offset, 0.5, 0.1],
            [boundary + offset, 0.1, 0.5],
        ])
        faces = np.asarray([[0, 1, 2]], dtype=int)
        area = np.asarray([0.08])
        receipt = aggregate_surface_field_on_physical_patches_3d(
            np.asarray([3.0]), area, vertices, faces,
            np.asarray([[1.0, 0.0, 0.0]]), np.asarray([1]), 0.5,
            mesh_length_unit_m=1.0)
        receipts.append(receipt)

        assert len(receipt.patch_key) == 1
        assert receipt.patch_area_m2[0] == pytest.approx(area[0], abs=2e-17)
        assert receipt.contribution_area_m2[0] == pytest.approx(area[0], abs=2e-17)
        assert receipt.integrated_field_area[0] == pytest.approx(3.0 * area[0])

    assert np.array_equal(receipts[0].patch_key, receipts[1].patch_key)
    assert np.array_equal(receipts[1].patch_key, receipts[2].patch_key)
    assert np.allclose(
        receipts[0].contribution_area_m2,
        receipts[2].contribution_area_m2,
        rtol=0.0, atol=2e-17)


def test_float32_mesh_quantization_does_not_create_a_microscopic_patch():
    boundary = np.float32(0.5)
    just_above = np.nextafter(boundary, np.float32(np.inf))
    vertices = np.asarray([
        [0.0, 0.0, 0.0],
        [just_above, 0.0, 0.0],
        [0.0, boundary, 0.0],
    ], dtype=np.float32)
    faces = np.asarray([[0, 1, 2]], dtype=int)
    area = np.asarray([
        0.5 * np.linalg.norm(np.cross(
            vertices[1].astype(float) - vertices[0],
            vertices[2].astype(float) - vertices[0]))])
    receipt = aggregate_surface_field_on_physical_patches_3d(
        np.asarray([1.0]), area, vertices, faces,
        np.asarray([[0.0, 0.0, 1.0]]), np.asarray([1]), 0.5,
        mesh_length_unit_m=1.0)

    assert len(receipt.patch_key) == 1
    assert receipt.patch_key[0, 3:].tolist() == [0, 0, 0]
    assert receipt.patch_area_m2[0] == pytest.approx(area[0], abs=2e-16)


def test_mixed_error_uses_absolute_scale_near_zero_and_reports_worst_face():
    reference = aggregate_surface_field_on_physical_patches_3d(
        patch_scale_m=0.5, **_surface(2, 0.0))
    candidate = aggregate_surface_field_on_physical_patches_3d(
        patch_scale_m=0.5, **_surface(2, 1e-10))
    score = compare_physical_patch_fields_3d(
        reference, candidate, absolute_tolerance=1e-9,
        relative_tolerance=0.05)

    assert score.mean_mixed_normalized_linf == pytest.approx(0.1 / 1.005)
    assert score.integrated_mixed_normalized_linf == pytest.approx(0.1 / 1.005)
    assert score.reference_maximum_absolute_face_value == 0.0
    assert score.candidate_maximum_absolute_face_value == pytest.approx(1e-10)


def test_patch_keys_are_deterministic_and_separate_oriented_surfaces():
    centroid = np.asarray([[0.25, 0.25, 0.25], [0.25, 0.25, 0.25]])
    normal = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    material = np.ones(2, dtype=int)
    key, group = physical_surface_patch_keys_3d(
        centroid, normal, material, 1.0, mesh_length_unit_m=1.0)
    replay = physical_surface_patch_groups_3d(
        centroid, normal, material, 1.0, mesh_length_unit_m=1.0)

    assert len(key) == 2
    assert group[0] != group[1]
    assert np.array_equal(group, replay)
    assert np.array_equal(key[:, :3], [[1, 0, 1], [1, 2, 1]])


def test_comparison_refuses_patch_coverage_or_contract_mismatch():
    reference = aggregate_surface_field_on_physical_patches_3d(
        patch_scale_m=0.5, **_surface(1, 1.0))
    shifted = _surface(1, 1.0)
    shifted["mesh_origin_m"] = (1.0, 0.0, 0.0)
    candidate = aggregate_surface_field_on_physical_patches_3d(
        patch_scale_m=0.5, **shifted)

    with pytest.raises(ValueError, match="patch keys differ"):
        compare_physical_patch_fields_3d(
            reference, candidate, absolute_tolerance=1e-9,
            relative_tolerance=0.05)
    with pytest.raises(ValueError, match="at least two"):
        score_surface_field_refinement_at_physical_scales_3d(
            _surface(1, 1.0), _surface(2, 1.0), (0.5,),
            absolute_tolerance=1e-9, relative_tolerance=0.05)
    with pytest.raises(ValueError, match="positive absolute"):
        compare_physical_patch_fields_3d(
            reference, reference, absolute_tolerance=0.0,
            relative_tolerance=0.05)


def test_replicated_patch_score_separates_interval_from_mean_operator_bias():
    vertices, faces, area, normal = _plane(2)
    material = np.ones(len(faces), dtype=int)
    authority = np.full(len(faces), 2.0)
    replicates = np.stack([
        np.full(len(faces), value) for value in (1.9, 2.1, 1.95, 2.05)])
    score = score_replicated_surface_field_on_physical_patches_3d(
        authority, replicates, area, vertices, faces, normal, material, 0.5,
        absolute_tolerance=0.5, relative_tolerance=0.0,
        mesh_length_unit_m=1.0)

    assert score.replicate_count == 4
    assert score.minimum_mean_support_fraction == pytest.approx(0.1)
    assert score.maximum_mean_confidence_mixed_normalized > 0.0
    assert score.maximum_mean_authority_bias_mixed_normalized < 1e-14
    assert score.maximum_mean_combined_mixed_normalized < 1.0
    assert score.all_mixed_tolerances_pass
    assert not score.mean_confidence_half_width.flags.writeable

    biased = score_replicated_surface_field_on_physical_patches_3d(
        authority, np.full((4, len(faces)), 4.0), area, vertices, faces,
        normal, material, 0.5, absolute_tolerance=0.1,
        relative_tolerance=0.0, mesh_length_unit_m=1.0)
    assert biased.maximum_mean_confidence_mixed_normalized == 0.0
    assert biased.maximum_mean_authority_bias_mixed_normalized == pytest.approx(20.0)
    assert not biased.all_mixed_tolerances_pass


def test_replicated_patch_score_requires_independence_and_two_physical_scales():
    vertices, faces, area, normal = _plane(1)
    common = dict(
        authority_face_field=np.ones(len(faces)),
        face_area_m2=area, verts=vertices, faces=faces,
        face_gas_normals=normal, face_material_id=np.ones(len(faces), dtype=int),
        absolute_tolerance=0.1, relative_tolerance=0.05,
        mesh_length_unit_m=1.0)
    with pytest.raises(ValueError, match="at least four"):
        score_replicated_surface_field_on_physical_patches_3d(
            replicate_face_fields=np.ones((3, len(faces))), patch_scale_m=0.5,
            **common)
    with pytest.raises(ValueError, match="at least two"):
        score_replicated_surface_field_at_physical_scales_3d(
            replicate_face_fields=np.ones((4, len(faces))), patch_scales_m=(0.5,),
            **common)
