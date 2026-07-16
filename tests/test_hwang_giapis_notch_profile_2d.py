import numpy as np
import pytest

from petch.boundary_transport_3d import _trace_field_events_float64_3d
from petch.hwang_giapis_notch_profile_2d import (
    HwangGiapisLocalIonEntry2D,
    HwangGiapisLocalNotchCheckpoint2D,
    _detached_poly_cells,
    _five_cell_surface_reaction_normals,
    _local_entry_origins_3d,
    _local_surface_mesh_2d,
    evolve_hwang_giapis_local_notch_2d,
    evolve_hwang_giapis_local_notch_event_driven_2d,
    hwang_giapis_exposed_oxide_potential_v,
    hwang_giapis_local_boundary_from_edge_array_result,
    solve_hwang_giapis_local_laplace_2d,
)


def test_detached_fragment_mask_requires_finite_area_connection_to_bulk():
    poly = np.ones((6, 6), dtype=bool)
    poly[0:4, 1] = False
    poly[1:4, 0] = False
    poly[0, 0] = False
    poly[1, 2] = False
    poly[1, 1] = True

    detached = _detached_poly_cells(poly)

    assert np.argwhere(detached).tolist() == [[1, 1]]
    assert not detached[4, 1]
    assert not np.any(detached[-1])
    assert not np.any(detached[:, -1])


def test_checkpoint_refuses_floating_material_and_reactive_ledger_loss():
    poly = np.ones((6, 6), dtype=bool)
    poly[0:4, 1] = False
    poly[0:4, 0] = False
    poly[1, 2] = False
    poly[1, 1] = True
    inventory = np.zeros_like(poly, dtype=float)
    oxide = np.zeros(poly.shape[0])
    with pytest.raises(ValueError, match="invalid event-driven"):
        HwangGiapisLocalNotchCheckpoint2D(
            poly, inventory, oxide, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0,
            0, 0, 0.0, 0.0, 0)

    connected = np.ones((6, 6), dtype=bool)
    with pytest.raises(ValueError, match="invalid event-driven"):
        HwangGiapisLocalNotchCheckpoint2D(
            connected, inventory, oxide, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0,
            0, 0, 0.0, 0.0, 0)


def test_five_cell_reaction_normal_recovers_vertical_and_linear_fronts():
    vertical = np.ones((10, 7), dtype=bool)
    assert np.allclose(
        _five_cell_surface_reaction_normals(vertical, 0.005),
        np.tile([-1.0, 0.0, 0.0], (7, 1)))

    diagonal = np.ones((10, 7), dtype=bool)
    for height in range(7):
        diagonal[:height, height] = False
    expected = np.asarray([-1.0, 0.0, 1.0]) / np.sqrt(2.0)
    assert np.allclose(
        _five_cell_surface_reaction_normals(diagonal, 0.005),
        np.tile(expected, (7, 1)))


def test_local_entries_start_in_gas_and_intact_sidewall_blocks_tunneling():
    dx = 0.005
    poly = np.ones((8, 6), dtype=bool)
    poly[0, 1] = False
    local = solve_hwang_giapis_local_laplace_2d(
        poly, np.zeros(6), np.zeros(8), poly_potential_v=0.0,
        cell_size_um=dx)
    mesh = _local_surface_mesh_2d(poly, dx)
    origin = _local_entry_origins_3d(
        np.array([0.5 * dx, 1.5 * dx]), dx, 6 * dx)
    velocity = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

    hit_face, _cosine, _energy, termination, *_ = (
        _trace_field_events_float64_3d(
            origin, velocity, 1.0, local.potential_v,
            np.asarray(local.origin_um), np.asarray(local.spacing_um),
            mesh.vertices_um, mesh.faces, 2.5e-4, 8192, False))

    assert np.all(origin[:, 0] < 0.0)
    assert np.array_equal(termination, [1, 1])
    # The intact lower cell stops at x=0.  The opened height crosses that
    # plane and stops on the next poly-Si cell at x=dx.
    assert mesh.cell_x[hit_face[0]] == 0
    assert mesh.cell_x[hit_face[1]] == 1


def test_hwang_giapis_eq_4_4_replays_printed_boundary_condition():
    result = hwang_giapis_exposed_oxide_potential_v(
        [100.0, 50.0, 25.0], poly_potential_v=8.0,
        initial_sidewall_floor_potential_v=20.0)

    # V0=(20+8)/2=14 V; V=Vp+N/N1*(V0-Vp).
    assert np.allclose(result, [14.0, 11.0, 9.5])
    assert np.allclose(
        hwang_giapis_exposed_oxide_potential_v(
            [0.0, 0.0], poly_potential_v=8.0,
            initial_sidewall_floor_potential_v=20.0),
        8.0)


def test_local_laplace_preserves_constant_dirichlet_solution():
    poly = np.ones((8, 6), dtype=bool)
    poly[:4, :4] = False
    result = solve_hwang_giapis_local_laplace_2d(
        poly, np.full(6, 7.5), np.ones(8),
        poly_potential_v=7.5, cell_size_um=0.005)

    assert np.allclose(result.potential_v, 7.5)
    assert np.allclose(result.gas_cell_potential_v, 7.5)
    assert np.allclose(result.oxide_potential_v, 7.5)


def test_local_laplace_refuses_refinement_boundary_breakout():
    poly = np.ones((8, 6), dtype=bool)
    poly[-1, 0] = False
    with pytest.raises(ValueError, match="refinement boundary"):
        solve_hwang_giapis_local_laplace_2d(
            poly, np.zeros(6), np.zeros(8), poly_potential_v=0.0)


def test_upper_undercut_exposes_an_inert_photoresist_ceiling():
    poly = np.ones((8, 6), dtype=bool)
    poly[0, -1] = False

    local = solve_hwang_giapis_local_laplace_2d(
        poly, np.full(6, 7.5), np.ones(8),
        poly_potential_v=7.5, cell_size_um=0.005)
    mesh = _local_surface_mesh_2d(poly, 0.005)
    photoresist = mesh.material_id == 2

    assert np.allclose(local.potential_v, 7.5)
    assert np.count_nonzero(photoresist) == 2
    assert np.allclose(mesh.gas_normal[photoresist], [0.0, 0.0, -1.0])
    assert np.all(mesh.cell_x[photoresist] == 0)
    assert np.all(mesh.cell_z[photoresist] == -2)


def test_event_driven_upper_undercut_continues_and_ledgers_photoresist_hits():
    entries = HwangGiapisLocalIonEntry2D(
        height_um=np.array([0.0275]),
        velocity_xz_sqrt_eV=np.array([
            [np.sqrt(177.0), np.sqrt(177.0)]]),
        expected_count=np.array([100.0]))

    result = evolve_hwang_giapis_local_notch_event_driven_2d(
        entries, np.zeros(6), poly_potential_v=0.0,
        cell_size_um=0.005, line_width_um=0.04,
        poly_thickness_um=0.03, include_forward_scatter=False)

    assert result.notch_depth_by_height_um[-1] == pytest.approx(0.005)
    assert result.landed_photoresist_count > 0.0
    assert result.landed_poly_count + result.landed_photoresist_count == (
        pytest.approx(result.launched_count))


def test_edge_array_lineage_reduces_to_weighted_local_boundary_without_profile_fit():
    nz = 20
    potential = np.zeros((10, nz))
    potential[5] = np.linspace(12.0, 30.0, nz)
    vx = np.array([-3.0, -4.0, 0.0, 0.0])
    vz = np.array([4.0, 3.0, 0.0, 0.0])
    result = {
        "V": potential,
        "V_poly_edge": 7.8,
        "geom": {"nz": nz, "trench0": 5, "poly_cells": 15},
        "final_ion_lineage": {
            "hit_type": np.array([4, 4, 1, 0]),
            "hit_z_grid": np.array([18.0, 17.0, 19.0, 0.0]),
            "hit_vx_sqrt_eV": vx,
            "hit_vz_sqrt_eV": vz,
            "impact_energy_eV": vx * vx + vz * vz,
            "source_particle_count": 4,
            "cell_size_um": 0.02,
            "source_width_um": 1.0,
            "edge_inner_poly_hit_type": 4,
        },
    }

    boundary = hwang_giapis_local_boundary_from_edge_array_result(result)

    assert np.allclose(boundary.entries.height_um, [0.02, 0.04])
    assert np.allclose(
        boundary.entries.velocity_xz_sqrt_eV, [[3.0, -4.0], [4.0, -3.0]])
    assert np.allclose(boundary.entries.expected_count, 9.35e6)
    assert boundary.poly_potential_v == 7.8
    assert boundary.sidewall_potential_v.shape == (60,)
    assert boundary.target_event_count == 2


def test_local_notch_removes_only_bombarded_height_band_and_closes_count_ledger():
    entries = HwangGiapisLocalIonEntry2D(
        height_um=np.array([0.04, 0.05]),
        velocity_xz_sqrt_eV=np.array([[np.sqrt(177.0), 0.0],
                                      [np.sqrt(177.0), 0.0]]),
        expected_count=np.array([70.0, 70.0]))
    result = evolve_hwang_giapis_local_notch_2d(
        entries, np.zeros(60), poly_potential_v=0.0,
        batches=8, include_forward_scatter=False)

    assert result.removed_cell_count > 0
    assert result.removed_cell_count == (
        result.threshold_removed_cell_count + result.detached_cell_count)
    assert result.maximum_notch_depth_um >= 0.005
    assert np.all(result.notch_depth_by_height_um[20:] == 0.0)
    assert result.scattered_reactive_collisions == 0.0
    assert (
        result.landed_poly_count + result.landed_oxide_count
        + result.landed_photoresist_count + result.escaped_count
        == pytest.approx(result.launched_count))
    assert (
        result.threshold_removed_reactive_collisions
        + result.detached_reactive_collisions
        + np.sum(result.reactive_collision_inventory)
        == pytest.approx(
            result.direct_reactive_collisions
            + result.scattered_reactive_collisions))


def test_event_driven_local_notch_uses_the_full_fluence_without_batch_clock():
    entries = HwangGiapisLocalIonEntry2D(
        height_um=np.array([0.04, 0.05]),
        velocity_xz_sqrt_eV=np.array([[np.sqrt(177.0), 0.0],
                                      [np.sqrt(177.0), 0.0]]),
        expected_count=np.array([70.0, 70.0]))

    result = evolve_hwang_giapis_local_notch_event_driven_2d(
        entries, np.zeros(60), poly_potential_v=0.0,
        include_forward_scatter=False)

    assert result.removed_cell_count > 0
    assert result.removed_cell_count == (
        result.threshold_removed_cell_count + result.detached_cell_count)
    assert result.maximum_notch_depth_um >= 0.005
    assert np.all(result.notch_depth_by_height_um[20:] == 0.0)
    assert (
        result.landed_poly_count + result.landed_oxide_count
        + result.landed_photoresist_count + result.escaped_count
        == pytest.approx(result.launched_count))
    assert (
        result.threshold_removed_reactive_collisions
        + result.detached_reactive_collisions
        + np.sum(result.reactive_collision_inventory)
        == pytest.approx(
            result.direct_reactive_collisions
            + result.scattered_reactive_collisions))
    assert "event-driven" in result.provenance["front_integrator"]
    assert "five-cell least-squares" in result.provenance["reaction_angle"]
    assert result.provenance["exposed_oxide_charging"]

    control = evolve_hwang_giapis_local_notch_event_driven_2d(
        entries, np.zeros(60), poly_potential_v=0.0,
        include_forward_scatter=False,
        include_exposed_oxide_charging=False)
    assert not control.provenance["exposed_oxide_charging"]
    assert (
        control.landed_poly_count + control.landed_oxide_count
        + control.landed_photoresist_count + control.escaped_count
        == pytest.approx(control.launched_count))


def test_event_driven_local_notch_resumes_exactly_from_a_front_checkpoint():
    entries = HwangGiapisLocalIonEntry2D(
        height_um=np.array([0.04, 0.05]),
        velocity_xz_sqrt_eV=np.array([[np.sqrt(177.0), 0.0],
                                      [np.sqrt(177.0), 0.0]]),
        expected_count=np.array([140.0, 140.0]))
    checkpoints = []
    full = evolve_hwang_giapis_local_notch_event_driven_2d(
        entries, np.zeros(60), poly_potential_v=0.0,
        include_forward_scatter=False,
        checkpoint_callback=checkpoints.append)

    checkpoint = next(
        item for item in checkpoints
        if isinstance(item, HwangGiapisLocalNotchCheckpoint2D)
        and 0.0 < item.remaining_campaign_fraction < 1.0)
    resumed = evolve_hwang_giapis_local_notch_event_driven_2d(
        entries, np.zeros(60), poly_potential_v=0.0,
        include_forward_scatter=False, initial_checkpoint=checkpoint)

    assert np.array_equal(resumed.poly_cell, full.poly_cell)
    assert np.array_equal(
        resumed.reactive_collision_inventory,
        full.reactive_collision_inventory)
    assert np.array_equal(
        resumed.cumulative_oxide_ion_count,
        full.cumulative_oxide_ion_count)
    assert resumed.direct_reactive_collisions == full.direct_reactive_collisions
    assert resumed.landed_poly_count == full.landed_poly_count
    assert (
        resumed.threshold_removed_cell_count
        == full.threshold_removed_cell_count)
    assert resumed.detached_cell_count == full.detached_cell_count
    assert (
        resumed.threshold_removed_reactive_collisions
        == full.threshold_removed_reactive_collisions)
    assert (
        resumed.detached_reactive_collisions
        == full.detached_reactive_collisions)
