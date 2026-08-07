import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from petch.feature_step_3d import (
    FeatureGeometry3D,
    make_rectangular_trench_geometry_3d,
)


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "krueger_2024_opening_certification",
    ROOT / "scripts" / "krueger_2024_opening_certification.py")
OPENING = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(OPENING)
PILOT_SPEC = importlib.util.spec_from_file_location(
    "krueger_2024_trench_pilot_opening_test",
    ROOT / "scripts" / "krueger_2024_trench_pilot.py")
PILOT = importlib.util.module_from_spec(PILOT_SPEC)
PILOT_SPEC.loader.exec_module(PILOT)


def _checkpoint(path, *, closed=False):
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.4, cell_length=0.1, domain_height=0.8,
        dx=0.05, opening_width=0.2, mask_thickness=0.3,
        substrate_top=0.3, etched_depth=0.2)
    mask = np.asarray(geometry.material_levelsets[2]).copy()
    if closed:
        x, y, z = geometry.coordinate_arrays
        x_grid, _, z_grid = np.meshgrid(x, y, z, indexing="ij")
        cap = np.minimum.reduce((
            z_grid - 0.425,
            0.525 - z_grid,
            0.1 - np.abs(x_grid - 0.2),
        ))
        mask = np.maximum(mask, cap)
    phi = np.maximum(geometry.material_levelsets[1], mask)
    metadata = {
        "dx": geometry.dx,
        "mesh_length_unit_m": geometry.mesh_length_unit_m,
        "physical_time_s": 7.0,
    }
    np.savez_compressed(
        path,
        phi=phi,
        material_id=geometry.material_id,
        material_levelset_1=geometry.material_levelsets[1],
        material_levelset_2=mask,
        metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
    )


def _audit(path, *, closed=False):
    event = ({
        "accepted": True,
        "kind": "gas_cavity_enclosed" if closed else "gas_cavity_opened",
    } if closed else None)
    history = [{
        "step": 1,
        "physical_time_s": 7.0,
        "metrics": {"mask_opening_nm": 200.0},
        "topology_event": event,
    }]
    path.write_text(json.dumps({
        "history": history,
        "final_metrics": {"mask_opening_nm": 200.0},
    }), encoding="utf-8")


def test_open_checkpoint_certifies_mask_material_width(tmp_path):
    checkpoint = tmp_path / "checkpoint.npz"
    audit = tmp_path / "audit.json"
    _checkpoint(checkpoint)
    _audit(audit)

    report = OPENING.build_report(
        audit, checkpoint, substrate_top_um=0.3,
        opening_center_um=0.2, opening_width_um=0.2,
        dense_samples_per_cell=8)

    assert report["definition_certified"]
    assert report["connectivity"]["exterior_to_mask_exit_open"]
    assert np.isclose(report["paper_connectivity_qualified_opening_nm"], 200.0)
    assert not report["grid_authority_certified"]


def test_sealed_pocket_width_is_not_reported_as_a_mask_opening(tmp_path):
    checkpoint = tmp_path / "checkpoint.npz"
    audit = tmp_path / "audit.json"
    _checkpoint(checkpoint, closed=True)
    _audit(audit, closed=True)

    report = OPENING.build_report(
        audit, checkpoint, substrate_top_um=0.3,
        opening_center_um=0.2, opening_width_um=0.2,
        dense_samples_per_cell=8)

    assert not report["definition_certified"]
    assert not report["connectivity"]["exterior_to_mask_exit_open"]
    assert report["mask_material_levelset_measurement"][
        "dense_bilinear_scan"]["width_nm"] > 0.0
    assert report["paper_connectivity_qualified_opening_nm"] == 0.0
    assert report["timeline"]["legacy_nonzero_while_closed_count"] == 1


def test_timeline_prices_one_step_flicker_without_filtering_it():
    widths = (20.0, 19.0, 18.0, 17.0, 20.0, 10.0, 20.0)
    history = []
    for step, width in enumerate(widths, start=1):
        event = None
        if step == 2:
            event = {"accepted": True, "kind": "gas_cavity_enclosed"}
        elif step == 4:
            event = {"accepted": True, "kind": "gas_cavity_opened"}
        history.append({
            "step": step,
            "physical_time_s": float(step),
            "metrics": {"mask_opening_nm": width},
            "topology_event": event,
        })

    result = OPENING._timeline({"history": history}, dx_nm=5.0)

    assert result["closed_history_count"] == 2
    assert result["legacy_nonzero_while_closed_count"] == 2
    assert result["maximum_open_state_one_step_flicker"]["step"] == 6
    assert result["maximum_open_state_one_step_flicker"][
        "absolute_residual_nm"] == 10.0


def test_missing_substrate_reference_is_refused(tmp_path):
    checkpoint = tmp_path / "checkpoint.npz"
    audit = tmp_path / "audit.json"
    _checkpoint(checkpoint)
    _audit(audit)

    try:
        OPENING.build_report(audit, checkpoint, dense_samples_per_cell=8)
    except ValueError as error:
        assert "substrate top" in str(error)
    else:
        raise AssertionError("missing benchmark geometry was silently defaulted")

    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.4, cell_length=0.1, domain_height=0.8,
        dx=0.05, opening_width=0.2, mask_thickness=0.3,
        substrate_top=0.3, etched_depth=0.2)
    with pytest.raises(TypeError, match="substrate_top_um"):
        PILOT.measure_krueger_metrics(geometry)


def test_explicit_substrate_reference_mismatch_is_refused(tmp_path):
    checkpoint = tmp_path / "checkpoint.npz"
    audit = tmp_path / "audit.json"
    _checkpoint(checkpoint)
    _audit(audit)
    payload = json.loads(audit.read_text(encoding="utf-8"))
    payload["configuration"] = {"geometry": {
        "substrate_top_um": 0.3,
        "opening_width_um": 0.2,
        "cell_width_um": 0.4,
    }}
    audit.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="disagrees with audit geometry"):
        OPENING.build_report(
            audit, checkpoint, substrate_top_um=0.1,
            dense_samples_per_cell=8)


@pytest.mark.parametrize("dx", (0.005, 0.01))
def test_shallow_subcell_etch_is_still_an_open_mask_throat(dx):
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.13, cell_length=0.02, domain_height=2.8,
        dx=dx, opening_width=0.09, mask_thickness=0.85,
        substrate_top=1.8, etched_depth=0.00683)

    connectivity = OPENING._opening_connectivity(
        geometry.phi, dx=dx, substrate_top=1.8,
        opening_center=0.065, opening_width=0.09)

    assert connectivity["exterior_to_mask_exit_open"]

    metrics = PILOT.measure_krueger_metrics(
        geometry, substrate_top_um=1.8)
    assert metrics["mask_opening_connected_to_exterior"]
    assert metrics["mask_opening_nm"] == metrics["mask_pocket_width_nm"]
    assert metrics["mask_opening_nm"] > 0.0


def test_pilot_reports_zero_but_retains_pocket_width_for_sealed_cap():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.4, cell_length=0.1, domain_height=0.8,
        dx=0.05, opening_width=0.2, mask_thickness=0.3,
        substrate_top=0.3, etched_depth=0.2)
    x, y, z = geometry.coordinate_arrays
    x_grid, _, z_grid = np.meshgrid(x, y, z, indexing="ij")
    cap = np.minimum.reduce((
        z_grid - 0.425,
        0.525 - z_grid,
        0.1 - np.abs(x_grid - 0.2),
    ))
    substrate = np.asarray(geometry.material_levelsets[1])
    mask = np.maximum(geometry.material_levelsets[2], cap)
    phi = np.maximum(substrate, mask)
    material = np.zeros(phi.shape, dtype=int)
    material[(phi >= 0.0) & (substrate >= mask)] = 1
    material[(phi >= 0.0) & (mask > substrate)] = 2
    sealed = FeatureGeometry3D(
        phi, material, geometry.dx, geometry.mesh_length_unit_m,
        material_levelsets={1: substrate, 2: mask})

    metrics = PILOT.measure_krueger_metrics(
        sealed, substrate_top_um=0.3, opening_center_um=0.2,
        opening_width_um=0.2)

    assert not metrics["mask_opening_connected_to_exterior"]
    assert metrics["mask_opening_nm"] == 0.0
    assert metrics["mask_pocket_width_nm"] > 0.0


def test_mask_metric_is_independent_of_substrate_union_shape():
    values = []
    for etched_depth in (0.1, 0.2):
        geometry = make_rectangular_trench_geometry_3d(
            cell_width=0.4, cell_length=0.1, domain_height=0.8,
            dx=0.05, opening_width=0.2, mask_thickness=0.3,
            substrate_top=0.3, etched_depth=etched_depth)
        values.append(PILOT.measure_krueger_metrics(
            geometry, substrate_top_um=0.3, opening_center_um=0.2,
            opening_width_um=0.2)["mask_opening_nm"])

    assert np.isclose(values[0], values[1])


def test_resolved_mirror_symmetry_is_distinct_from_subcell_offset():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.4, cell_length=0.1, domain_height=0.8,
        dx=0.05, opening_width=0.2, mask_thickness=0.3,
        substrate_top=0.3, etched_depth=0.2)
    phi = np.asarray(geometry.phi).copy()
    substrate = np.asarray(geometry.material_levelsets[1]).copy()
    mask = np.asarray(geometry.material_levelsets[2]).copy()
    # A small positive-side signed-distance perturbation changes no resolved
    # phase or material label, but can move an interpolated zero crossing.
    positive = (phi > 0.0) & (np.arange(phi.shape[0])[:, None, None]
                             < phi.shape[0] // 2)
    phi[positive] += 1.0e-5
    substrate[positive] += 1.0e-5
    mask[positive] += 1.0e-5
    perturbed = FeatureGeometry3D(
        phi, geometry.material_id, geometry.dx, geometry.mesh_length_unit_m,
        material_levelsets={1: substrate, 2: mask})

    symmetry = PILOT._mirror_symmetry_diagnostics(perturbed)

    assert symmetry["asymmetry_cell_count"] == 0
    assert symmetry["mirrored_node_sign_mismatch_pair_count"] == 0
    assert symmetry["mirrored_material_label_mismatch_pair_count"] == 0


def test_resolved_phase_asymmetry_counts_mirror_pair():
    geometry = make_rectangular_trench_geometry_3d(
        cell_width=0.4, cell_length=0.1, domain_height=0.8,
        dx=0.05, opening_width=0.2, mask_thickness=0.3,
        substrate_top=0.3, etched_depth=0.2)
    phi = np.asarray(geometry.phi).copy()
    material = np.asarray(geometry.material_id).copy()
    substrate = np.asarray(geometry.material_levelsets[1]).copy()
    mask = np.asarray(geometry.material_levelsets[2]).copy()
    cell = np.s_[1:3, 0:2, 0:2]
    phi[cell] = -geometry.dx
    substrate[cell] = -geometry.dx
    mask[cell] = -2.0 * geometry.dx
    material[cell] = 0
    asymmetric = FeatureGeometry3D(
        phi, material, geometry.dx, geometry.mesh_length_unit_m,
        material_levelsets={1: substrate, 2: mask})

    symmetry = PILOT._mirror_symmetry_diagnostics(asymmetric)

    assert symmetry["mirrored_node_sign_mismatch_pair_count"] > 0
    assert symmetry["asymmetry_cell_count"] > 0


def test_local_real_checkpoints_retain_open_throat_and_r18_width():
    root = ROOT / "results"
    shallow = root / "krueger_2024_multiresolution_audit" / (
        "initial_gpu_r18_20260717")
    paths = [
        shallow / "initial_5nm" / "checkpoint.npz",
        shallow / "initial_10nm" / "checkpoint.npz",
    ]
    r18 = root / "krueger_2024_base_calibration_r18" / (
        "mixed_operator_topology_continuation") / "checkpoint.npz"
    if not all(path.exists() for path in (*paths, r18)):
        pytest.skip("local checksum-bound Krueger checkpoints are unavailable")

    for path in paths:
        geometry, _, _, _ = PILOT._load_checkpoint(path)
        metrics = PILOT.measure_krueger_metrics(
            geometry, substrate_top_um=1.8)
        assert metrics["mask_opening_connected_to_exterior"]
        assert metrics["mask_opening_nm"] > 87.0

    geometry, _, _, _ = PILOT._load_checkpoint(r18)
    metrics = PILOT.measure_krueger_metrics(
        geometry, substrate_top_um=1.8)
    assert metrics["mask_opening_connected_to_exterior"]
    assert np.isclose(metrics["mask_opening_nm"], 39.469977232628885)


def test_pilot_and_certifier_share_periodic_seam_connectivity():
    phi = np.ones((5, 3, 5), dtype=float)
    # Exterior path approaches one side of the periodic x seam; the exit path
    # approaches the other. They meet only after periodic component merging.
    phi[3, 0:2, 2:5] = -1.0
    phi[0, 0:2, 1:3] = -1.0
    phi[-1, :, :] = phi[0, :, :]
    phi[:, -1, :] = phi[:, 0, :]
    coordinates = tuple(
        np.arange(size, dtype=float) for size in phi.shape)
    geometry = SimpleNamespace(
        phi=phi, dx=1.0, coordinate_arrays=coordinates)

    pilot = PILOT._mask_throat_connectivity(
        geometry, substrate_top_um=0.0,
        opening_center_um=0.0, opening_width_um=0.5)
    certifier = OPENING._opening_connectivity(
        phi, dx=1.0, substrate_top=0.0,
        opening_center=0.0, opening_width=0.5)

    assert pilot["exterior_to_mask_exit_open"]
    assert certifier["exterior_to_mask_exit_open"]
    assert pilot["shared_component_count"] == certifier["shared_component_count"]
