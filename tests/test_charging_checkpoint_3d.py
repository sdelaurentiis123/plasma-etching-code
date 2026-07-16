import numpy as np

from petch.charging_checkpoint_3d import (
    CHARGING_CHECKPOINT_SCHEMA,
    PhysicalChargingCheckpoint3D,
)
from petch.feature_step_3d import FeatureGeometry3D
from petch.fluorocarbon_lamagna import LaMagnaFluorocarbonState
from petch.material_mechanism_3d import MaterialSurfaceState3D
from petch.physical_sputtering import PhysicalSputterState


def test_charging_checkpoint_safe_npz_round_trip(tmp_path):
    phi = np.ones((2, 2, 2))
    phi[:, :, 1] = -1.0
    geometry = FeatureGeometry3D(phi, np.where(phi > 0.0, 1, 0), 0.5, 1e-6)
    checkpoint = PhysicalChargingCheckpoint3D(
        geometry=geometry,
        sigma_c_per_m2=np.array([-2e-6, 3e-6]),
        surface_state_type="petch.PhysicalSputterState",
        surface_state_fields={"removed_material_units_m2": np.array([2.0, 4.0])},
        surface_state_mesh_fingerprint="mesh-sha256",
        completed_duration_s=0.25,
        completed_steps=3,
        source_manifest_sha256="a" * 64)
    path = tmp_path / "checkpoint.npz"

    checkpoint.save(path)
    restored = PhysicalChargingCheckpoint3D.load(path)

    assert restored.schema == CHARGING_CHECKPOINT_SCHEMA
    assert np.array_equal(restored.geometry.phi, checkpoint.geometry.phi)
    assert np.array_equal(restored.geometry.material_id, checkpoint.geometry.material_id)
    assert np.array_equal(restored.sigma_c_per_m2, checkpoint.sigma_c_per_m2)
    state = restored.restore_surface_state()
    assert isinstance(state, PhysicalSputterState)
    assert np.array_equal(state.removed_material_units_m2, [2.0, 4.0])


def test_charging_checkpoint_restores_namespaced_multi_material_state(tmp_path):
    phi = np.ones((2, 2, 2)); phi[:, :, 1] = -1.0
    geometry = FeatureGeometry3D(phi, np.where(phi > 0.0, 1, 0), 0.5, 1e-6)
    checkpoint = PhysicalChargingCheckpoint3D(
        geometry=geometry, sigma_c_per_m2=np.array([0.0, 0.0]),
        surface_state_type="petch.MaterialSurfaceState3D",
        surface_state_fields={
            "m1__removed_material_units_m2": np.array([2.0, 0.0]),
            "m2__removed_material_units_m2": np.array([0.0, 4.0])},
        surface_state_upper_bounds={
            "m1__removed_material_units_m2": None,
            "m2__removed_material_units_m2": None},
        surface_state_remap_modes={
            "m1__removed_material_units_m2": "intensive",
            "m2__removed_material_units_m2": "conservative"},
        surface_state_mesh_fingerprint="mesh-sha256", completed_duration_s=0.25,
        completed_steps=3, source_manifest_sha256="b" * 64)
    path = tmp_path / "multi_material_checkpoint.npz"

    checkpoint.save(path)
    state = PhysicalChargingCheckpoint3D.load(path).restore_surface_state()

    assert isinstance(state, MaterialSurfaceState3D)
    assert np.array_equal(
        state.fields["m1__removed_material_units_m2"], [2.0, 0.0])
    assert np.array_equal(
        state.fields["m2__removed_material_units_m2"], [0.0, 4.0])
    assert state.surface_field_remap_modes() == {
        "m1__removed_material_units_m2": "intensive",
        "m2__removed_material_units_m2": "conservative"}


def test_charging_checkpoint_restores_lamagna_coverages_and_film_inventory(tmp_path):
    phi = np.ones((2, 2, 2)); phi[:, :, 1] = -1.0
    geometry = FeatureGeometry3D(phi, np.where(phi > 0.0, 1, 0), 0.5, 1e-6)
    state = LaMagnaFluorocarbonState(
        [0.1, 0.2], [0.3, 0.4], [0.5, 0.6], [2.0e18, 3.0e18], [4.0, 5.0])
    checkpoint = PhysicalChargingCheckpoint3D(
        geometry=geometry, sigma_c_per_m2=np.array([0.0, 0.0]),
        surface_state_type="petch.LaMagnaFluorocarbonState",
        surface_state_fields=state.conservative_surface_fields(),
        surface_state_mesh_fingerprint="mesh-sha256", completed_duration_s=1.0,
        completed_steps=2, source_manifest_sha256="c" * 64)
    path = tmp_path / "lamagna_checkpoint.npz"

    checkpoint.save(path)
    restored = PhysicalChargingCheckpoint3D.load(path).restore_surface_state()

    assert isinstance(restored, LaMagnaFluorocarbonState)
    for name, value in state.conservative_surface_fields().items():
        assert np.array_equal(restored.conservative_surface_fields()[name], value)
    assert restored.surface_field_remap_modes()["polymer_coverage"] == "intensive"
