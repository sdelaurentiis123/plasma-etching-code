import importlib.util
from hashlib import sha256
from pathlib import Path
import sys

import numpy as np

from petch.feature_step_3d import FeatureGeometry3D
from petch.material_mechanism_3d import MaterialSurfaceState3D


ROOT = Path(__file__).parents[1]


def _pilot_module():
    path = ROOT / "scripts" / "krueger_2024_trench_pilot.py"
    spec = importlib.util.spec_from_file_location("krueger_2024_trench_pilot", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PILOT = _pilot_module()


def _resume_configuration(policy=None, remap_backend=None):
    result = {
        "minimum_step_duration_s": 0.001,
        "maximum_accepted_steps": 10000,
        "operator": "manufactured-fixed",
    }
    if policy is not None:
        result["topology_change_policy"] = policy
    if remap_backend is not None:
        result["surface_state_remap_backend"] = remap_backend
    return result


def test_resume_may_only_promote_refusal_to_declared_gas_cavity_continuation():
    previous = _resume_configuration()
    continued = _resume_configuration("continue_gas_cavity")

    compatible, changes = PILOT._monotone_resume_refinement(
        previous, continued)

    assert compatible
    assert changes == {
        "topology_change_policy": {
            "old": "refuse",
            "new": "continue_gas_cavity",
            "classification": "explicit_gas_cavity_continuation_scope",
        }}
    incompatible, _ = PILOT._monotone_resume_refinement(
        continued, _resume_configuration("refuse"))
    assert not incompatible


def test_resume_makes_historical_legacy_remap_explicit_but_never_changes_backend():
    compatible, changes = PILOT._monotone_resume_refinement(
        _resume_configuration(),
        _resume_configuration(remap_backend="legacy_knn"))

    assert compatible
    assert changes == {
        "surface_state_remap_backend_declaration": {
            "old": "implicit legacy_knn",
            "new": "explicit legacy_knn",
            "classification": "provenance_only_operator_declaration",
        }}
    incompatible, _ = PILOT._monotone_resume_refinement(
        _resume_configuration(remap_backend="indexed_knn"),
        _resume_configuration(remap_backend="common_refinement"))
    assert not incompatible


def test_pilot_cli_accepts_explicit_surface_state_remap_backend(monkeypatch):
    monkeypatch.setattr(sys, "argv", [
        "krueger_2024_trench_pilot.py",
        "--surface-state-remap-backend", "common_refinement",
    ])

    args = PILOT.parse_args()

    assert args.surface_state_remap_backend == "common_refinement"


def test_visibility_history_summary_preserves_shared_recovery_canary():
    diagnostics = {
        "CF2": {
            "visibility_float64_evaluated_count": 2,
            "visibility_recovered_hit_count": 1,
            "visibility_derived_horizon_extension_count": 1,
            "visibility_maximum_wrap_count": 1110,
            "visibility_final_maximum_wraps": 12556,
            "visibility_source_support_face_count": 29,
            "visibility_source_support_area_fraction": 8.4e-6,
            "visibility_maximum_source_support_distance": 0.0034,
            "visibility_source_relaunch_count": 3,
            "visibility_maximum_source_relaunch_distance": 0.0021,
            "visibility_overlap_skip_count": 5,
            "visibility_maximum_overlap_skip_depth": 0.0016,
        },
        # The sampled form-factor operator is shared by species.  Maxima preserve its one receipt;
        # summing would falsely report two independent recoveries.
        "CF": {
            "visibility_float64_evaluated_count": 2,
            "visibility_recovered_hit_count": 1,
            "visibility_derived_horizon_extension_count": 1,
            "visibility_maximum_wrap_count": 1110,
            "visibility_final_maximum_wraps": 12556,
            "visibility_source_support_face_count": 29,
            "visibility_source_support_area_fraction": 8.4e-6,
            "visibility_maximum_source_support_distance": 0.0034,
            "visibility_source_relaunch_count": 3,
            "visibility_maximum_source_relaunch_distance": 0.0021,
            "visibility_overlap_skip_count": 5,
            "visibility_maximum_overlap_skip_depth": 0.0016,
        },
    }

    assert PILOT._visibility_history_summary(diagnostics) == {
        "maximum_visibility_float64_evaluated_count": 2,
        "maximum_visibility_recovered_hit_count": 1,
        "maximum_visibility_derived_horizon_extension_count": 1,
        "maximum_visibility_wrap_count": 1110,
        "maximum_visibility_final_horizon_wraps": 12556,
        "maximum_visibility_source_support_face_count": 29,
        "maximum_visibility_source_support_area_fraction": 8.4e-6,
        "maximum_visibility_source_support_distance_um": 0.0034,
        "maximum_visibility_source_relaunch_count": 3,
        "maximum_visibility_source_relaunch_distance_um": 0.0021,
        "maximum_visibility_overlap_skip_count": 5,
        "maximum_visibility_overlap_skip_depth_um": 0.0016,
    }


def test_resume_retains_prior_refusal_and_accepted_event_time_provenance_once():
    terminal = {
        "geometry_event_kind": "gas_cavity_enclosed",
        "physical_time_lower_s": 56.9,
        "physical_time_upper_s": 56.91,
    }
    events = PILOT._resume_topology_events({"terminal_event": terminal})
    replayed = PILOT._resume_topology_events({
        "terminal_event": terminal, "topology_events": events})

    assert events == replayed
    assert events[0]["accepted"] is False
    assert events[0]["source"] == "prior_terminal_refusal"

    accepted = PILOT._accepted_topology_event(
        {"accepted": True, "kind": "gas_cavity_enclosed",
         "policy": "continue_gas_cavity"},
        physical_time_s=56.92, step_duration_s=0.01, step=1815)
    assert accepted["physical_time_lower_s"] == 56.910000000000004
    assert accepted["physical_time_upper_s"] == 56.92
    assert accepted["accepted_step"] == 1815
    assert accepted["source"] == "accepted_feature_step"


def test_policy_promotion_does_not_mutate_exact_checkpoint_state(tmp_path):
    phi = np.broadcast_to(
        np.array([1.0, 0.5, -0.5, -1.0]), (3, 3, 4)).copy()
    geometry = FeatureGeometry3D(
        phi, np.where(phi > 0.0, 1, 0), 0.1, 1e-6,
        material_levelsets={1: phi})
    state = MaterialSurfaceState3D(
        {"m1__coverage": np.array([0.1, 0.2])},
        {"m1__coverage": 1.0}, {"m1__coverage": "intensive"})
    checkpoint = tmp_path / "checkpoint.npz"
    PILOT._checkpoint(
        checkpoint, geometry, state, "fingerprint", 1814, 56.9, 0.01)
    before = sha256(checkpoint.read_bytes()).hexdigest()

    compatible, _ = PILOT._monotone_resume_refinement(
        _resume_configuration("refuse"),
        _resume_configuration("continue_gas_cavity"))
    restored_geometry, restored_state, fingerprint, metadata = (
        PILOT._load_checkpoint(checkpoint))

    assert compatible
    assert sha256(checkpoint.read_bytes()).hexdigest() == before
    assert np.array_equal(restored_geometry.phi, geometry.phi)
    assert np.array_equal(
        restored_state.fields["m1__coverage"],
        state.fields["m1__coverage"])
    assert fingerprint == "fingerprint"
    assert metadata["step"] == 1814
    assert metadata["physical_time_s"] == 56.9
