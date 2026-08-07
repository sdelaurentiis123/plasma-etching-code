import json
from pathlib import Path

import pytest


ROOT = Path(__file__).parents[1]
AUDIT = (
    ROOT
    / "results"
    / "curated"
    / "krueger_qmc_receiver_singularity"
    / "audit.json"
)


def _audit():
    return json.loads(AUDIT.read_text(encoding="utf-8"))


def test_paired_control_is_same_face_and_zero_neutral_flux():
    report = _audit()
    qmc = report["with_radiosity"]["target_face"]
    control = report["without_radiosity"]["same_target_face"]
    assert qmc["global_face_index"] == control["global_face_index"] == 652
    assert qmc["centroid_mesh_units"] == control["centroid_mesh_units"]
    assert set(control["neutral_flux_m2_s"].values()) == {0.0}
    assert control["etch_velocity_m_s"] == 0.0
    assert control["growth_velocity_m_s"] == 0.0


def test_receiver_area_explains_qmc_speed_spike_without_chemistry_change():
    report = _audit()
    face = report["with_radiosity"]["target_face"]
    assert face["area_m2"] == pytest.approx(1.0586623602343975e-21)
    assert face["area_over_active_median"] < 1.0e-4
    assert face["energetic_flux_m2_s"]["ions"] == 0.0
    assert face["neutral_flux_m2_s"]["C3F4"] > 6.0e24
    assert report["speed_ratio_qmc_over_control"] > 100.0
    assert report["all_diagnostic_gates_passed"]


def test_diagnosis_does_not_authorize_smoothing_or_physics_fit():
    decision = _audit()["decision"]
    assert not decision["qmc_moving_profile_authority"]
    assert not decision["chemistry_change_authorized"]
    assert not decision["surface_smoothing_authorized"]
    assert decision["krueger_authority_candidate"] == (
        "deterministic_extruded_2d"
    )
