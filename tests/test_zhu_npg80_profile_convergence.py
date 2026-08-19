from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT / "results" / "curated" / "zhu_npg80_profile_convergence_v1"
    / "audit.json"
)


def _load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def test_profile_convergence_is_target_free_and_checksum_bound():
    audit = _load(AUDIT)

    assert audit["target_sem_used"] is False
    assert audit["target_depth_used"] is False
    assert audit["supports_absolute_target_profile_prediction"] is False
    assert audit["supports_atomic_accuracy"] is False
    preregistration = audit["inputs"]["convergence_preregistration"]
    preregistration_path = ROOT / preregistration["path"]
    assert sha256(preregistration_path.read_bytes()).hexdigest() == (
        preregistration["sha256"]
    )
    for receipt in audit["inputs"]["cases"].values():
        path = ROOT / receipt["path"]
        assert sha256(path.read_bytes()).hexdigest() == receipt["sha256"]


def test_profile_convergence_preserves_the_frozen_failure():
    audit = _load(AUDIT)
    gates = audit["gate_results"]

    assert gates == {
        "fine_grid_xy_symmetry": True,
        "grid_cd": False,
        "grid_depth": True,
        "particle_balance": True,
        "state_remap_conservation": True,
        "timestep_cd": True,
        "timestep_depth": True,
    }
    assert audit["all_numerical_gates_pass"] is False
    assert audit["timestep_comparison"]["depth_relative_change"] < 0.006
    assert audit["timestep_comparison"][
        "maximum_cd_absolute_change_nm"] < 4.0
    assert audit["grid_comparison"]["depth_relative_change"] < 0.009
    assert 24.0 < audit["grid_comparison"][
        "maximum_cd_absolute_change_nm"] < 24.3
    assert audit["fine_grid_maximum_xy_asymmetry_nm"] < 0.001


def test_profile_convergence_case_specs_separate_time_and_space_error():
    audit = _load(AUDIT)
    case_receipts = audit["inputs"]["cases"]
    cases = {
        name: _load(ROOT / receipt["path"])
        for name, receipt in case_receipts.items()
    }

    assert cases["coarse_dt8"]["case_specification"] == {
        "mesh_spacing_nm": 20.0,
        "maximum_step_s": 8.0,
    }
    assert cases["coarse_dt4"]["case_specification"] == {
        "mesh_spacing_nm": 20.0,
        "maximum_step_s": 4.0,
    }
    assert cases["fine_dt4"]["case_specification"] == {
        "mesh_spacing_nm": 10.0,
        "maximum_step_s": 4.0,
    }
    for case in cases.values():
        assert case["target_sem_used"] is False
        assert case["target_depth_used"] is False
        assert case["profile"]["tio2_clearance_detected"] is False
        assert case["profile"]["validity"][
            "parameter_evidence_supports_prediction"] is False
