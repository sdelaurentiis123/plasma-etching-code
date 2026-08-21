from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = ROOT / "tests" / "data" / "zhu_v5_prefailure_8ccd278.npz"
CHECKPOINT_SHA256 = (
    "95d43ec8df8372ef14df57cd28df93ac90fe5c2e72da84cceaf66965f7c50998"
)
V6_ACCEPTANCE = (
    ROOT / "results" / "curated" / "zhu_npg80_moving_cr_profiles_v1"
    / "trajectories"
    / "w200_s18.017_ion_high_tail_0p0_c904f1b4f4500a0e.json"
)
V6_ACCEPTANCE_SHA256 = (
    "6ee2f2dbb6cf1ef332fb9da182297a04a9424e4a8c80aa6f81f31deb30d6e23a"
)


def test_zhu_v5_prefailure_checkpoint_is_pinned_and_self_describing():
    assert sha256(CHECKPOINT.read_bytes()).hexdigest() == CHECKPOINT_SHA256

    with np.load(CHECKPOINT, allow_pickle=False) as stored:
        metadata = json.loads(str(stored["metadata_json"]))
        assert stored["phi"].shape == (41, 41, 106)
        assert stored["material_id"].shape == stored["phi"].shape
        assert metadata["material_levelset_ids"] == [1, 2, 3]
        assert metadata["elapsed_s"] == 719.8619631901854
        assert metadata["accepted_steps"] == 149
        assert metadata["step_duration_s"] == 4.831288343558282
        assert metadata["surface_state_names"] == [
            "m1__removed_material_units_m2",
            "m2__removed_material_units_m2",
        ]


def test_zhu_v5_prefailure_checkpoint_advances_one_certified_step():
    command = [
        sys.executable,
        "-u",
        "scripts/reproduce_zhu_npg80_moving_cr_cell.py",
        "--width-nm", "200",
        "--scenario", "ion_high_tail_0p0",
        "--selectivity", "18.016664028610727",
        "--duration-s", "1200",
        "--dx-nm", "10",
        "--transport-device", "cpu",
        "--resume-checkpoint", str(CHECKPOINT),
    ]
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(ROOT / "src")
    completed = subprocess.run(
        command, cwd=ROOT, env=environment, text=True,
        capture_output=True, timeout=60, check=False)
    assert completed.returncode == 0, completed.stderr
    json_start = completed.stdout.rfind("\n{")
    assert json_start >= 0, completed.stdout
    payload = json.loads(completed.stdout[json_start + 1:])
    assert payload["status"] == "checkpoint step passed"
    assert payload["reassigned_unresolved_material_nodes"] == 393
    assert payload["elapsed_s_before_step"] == 719.8619631901854


def test_zhu_v6_clean_acceptance_trajectory_is_pinned_and_conservative():
    assert sha256(V6_ACCEPTANCE.read_bytes()).hexdigest() == V6_ACCEPTANCE_SHA256
    payload = json.loads(V6_ACCEPTANCE.read_text())
    spec = payload["job_spec"]
    assert spec["model_revision"] == (
        "two-material-moving-tio2-cr-owner-projection-v6")
    assert spec["width_nm"] == 200.0
    assert spec["scenario"]["name"] == "ion_high_tail_0p0"
    assert spec["selectivity"] == 18.016664028610727
    assert len(payload["profiles"]) == 2
    assert {item["terminal_reason"] for item in payload["profiles"]} == {
        "requested_duration", "domain_gas_breakthrough"}
    for item in payload["profiles"]:
        assert item["maximum_transport_relative_particle_balance_error"] == 0.0
        assert item[
            "maximum_state_remap_relative_conservation_residual"] < 1.0e-14
        assert item["parameter_evidence_supports_prediction"] is False
        assert item["profile"]["etched_depth_nm"] == 679.4077751001814
