from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = (
    ROOT / "results" / "quarantine" / "zhu_npg80_v6_failure_20260821"
)
CHECKPOINT = EVIDENCE / "zhu_v6_w120_failure_checkpoint.npz"
CHECKPOINT_SHA256 = (
    "1fbd52695a1df0379c67d22d9743d8e74bee29f2322142fd41737cafd7881dd6"
)
DIAGNOSTIC_LOG = EVIDENCE / "zhu_v6_w120_diagnostic.log"
DIAGNOSTIC_LOG_SHA256 = (
    "6683045c63d1efd0978a46d8940325043d758bb8c5a5aaac446c5bf092fc9a16"
)


def test_zhu_v6_sparse_remap_failure_checkpoint_is_pinned():
    assert sha256(CHECKPOINT.read_bytes()).hexdigest() == CHECKPOINT_SHA256
    assert sha256(DIAGNOSTIC_LOG.read_bytes()).hexdigest() == DIAGNOSTIC_LOG_SHA256
    diagnostic = DIAGNOSTIC_LOG.read_text(encoding="utf-8")
    assert "negative_weights count=1, first_entry=24516" in diagnostic
    assert "minimum=-2.2204460492503131e-16" in diagnostic
    with np.load(CHECKPOINT, allow_pickle=False) as stored:
        metadata = json.loads(str(stored["metadata_json"]))
        assert stored["phi"].shape == (41, 41, 106)
        assert stored["material_id"].shape == stored["phi"].shape
        assert metadata["material_levelset_ids"] == [1, 2, 3]
        assert metadata["surface_state_names"] == [
            "m1__removed_material_units_m2",
            "m2__removed_material_units_m2",
        ]
        assert metadata["accepted_steps"] == 133
        assert metadata["elapsed_s"] == 642.5613496932526
        assert metadata["step_duration_s"] == 4.831288343558282


def test_zhu_v7_sparse_remap_checkpoint_advances_without_negative_weight():
    command = [
        sys.executable,
        "-u",
        "scripts/reproduce_zhu_npg80_moving_cr_cell.py",
        "--width-nm", "120",
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
    assert payload["accepted_steps_before_step"] == 133
    assert payload["elapsed_s_before_step"] == 642.5613496932526
    assert payload["step_duration_s"] == 4.831288343558282
    assert len(payload["surface_state_transfer_fingerprint"]) == 64
    assert payload[
        "maximum_state_remap_relative_conservation_residual"] < 1.0e-14
