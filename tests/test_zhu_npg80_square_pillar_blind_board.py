from hashlib import sha256
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT / "results" / "curated" / "zhu_npg80_square_pillar_blind_v1"
    / "audit.json"
)


def _load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def test_blind_board_is_target_free_and_input_checksum_bound():
    audit = _load(AUDIT)
    assert audit["target_sem_used"] is False
    assert audit["target_depth_used"] is False
    assert audit["coefficient_selected_from_target"] is None
    assert audit["certification"]["target_outcome_used"] is False
    for receipt in audit["inputs"].values():
        path = ROOT / receipt["path"]
        assert sha256(path.read_bytes()).hexdigest() == receipt["sha256"]


def test_transport_board_conserves_and_is_grid_converged():
    audit = _load(AUDIT)
    certification = audit["certification"]
    assert certification["periodic_3d_transport_is_deterministic"]
    assert certification["incident_measure_conserved_below_2e_12"]
    assert certification["transport_convergence_passed_5_percent"]
    assert certification["maximum_incident_measure_relative_residual"] < 2e-12
    assert audit["transport_convergence"][
        "maximum_floor_transmission_relative_change"] < 0.01
    assert len(audit["transport_snapshot_board"]) == 7 * 6


def test_tightest_square_retains_bounded_floor_dose_but_not_unique_depth():
    audit = _load(AUDIT)
    tight = next(
        row for row in audit["width_summary"] if row["width_nm"] == 320.0)
    analog = next(
        row for row in audit["published_cross_machine_analog_slice"]
        ["width_summary"] if row["width_nm"] == 320.0)
    depth_660 = next(
        row for row in audit["transport_snapshot_board"]
        if row["width_nm"] == 320.0 and row["etched_depth_nm"] == 660.0)
    floor = [
        row["floor_local_transmission"]
        for row in depth_660["scenario"].values()
    ]
    assert min(floor) > 0.75
    assert max(floor) <= 1.0
    assert tight["mask_pinned_depth_envelope_nm"][0] < 250.0
    assert tight["mask_pinned_depth_envelope_nm"][1] == 700.0
    assert 600.0 < analog["mask_pinned_depth_envelope_nm"][0] < 650.0
    assert analog["mask_pinned_depth_envelope_nm"][1] == 700.0
    assert audit["certification"]["supports_unique_target_depth"] is False
    assert audit["certification"]["supports_unique_target_sem"] is False
