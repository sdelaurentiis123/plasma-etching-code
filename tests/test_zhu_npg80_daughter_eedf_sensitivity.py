import json
from pathlib import Path

import pytest

from scripts.audit_zhu_npg80_daughter_eedf_sensitivity import _check


ROOT = Path(__file__).resolve().parents[1]
RECEIPT = (
    ROOT / "results" / "curated"
    / "zhu_npg80_daughter_eedf_sensitivity_v1" / "audit.json"
)


def test_committed_daughter_eedf_board_is_current_and_target_blind():
    _check(RECEIPT)
    payload = json.loads(RECEIPT.read_text(encoding="utf-8"))
    assert payload["sem_or_depth_target_used"] is False
    assert payload["source_inputs"]["raw_lxcat_bytes_committed"] is False
    assert payload["method"]["heavy_particle_composition_frozen"] is True
    assert payload["method"]["nonlinear_reactor_reclosed"] is False
    assert not payload["derived_inputs"]["complete_hf_eedf"]


def test_hf_materially_changes_all_frozen_power_nodes_and_forces_reclose():
    payload = json.loads(RECEIPT.read_text(encoding="utf-8"))
    for row in payload["power_board"]:
        variants = row["variants"]
        parent = variants["parent_only_replay"]
        hf = variants["partial_hf"]
        hf_f2 = variants["partial_hf_plus_f2"]
        assert parent["mean_electron_energy_eV"] == pytest.approx(
            row["stored_parent_only_mean_energy_eV"], rel=2.0e-13
        )
        maximum_ratio = 0.73 if row["absorbed_power_W"] == 60 else 0.75
        assert hf["mean_electron_energy_eV"] < maximum_ratio * parent[
            "mean_electron_energy_eV"
        ]
        assert parent["net_growth_rate_coefficient_m3_s"] > 0.0
        assert hf["net_growth_rate_coefficient_m3_s"] < 0.0
        assert hf_f2["net_growth_rate_coefficient_m3_s"] < 0.0
        assert abs(
            hf_f2["mean_electron_energy_eV"]
            / hf["mean_electron_energy_eV"] - 1.0
        ) < 0.007
    assert payload["finding"]["nonlinear_reclose_required"] is True
    assert not payload["certification"]["supports_feature_depth_prediction"]
