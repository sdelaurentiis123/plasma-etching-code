import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT
    / "results"
    / "curated"
    / "vella_hao_ale_depth"
    / "audit.json"
)


def test_vella_hao_board_is_no_fit_atom_balanced_and_passes_nominal_gate():
    audit = json.loads(AUDIT.read_text())
    board = audit["absolute_depth_board"]

    assert "no-depth-fit" in audit["claim"]
    assert "not a blind" in audit["claim"]
    assert not board["fit_to_depth_used"]
    assert board["nominal_gate_passed"]
    assert board["maximum_absolute_relative_error"] < 0.13
    assert np.allclose(
        board["predicted_depth_nm"],
        [0.48421003142705055, 0.945257376794254, 1.4006484721614578],
    )
    assert all(
        abs(item["atom_balance_residual_cm2"]) < 1.0
        for item in board["point_ledgers"]
    )


def test_vella_hao_board_keeps_boundary_and_uncertainty_limits_explicit():
    audit = json.loads(AUDIT.read_text())
    boundary = audit["boundary_and_units"]
    validity = audit["uncertainty_and_validity"]

    assert boundary["measured_positive_ion_flux_cm2_s"] == 3.7e16
    assert boundary["ion_bombardment_duration_s"] == 3.0
    assert not boundary["monolayer_alias_used_for_fluence_conversion"]
    assert not validity["combined_uncertainty_claimed"]
    assert "experimental IEDF and ion angular distribution" in (
        validity["incomplete_terms"])
    assert "inferred mean-energy" in validity["evidence_tier"]


def test_vella_hao_source_markers_pass_checksum_and_vision_gate():
    audit = json.loads(AUDIT.read_text())
    vision = audit["vision_digitization_audit"]

    assert vision["status"] == "passed"
    assert vision["paper_sha256"] == (
        "789bf50302fc2ed9175403c47f895e1fb8db5481be5fb980739cad97f33c2218")
    assert vision["all_markers_match_within_1p25_pixels"]
    assert vision["maximum_marker_center_error_pixels"] < 0.2
    assert vision["plot_frame_pixels"] == {
        "left": 657.0,
        "right": 1450.0,
        "top": 268.0,
        "bottom": 827.0,
    }


def test_legacy_vella_graves_rom_is_not_mislabeled_atom_conservative():
    audit = json.loads(AUDIT.read_text())
    rom = audit["legacy_transient_rom_atom_balance"]

    assert rom["status"] == "failed_elemental_conservation"
    assert rom["chlorine_created_per_ar_at_theta_top_0p5_theta_mixed_0p5"] > 0.0
    assert rom["printed_rom_depth_nm"] > (
        2.0 * rom["atom_conservative_theta_squared_depth_nm"])
    assert max(rom["time_step_refinement_difference_nm"].values()) < 1.0e-5
