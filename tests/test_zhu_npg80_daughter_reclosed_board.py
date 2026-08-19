from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

import pytest

from scripts.run_zhu_npg80_daughter_reclosed_board import (
    SELF_BIAS_MAGNITUDE_V,
    _plasma_potential_V,
    _run_args,
)
from scripts.audit_zhu_npg80_daughter_wafer_dose_board import (
    CANDIDATE_FEATURE_TRANSMISSIONS,
    CANDIDATE_SURFACE_YIELDS,
    _depth_sensitivity,
)


ROOT = Path(__file__).resolve().parents[1]


def test_reclosed_board_args_preserve_self_bias_and_full_daughter_inputs(tmp_path):
    args = _run_args(
        source_workbook=tmp_path / "o2.xlsx",
        hcl_lxcat=tmp_path / "hcl.txt",
        f2_lxcat=tmp_path / "f2.txt",
        initial_state=tmp_path / "prior.json",
        absorbed_power_W=90,
        grounded_sheath_V=19.5,
        maximum_evaluations=500,
        nonlinear_verbose=0,
    )
    assert (
        args.powered_electrode_sheath_drop_V
        - args.grounded_surface_sheath_drop_V
        == SELF_BIAS_MAGNITUDE_V
    )
    assert args.hcl_lxcat.name == "hcl.txt"
    assert args.f2_lxcat.name == "f2.txt"
    assert args.minimum_reduced_field_Td == 40.0
    assert args.maximum_reduced_field_Td == 900.0


def test_plasma_potential_uses_multi_ion_inventory_not_flux_target():
    payload = {
        "state": {
            "mean_electron_energy_eV": 6.0,
            "densities_m3": {"F+": 2.0, "CF3+": 1.0},
            "axial_positive_ion_flux_m2_s": {"F+": 7.0, "CF3+": 9.0},
        }
    }
    original = _plasma_potential_V(payload)
    payload = json.loads(json.dumps(payload))
    payload["state"]["axial_positive_ion_flux_m2_s"] = {
        "F+": 700.0,
        "CF3+": 0.09,
    }
    assert _plasma_potential_V(payload) == pytest.approx(original)


def test_wafer_dose_sensitivity_keeps_surface_and_transport_unfitted():
    board = _depth_sensitivity(
        ion_flux_m2_s=1.0e19,
        duration_s=1200.0,
        film_thickness_nm=700.0,
    )
    assert board["surface_yield_candidates_are_fitted"] is False
    assert board["feature_transmission_candidates_are_fitted"] is False
    assert len(board["rows"]) == (
        len(CANDIDATE_SURFACE_YIELDS)
        * len(CANDIDATE_FEATURE_TRANSMISSIONS)
    )
    lowest = board["rows"][0]["film_capped_depth_nm_by_density_endpoint"]
    highest = board["rows"][-1]["film_capped_depth_nm_by_density_endpoint"]
    assert lowest[0] > lowest[1]
    assert highest[0] == highest[1] == 700.0


def test_wafer_dose_script_is_directly_executable():
    completed = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" /
                "audit_zhu_npg80_daughter_wafer_dose_board.py"),
            "--help",
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    assert "--reactor-board" in completed.stdout


def test_committed_reclosed_and_wafer_dose_boards_fail_closed_on_authority():
    reactor_path = (
        ROOT / "results" / "curated" / "zhu_npg80_daughter_reclosed_v1"
        / "audit.json"
    )
    wafer_path = (
        ROOT / "results" / "curated"
        / "zhu_npg80_daughter_wafer_dose_v1" / "audit.json"
    )
    reactor = json.loads(reactor_path.read_text(encoding="utf-8"))
    wafer = json.loads(wafer_path.read_text(encoding="utf-8"))
    assert reactor["target_outcome_used"] is False
    assert reactor["certification"][
        "daughter_collision_nonlinear_reclose_completed"
    ] is True
    assert reactor["certification"]["supports_unique_sem_profile"] is False
    assert [row["absorbed_power_W"] for row in reactor["state_board"]] == [
        60, 90, 105, 120,
    ]
    for row in reactor["state_board"]:
        state_path = ROOT / row["state_path"]
        assert sha256(state_path.read_bytes()).hexdigest() == row[
            "state_sha256"
        ]
        assert row["maximum_normalized_reactor_residual"] < 2.0e-6
        assert abs(row["sheath_fixed_point_residual_V"]) < 0.01

    certification = wafer["certification"]
    assert certification[
        "all_axisymmetric_lifts_conserved_and_grid_converged"
    ] is True
    assert certification["supports_conditional_atom_counted_depths"] is True
    assert certification["supports_unique_sem_profile"] is False
    assert certification["species_resolved_tio2_surface_law_validated"] is False
    rows = wafer["power_board"]
    assert min(
        row["electron_collision_basis_neutral_fraction"] for row in rows
    ) > 0.75
    assert max(row["global_positive_ion_flux_m2_s"] for row in rows) / min(
        row["global_positive_ion_flux_m2_s"] for row in rows
    ) < 1.04
    assert rows[-1]["global_neutral_F_thermal_flux_m2_s"] / rows[0][
        "global_neutral_F_thermal_flux_m2_s"
    ] > 8.0
