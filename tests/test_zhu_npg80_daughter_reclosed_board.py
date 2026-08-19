import json

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
