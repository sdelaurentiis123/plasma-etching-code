import json

import pytest

from scripts.run_zhu_npg80_daughter_reclosed_board import (
    SELF_BIAS_MAGNITUDE_V,
    _plasma_potential_V,
    _run_args,
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
