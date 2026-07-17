import importlib.util
import json
from pathlib import Path

import pytest


_SCRIPT = (
    Path(__file__).parents[1] / "scripts" /
    "krueger_2024_cuda_profile_report.py")
_SPEC = importlib.util.spec_from_file_location(
    "krueger_2024_cuda_profile_report", _SCRIPT)
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)


def _profile(*, unified):
    end = {
        "etch_depth_nm": 1.0,
        "floor_z_um": 1.799,
        "mask_opening_nm": 89.0,
        "mask_opening_per_y_widths_nm": [88.9, 89.0, 89.1],
        "mask_top_z_um": 2.65,
        "maximum_feature_width_nm": 88.0,
        "remaining_mask_thickness_nm": 850.0,
        "top_feature_width_nm": 88.0,
    }
    if unified:
        end["mask_opening_nm"] += 1e-8
    return {
        "schema": _MODULE.PROFILE_SCHEMA,
        "status": "complete",
        "held_out_profile_data_read": False,
        "configuration": {
            "dx_nm": 5.0,
            "step_duration_s": 0.025,
            "seed": 241,
            "total_profiled_physical_time_s": 0.05 if not unified else 0.025,
            **({"total_executed_physical_time_s": 0.05} if unified else {}),
            "calibration_parameters": {"f": 0.9, "s": 0.5},
            "operator": {"transport": "face_gather"},
        },
        "hardware": ({
            "unified_device_selection": True,
            "transport_device": "cuda:0",
            "level_set_device": "cuda:0",
        } if unified else {"device": "cuda:0"}),
        "end_metrics": end,
        "runtime": {
            "operator_initialization_wall_time_s": 10.0,
            "positive_warmup_step_wall_time_s": [12.0],
            "mean_profiled_step_wall_time_s": 8.0,
        },
        "step_receipts": [{
            "maximum_material_ledger_residual_units_m2": 0.0,
            "maximum_radiosity_relative_balance_error": 1e-12,
            "topology_event": None,
        }],
        "top_cumulative_functions": [
            {"path": "src/petch/feature_step_3d.py", "function": "advance_feature_step_3d",
             "cumulative_time_s": 7.8},
            {"path": "src/petch/material_mechanism_3d.py", "function": "advance_by_material",
             "cumulative_time_s": 2.3},
            {"path": "src/petch/boundary_transport_3d.py",
             "function": "gather_boundary_state_ballistic_3d", "cumulative_time_s": 2.2},
            {"path": "src/petch/feature_step_3d.py",
             "function": "_apply_diffuse_neutral_transport", "cumulative_time_s": 1.2},
            {"path": "src/petch/feature_step_3d.py",
             "function": "conservative_remap_surface_state", "cumulative_time_s": 1.0},
            {"path": "src/petch/feature_step_3d.py", "function": "_redistance_feature_field",
             "cumulative_time_s": 0.05},
        ],
    }


def _multiresolution():
    values = {
        "etch_depth_nm_s": -0.00073,
        "mask_opening_nm_s": -0.0143,
        "remaining_mask_thickness_nm_s": 0.0103,
        "maximum_feature_width_nm_s": 0.712,
        "top_feature_width_nm_s": 0.712,
    }
    return {
        "paired_10nm_vs_5nm": {"initial": {
            name: {"coarse_10nm": 1.0 + relative, "fine_5nm": 1.0,
                   "relative_to_fine": relative}
            for name, relative in values.items()
        }},
        "cases": [
            {"dx_nm": 5.0, "wall_time_s": 70.0},
            {"dx_nm": 10.0, "wall_time_s": 17.0},
        ],
    }


def _write(tmp_path, name, payload):
    path = tmp_path / name
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_summary_separates_parity_timing_and_refinement_claims(tmp_path):
    mixed = _write(tmp_path, "mixed.json", _profile(unified=False))
    unified = _write(tmp_path, "unified.json", _profile(unified=True))
    multires = _write(tmp_path, "multires.json", _multiresolution())

    result = _MODULE.summarize(mixed, unified, multires)

    assert result["paired_cuda_parity"]["within_1e_6_nm"] is True
    assert result["unified_cuda_runtime"]["steady_profiled_step_wall_time_s"] == 8.0
    assert result["unified_cuda_runtime"][
        "naive_fixed_0p025s_step_projection_for_60s"]["wall_time_h"] == pytest.approx(
            5.333333333333333)
    budget = result["unified_cuda_runtime"]["timing_budget"]
    assert sum(row["time_s"] for row in budget) == pytest.approx(7.8)
    assert result["paired_10nm_5nm_initial_refinement"]["metrics"][
        "maximum_feature_width_nm_s"]["absolute_relative_percent"] == pytest.approx(71.2)
    assert result["held_out_profile_data_read"] is False
    assert result["calibration_performed"] is False


def test_summary_refuses_mixed_device_timing_authority(tmp_path):
    mixed = _write(tmp_path, "mixed.json", _profile(unified=False))
    bad = _profile(unified=True)
    bad["hardware"]["level_set_device"] = "cpu"
    unified = _write(tmp_path, "unified.json", bad)
    multires = _write(tmp_path, "multires.json", _multiresolution())

    with pytest.raises(ValueError, match="one declared CUDA device"):
        _MODULE.summarize(mixed, unified, multires)


def test_plot_summary_writes_nonempty_png(tmp_path):
    mixed = _write(tmp_path, "mixed.json", _profile(unified=False))
    unified = _write(tmp_path, "unified.json", _profile(unified=True))
    multires = _write(tmp_path, "multires.json", _multiresolution())
    result = _MODULE.summarize(mixed, unified, multires)

    output = tmp_path / "summary.png"
    _MODULE.plot_summary(result, output)

    assert output.read_bytes().startswith(b"\x89PNG")
    assert output.stat().st_size > 10_000
