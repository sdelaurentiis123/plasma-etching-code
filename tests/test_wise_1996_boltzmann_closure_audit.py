import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/audit_wise_1996_boltzmann_closure.py"


def _module():
    spec = importlib.util.spec_from_file_location(
        "audit_wise_1996_boltzmann_closure", SCRIPT
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_wise_boltzmann_closure_replays_direct_marker_board():
    module = _module()
    result = module.audit(
        module.DEFAULT_DATA, module.DEFAULT_MANIFEST, module.DEFAULT_GEOMETRY
    )
    local = result["pressure_gradient_with_measured_local_temperature"]
    assert local["axis_to_outer_drop_V"] == pytest.approx(
        -7.67978375687332, rel=2.0e-12
    )
    assert local["unweighted_mape_percent_excluding_gauge"] == pytest.approx(
        5.066249972222311, rel=2.0e-12
    )
    assert result["density_width"][
        "digitized_full_width_half_maximum_m"
    ] == pytest.approx(0.07533409080957697, rel=2.0e-12)
    assert result["density_width"]["inside_independent_interval"] is True
    assert result["certification"]["feature_depth_used"] is False
    assert result["certification"]["formal_uncertainty_weighted_pass"] is False
    assert result["certification"][
        "supports_spatial_reactor_state_prediction"
    ] is False


def test_wise_boltzmann_audit_rejects_corrupted_csv(tmp_path):
    module = _module()
    corrupt = tmp_path / "wise.csv"
    corrupt.write_bytes(module.DEFAULT_DATA.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="checksum"):
        module.audit(corrupt, module.DEFAULT_MANIFEST, module.DEFAULT_GEOMETRY)
