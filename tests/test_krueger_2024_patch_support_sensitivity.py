import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "krueger_2024_patch_support_sensitivity",
    ROOT / "scripts" / "krueger_2024_patch_support_sensitivity.py",
)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def test_geometry_only_table_normalizes_40nm_patch_over_20nm_periodic_y():
    vertices = np.asarray([
        [0.0, 0.0, 0.0], [40.0, 0.0, 0.0],
        [40.0, 20.0, 0.0], [0.0, 20.0, 0.0],
    ])
    faces = np.asarray([[0, 1, 2], [0, 2, 3]], dtype=int)
    direct = {
        "verts": vertices,
        "faces": faces,
        "face_area_m2": np.full(2, 400.0e-18),
        "gas_normals": np.tile([0.0, 0.0, 1.0], (2, 1)),
        "face_material_id": np.ones(2, dtype=int),
    }
    source = {"geometry": SimpleNamespace(
        phi=np.zeros((5, 3, 2)), dx=10.0, mesh_length_unit_m=1.0e-9,
        mesh_origin_m=np.zeros(3))}

    table = AUDIT._support_table(direct, source)

    forty = next(item for item in table if item["patch_scale_m"] == 40.0e-9)
    assert forty["patch_count"] == 1
    assert forty["represented_nominal_projected_area_m2"]["minimum"] == pytest.approx(
        40.0e-9 * 20.0e-9)
    assert forty["projected_support_fraction"]["minimum"] == pytest.approx(1.0)
    assert forty["surface_area_inventory_relative_error"] == pytest.approx(0.0)
    assert [item["minimum_mean_support_fraction"] for item in forty["thresholds"]] == list(
        AUDIT.closure.PATCH_SUPPORT_SENSITIVITY_THRESHOLDS)
    assert all(item["eligible_mean_patch_count"] == 1 for item in forty["thresholds"])
    assert all(item["excluded_surface_area_fraction"] == 0.0
               for item in forty["thresholds"])


def test_primary_threshold_is_predeclared_ten_percent_not_sensitivity_selected():
    assert AUDIT.closure.DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION == pytest.approx(0.10)
    assert AUDIT.closure.PATCH_SUPPORT_SENSITIVITY_THRESHOLDS == (
        0.05, 0.075, 0.10, 0.25, 0.50, 0.75)
