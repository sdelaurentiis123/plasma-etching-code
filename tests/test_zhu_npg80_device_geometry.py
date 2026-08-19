import json
from pathlib import Path

import numpy as np

from petch.feature_step_3d import make_square_pillar_mask_geometry_3d


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = (
    ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
    / "device_geometry_evidence.json"
)


def test_same_group_geometry_prior_is_not_promoted_to_target_layout():
    evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
    prior = evidence["same_group_public_geometry_prior"]
    guard = evidence["process_nontransfer_guard"]

    assert prior["pitch_nm"] == 400.0
    assert prior["pillar_height_nm"] == 700.0
    assert prior["minimum_reliable_geometry_nm"] == 80.0
    assert prior["maximum_lateral_dimension_nm"] == 320.0
    assert guard["same_process"] is False
    assert guard["pitch_is_target_measured"] is False
    assert guard["width_scan_is_target_layout"] is False
    assert evidence["sem_target_used"] is False
    assert evidence["measured_profile_target_used"] is False


def test_square_geometry_has_film_mask_base_and_open_gas():
    geometry = make_square_pillar_mask_geometry_3d(
        pitch=0.4,
        domain_height=1.05,
        dx=0.025,
        pillar_width=0.2,
        film_thickness=0.7,
        mask_thickness=0.05,
        base_top=0.1,
        film_material_id=11,
        mask_material_id=12,
        base_material_id=13,
    )
    x = np.arange(geometry.phi.shape[0]) * geometry.dx
    z = np.arange(geometry.phi.shape[2]) * geometry.dx
    center = int(np.argmin(np.abs(x - 0.2)))
    edge = int(np.argmin(np.abs(x - 0.0)))
    base = int(np.argmin(np.abs(z - 0.05)))
    film = int(np.argmin(np.abs(z - 0.4)))
    mask = int(np.argmin(np.abs(z - 0.825)))
    gas = int(np.argmin(np.abs(z - 0.95)))

    assert geometry.material_id[center, center, base] == 13
    assert geometry.material_id[center, center, film] == 11
    assert geometry.material_id[center, center, mask] == 12
    assert geometry.material_id[edge, edge, mask] == 0
    assert geometry.material_id[center, center, gas] == 0
    assert set(geometry.material_levelsets) == {11, 12, 13}


def test_square_geometry_etched_snapshot_lowers_only_exposed_film():
    geometry = make_square_pillar_mask_geometry_3d(
        pitch=0.4,
        domain_height=1.05,
        dx=0.025,
        pillar_width=0.2,
        film_thickness=0.7,
        mask_thickness=0.05,
        base_top=0.1,
        etched_depth=0.4,
        film_material_id=11,
        mask_material_id=12,
        base_material_id=13,
    )
    x = np.arange(geometry.phi.shape[0]) * geometry.dx
    z = np.arange(geometry.phi.shape[2]) * geometry.dx
    center = int(np.argmin(np.abs(x - 0.2)))
    edge = int(np.argmin(np.abs(x - 0.0)))
    below_floor = int(np.argmin(np.abs(z - 0.35)))
    above_floor = int(np.argmin(np.abs(z - 0.55)))
    protected = int(np.argmin(np.abs(z - 0.7)))

    assert geometry.material_id[edge, edge, below_floor] == 11
    assert geometry.material_id[edge, edge, above_floor] == 0
    assert geometry.material_id[center, center, protected] == 11


def test_square_geometry_refuses_depth_beyond_film():
    with np.testing.assert_raises_regex(ValueError, "invalid square-pillar"):
        make_square_pillar_mask_geometry_3d(
            pitch=0.4,
            domain_height=1.05,
            dx=0.025,
            pillar_width=0.2,
            film_thickness=0.7,
            mask_thickness=0.05,
            base_top=0.1,
            etched_depth=0.701,
        )


def test_geometry_evidence_square_cell_arithmetic_is_exact():
    evidence = json.loads(EVIDENCE.read_text(encoding="utf-8"))
    pitch = evidence["same_group_public_geometry_prior"]["pitch_nm"]
    for row in evidence["derived_square_unit_cells"]:
        width = row["width_nm"]
        mask_fraction = (width / pitch) ** 2
        assert row["minimum_neighbor_gap_nm"] == pitch - width
        assert np.isclose(row["mask_area_fraction"], mask_fraction)
        assert np.isclose(row["exposed_area_fraction"], 1.0 - mask_fraction)
