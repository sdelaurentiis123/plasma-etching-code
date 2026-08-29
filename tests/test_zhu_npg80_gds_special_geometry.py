import pytest

from scripts.audit_zhu_npg80_gds_special_geometry import build


def test_all_non_square_gds_cells_build_as_periodic_layered_geometry():
    audit = build()
    rows = {row["cell_name"]: row for row in audit["rows"]}

    assert set(rows) == {
        "CROSS_250x105", "INVHOLE_105", "INVHOLE_250", "RECT_250x105",
    }
    assert rows["RECT_250x105"]["analytic_mask_area_fraction"] == pytest.approx(
        250.0 * 105.0 / 350.0 ** 2)
    assert rows["CROSS_250x105"]["analytic_mask_area_fraction"] == pytest.approx(
        (2.0 * 250.0 * 105.0 - 105.0 ** 2) / 350.0 ** 2)
    assert rows["INVHOLE_105"]["analytic_mask_area_fraction"] == pytest.approx(
        1.0 - 105.0 ** 2 / 350.0 ** 2)
    assert rows["INVHOLE_250"]["analytic_mask_area_fraction"] == pytest.approx(
        1.0 - 250.0 ** 2 / 350.0 ** 2)
    assert rows["RECT_250x105"]["center_is_mask"] is True
    assert rows["CROSS_250x105"]["center_is_mask"] is True
    assert rows["INVHOLE_105"]["center_is_mask"] is False
    assert rows["INVHOLE_250"]["center_is_mask"] is False
    assert all(row["footprint_periodic_x"] for row in rows.values())
    assert all(row["footprint_periodic_y"] for row in rows.values())
    assert all(row["finite_levelsets"] for row in rows.values())
    assert all(row["material_ids"] == [1, 2, 3] for row in rows.values())
