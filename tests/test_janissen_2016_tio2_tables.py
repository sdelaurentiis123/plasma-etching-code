import csv

from scripts.extract_janissen_2016_tio2_tables import (
    MANIFEST_PATH,
    TABLE_S31_CSV,
    TABLE_S32_CSV,
    TABLE_S33_CSV,
    manifest_text,
    outputs,
)


def _rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_committed_janissen_tables_replay_exactly():
    payloads = outputs()
    for path, expected in payloads.items():
        assert path.read_text(encoding="utf-8") == expected
    assert MANIFEST_PATH.read_text(encoding="utf-8") == manifest_text(payloads)


def test_power_and_pressure_sweeps_preserve_printed_batch_pairing():
    rows = {row["sample"]: row for row in _rows(TABLE_S31_CSV)}
    assert [float(rows[key]["tio2_etch_rate_nm_min"])
            for key in ("Ra3", "Ra2", "Ra1")] == [30.0, 58.0, 68.0]
    assert [float(rows[key]["rf_power_W"])
            for key in ("Ra3", "Ra2", "Ra1")] == [100.0, 165.0, 200.0]
    assert float(rows["Ra4"]["pressure_ubar"]) == 10.0
    assert float(rows["Ra4"]["tio2_to_cr_selectivity"]) == 4.0
    assert float(rows["Ra1"]["tio2_to_cr_selectivity"]) == 18.0


def test_closest_stack_witness_keeps_machine_and_material_boundaries_visible():
    rows = {row["figure"]: row for row in _rows(TABLE_S32_CSV)}
    witness = rows["3.2a"]
    assert witness["system_model"] == "Fluor Z401S"
    assert float(witness["mask_height_nm"]) == 45.0
    assert float(witness["mask_diameter_nm"]) == 175.0
    assert float(witness["rf_power_W"]) == 200.0
    assert float(witness["pressure_ubar"]) == 50.0
    assert float(witness["dc_bias_V_signed"]) == -950.0
    assert float(witness["etch_time_min"]) == 11.0


def test_feature_depth_rows_are_measured_outputs_not_surface_coefficients():
    rows = {row["figure"]: row for row in _rows(TABLE_S33_CSV)}
    assert float(rows["3.3"]["average_height_nm"]) == 652.0
    assert float(rows["3.3"]["height_global_rsd_percent"]) == 1.4
    assert float(rows["S3.4"]["average_height_nm"]) == 273.0
    assert float(rows["S3.4"]["height_global_rsd_percent"]) == 3.1
