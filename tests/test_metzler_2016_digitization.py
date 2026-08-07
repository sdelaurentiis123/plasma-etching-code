import csv
from pathlib import Path

import pytest

from scripts.digitize_metzler_2016_fc_ale import (
    DEPTH_POINTS,
    YIELD_POINTS,
    XPS_POINTS,
    expected_files,
)


ROOT = Path(__file__).parents[1]


def _xps(thickness, observable):
    return [
        point for point in XPS_POINTS
        if point.film_thickness_A == thickness
        and point.observable == observable
    ]


def test_complete_depth_yield_and_xps_marker_boards_are_retained():
    assert len(DEPTH_POINTS) == 42
    assert len(YIELD_POINTS) == 25
    assert len(XPS_POINTS) == 32

    assert {
        (point.energy_eV, point.substrate)
        for point in DEPTH_POINTS
    } == {
        (20, "Si"), (20, "SiO2"),
        (25, "Si"), (25, "SiO2"),
        (30, "Si"), (30, "SiO2"),
    }
    assert {
        (point.film_thickness_A, point.observable)
        for point in XPS_POINTS
    } == {
        (5, "film_F_over_C_from_C1s"),
        (5, "delta_F_over_C_substrate_proxy"),
        (5, "CFx_C1s_intensity_kcps"),
        (5, "F1s_intensity_kcps"),
        (11, "film_F_over_C_from_C1s"),
        (11, "delta_F_over_C_substrate_proxy"),
        (11, "CFx_C1s_intensity_kcps"),
        (11, "F1s_intensity_kcps"),
    }


def test_cycle_normalized_board_preserves_units_and_supply_capacity_split():
    sio2_25 = [
        point for point in YIELD_POINTS
        if point.energy_eV == 25 and point.substrate == "SiO2"
    ]
    si_25 = [
        point for point in YIELD_POINTS
        if point.energy_eV == 25 and point.substrate == "Si"
    ]

    # SiO2 is ion-capacity limited at 25 eV: a roughly fourfold change in
    # available F/ion produces less than a factor-three yield span.
    sio2_supply_span = (
        max(point.fluorine_per_incident_ion for point in sio2_25)
        / min(point.fluorine_per_incident_ion for point in sio2_25)
    )
    sio2_yield_span = (
        max(point.substrate_units_per_incident_ion for point in sio2_25)
        / min(point.substrate_units_per_incident_ion for point in sio2_25)
    )
    assert sio2_supply_span > 3.0
    assert sio2_yield_span < 2.5

    # Si is more supply-sensitive under the same nominal ion energy.
    si_yield_span = (
        max(point.substrate_units_per_incident_ion for point in si_25)
        / min(point.substrate_units_per_incident_ion for point in si_25)
    )
    assert si_yield_span > 5.0


def test_cycle_normalized_duplicate_open_si_marker_is_preserved():
    duplicates = [
        point for point in YIELD_POINTS
        if point.panel == "6.9b"
        and point.energy_eV == 30
        and point.substrate == "Si"
        and point.etch_step_s == 40
    ]
    assert len(duplicates) == 2
    assert {point.replicate for point in duplicates} == {1, 2}


def test_published_duplicate_open_si_marker_is_preserved_not_averaged():
    duplicates = [
        point for point in DEPTH_POINTS
        if point.panel == "6.6"
        and point.energy_eV == 30
        and point.substrate == "Si"
        and point.etch_step_s == 40
        and point.marker_fill == "open"
    ]
    assert len(duplicates) == 2
    assert {point.replicate for point in duplicates} == {1, 2}
    assert abs(
        duplicates[0].etch_depth_A_per_cycle
        - duplicates[1].etch_depth_A_per_cycle
    ) > 0.3


def test_thin_and_thick_films_require_finite_transfer_depth():
    thin_delta = _xps(5, "delta_F_over_C_substrate_proxy")
    thick_delta = _xps(11, "delta_F_over_C_substrate_proxy")
    thin_fc = _xps(5, "film_F_over_C_from_C1s")
    thick_fc = _xps(11, "film_F_over_C_from_C1s")

    assert [point.time_s for point in thin_delta] == [0, 5, 15, 40]
    assert thin_delta[-1].value > 1.0
    assert max(point.value for point in thick_delta) < 0.07
    assert thin_delta[-1].value > (
        10.0 * max(point.value for point in thick_delta)
    )

    assert thin_fc[-1].value < 0.4 * thin_fc[0].value
    assert thick_fc[-1].value < 0.5 * thick_fc[0].value


def test_five_angstrom_si_depth_is_nonmonotone_with_ion_step_duration():
    points = [
        point for point in DEPTH_POINTS
        if point.panel == "6.6"
        and point.energy_eV == 25
        and point.substrate == "Si"
    ]
    by_time = {point.etch_step_s: point.etch_depth_A_per_cycle for point in points}
    assert by_time[40] > by_time[20] > by_time[60]


def test_csv_refuses_absolute_atom_inventory_semantics():
    xps_path = (
        ROOT / "data/experimental/metzler_2016/"
        "figures6_14_6_15_xps_cycle_state.csv"
    )
    with xps_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    assert rows
    assert {row["supports_absolute_areal_atom_inventory"] for row in rows} == {
        "false"
    }
    assert all(
        "absolute" in row["quantity_semantics"]
        or "background_sensitive" in row["quantity_semantics"]
        for row in rows
    )


def test_committed_digitization_is_exactly_replayable():
    for path, expected in expected_files().items():
        assert path.is_relative_to(ROOT)
        assert path.read_text(encoding="utf-8") == expected


def test_depth_values_are_angstrom_scale_not_mislabeled_atomic_yields():
    values = [point.etch_depth_A_per_cycle for point in DEPTH_POINTS]
    assert min(values) == pytest.approx(0.0, abs=0.01)
    assert 6.7 < max(values) < 7.1
