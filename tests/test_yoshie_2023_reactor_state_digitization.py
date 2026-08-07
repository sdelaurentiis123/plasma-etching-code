import csv
from pathlib import Path

import pytest

from scripts.digitize_yoshie_2023_reactor_state import (
    BIAS_WINDOWS,
    OES_POINTS,
    expected_files,
)


ROOT = Path(__file__).parents[1]


def _window(cycle, timing):
    return next(
        item for item in BIAS_WINDOWS
        if item.cycle_duration_s == cycle and item.timing == timing
    )


def _oes(signal):
    return [item for item in OES_POINTS if item.signal == signal]


def test_all_published_bias_windows_are_retained_without_calling_ne_ion_flux():
    assert len(BIAS_WINDOWS) == 7
    assert {
        (item.cycle_duration_s, item.timing) for item in BIAS_WINDOWS
    } == {
        (4, "I"), (4, "II"), (4, "III"),
        (8, "I"), (8, "II"), (8, "III"), (8, "IV"),
    }

    payload = expected_files()[
        ROOT / "data/experimental/yoshie_2023/"
        "figure12_bias_window_electron_density.csv"
    ]
    assert "bulk_electron_density_window_not_positive_ion_flux" in payload
    assert ",false,B_facility_diagnostic," in payload


def test_fast_eight_second_timing_ii_transition_is_integrated_not_sampled_once():
    timing_ii = _window(8, "II")
    start, midpoint, end = timing_ii.sample_density_1e10_cm3

    assert start > midpoint > end
    assert start / end > 4.0
    assert timing_ii.simpson_average_1e10_cm3 == pytest.approx(1.0737, abs=5e-5)
    assert timing_ii.window_uncertainty_1e10_cm3 == 0.15


def test_oes_board_preserves_measured_phase_but_is_not_a_flux_deck():
    assert len(OES_POINTS) == 60
    assert len(_oes("CF")) == len(_oes("CF2")) == len(_oes("F")) == 20

    cf_peak = max(_oes("CF"), key=lambda item: item.ratio_to_ar)
    cf2_peak = max(_oes("CF2"), key=lambda item: item.ratio_to_ar)
    f_peak = max(_oes("F"), key=lambda item: item.ratio_to_ar)

    assert cf_peak.time_s == pytest.approx(-1.5)
    assert cf2_peak.time_s == pytest.approx(-1.5)
    assert f_peak.time_s == pytest.approx(0.7)

    payload = expected_files()[
        ROOT / "data/experimental/yoshie_2023/"
        "figure14_phase_resolved_oes.csv"
    ]
    assert "optical_emission_ratio_not_ground_state_density_or_flux" in payload
    assert "supports_absolute_ground_state_flux" in payload


def test_committed_reactor_state_digitization_is_exactly_replayable():
    for path, expected in expected_files().items():
        assert path.is_relative_to(ROOT)
        assert path.read_text(encoding="utf-8") == expected


def test_published_xps_table_requires_persistent_c_f_s_surface_state():
    path = (
        ROOT / "data/experimental/yoshie_2023/"
        "table1_xps_surface_composition.csv"
    )
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    assert len(rows) == 19
    elemental = {
        row["chemical_state"]: (
            float(row["timing_II_percent"]),
            float(row["timing_III_percent"]),
        )
        for row in rows
        if row["quantity_group"] == "elemental_composition"
    }
    assert elemental == {
        "S": (32.2, 18.5),
        "C": (39.5, 45.8),
        "F": (28.3, 35.7),
    }
    assert sum(value[0] for value in elemental.values()) == pytest.approx(100.0)
    assert sum(value[1] for value in elemental.values()) == pytest.approx(100.0)

    states = [
        row for row in rows if row["quantity_group"] == "chemical_state_fraction"
    ]
    for orbital in {"S 2p", "C 1s", "F 1s", "Si 2p"}:
        partition = [row for row in states if row["orbital"] == orbital]
        assert sum(float(row["timing_II_percent"]) for row in partition) == (
            pytest.approx(100.0))
        assert sum(float(row["timing_III_percent"]) for row in partition) == (
            pytest.approx(100.0))

    assert any(
        row["chemical_state"] == "S-CF"
        and float(row["timing_III_percent"])
        > float(row["timing_II_percent"])
        for row in states
    )
