import hashlib
from pathlib import Path

import pytest

from petch.reactor_global import (
    MALYSHEV_1998_ELECTRON_TEMPERATURE_CSV_SHA256,
    MalyshevMeasuredElectronTemperatureProvider,
    REACTOR_SCALAR_EVIDENCE_KINDS,
)


ROOT = Path(__file__).resolve().parents[1]
EXPERIMENTAL_CSV = (
    ROOT / "data" / "experimental" / "malyshev_1998_lam"
    / "figure3_electron_temperature.csv"
)
PACKAGE_CSV = (
    ROOT / "src" / "petch" / "reactor_global" / "data"
    / "malyshev_1998_lam_electron_temperature.csv"
)


def _provider():
    return MalyshevMeasuredElectronTemperatureProvider.from_package_data()


def test_packaged_lam_temperature_board_is_the_audited_data_exactly():
    assert PACKAGE_CSV.read_bytes() == EXPERIMENTAL_CSV.read_bytes()
    assert hashlib.sha256(PACKAGE_CSV.read_bytes()).hexdigest() == (
        MALYSHEV_1998_ELECTRON_TEMPERATURE_CSV_SHA256
    )
    assert len(_provider().markers) == 62
    assert "interpolated_measurement" in REACTOR_SCALAR_EVIDENCE_KINDS


def test_exact_marker_stays_measured_but_not_predictive_without_uncertainty():
    state = _provider().evaluate(
        window_to_wafer_gap_cm=11.0,
        pressure_mTorr=0.5,
        tcp_source_power_W=300.0,
    )

    assert state.method == "exact_marker"
    assert state.electron_temperature.value == pytest.approx(2.82038)
    assert state.electron_temperature.evidence_kind == "measured"
    assert state.electron_temperature.relative_uncertainty is None
    assert not state.electron_temperature.supports_prediction
    assert not state.supports_prediction


def test_interpolation_is_local_explicit_and_not_promoted_to_measurement():
    provider = _provider()
    state = provider.evaluate(
        window_to_wafer_gap_cm=11.0,
        pressure_mTorr=0.5,
        tcp_source_power_W=400.0,
    )
    left, right = state.support_markers
    fraction = (
        (400.0 - left.tcp_source_power_W)
        / (right.tcp_source_power_W - left.tcp_source_power_W)
    )
    expected = (
        left.electron_temperature_eV
        + fraction
        * (right.electron_temperature_eV - left.electron_temperature_eV)
    )

    assert state.method == "linear_interpolation"
    assert state.electron_temperature.value == pytest.approx(expected)
    assert state.electron_temperature.evidence_kind == (
        "interpolated_measurement"
    )
    assert not state.electron_temperature.supports_prediction
    assert "linear_interpolation" in state.electron_temperature.source


def test_provider_refuses_gap_pressure_power_extrapolation_and_ambiguity():
    provider = _provider()
    with pytest.raises(ValueError, match="gap/pressure"):
        provider.evaluate(
            window_to_wafer_gap_cm=8.0,
            pressure_mTorr=1.0,
            tcp_source_power_W=300.0,
        )
    with pytest.raises(ValueError, match="outside"):
        provider.evaluate(
            window_to_wafer_gap_cm=11.0,
            pressure_mTorr=0.5,
            tcp_source_power_W=10.0,
        )
    with pytest.raises(ValueError, match="multiple"):
        provider.evaluate(
            window_to_wafer_gap_cm=6.5,
            pressure_mTorr=1.0,
            tcp_source_power_W=900.0,
        )
    with pytest.raises(ValueError, match="ambiguous marker cluster"):
        provider.evaluate(
            window_to_wafer_gap_cm=6.5,
            pressure_mTorr=1.0,
            tcp_source_power_W=850.0,
        )


def test_exact_only_mode_refuses_to_hide_an_interpolation():
    with pytest.raises(ValueError, match="not an exact"):
        _provider().evaluate(
            window_to_wafer_gap_cm=11.0,
            pressure_mTorr=2.0,
            tcp_source_power_W=400.0,
            allow_linear_interpolation=False,
        )
