import hashlib
from pathlib import Path

import pytest

from petch.reactor_global import (
    MALYSHEV_1998_ELECTRON_TEMPERATURE_CSV_SHA256,
    MALYSHEV_1998_LAM_CONTROL_VOLUME_M3,
    MALYSHEV_1998_LAM_RADIUS_M,
    MalyshevMeasuredElectronTemperatureProvider,
    REACTOR_SCALAR_EVIDENCE_KINDS,
    malyshev_1998_lam_geometry,
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


def test_reported_lam_active_and_control_volumes_are_not_collapsed():
    large = malyshev_1998_lam_geometry(11.0)
    small = malyshev_1998_lam_geometry(6.5)

    assert large.active_geometry.radius_m == MALYSHEV_1998_LAM_RADIUS_M
    assert large.active_geometry.volume_m3 == pytest.approx(
        16.0e-3, rel=2.0e-3)
    assert large.neutral_control_volume.value == (
        MALYSHEV_1998_LAM_CONTROL_VOLUME_M3)
    assert large.neutral_control_volume.evidence_kind == "reported_equipment"
    assert large.calculated_effective_length_m == pytest.approx(
        large.reported_effective_length_m, abs=4.0e-4)
    assert small.calculated_effective_length_m == pytest.approx(
        small.reported_effective_length_m, abs=5.0e-5)
    assert small.active_volume_fraction < large.active_volume_fraction < 1.0
    assert not large.supports_prediction
    with pytest.raises(ValueError, match="two reported"):
        malyshev_1998_lam_geometry(8.0)


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
