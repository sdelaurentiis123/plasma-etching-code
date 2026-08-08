"""Measured Lam Alliance electron-state conditioning for chlorine models.

Malyshev et al. measured electron temperature by OES versus forward TCP
power, pressure, and window-to-wafer gap, and reported volume-average electron
density derived from Langmuir-probe analysis.  This module exposes those
markers without turning forward power into absorbed power, a volume average
into a local sheath density, or either observable into a solved state.  Exact
markers remain measurements; a linear value between two unambiguous markers
is explicitly labeled an interpolated measurement.  Pressure and gap are
never interpolated, and overlapping/duplicate marker clusters fail closed.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
from importlib.resources import files
import io
import math

from .chlorine_particle_model import ReactorScalarInput
from .geometry import CylindricalReactor


MALYSHEV_1998_ELECTRON_TEMPERATURE_CSV_SHA256 = (
    "cd66b1bce25739be4ab555203d78f8b93f1af2f690050ee7d008e567bd3acb7a"
)
MALYSHEV_1998_ELECTRON_DENSITY_CSV_SHA256 = (
    "0fe210a294866b5a8cc9f0f09011ab099458162e4c346db7d08c1f756cc5bd91"
)
_TEMPERATURE_PACKAGE_DATA_NAME = (
    "malyshev_1998_lam_electron_temperature.csv")
_DENSITY_PACKAGE_DATA_NAME = "malyshev_1998_lam_electron_density.csv"
_VALID_METHODS = frozenset({"exact_marker", "linear_interpolation"})
MALYSHEV_1998_LAM_RADIUS_M = 0.215
MALYSHEV_1998_LAM_CONTROL_VOLUME_M3 = 0.043
_REPORTED_EFFECTIVE_LENGTH_M = {11.0: 0.036, 6.5: 0.025}


@dataclass(frozen=True)
class MalyshevLamGeometryState:
    """Reported Lam chamber inventory and active cylindrical plasma region."""

    active_geometry: CylindricalReactor
    neutral_control_volume: ReactorScalarInput
    window_to_wafer_gap_cm: float
    reported_effective_length_m: float
    source: str = "malyshev-1998-lam-cl2 apparatus and Eqs. 6-10"

    def __post_init__(self):
        values = (
            self.window_to_wafer_gap_cm,
            self.reported_effective_length_m,
        )
        if (
            not isinstance(self.active_geometry, CylindricalReactor)
            or not isinstance(self.neutral_control_volume, ReactorScalarInput)
            or self.neutral_control_volume.unit != "m3"
            or self.neutral_control_volume.value < self.active_geometry.volume_m3
            or any(not math.isfinite(float(value)) for value in values)
            or any(float(value) <= 0.0 for value in values)
            or not str(self.source).strip()
        ):
            raise ValueError("invalid Malyshev Lam geometry state")

    @property
    def calculated_effective_length_m(self) -> float:
        return float(
            self.active_geometry.volume_m3
            / self.active_geometry.physical_area_m2
        )

    @property
    def active_volume_fraction(self) -> float:
        return float(
            self.active_geometry.volume_m3
            / self.neutral_control_volume.value
        )

    @property
    def supports_prediction(self) -> bool:
        """The source reports no dimensional uncertainties."""
        return False


def malyshev_1998_lam_geometry(
    window_to_wafer_gap_cm: float,
) -> MalyshevLamGeometryState:
    """Return one of the two reported Lam Alliance active geometries."""
    gap = float(window_to_wafer_gap_cm)
    if gap not in _REPORTED_EFFECTIVE_LENGTH_M:
        raise ValueError("Lam gap must be one of the two reported geometries")
    geometry = CylindricalReactor(
        radius_m=MALYSHEV_1998_LAM_RADIUS_M,
        length_m=gap / 100.0,
    )
    return MalyshevLamGeometryState(
        active_geometry=geometry,
        neutral_control_volume=ReactorScalarInput(
            value=MALYSHEV_1998_LAM_CONTROL_VOLUME_M3,
            unit="m3",
            source=(
                "malyshev-1998-lam-cl2 Eq. 6 text: "
                "reported chamber volume 43000 cm3"
            ),
            evidence_kind="reported_equipment",
            relative_uncertainty=None,
        ),
        window_to_wafer_gap_cm=gap,
        reported_effective_length_m=_REPORTED_EFFECTIVE_LENGTH_M[gap],
    )


@dataclass(frozen=True)
class MalyshevElectronTemperatureMarker:
    """One independently resolved Figure-3 OES marker."""

    window_to_wafer_gap_cm: float
    pressure_mTorr: float
    tcp_source_power_W: float
    electron_temperature_eV: float
    marker: str
    digitization_power_uncertainty_W: float
    digitization_temperature_uncertainty_eV: float

    def __post_init__(self):
        values = (
            self.window_to_wafer_gap_cm,
            self.pressure_mTorr,
            self.tcp_source_power_W,
            self.electron_temperature_eV,
            self.digitization_power_uncertainty_W,
            self.digitization_temperature_uncertainty_eV,
        )
        if (
            any(not math.isfinite(float(value)) for value in values)
            or any(float(value) <= 0.0 for value in values)
            or not str(self.marker).strip()
        ):
            raise ValueError("invalid Malyshev electron-temperature marker")


@dataclass(frozen=True)
class ElectronTemperatureConditioningState:
    """A measured or transparently interpolated electron-temperature input."""

    electron_temperature: ReactorScalarInput
    requested_gap_cm: float
    requested_pressure_mTorr: float
    requested_tcp_source_power_W: float
    method: str
    support_markers: tuple[MalyshevElectronTemperatureMarker, ...]
    digitization_power_uncertainty_W: float
    digitization_temperature_uncertainty_eV: float
    reported_measurement_uncertainty: float | None = None

    def __post_init__(self):
        values = (
            self.requested_gap_cm,
            self.requested_pressure_mTorr,
            self.requested_tcp_source_power_W,
            self.digitization_power_uncertainty_W,
            self.digitization_temperature_uncertainty_eV,
        )
        expected_support = 1 if self.method == "exact_marker" else 2
        if (
            not isinstance(self.electron_temperature, ReactorScalarInput)
            or self.electron_temperature.unit != "eV"
            or self.method not in _VALID_METHODS
            or len(self.support_markers) != expected_support
            or any(
                not isinstance(marker, MalyshevElectronTemperatureMarker)
                for marker in self.support_markers
            )
            or any(not math.isfinite(float(value)) for value in values)
            or any(float(value) <= 0.0 for value in values)
            or self.reported_measurement_uncertainty is not None
        ):
            raise ValueError("invalid electron-temperature conditioning state")

    @property
    def supports_prediction(self) -> bool:
        """Figure 3 reports no Te measurement uncertainty."""
        return False


@dataclass(frozen=True)
class MalyshevMeasuredElectronTemperatureProvider:
    """Strict evaluator over the native-pixel-audited Figure-3 marker board."""

    markers: tuple[MalyshevElectronTemperatureMarker, ...]
    name: str = "malyshev_1998_measured_electron_temperature"
    version: str = "1"

    def __post_init__(self):
        markers = tuple(self.markers)
        if (
            not markers
            or any(
                not isinstance(marker, MalyshevElectronTemperatureMarker)
                for marker in markers
            )
            or not str(self.name).strip()
            or not str(self.version).strip()
        ):
            raise ValueError("invalid measured electron-temperature provider")
        object.__setattr__(self, "markers", markers)

    @classmethod
    def from_package_data(cls) -> "MalyshevMeasuredElectronTemperatureProvider":
        payload = files(__package__).joinpath(
            "data", _TEMPERATURE_PACKAGE_DATA_NAME).read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        if digest != MALYSHEV_1998_ELECTRON_TEMPERATURE_CSV_SHA256:
            raise RuntimeError(
                "packaged Malyshev electron-temperature data hash mismatch"
            )
        records = list(csv.DictReader(io.StringIO(payload.decode("utf-8"))))
        if len(records) != 62:
            raise RuntimeError("incomplete packaged Malyshev marker board")
        markers = []
        for row in records:
            if (
                row["source_figure"] != "Figure 3"
                or row["measurement_method"]
                != "optical_emission_spectroscopy"
                or row["reported_measurement_uncertainty"]
                != "not_reported_in_article"
                or row["tcp_power_semantics"]
                != "power_into_matching_network_not_absorbed_power"
                or row["supports_absorbed_power"] != "false"
                or row["supports_wafer_flux"] != "false"
                or row["validation_role"]
                != "measured_electron_state_conditioning_input"
            ):
                raise RuntimeError("Malyshev marker use boundary is corrupted")
            markers.append(MalyshevElectronTemperatureMarker(
                window_to_wafer_gap_cm=float(
                    row["window_to_wafer_gap_cm"]),
                pressure_mTorr=float(row["pressure_mTorr"]),
                tcp_source_power_W=float(row["tcp_source_power_W"]),
                electron_temperature_eV=float(
                    row["electron_temperature_eV"]),
                marker=row["marker"],
                digitization_power_uncertainty_W=float(
                    row["digitization_power_uncertainty_W"]),
                digitization_temperature_uncertainty_eV=float(
                    row["digitization_temperature_uncertainty_eV"]),
            ))
        return cls(tuple(markers))

    def evaluate(
        self,
        *,
        window_to_wafer_gap_cm: float,
        pressure_mTorr: float,
        tcp_source_power_W: float,
        allow_linear_interpolation: bool = True,
    ) -> ElectronTemperatureConditioningState:
        """Evaluate one supported condition without pressure/gap extrapolation."""
        gap = float(window_to_wafer_gap_cm)
        pressure = float(pressure_mTorr)
        power = float(tcp_source_power_W)
        if (
            not all(math.isfinite(value) for value in (gap, pressure, power))
            or min(gap, pressure, power) <= 0.0
        ):
            raise ValueError("Lam electron-state query must be positive")
        series = tuple(
            marker for marker in self.markers
            if math.isclose(
                marker.window_to_wafer_gap_cm, gap, rel_tol=0.0, abs_tol=1e-12)
            and math.isclose(
                marker.pressure_mTorr, pressure, rel_tol=0.0, abs_tol=1e-12)
        )
        if not series:
            raise ValueError(
                "gap/pressure is outside the measured Figure-3 series"
            )

        exact = tuple(
            marker for marker in series
            if abs(marker.tcp_source_power_W - power)
            <= marker.digitization_power_uncertainty_W
        )
        if len(exact) > 1:
            raise ValueError(
                "multiple Figure-3 markers support this power; "
                "electron temperature is ambiguous"
            )
        if len(exact) == 1:
            marker = exact[0]
            return self._state(
                gap=gap,
                pressure=pressure,
                power=power,
                temperature_eV=marker.electron_temperature_eV,
                method="exact_marker",
                markers=(marker,),
            )
        if not allow_linear_interpolation:
            raise ValueError("requested power is not an exact Figure-3 marker")

        clusters = self._power_clusters(series)
        if power < min(item.tcp_source_power_W for item in series) or power > max(
            item.tcp_source_power_W for item in series
        ):
            raise ValueError("TCP power is outside the measured Figure-3 series")
        bracket = None
        for left, right in zip(clusters[:-1], clusters[1:]):
            left_power = sum(item.tcp_source_power_W for item in left) / len(left)
            right_power = sum(item.tcp_source_power_W for item in right) / len(right)
            if left_power < power < right_power:
                bracket = (left, right)
                break
        if bracket is None:
            raise ValueError("no unambiguous Figure-3 interpolation bracket")
        if len(bracket[0]) != 1 or len(bracket[1]) != 1:
            raise ValueError(
                "Figure-3 interpolation touches an ambiguous marker cluster"
            )
        left, right = bracket[0][0], bracket[1][0]
        fraction = (
            (power - left.tcp_source_power_W)
            / (right.tcp_source_power_W - left.tcp_source_power_W)
        )
        temperature = (
            left.electron_temperature_eV
            + fraction
            * (right.electron_temperature_eV - left.electron_temperature_eV)
        )
        return self._state(
            gap=gap,
            pressure=pressure,
            power=power,
            temperature_eV=temperature,
            method="linear_interpolation",
            markers=(left, right),
        )

    @staticmethod
    def _power_clusters(
        series: tuple[MalyshevElectronTemperatureMarker, ...],
    ) -> tuple[tuple[MalyshevElectronTemperatureMarker, ...], ...]:
        ordered = sorted(series, key=lambda marker: marker.tcp_source_power_W)
        clusters: list[list[MalyshevElectronTemperatureMarker]] = []
        for marker in ordered:
            if not clusters:
                clusters.append([marker])
                continue
            previous = clusters[-1][-1]
            overlap = max(
                previous.digitization_power_uncertainty_W,
                marker.digitization_power_uncertainty_W,
            )
            if marker.tcp_source_power_W - previous.tcp_source_power_W <= overlap:
                clusters[-1].append(marker)
            else:
                clusters.append([marker])
        return tuple(tuple(cluster) for cluster in clusters)

    def _state(
        self,
        *,
        gap: float,
        pressure: float,
        power: float,
        temperature_eV: float,
        method: str,
        markers: tuple[MalyshevElectronTemperatureMarker, ...],
    ) -> ElectronTemperatureConditioningState:
        support = ", ".join(
            f"({marker.tcp_source_power_W:.3f} W, "
            f"{marker.electron_temperature_eV:.5f} eV)"
            for marker in markers
        )
        evidence_kind = (
            "measured" if method == "exact_marker"
            else "interpolated_measurement"
        )
        return ElectronTemperatureConditioningState(
            electron_temperature=ReactorScalarInput(
                value=temperature_eV,
                unit="eV",
                source=(
                    "malyshev-1998-lam-cl2 Figure 3 OES; "
                    f"{method}; support {support}; no reported Te uncertainty"
                ),
                evidence_kind=evidence_kind,
                relative_uncertainty=None,
            ),
            requested_gap_cm=gap,
            requested_pressure_mTorr=pressure,
            requested_tcp_source_power_W=power,
            method=method,
            support_markers=markers,
            digitization_power_uncertainty_W=max(
                marker.digitization_power_uncertainty_W for marker in markers
            ),
            digitization_temperature_uncertainty_eV=max(
                marker.digitization_temperature_uncertainty_eV
                for marker in markers
            ),
            reported_measurement_uncertainty=None,
        )


@dataclass(frozen=True)
class MalyshevElectronDensityMarker:
    """One independently resolved Figure-11 volume-average marker."""

    window_to_wafer_gap_cm: float
    pressure_mTorr: float
    tcp_source_power_W: float
    volume_average_electron_density_cm3: float
    marker: str
    digitization_power_uncertainty_W: float
    digitization_electron_density_uncertainty_cm3: float

    def __post_init__(self):
        values = (
            self.window_to_wafer_gap_cm,
            self.pressure_mTorr,
            self.tcp_source_power_W,
            self.volume_average_electron_density_cm3,
            self.digitization_power_uncertainty_W,
            self.digitization_electron_density_uncertainty_cm3,
        )
        if (
            any(not math.isfinite(float(value)) for value in values)
            or any(float(value) <= 0.0 for value in values)
            or not str(self.marker).strip()
        ):
            raise ValueError("invalid Malyshev electron-density marker")

    @property
    def volume_average_electron_density_m3(self) -> float:
        return float(self.volume_average_electron_density_cm3 * 1.0e6)


@dataclass(frozen=True)
class ElectronDensityConditioningState:
    """A measured or interpolated volume-average electron-density input."""

    volume_average_electron_density: ReactorScalarInput
    requested_gap_cm: float
    requested_pressure_mTorr: float
    requested_tcp_source_power_W: float
    method: str
    support_markers: tuple[MalyshevElectronDensityMarker, ...]
    digitization_power_uncertainty_W: float
    digitization_electron_density_uncertainty_m3: float
    reported_measurement_uncertainty: float | None = None

    def __post_init__(self):
        values = (
            self.requested_gap_cm,
            self.requested_pressure_mTorr,
            self.requested_tcp_source_power_W,
            self.digitization_power_uncertainty_W,
            self.digitization_electron_density_uncertainty_m3,
        )
        expected_support = 1 if self.method == "exact_marker" else 2
        if (
            not isinstance(
                self.volume_average_electron_density, ReactorScalarInput)
            or self.volume_average_electron_density.unit != "m^-3"
            or self.method not in _VALID_METHODS
            or len(self.support_markers) != expected_support
            or any(
                not isinstance(marker, MalyshevElectronDensityMarker)
                for marker in self.support_markers
            )
            or any(not math.isfinite(float(value)) for value in values)
            or any(float(value) <= 0.0 for value in values)
            or self.reported_measurement_uncertainty is not None
        ):
            raise ValueError("invalid electron-density conditioning state")

    @property
    def supports_prediction(self) -> bool:
        """Figure 11 reports no electron-density measurement uncertainty."""
        return False

    @property
    def supports_local_wafer_electron_density(self) -> bool:
        """The source reports a reconstructed volume average, not a local."""
        return False

    @property
    def supports_wafer_flux(self) -> bool:
        return False


@dataclass(frozen=True)
class MalyshevMeasuredElectronDensityProvider:
    """Strict evaluator over the audited Figure-11 volume-average board."""

    markers: tuple[MalyshevElectronDensityMarker, ...]
    name: str = "malyshev_1998_measured_volume_average_electron_density"
    version: str = "1"

    def __post_init__(self):
        markers = tuple(self.markers)
        if (
            not markers
            or any(
                not isinstance(marker, MalyshevElectronDensityMarker)
                for marker in markers
            )
            or not str(self.name).strip()
            or not str(self.version).strip()
        ):
            raise ValueError("invalid measured electron-density provider")
        object.__setattr__(self, "markers", markers)

    @classmethod
    def from_package_data(cls) -> "MalyshevMeasuredElectronDensityProvider":
        payload = files(__package__).joinpath(
            "data", _DENSITY_PACKAGE_DATA_NAME).read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        if digest != MALYSHEV_1998_ELECTRON_DENSITY_CSV_SHA256:
            raise RuntimeError(
                "packaged Malyshev electron-density data hash mismatch"
            )
        records = list(csv.DictReader(io.StringIO(payload.decode("utf-8"))))
        if len(records) != 27:
            raise RuntimeError("incomplete packaged Malyshev density board")
        markers = []
        for row in records:
            if (
                row["source_figure"] != "Figure 11"
                or row["measurement_method"]
                != "Langmuir_probe_analysis_reported_elsewhere"
                or row["volume_average_conversion"]
                != "radial_symmetry_and_axial_sin_pi_h_over_gap"
                or row["reported_measurement_uncertainty"]
                != "not_reported_in_article"
                or row["tcp_power_semantics"]
                != "power_into_matching_network_not_absorbed_power"
                or row["supports_local_wafer_electron_density"] != "false"
                or row["supports_wafer_flux"] != "false"
                or row["validation_role"]
                != "measured_volume_average_electron_state_conditioning_input"
            ):
                raise RuntimeError(
                    "Malyshev density-marker use boundary is corrupted")
            markers.append(MalyshevElectronDensityMarker(
                window_to_wafer_gap_cm=float(
                    row["window_to_wafer_gap_cm"]),
                pressure_mTorr=float(row["pressure_mTorr"]),
                tcp_source_power_W=float(row["tcp_source_power_W"]),
                volume_average_electron_density_cm3=float(
                    row["volume_average_electron_density_cm3"]),
                marker=row["marker"],
                digitization_power_uncertainty_W=float(
                    row["digitization_power_uncertainty_W"]),
                digitization_electron_density_uncertainty_cm3=float(
                    row[
                        "digitization_electron_density_uncertainty_cm3"]),
            ))
        return cls(tuple(markers))

    def evaluate(
        self,
        *,
        window_to_wafer_gap_cm: float,
        pressure_mTorr: float,
        tcp_source_power_W: float,
        allow_linear_interpolation: bool = True,
    ) -> ElectronDensityConditioningState:
        """Evaluate only a fixed measured gap/pressure series."""
        gap = float(window_to_wafer_gap_cm)
        pressure = float(pressure_mTorr)
        power = float(tcp_source_power_W)
        if (
            not all(math.isfinite(value) for value in (gap, pressure, power))
            or min(gap, pressure, power) <= 0.0
        ):
            raise ValueError("Lam electron-density query must be positive")
        series = tuple(
            marker for marker in self.markers
            if math.isclose(
                marker.window_to_wafer_gap_cm, gap,
                rel_tol=0.0, abs_tol=1e-12)
            and math.isclose(
                marker.pressure_mTorr, pressure,
                rel_tol=0.0, abs_tol=1e-12)
        )
        if not series:
            raise ValueError(
                "gap/pressure is outside the measured Figure-11 series")

        exact = tuple(
            marker for marker in series
            if abs(marker.tcp_source_power_W - power)
            <= marker.digitization_power_uncertainty_W
        )
        if len(exact) > 1:
            raise ValueError(
                "multiple Figure-11 markers support this power; "
                "electron density is ambiguous"
            )
        if len(exact) == 1:
            marker = exact[0]
            return self._state(
                gap=gap,
                pressure=pressure,
                power=power,
                density_m3=marker.volume_average_electron_density_m3,
                method="exact_marker",
                markers=(marker,),
            )
        if not allow_linear_interpolation:
            raise ValueError("requested power is not an exact Figure-11 marker")

        ordered = sorted(series, key=lambda marker: marker.tcp_source_power_W)
        if (
            power < ordered[0].tcp_source_power_W
            or power > ordered[-1].tcp_source_power_W
        ):
            raise ValueError(
                "TCP power is outside the measured Figure-11 series")
        bracket = next(
            (
                (left, right)
                for left, right in zip(ordered[:-1], ordered[1:])
                if left.tcp_source_power_W < power < right.tcp_source_power_W
            ),
            None,
        )
        if bracket is None:
            raise ValueError("no Figure-11 interpolation bracket")
        left, right = bracket
        fraction = (
            (power - left.tcp_source_power_W)
            / (right.tcp_source_power_W - left.tcp_source_power_W)
        )
        density_m3 = (
            left.volume_average_electron_density_m3
            + fraction
            * (
                right.volume_average_electron_density_m3
                - left.volume_average_electron_density_m3
            )
        )
        return self._state(
            gap=gap,
            pressure=pressure,
            power=power,
            density_m3=density_m3,
            method="linear_interpolation",
            markers=(left, right),
        )

    @staticmethod
    def _state(
        *,
        gap: float,
        pressure: float,
        power: float,
        density_m3: float,
        method: str,
        markers: tuple[MalyshevElectronDensityMarker, ...],
    ) -> ElectronDensityConditioningState:
        support = ", ".join(
            f"({marker.tcp_source_power_W:.3f} W, "
            f"{marker.volume_average_electron_density_cm3:.6e} cm^-3)"
            for marker in markers
        )
        evidence_kind = (
            "measured" if method == "exact_marker"
            else "interpolated_measurement"
        )
        return ElectronDensityConditioningState(
            volume_average_electron_density=ReactorScalarInput(
                value=density_m3,
                unit="m^-3",
                source=(
                    "malyshev-1998-lam-cl2 Figure 11 Langmuir analysis; "
                    "volume average from radial symmetry and axial "
                    f"sin(pi*h/gap); {method}; support {support}; "
                    "no reported ne uncertainty"
                ),
                evidence_kind=evidence_kind,
                relative_uncertainty=None,
            ),
            requested_gap_cm=gap,
            requested_pressure_mTorr=pressure,
            requested_tcp_source_power_W=power,
            method=method,
            support_markers=markers,
            digitization_power_uncertainty_W=max(
                marker.digitization_power_uncertainty_W
                for marker in markers
            ),
            digitization_electron_density_uncertainty_m3=max(
                marker.digitization_electron_density_uncertainty_cm3
                for marker in markers
            ) * 1.0e6,
            reported_measurement_uncertainty=None,
        )
