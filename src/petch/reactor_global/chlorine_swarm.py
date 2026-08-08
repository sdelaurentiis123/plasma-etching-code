"""Hash-locked direct electron-swarm measurements in pure chlorine.

These observations validate an electron-collision cross-section set. They do
not measure an operating plasma reactor, absorbed power, species-resolved
wafer flux, or feature depth.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
from importlib.resources import files
import io
import math


GONZALEZ_MAGANA_2018_PURE_CL2_SWARM_CSV_SHA256 = (
    "a2bd033aeacaecb302a3fa7c3d3892b478a1ce119e0d7890dfaa4fd66b783d5c"
)
_PACKAGE_DATA_NAME = "gonzalez_magana_2018_pure_cl2_swarm.csv"
_UNITS = {
    "electron_drift_velocity": "m s^-1",
    "effective_ionization_coefficient": "m^2",
    "density_normalized_longitudinal_diffusion": "m^-1 s^-1",
}


@dataclass(frozen=True)
class ChlorineSwarmMeasurement:
    """One pixel-audited, printed pure-Cl2 swarm marker."""

    observation_id: str
    observable: str
    reduced_field_Td: float
    value_si: float
    si_unit: str
    relative_uncertainty_min: float
    relative_uncertainty_max: float
    source_table: str
    source_pdf_page: int
    source_print_page: int
    measurement_method: str
    gas_temperature_K_min: float
    gas_temperature_K_max: float
    pressure_Torr_min: float
    pressure_Torr_max: float

    def __post_init__(self):
        finite = (
            self.reduced_field_Td,
            self.value_si,
            self.relative_uncertainty_min,
            self.relative_uncertainty_max,
            self.gas_temperature_K_min,
            self.gas_temperature_K_max,
            self.pressure_Torr_min,
            self.pressure_Torr_max,
        )
        if (
            not str(self.observation_id).strip()
            or self.observable not in _UNITS
            or self.si_unit != _UNITS.get(self.observable)
            or any(not math.isfinite(value) for value in finite)
            or self.reduced_field_Td <= 0.0
            or (
                self.observable != "effective_ionization_coefficient"
                and self.value_si <= 0.0
            )
            or not (
                0.0 < self.relative_uncertainty_min
                <= self.relative_uncertainty_max < 1.0
            )
            or self.gas_temperature_K_min <= 0.0
            or self.gas_temperature_K_max < self.gas_temperature_K_min
            or self.pressure_Torr_min <= 0.0
            or self.pressure_Torr_max < self.pressure_Torr_min
            or self.source_pdf_page <= 0
            or self.source_print_page <= 0
            or not str(self.source_table).strip()
            or self.measurement_method != "pulsed_townsend_transient"
        ):
            raise ValueError("invalid pure-chlorine swarm measurement")

    @property
    def supports_cross_section_validation(self) -> bool:
        return True

    @property
    def supports_reactor_state_prediction(self) -> bool:
        return False

    @property
    def supports_wafer_flux(self) -> bool:
        return False

    @property
    def supports_feature_depth(self) -> bool:
        return False


@dataclass(frozen=True)
class GonzalezMaganaPureChlorineSwarmBoard:
    """Direct, hash-locked transport board for collision-set validation."""

    measurements: tuple[ChlorineSwarmMeasurement, ...]
    source_bibkey: str = "gonzalez-magana-de-urquijo-2018-cl2"
    version: str = "1"

    def __post_init__(self):
        measurements = tuple(self.measurements)
        if (
            len(measurements) != 52
            or any(
                not isinstance(item, ChlorineSwarmMeasurement)
                for item in measurements
            )
            or len({item.observation_id for item in measurements}) != 52
            or not str(self.source_bibkey).strip()
            or not str(self.version).strip()
        ):
            raise ValueError("invalid pure-chlorine swarm board")
        object.__setattr__(self, "measurements", measurements)

    @classmethod
    def from_package_data(cls) -> "GonzalezMaganaPureChlorineSwarmBoard":
        payload = files(__package__).joinpath(
            "data", _PACKAGE_DATA_NAME).read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        if digest != GONZALEZ_MAGANA_2018_PURE_CL2_SWARM_CSV_SHA256:
            raise RuntimeError("packaged chlorine-swarm data hash mismatch")
        records = list(csv.DictReader(io.StringIO(payload.decode("utf-8"))))
        measurements = []
        for row in records:
            if (
                row["supports_cross_section_validation"] != "true"
                or row["supports_reactor_state_prediction"] != "false"
                or row["supports_wafer_flux"] != "false"
                or row["supports_feature_depth"] != "false"
            ):
                raise RuntimeError("chlorine-swarm evidence boundary corrupted")
            measurements.append(ChlorineSwarmMeasurement(
                observation_id=row["observation_id"],
                observable=row["observable"],
                reduced_field_Td=float(row["reduced_field_Td"]),
                value_si=float(row["value_si"]),
                si_unit=row["si_unit"],
                relative_uncertainty_min=float(
                    row["relative_uncertainty_min"]),
                relative_uncertainty_max=float(
                    row["relative_uncertainty_max"]),
                source_table=row["source_table"],
                source_pdf_page=int(row["source_pdf_page"]),
                source_print_page=int(row["source_print_page"]),
                measurement_method=row["measurement_method"],
                gas_temperature_K_min=float(row["gas_temperature_K_min"]),
                gas_temperature_K_max=float(row["gas_temperature_K_max"]),
                pressure_Torr_min=float(row["pressure_Torr_min"]),
                pressure_Torr_max=float(row["pressure_Torr_max"]),
            ))
        return cls(tuple(measurements))

    def for_observable(
        self, observable: str,
    ) -> tuple[ChlorineSwarmMeasurement, ...]:
        if observable not in _UNITS:
            raise ValueError(f"unsupported chlorine-swarm observable {observable!r}")
        return tuple(
            item for item in self.measurements
            if item.observable == observable
        )

    @property
    def supports_cross_section_validation(self) -> bool:
        return True

    @property
    def supports_reactor_state_prediction(self) -> bool:
        return False

    @property
    def supports_wafer_flux(self) -> bool:
        return False

    @property
    def supports_feature_depth(self) -> bool:
        return False
