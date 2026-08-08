"""Measured-state inversion of Lam Alliance chlorine dissociation.

This module implements the fast-reaction limit printed as Eq. 7 by Malyshev
et al. It combines independently digitized Cl2 density, electron temperature,
and volume-average electron density with the Hamilton state-resolved neutral
dissociation rate and the retained Lee--Lieberman attachment rate. The result
is the first-order wall-return frequency required by that reduced equation.

It is a diagnostic inversion, not a fitted wall probability: Figure 11 is a
reconstructed volume average, Te and ne uncertainties are incomplete, and a
transport model is still required to map frequency to a surface probability.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
from importlib.resources import files
import io
import math

from .chlorine_lam import (
    ElectronDensityConditioningState,
    ElectronTemperatureConditioningState,
    MalyshevMeasuredElectronDensityProvider,
    MalyshevMeasuredElectronTemperatureProvider,
)
from .evaluated_chlorine import (
    build_hamilton_dissociation_chlorine_particle_network,
)
from .network import RateContext


MALYSHEV_1998_CHLORINE_DISSOCIATION_CSV_SHA256 = (
    "df312da3ca72f1424be84eef00a17488298ef09bd34bedfdc87c2c926fcd9540"
)
_PACKAGE_DATA_NAME = "malyshev_1998_lam_chlorine_dissociation.csv"
_VALIDATION_ROLES = frozenset({
    "reactor_dissociation_validation_candidate",
    "diagnostic_flow_check",
})


@dataclass(frozen=True)
class MalyshevChlorineDissociationMarker:
    """One audited Figure-7/8 relative-Cl2 measurement."""

    source_figure: str
    window_to_wafer_gap_cm: float
    pressure_mTorr: float
    tcp_source_power_W: float
    relative_cl2_density_percent: float
    cl2_dissociation_percent: float
    cl2_flow_sccm: float
    rare_gas_flow_sccm: float
    flow_condition: str
    marker: str
    digitization_power_uncertainty_W: float
    digitization_relative_cl2_uncertainty_percentage_point: float
    reported_absolute_density_relative_uncertainty_percent: float
    validation_role: str

    def __post_init__(self):
        values = (
            self.window_to_wafer_gap_cm,
            self.pressure_mTorr,
            self.tcp_source_power_W,
            self.relative_cl2_density_percent,
            self.cl2_flow_sccm,
            self.rare_gas_flow_sccm,
            self.digitization_power_uncertainty_W,
            self.digitization_relative_cl2_uncertainty_percentage_point,
            self.reported_absolute_density_relative_uncertainty_percent,
        )
        if (
            self.source_figure not in {"Figure 7", "Figure 8"}
            or any(not math.isfinite(float(value)) for value in values)
            or any(float(value) <= 0.0 for value in values)
            or not math.isfinite(float(self.cl2_dissociation_percent))
            or not math.isclose(
                self.cl2_dissociation_percent,
                100.0 - self.relative_cl2_density_percent,
                rel_tol=0.0,
                abs_tol=1.0e-4,
            )
            or not str(self.flow_condition).strip()
            or not str(self.marker).strip()
            or self.validation_role not in _VALIDATION_ROLES
        ):
            raise ValueError("invalid Malyshev chlorine-dissociation marker")

    @property
    def supports_eq7_inversion(self) -> bool:
        return (
            self.validation_role
            == "reactor_dissociation_validation_candidate"
            and 0.0 < self.relative_cl2_density_percent < 100.0
        )


@dataclass(frozen=True)
class MalyshevMeasuredChlorineDissociationProvider:
    """Hash-locked provider for the audited Figures 7--8 marker board."""

    markers: tuple[MalyshevChlorineDissociationMarker, ...]
    name: str = "malyshev_1998_measured_relative_chlorine_density"
    version: str = "1"

    def __post_init__(self):
        markers = tuple(self.markers)
        if (
            not markers
            or any(
                not isinstance(marker, MalyshevChlorineDissociationMarker)
                for marker in markers
            )
            or not str(self.name).strip()
            or not str(self.version).strip()
        ):
            raise ValueError("invalid measured chlorine-dissociation provider")
        object.__setattr__(self, "markers", markers)

    @classmethod
    def from_package_data(
        cls,
    ) -> "MalyshevMeasuredChlorineDissociationProvider":
        payload = files(__package__).joinpath(
            "data", _PACKAGE_DATA_NAME).read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        if digest != MALYSHEV_1998_CHLORINE_DISSOCIATION_CSV_SHA256:
            raise RuntimeError(
                "packaged Malyshev chlorine-dissociation data hash mismatch"
            )
        records = list(csv.DictReader(io.StringIO(payload.decode("utf-8"))))
        if len(records) != 38:
            raise RuntimeError("incomplete packaged Malyshev dissociation board")
        markers = []
        for row in records:
            if (
                row["error_bar_semantics"]
                != "range_between_Ar_and_Xe_reductions_not_sigma"
                or row["tcp_power_semantics"]
                != "power_into_matching_network_not_absorbed_power"
                or row["supports_absorbed_power"] != "false"
                or row["supports_wafer_flux"] != "false"
            ):
                raise RuntimeError(
                    "Malyshev dissociation-marker boundary is corrupted")
            markers.append(MalyshevChlorineDissociationMarker(
                source_figure=row["source_figure"],
                window_to_wafer_gap_cm=float(
                    row["window_to_wafer_gap_cm"]),
                pressure_mTorr=float(row["pressure_mTorr"]),
                tcp_source_power_W=float(row["tcp_source_power_W"]),
                relative_cl2_density_percent=float(
                    row["relative_cl2_density_percent"]),
                cl2_dissociation_percent=float(
                    row["cl2_dissociation_percent"]),
                cl2_flow_sccm=float(row["cl2_flow_sccm"]),
                rare_gas_flow_sccm=float(row["rare_gas_flow_sccm"]),
                flow_condition=row["flow_condition"],
                marker=row["marker"],
                digitization_power_uncertainty_W=float(
                    row["digitization_power_uncertainty_W"]),
                digitization_relative_cl2_uncertainty_percentage_point=float(
                    row[
                        "digitization_relative_cl2_uncertainty_percentage_point"]
                ),
                reported_absolute_density_relative_uncertainty_percent=float(
                    row[
                        "reported_absolute_density_relative_uncertainty_percent"]
                ),
                validation_role=row["validation_role"],
            ))
        return cls(tuple(markers))


@dataclass(frozen=True)
class MalyshevEq7WallReturnInversion:
    """Wall-return frequency implied by one measured-state Eq.-7 closure."""

    dissociation_marker: MalyshevChlorineDissociationMarker
    electron_temperature_state: ElectronTemperatureConditioningState
    electron_density_state: ElectronDensityConditioningState
    hamilton_neutral_dissociation_rate_m3_s: float
    lee_dissociative_attachment_rate_m3_s: float
    electron_driven_cl2_destruction_frequency_s_inv: float
    required_wall_return_frequency_s_inv: float
    cl_to_cl2_number_density_ratio: float
    reported_cl2_uncertainty_lower_frequency_s_inv: float
    reported_cl2_uncertainty_upper_frequency_s_inv: float | None
    method: str = "malyshev_1998_eq7_fast_reaction_limit"

    def __post_init__(self):
        required = (
            self.hamilton_neutral_dissociation_rate_m3_s,
            self.lee_dissociative_attachment_rate_m3_s,
            self.electron_driven_cl2_destruction_frequency_s_inv,
            self.required_wall_return_frequency_s_inv,
            self.cl_to_cl2_number_density_ratio,
            self.reported_cl2_uncertainty_lower_frequency_s_inv,
        )
        if (
            not isinstance(
                self.dissociation_marker,
                MalyshevChlorineDissociationMarker,
            )
            or not self.dissociation_marker.supports_eq7_inversion
            or not isinstance(
                self.electron_temperature_state,
                ElectronTemperatureConditioningState,
            )
            or not isinstance(
                self.electron_density_state,
                ElectronDensityConditioningState,
            )
            or any(not math.isfinite(float(value)) for value in required)
            or any(float(value) <= 0.0 for value in required)
            or (
                self.reported_cl2_uncertainty_upper_frequency_s_inv
                is not None
                and (
                    not math.isfinite(float(
                        self.reported_cl2_uncertainty_upper_frequency_s_inv))
                    or self.reported_cl2_uncertainty_upper_frequency_s_inv
                    <= self.required_wall_return_frequency_s_inv
                )
            )
            or not str(self.method).strip()
        ):
            raise ValueError("invalid Malyshev Eq.-7 wall-return inversion")
        if not math.isclose(
            self.reproduced_relative_cl2_density_percent,
            self.dissociation_marker.relative_cl2_density_percent,
            rel_tol=1.0e-13,
            abs_tol=1.0e-11,
        ):
            raise ValueError("Eq.-7 inversion does not reproduce its marker")

    @property
    def reproduced_relative_cl2_density_percent(self) -> float:
        ratio = (
            self.electron_driven_cl2_destruction_frequency_s_inv
            / (2.0 * self.required_wall_return_frequency_s_inv)
        )
        return float(100.0 / (1.0 + ratio))

    @property
    def supports_prediction(self) -> bool:
        return False

    @property
    def supports_wall_probability_inference(self) -> bool:
        return False

    @property
    def supports_wafer_flux(self) -> bool:
        return False

    @property
    def supports_feature_depth(self) -> bool:
        return False


def malyshev_1998_eq7_wall_return_inversion(
    marker: MalyshevChlorineDissociationMarker,
    *,
    electron_temperature_provider: (
        MalyshevMeasuredElectronTemperatureProvider | None) = None,
    electron_density_provider: (
        MalyshevMeasuredElectronDensityProvider | None) = None,
) -> MalyshevEq7WallReturnInversion:
    """Invert Eq. 7 without tuning a transport coefficient or feature depth."""
    if not isinstance(marker, MalyshevChlorineDissociationMarker):
        raise TypeError("a Malyshev chlorine-dissociation marker is required")
    if not marker.supports_eq7_inversion:
        raise ValueError("marker cannot support the physical Eq.-7 inversion")
    if electron_temperature_provider is None:
        electron_temperature_provider = (
            MalyshevMeasuredElectronTemperatureProvider.from_package_data())
    if electron_density_provider is None:
        electron_density_provider = (
            MalyshevMeasuredElectronDensityProvider.from_package_data())

    query = {
        "window_to_wafer_gap_cm": marker.window_to_wafer_gap_cm,
        "pressure_mTorr": marker.pressure_mTorr,
        "tcp_source_power_W": marker.tcp_source_power_W,
    }
    temperature = electron_temperature_provider.evaluate(**query)
    density = electron_density_provider.evaluate(**query)
    context = RateContext(temperature.electron_temperature.value)
    network = build_hamilton_dissociation_chlorine_particle_network()
    hamilton_rate = sum(
        reaction.rate_coefficient.coefficient_si(context)
        for reaction in network.reactions
        if reaction.name.startswith("e_Cl2_dissociation_")
    )
    attachment_rate = next(
        reaction.rate_coefficient.coefficient_si(context)
        for reaction in network.reactions
        if reaction.name == "e_Cl2_dissociative_attachment"
    )
    destruction_frequency = (
        (hamilton_rate + attachment_rate)
        * density.volume_average_electron_density.value
    )
    relative = marker.relative_cl2_density_percent / 100.0
    required_frequency = (
        destruction_frequency * relative / (2.0 * (1.0 - relative))
    )
    cl_to_cl2_ratio = 2.0 * (1.0 - relative) / relative

    uncertainty = (
        marker.reported_absolute_density_relative_uncertainty_percent / 100.0
    )
    relative_lower = relative * (1.0 - uncertainty)
    relative_upper = relative * (1.0 + uncertainty)
    lower_frequency = (
        destruction_frequency
        * relative_lower
        / (2.0 * (1.0 - relative_lower))
    )
    upper_frequency = None
    if relative_upper < 1.0:
        upper_frequency = (
            destruction_frequency
            * relative_upper
            / (2.0 * (1.0 - relative_upper))
        )

    return MalyshevEq7WallReturnInversion(
        dissociation_marker=marker,
        electron_temperature_state=temperature,
        electron_density_state=density,
        hamilton_neutral_dissociation_rate_m3_s=hamilton_rate,
        lee_dissociative_attachment_rate_m3_s=attachment_rate,
        electron_driven_cl2_destruction_frequency_s_inv=(
            destruction_frequency),
        required_wall_return_frequency_s_inv=required_frequency,
        cl_to_cl2_number_density_ratio=cl_to_cl2_ratio,
        reported_cl2_uncertainty_lower_frequency_s_inv=lower_frequency,
        reported_cl2_uncertainty_upper_frequency_s_inv=upper_frequency,
    )
