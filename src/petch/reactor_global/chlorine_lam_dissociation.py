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

from scipy.optimize import brentq

from .chlorine_lam import (
    ElectronDensityConditioningState,
    ElectronTemperatureConditioningState,
    MalyshevLamGeometryState,
    MalyshevMeasuredElectronDensityProvider,
    MalyshevMeasuredElectronTemperatureProvider,
    malyshev_1998_lam_geometry,
)
from .chlorine_transport import (
    malyshev_1998_chlorine_in_chlorine_diffusivity,
)
from .chlorine_wall import (
    BOLTZMANN_J_K,
    ChlorineIncidentVelocityState,
    thermalized_chlorine_incident_velocity_state,
)
from .evaluated_chlorine import (
    build_hamilton_dissociation_chlorine_particle_network,
)
from .model import PASCAL_PER_MTORR
from .network import RateContext
from .neutral_transport import (
    CylindricalNeutralWallLoss,
    NeutralDiffusivityState,
    solve_cylindrical_neutral_wall_loss,
)


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
            raise RuntimeError(
                "incomplete packaged Malyshev dissociation board")
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
                        "digitization_relative_cl2_uncertainty_"
                        "percentage_point"]
                ),
                reported_absolute_density_relative_uncertainty_percent=float(
                    row[
                        "reported_absolute_density_relative_"
                        "uncertainty_percent"]
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


@dataclass(frozen=True)
class MalyshevEq7TransportDiagnostic:
    """Exact-cylinder transport mapping at one declared gas temperature.

    The source reports an initial gas/wall temperature of 333 K and says the
    gas heats with power, but it does not publish the powered gas temperature.
    Consequently this object reports a model-conditioned effective wall
    probability at the caller-declared temperature, never a local measured
    probability or predictive boundary condition.
    """

    eq7_inversion: MalyshevEq7WallReturnInversion
    geometry_state: MalyshevLamGeometryState
    gas_temperature_K: float
    gas_temperature_basis: str
    gauge_pressure_Pa: float
    initial_cl2_particle_fraction: float
    particle_pressure_multiplier: float
    bulk_particle_pressure_Pa: float
    bulk_neutral_density_m3: float
    diffusivity_state: NeutralDiffusivityState
    incident_velocity_state: ChlorineIncidentVelocityState
    absorbing_wall_state: CylindricalNeutralWallLoss
    effective_wall_recombination_probability: float | None
    matched_wall_state: CylindricalNeutralWallLoss | None
    status: str

    def __post_init__(self):
        numeric = (
            self.gas_temperature_K,
            self.gauge_pressure_Pa,
            self.initial_cl2_particle_fraction,
            self.particle_pressure_multiplier,
            self.bulk_particle_pressure_Pa,
            self.bulk_neutral_density_m3,
        )
        probability = self.effective_wall_recombination_probability
        if probability is not None:
            probability = float(probability)
        if (
            not isinstance(
                self.eq7_inversion, MalyshevEq7WallReturnInversion)
            or not isinstance(self.geometry_state, MalyshevLamGeometryState)
            or not isinstance(self.diffusivity_state, NeutralDiffusivityState)
            or not isinstance(
                self.incident_velocity_state, ChlorineIncidentVelocityState)
            or not isinstance(
                self.absorbing_wall_state, CylindricalNeutralWallLoss)
            or any(not math.isfinite(float(value)) for value in numeric)
            or any(float(value) <= 0.0 for value in numeric)
            or not str(self.gas_temperature_basis).strip()
            or self.status not in {
                "model_conditioned_effective_probability",
                "target_exceeds_absorbing_wall_limit",
            }
        ):
            raise ValueError("invalid Malyshev Eq.-7 transport diagnostic")
        marker = self.eq7_inversion.dissociation_marker
        expected_initial_cl2_fraction = (
            marker.cl2_flow_sccm
            / (marker.cl2_flow_sccm + marker.rare_gas_flow_sccm)
        )
        if not math.isclose(
            self.initial_cl2_particle_fraction,
            expected_initial_cl2_fraction,
            rel_tol=1.0e-13,
            abs_tol=0.0,
        ):
            raise ValueError("initial Cl2 particle fraction is inconsistent")
        relative_cl2 = marker.relative_cl2_density_percent / 100.0
        expected_multiplier = (
            1.0
            + self.initial_cl2_particle_fraction * (1.0 - relative_cl2)
        )
        if not math.isclose(
            self.particle_pressure_multiplier,
            expected_multiplier,
            rel_tol=1.0e-13,
            abs_tol=0.0,
        ):
            raise ValueError("particle-pressure multiplier violates Eq. 11")
        if not math.isclose(
            self.bulk_particle_pressure_Pa,
            self.gauge_pressure_Pa * self.particle_pressure_multiplier,
            rel_tol=1.0e-13,
            abs_tol=0.0,
        ):
            raise ValueError("bulk particle pressure is inconsistent")
        if not math.isclose(
            self.bulk_neutral_density_m3,
            self.bulk_particle_pressure_Pa
            / (BOLTZMANN_J_K * self.gas_temperature_K),
            rel_tol=1.0e-13,
            abs_tol=0.0,
        ):
            raise ValueError("bulk particle density is inconsistent")
        if not math.isclose(
            self.diffusivity_state.total_neutral_density_m3,
            self.bulk_neutral_density_m3,
            rel_tol=1.0e-13,
            abs_tol=0.0,
        ):
            raise ValueError("transport diffusivity uses a different density")
        if not math.isclose(
            self.diffusivity_state.gas_temperature_K,
            self.gas_temperature_K,
            rel_tol=1.0e-13,
            abs_tol=0.0,
        ):
            raise ValueError(
                "transport diffusivity uses a different temperature")
        if probability is None:
            if (
                self.matched_wall_state is not None
                or self.status != "target_exceeds_absorbing_wall_limit"
                or self.required_wall_return_frequency_s_inv
                <= self.absorbing_wall_state.exact_loss_frequency_s_inv
            ):
                raise ValueError("invalid unattainable transport target")
        else:
            if (
                not 0.0 < probability <= 1.0
                or not isinstance(
                    self.matched_wall_state, CylindricalNeutralWallLoss)
                or self.status
                != "model_conditioned_effective_probability"
                or not math.isclose(
                    self.matched_wall_state.wall_reaction_probability,
                    probability,
                    rel_tol=0.0,
                    abs_tol=1.0e-14,
                )
                or not math.isclose(
                    self.matched_wall_state.exact_loss_frequency_s_inv,
                    self.required_wall_return_frequency_s_inv,
                    rel_tol=2.0e-12,
                    abs_tol=2.0e-9,
                )
            ):
                raise ValueError("transport inversion does not close")
        object.__setattr__(
            self,
            "effective_wall_recombination_probability",
            probability,
        )

    @property
    def required_wall_return_frequency_s_inv(self) -> float:
        return float(
            self.eq7_inversion.required_wall_return_frequency_s_inv)

    @property
    def target_is_transport_attainable(self) -> bool:
        return self.effective_wall_recombination_probability is not None

    @property
    def supports_prediction(self) -> bool:
        return False

    @property
    def supports_local_wall_probability_prediction(self) -> bool:
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


def malyshev_1998_eq7_transport_diagnostic(
    inversion: MalyshevEq7WallReturnInversion,
    *,
    gas_temperature_K: float,
    gas_temperature_basis: str,
) -> MalyshevEq7TransportDiagnostic:
    """Map an Eq.-7 frequency to effective gamma without tuning a target.

    Malyshev Eq. 11 implies that dissociation raises the in-reactor particle
    density above the plasma-off gauge-density basis. The source's explicit
    Cl2 and rare-gas flows retain the non-dissociating 5% actinometry
    inventory, giving ``1 + x_Cl2,0 * (1 - relative_Cl2)`` rather than silently
    applying the pure-Cl2 multiplier to total pressure. The
    source-parameterized Cl-in-Cl2 binary diffusivity, with Neufeld's evaluated
    collision integral, is evaluated at that particle density. The requested
    first-order wall loss is then inverted through the exact fundamental
    cylindrical Robin mode.

    If even a perfectly absorbing wall cannot supply the requested loss, the
    result fails closed with no probability.  Gas temperature is mandatory
    because the source did not measure its powered value.
    """
    if not isinstance(inversion, MalyshevEq7WallReturnInversion):
        raise TypeError("a Malyshev Eq.-7 inversion is required")
    temperature = float(gas_temperature_K)
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("gas temperature must be positive and finite")
    if not str(gas_temperature_basis).strip():
        raise ValueError("gas-temperature evidence basis is required")

    marker = inversion.dissociation_marker
    relative_cl2 = marker.relative_cl2_density_percent / 100.0
    gauge_pressure = marker.pressure_mTorr * PASCAL_PER_MTORR
    initial_cl2_fraction = (
        marker.cl2_flow_sccm
        / (marker.cl2_flow_sccm + marker.rare_gas_flow_sccm)
    )
    particle_multiplier = (
        1.0 + initial_cl2_fraction * (1.0 - relative_cl2))
    bulk_pressure = gauge_pressure * particle_multiplier
    bulk_density = bulk_pressure / (BOLTZMANN_J_K * temperature)
    diffusivity = (
        malyshev_1998_chlorine_in_chlorine_diffusivity().evaluate(
            total_neutral_density_m3=bulk_density,
            gas_temperature_K=temperature,
        )
    )
    velocity = thermalized_chlorine_incident_velocity_state(
        temperature,
        source=(
            "declared-temperature Maxwellian sensitivity for Malyshev "
            "Eq.-7 transport diagnostic"
        ),
        evidence_kind="assumed",
        relative_uncertainty=None,
        provenance={
            "gas_temperature_basis": gas_temperature_basis,
            "source_powered_gas_temperature_measured": False,
            "coefficient_selection_target": None,
            "feature_depth_target": None,
        },
    )
    geometry = malyshev_1998_lam_geometry(
        marker.window_to_wafer_gap_cm)
    solver_inputs = {
        "geometry": geometry.active_geometry,
        "diffusivity_m2_s": diffusivity.diffusivity_m2_s,
        "mean_thermal_speed_m_s": velocity.mean_speed_m_s,
    }
    absorbing = solve_cylindrical_neutral_wall_loss(
        **solver_inputs,
        wall_reaction_probability=1.0,
    )
    target = inversion.required_wall_return_frequency_s_inv
    if target > absorbing.exact_loss_frequency_s_inv:
        probability = None
        matched = None
        status = "target_exceeds_absorbing_wall_limit"
    else:
        probability = float(brentq(
            lambda gamma: (
                solve_cylindrical_neutral_wall_loss(
                    **solver_inputs,
                    wall_reaction_probability=gamma,
                ).exact_loss_frequency_s_inv
                - target
            ),
            0.0,
            1.0,
            xtol=1.0e-14,
            rtol=1.0e-14,
        ))
        matched = solve_cylindrical_neutral_wall_loss(
            **solver_inputs,
            wall_reaction_probability=probability,
        )
        status = "model_conditioned_effective_probability"

    return MalyshevEq7TransportDiagnostic(
        eq7_inversion=inversion,
        geometry_state=geometry,
        gas_temperature_K=temperature,
        gas_temperature_basis=gas_temperature_basis,
        gauge_pressure_Pa=gauge_pressure,
        initial_cl2_particle_fraction=initial_cl2_fraction,
        particle_pressure_multiplier=particle_multiplier,
        bulk_particle_pressure_Pa=bulk_pressure,
        bulk_neutral_density_m3=bulk_density,
        diffusivity_state=diffusivity,
        incident_velocity_state=velocity,
        absorbing_wall_state=absorbing,
        effective_wall_recombination_probability=probability,
        matched_wall_state=matched,
        status=status,
    )
