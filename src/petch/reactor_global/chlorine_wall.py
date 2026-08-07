"""Local, atom-conserving chlorine wall-recombination physics.

This module deliberately stops at the plasma-facing wall.  It converts a
local Cl number density and an explicit incident-velocity state into an
isotropic impingement flux and applies a condition-scoped recombination
probability:

    2 Cl(wall) -> Cl2(gas)

The velocity state is mandatory.  Guha et al. (2008) show that the common
300 K Maxwellian assumption can fail when freshly dissociated Cl reaches the
wall before thermalizing.  Keeping that state explicit prevents a gas
temperature from silently certifying the incident flux.

This module does not turn the local wall flux into a volume-averaged loss
frequency.  That step requires neutral diffusion/Knudsen transport and the
distribution of conditioned surface states.  Omitting it here prevents a
ballistic well-mixed assumption from silently becoming reactor closure.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .chlorine import CHLORINE_ATOM_MASS_AMU
from .transport import ATOMIC_MASS_UNIT_KG

BOLTZMANN_J_K = 1.380649e-23

WALL_RECOMBINATION_EVIDENCE_KINDS = frozenset({
    "measured",
    "validated_model",
    "regressed",
    "published_range_member",
    "assumed",
    "sensitivity",
})
_PREDICTIVE_EVIDENCE_KINDS = frozenset({"measured", "validated_model"})

CHLORINE_INCIDENT_VELOCITY_DISTRIBUTIONS = frozenset({
    "thermalized_maxwellian",
    "measured_isotropic",
    "modeled_isotropic",
    "sensitivity_isotropic",
})
CHLORINE_INCIDENT_VELOCITY_EVIDENCE_KINDS = frozenset({
    "measured",
    "validated_model",
    "assumed",
    "sensitivity",
})


def chlorine_atom_mean_thermal_speed_m_s(
        gas_temperature_K: float) -> float:
    """Return the Maxwellian mean speed ``sqrt(8 kT / (pi m_Cl))``."""
    temperature = float(gas_temperature_K)
    if not np.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("gas temperature must be positive and finite")
    mass_kg = CHLORINE_ATOM_MASS_AMU * ATOMIC_MASS_UNIT_KG
    return float(np.sqrt(
        8.0 * BOLTZMANN_J_K * temperature / (np.pi * mass_kg)))


@dataclass(frozen=True)
class ChlorineIncidentVelocityState:
    """Isotropic Cl velocity moment used by the kinetic wall boundary.

    ``mean_speed_m_s`` is the full speed moment ``<|v|>``.  Isotropy is a
    declared part of every allowed distribution kind, so the incident number
    flux is ``n <|v|> / 4``.  An anisotropic distribution requires its inward
    normal velocity moment and is intentionally outside this scalar closure.
    """

    mean_speed_m_s: float
    distribution_kind: str
    source: str
    evidence_kind: str
    relative_uncertainty: float | None
    reference_temperature_K: float | None = None
    provenance: Mapping[str, object] | None = None

    def __post_init__(self):
        mean_speed = float(self.mean_speed_m_s)
        uncertainty = self.relative_uncertainty
        if uncertainty is not None:
            uncertainty = float(uncertainty)
        reference_temperature = self.reference_temperature_K
        if reference_temperature is not None:
            reference_temperature = float(reference_temperature)
        if (
            not np.isfinite(mean_speed)
            or mean_speed <= 0.0
            or self.distribution_kind
            not in CHLORINE_INCIDENT_VELOCITY_DISTRIBUTIONS
            or not str(self.source).strip()
            or self.evidence_kind
            not in CHLORINE_INCIDENT_VELOCITY_EVIDENCE_KINDS
            or (
                uncertainty is not None
                and (
                    not np.isfinite(uncertainty)
                    or not 0.0 <= uncertainty < 1.0
                )
            )
            or (
                reference_temperature is not None
                and (
                    not np.isfinite(reference_temperature)
                    or reference_temperature <= 0.0
                )
            )
        ):
            raise ValueError("invalid chlorine incident-velocity state")
        if (
            self.distribution_kind == "thermalized_maxwellian"
            and reference_temperature is None
        ):
            raise ValueError(
                "thermalized Maxwellian state requires a reference temperature")
        if (
            self.distribution_kind != "thermalized_maxwellian"
            and reference_temperature is not None
        ):
            raise ValueError(
                "reference temperature is reserved for a Maxwellian state")
        object.__setattr__(self, "mean_speed_m_s", mean_speed)
        object.__setattr__(self, "relative_uncertainty", uncertainty)
        object.__setattr__(
            self, "reference_temperature_K", reference_temperature)
        object.__setattr__(
            self,
            "provenance",
            MappingProxyType(
                {} if self.provenance is None else dict(self.provenance)),
        )

    @property
    def supports_prediction(self) -> bool:
        return (
            self.evidence_kind in _PREDICTIVE_EVIDENCE_KINDS
            and self.relative_uncertainty is not None
        )

    def require_compatible_temperature(self, gas_temperature_K: float) -> None:
        """Reject a thermal state evaluated at a different gas temperature."""
        temperature = float(gas_temperature_K)
        if not np.isfinite(temperature) or temperature <= 0.0:
            raise ValueError("gas temperature must be positive and finite")
        if (
            self.distribution_kind == "thermalized_maxwellian"
            and not np.isclose(
                self.reference_temperature_K,
                temperature,
                rtol=1.0e-12,
                atol=0.0,
            )
        ):
            raise ValueError(
                "incident Maxwellian reference temperature does not match "
                "the reactor gas temperature")


def thermalized_chlorine_incident_velocity_state(
    gas_temperature_K: float,
    *,
    source: str,
    evidence_kind: str,
    relative_uncertainty: float | None,
    provenance: Mapping[str, object] | None = None,
) -> ChlorineIncidentVelocityState:
    """Construct an explicit Maxwellian Cl incident-velocity state."""
    temperature = float(gas_temperature_K)
    return ChlorineIncidentVelocityState(
        mean_speed_m_s=chlorine_atom_mean_thermal_speed_m_s(temperature),
        distribution_kind="thermalized_maxwellian",
        source=source,
        evidence_kind=evidence_kind,
        relative_uncertainty=relative_uncertainty,
        reference_temperature_K=temperature,
        provenance=provenance,
    )


@dataclass(frozen=True)
class ChlorineWallRecombinationBoundary:
    """One condition-scoped local wall-recombination boundary.

    The validity ranges are mandatory because Stafford et al. show that
    ``gamma_Cl`` changes with Cl/Cl2 ratio, pressure, wall conditioning, and
    material.  A scalar with no condition domain is therefore invalid.
    """

    recombination_probability: float
    surface_state: str
    source: str
    evidence_kind: str
    valid_cl_to_cl2_ratio: tuple[float, float]
    valid_pressure_Pa: tuple[float, float]
    valid_icp_power_W: tuple[float, float]
    valid_gas_temperature_K: tuple[float, float]
    relative_measurement_uncertainty: float | None = None
    provenance: Mapping[str, object] | None = None

    def __post_init__(self):
        probability = float(self.recombination_probability)
        ratio_domain = _bounded_domain(
            self.valid_cl_to_cl2_ratio,
            name="Cl/Cl2 ratio",
            allow_zero=True,
        )
        pressure_domain = _bounded_domain(
            self.valid_pressure_Pa,
            name="pressure",
            allow_zero=False,
        )
        power_domain = _bounded_domain(
            self.valid_icp_power_W,
            name="ICP power",
            allow_zero=False,
        )
        temperature_domain = _bounded_domain(
            self.valid_gas_temperature_K,
            name="gas temperature",
            allow_zero=False,
        )
        uncertainty = self.relative_measurement_uncertainty
        if uncertainty is not None:
            uncertainty = float(uncertainty)
        if (
            not np.isfinite(probability)
            or not 0.0 <= probability <= 1.0
            or not str(self.surface_state).strip()
            or not str(self.source).strip()
            or self.evidence_kind not in WALL_RECOMBINATION_EVIDENCE_KINDS
            or (
                uncertainty is not None
                and (
                    not np.isfinite(uncertainty)
                    or not 0.0 <= uncertainty < 1.0
                )
            )
        ):
            raise ValueError("invalid chlorine wall-recombination boundary")
        object.__setattr__(self, "recombination_probability", probability)
        object.__setattr__(self, "valid_cl_to_cl2_ratio", ratio_domain)
        object.__setattr__(self, "valid_pressure_Pa", pressure_domain)
        object.__setattr__(self, "valid_icp_power_W", power_domain)
        object.__setattr__(
            self, "valid_gas_temperature_K", temperature_domain)
        object.__setattr__(
            self, "relative_measurement_uncertainty", uncertainty)
        object.__setattr__(
            self,
            "provenance",
            MappingProxyType(
                {} if self.provenance is None else dict(self.provenance)),
        )

    @property
    def supports_local_prediction(self) -> bool:
        """Whether the local gamma carries predictive-grade evidence.

        This does not certify a volume-averaged reactor prediction.  Missing
        neutral transport and wall-state area fractions remain separate.
        """
        return (
            self.evidence_kind in _PREDICTIVE_EVIDENCE_KINDS
            and self.relative_measurement_uncertainty is not None
        )

    def require_applicable(
        self,
        *,
        cl_to_cl2_ratio: float,
        pressure_Pa: float,
        icp_power_W: float,
        gas_temperature_K: float,
    ) -> None:
        values = {
            "Cl/Cl2 ratio": (
                float(cl_to_cl2_ratio), self.valid_cl_to_cl2_ratio),
            "pressure": (float(pressure_Pa), self.valid_pressure_Pa),
            "ICP power": (float(icp_power_W), self.valid_icp_power_W),
            "gas temperature": (
                float(gas_temperature_K), self.valid_gas_temperature_K),
        }
        for name, (value, domain) in values.items():
            if not np.isfinite(value) or not domain[0] <= value <= domain[1]:
                raise ValueError(
                    f"{name} is outside the wall-boundary evidence domain")

    def evaluate(
        self,
        *,
        chlorine_atom_density_m3: float,
        incident_velocity_state: ChlorineIncidentVelocityState,
        gas_temperature_K: float,
        cl_to_cl2_ratio: float,
        pressure_Pa: float,
        icp_power_W: float,
    ) -> "ChlorineWallFlux":
        self.require_applicable(
            cl_to_cl2_ratio=cl_to_cl2_ratio,
            pressure_Pa=pressure_Pa,
            icp_power_W=icp_power_W,
            gas_temperature_K=gas_temperature_K,
        )
        density = float(chlorine_atom_density_m3)
        if not isinstance(
            incident_velocity_state, ChlorineIncidentVelocityState
        ):
            raise TypeError("chlorine incident-velocity state is required")
        incident_velocity_state.require_compatible_temperature(
            gas_temperature_K)
        if not np.isfinite(density) or density < 0.0:
            raise ValueError(
                "chlorine atom density must be finite and nonnegative")
        mean_speed = incident_velocity_state.mean_speed_m_s
        incident = 0.25 * density * mean_speed
        recombined_atoms = self.recombination_probability * incident
        returned_molecules = 0.5 * recombined_atoms
        return ChlorineWallFlux(
            incident_cl_atom_flux_m2_s=incident,
            recombined_cl_atom_flux_m2_s=recombined_atoms,
            returned_cl2_molecule_flux_m2_s=returned_molecules,
            mean_cl_speed_m_s=mean_speed,
            source=self.source,
            evidence_kind=self.evidence_kind,
            incident_velocity_state=incident_velocity_state,
            wall_boundary=self,
        )


@dataclass(frozen=True)
class LogLinearChlorineWallRecombinationProvider:
    """In-domain ``log10(gamma) = intercept + slope * nCl/nCl2`` law.

    This is a regression of direct wall measurements, not a site-kinetics
    mechanism. It therefore exposes fit residuals and returns no experimental
    uncertainty when the source figure provides none.
    """

    slope_per_ratio: float
    intercept_log10: float
    surface_state: str
    source: str
    valid_cl_to_cl2_ratio: tuple[float, float]
    valid_pressure_Pa: tuple[float, float]
    valid_icp_power_W: tuple[float, float]
    valid_gas_temperature_K: tuple[float, float]
    marker_count: int
    fit_rmse_log10: float
    fit_maximum_absolute_residual_log10: float
    leave_one_out_rmse_log10: float
    leave_one_out_maximum_absolute_residual_log10: float
    provenance: Mapping[str, object] | None = None
    name: str = "log_linear_chlorine_wall_recombination"
    version: str = "1"

    def __post_init__(self):
        scalar_values = np.asarray([
            self.slope_per_ratio,
            self.intercept_log10,
            self.fit_rmse_log10,
            self.fit_maximum_absolute_residual_log10,
            self.leave_one_out_rmse_log10,
            self.leave_one_out_maximum_absolute_residual_log10,
        ], dtype=float)
        if (
            np.any(~np.isfinite(scalar_values))
            or self.slope_per_ratio <= 0.0
            or np.any(scalar_values[2:] < 0.0)
            or int(self.marker_count) != self.marker_count
            or int(self.marker_count) < 3
            or not str(self.surface_state).strip()
            or not str(self.source).strip()
            or not str(self.name).strip()
            or not str(self.version).strip()
        ):
            raise ValueError("invalid chlorine wall regression provider")
        ratio_domain = _bounded_domain(
            self.valid_cl_to_cl2_ratio,
            name="Cl/Cl2 ratio",
            allow_zero=True,
        )
        pressure_domain = _bounded_domain(
            self.valid_pressure_Pa,
            name="pressure",
            allow_zero=False,
        )
        power_domain = _bounded_domain(
            self.valid_icp_power_W,
            name="ICP power",
            allow_zero=False,
        )
        temperature_domain = _bounded_domain(
            self.valid_gas_temperature_K,
            name="gas temperature",
            allow_zero=False,
        )
        endpoint_probabilities = 10.0 ** np.asarray([
            self.intercept_log10 + self.slope_per_ratio * ratio_domain[0],
            self.intercept_log10 + self.slope_per_ratio * ratio_domain[1],
        ])
        if (
            np.any(~np.isfinite(endpoint_probabilities))
            or np.any(endpoint_probabilities <= 0.0)
            or np.any(endpoint_probabilities > 1.0)
        ):
            raise ValueError("wall regression leaves the probability domain")
        object.__setattr__(self, "slope_per_ratio", scalar_values[0])
        object.__setattr__(self, "intercept_log10", scalar_values[1])
        object.__setattr__(self, "fit_rmse_log10", scalar_values[2])
        object.__setattr__(
            self,
            "fit_maximum_absolute_residual_log10",
            scalar_values[3],
        )
        object.__setattr__(
            self, "leave_one_out_rmse_log10", scalar_values[4])
        object.__setattr__(
            self,
            "leave_one_out_maximum_absolute_residual_log10",
            scalar_values[5],
        )
        object.__setattr__(self, "marker_count", int(self.marker_count))
        object.__setattr__(self, "valid_cl_to_cl2_ratio", ratio_domain)
        object.__setattr__(self, "valid_pressure_Pa", pressure_domain)
        object.__setattr__(self, "valid_icp_power_W", power_domain)
        object.__setattr__(
            self, "valid_gas_temperature_K", temperature_domain)
        object.__setattr__(
            self,
            "provenance",
            MappingProxyType(
                {} if self.provenance is None else dict(self.provenance)),
        )

    @property
    def supports_prediction(self) -> bool:
        """Direct-data fit lacks source measurement uncertainty/site state."""
        return False

    def predict(
        self,
        *,
        cl_to_cl2_ratio: float,
        pressure_Pa: float,
        icp_power_W: float,
        gas_temperature_K: float,
    ) -> ChlorineWallRecombinationBoundary:
        ratio = float(cl_to_cl2_ratio)
        values = {
            "Cl/Cl2 ratio": (ratio, self.valid_cl_to_cl2_ratio),
            "pressure": (float(pressure_Pa), self.valid_pressure_Pa),
            "ICP power": (float(icp_power_W), self.valid_icp_power_W),
            "gas temperature": (
                float(gas_temperature_K),
                self.valid_gas_temperature_K,
            ),
        }
        for quantity, (value, domain) in values.items():
            if not np.isfinite(value) or not domain[0] <= value <= domain[1]:
                raise ValueError(
                    f"{quantity} is outside the wall-provider evidence domain")
        boundary = ChlorineWallRecombinationBoundary(
            recombination_probability=(
                10.0 ** (
                    self.intercept_log10 + self.slope_per_ratio * ratio)
            ),
            surface_state=self.surface_state,
            source=self.source,
            evidence_kind="regressed",
            valid_cl_to_cl2_ratio=self.valid_cl_to_cl2_ratio,
            valid_pressure_Pa=self.valid_pressure_Pa,
            valid_icp_power_W=self.valid_icp_power_W,
            valid_gas_temperature_K=self.valid_gas_temperature_K,
            relative_measurement_uncertainty=None,
            provenance={
                **self.provenance,
                "fit_form": "log10(gamma) = intercept + slope * nCl/nCl2",
                "slope_per_ratio": self.slope_per_ratio,
                "intercept_log10": self.intercept_log10,
                "marker_count": self.marker_count,
                "fit_rmse_log10": self.fit_rmse_log10,
                "fit_maximum_absolute_residual_log10": (
                    self.fit_maximum_absolute_residual_log10),
                "leave_one_out_rmse_log10": (
                    self.leave_one_out_rmse_log10),
                "leave_one_out_maximum_absolute_residual_log10": (
                    self.leave_one_out_maximum_absolute_residual_log10),
                "experimental_uncertainty": (
                    "not reported in source Figure 8; digitization error "
                    "is not substituted for measurement uncertainty"
                ),
                "coefficient_selection_target": (
                    "direct wall gamma markers only; no reactor, feature, "
                    "or depth observable"
                ),
            },
        )
        boundary.require_applicable(
            cl_to_cl2_ratio=ratio,
            pressure_Pa=pressure_Pa,
            icp_power_W=icp_power_W,
            gas_temperature_K=gas_temperature_K,
        )
        return boundary


def stafford_2010_conditioned_wall_recombination_provider(
    material: str,
) -> LogLinearChlorineWallRecombinationProvider:
    """Return the exact unweighted Figure-8 first-order data regression."""
    material_key = str(material).strip().lower()
    common = {
        "valid_pressure_Pa": (
            1.25 * 0.1333223684,
            20.0 * 0.1333223684,
        ),
        "valid_icp_power_W": (100.0, 600.0),
        "valid_gas_temperature_K": (300.0, 300.0),
        "source": (
            "stafford-2010-cl-wall Figure 8 direct spinning-wall markers; "
            "unweighted first-order least-squares fit in log10(gamma)"
        ),
        "provenance": {
            "doi": "10.1351/PAC-CON-09-11-02",
            "digitized_dataset": (
                "data/experimental/stafford_2010/"
                "figure8_chlorine_wall_recombination.csv"
            ),
            "digitization_ratio_uncertainty": 0.0027924294135898233,
            "digitization_log10_gamma_uncertainty": 0.012326656394453005,
            "individual_power_per_marker": "not reported in Figure 8",
        },
        "name": "stafford_2010_conditioned_wall_log_linear",
        "version": "1",
    }
    if material_key == "anodized_aluminum":
        return LogLinearChlorineWallRecombinationProvider(
            slope_per_ratio=1.335091971549131,
            intercept_log10=-2.1773737439464496,
            surface_state=(
                "plasma-conditioned anodized aluminum with Si-oxychloride"
            ),
            valid_cl_to_cl2_ratio=(0.105610, 0.779646),
            marker_count=23,
            fit_rmse_log10=0.13803939372878862,
            fit_maximum_absolute_residual_log10=0.3756450974608243,
            leave_one_out_rmse_log10=0.1529041548822586,
            leave_one_out_maximum_absolute_residual_log10=(
                0.4026556294460728),
            **common,
        )
    if material_key == "stainless_steel":
        return LogLinearChlorineWallRecombinationProvider(
            slope_per_ratio=1.0928756558958963,
            intercept_log10=-2.4339733013961156,
            surface_state=(
                "plasma-conditioned stainless steel with Si-oxychloride"
            ),
            valid_cl_to_cl2_ratio=(0.105721, 0.779088),
            marker_count=16,
            fit_rmse_log10=0.11656559337098579,
            fit_maximum_absolute_residual_log10=0.29980451898009086,
            leave_one_out_rmse_log10=0.13710284271067782,
            leave_one_out_maximum_absolute_residual_log10=(
                0.3430112768690421),
            **common,
        )
    raise ValueError(
        "material must be 'anodized_aluminum' or 'stainless_steel'")


@dataclass(frozen=True)
class ChlorineWallFlux:
    """Species fluxes at one plasma-facing wall location."""

    incident_cl_atom_flux_m2_s: float
    recombined_cl_atom_flux_m2_s: float
    returned_cl2_molecule_flux_m2_s: float
    mean_cl_speed_m_s: float
    source: str
    evidence_kind: str
    incident_velocity_state: ChlorineIncidentVelocityState
    wall_boundary: ChlorineWallRecombinationBoundary

    def __post_init__(self):
        values = np.asarray([
            self.incident_cl_atom_flux_m2_s,
            self.recombined_cl_atom_flux_m2_s,
            self.returned_cl2_molecule_flux_m2_s,
            self.mean_cl_speed_m_s,
        ], dtype=float)
        if (
            np.any(~np.isfinite(values))
            or np.any(values < 0.0)
            or self.mean_cl_speed_m_s <= 0.0
            or not str(self.source).strip()
            or self.evidence_kind not in WALL_RECOMBINATION_EVIDENCE_KINDS
            or not isinstance(
                self.incident_velocity_state,
                ChlorineIncidentVelocityState,
            )
            or not isinstance(
                self.wall_boundary,
                ChlorineWallRecombinationBoundary,
            )
        ):
            raise ValueError("invalid chlorine wall flux")
        for name, value in (
            ("incident_cl_atom_flux_m2_s", values[0]),
            ("recombined_cl_atom_flux_m2_s", values[1]),
            ("returned_cl2_molecule_flux_m2_s", values[2]),
            ("mean_cl_speed_m_s", values[3]),
        ):
            object.__setattr__(self, name, float(value))
        if self.recombined_cl_atom_flux_m2_s > (
                self.incident_cl_atom_flux_m2_s):
            raise ValueError("recombined flux cannot exceed incident flux")
        if not np.isclose(
            self.recombined_cl_atom_flux_m2_s,
            2.0 * self.returned_cl2_molecule_flux_m2_s,
            rtol=1.0e-14,
            atol=0.0,
        ):
            raise ValueError("chlorine wall flux does not conserve atoms")
        if not np.isclose(
            self.mean_cl_speed_m_s,
            self.incident_velocity_state.mean_speed_m_s,
            rtol=0.0,
            atol=0.0,
        ):
            raise ValueError(
                "wall flux does not use the supplied incident velocity")
        if (
            self.source != self.wall_boundary.source
            or self.evidence_kind != self.wall_boundary.evidence_kind
        ):
            raise ValueError("wall flux evidence does not match its boundary")
        if self.incident_cl_atom_flux_m2_s > 0.0 and not np.isclose(
            self.recombined_cl_atom_flux_m2_s
            / self.incident_cl_atom_flux_m2_s,
            self.wall_boundary.recombination_probability,
            rtol=1.0e-14,
            atol=0.0,
        ):
            raise ValueError("wall flux does not use its boundary probability")

    @property
    def supports_local_prediction(self) -> bool:
        """Whether wall response and incident velocity are both evidenced."""
        return (
            self.wall_boundary.supports_local_prediction
            and self.incident_velocity_state.supports_prediction
        )

    @property
    def chlorine_atom_inventory_residual_m2_s(self) -> float:
        """Return the exact elemental balance for the closed wall return."""
        return float(
            -self.recombined_cl_atom_flux_m2_s
            + 2.0 * self.returned_cl2_molecule_flux_m2_s)

    @property
    def reactor_volume_closure_ready(self) -> bool:
        """Local wall kinetics alone never closes chamber transport."""
        return False


def _bounded_domain(
    values: tuple[float, float],
    *,
    name: str,
    allow_zero: bool,
) -> tuple[float, float]:
    try:
        lower, upper = (float(value) for value in values)
    except (TypeError, ValueError):
        raise ValueError(f"{name} domain must contain two numbers") from None
    threshold_ok = lower >= 0.0 if allow_zero else lower > 0.0
    if (
        not np.isfinite(lower)
        or not np.isfinite(upper)
        or not threshold_ok
        or upper < lower
    ):
        raise ValueError(f"invalid {name} domain")
    return lower, upper
