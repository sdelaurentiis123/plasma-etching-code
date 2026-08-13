"""Deterministic collision-order transport through a reduced RF sheath.

The collisionless finite-transit sheath supplies one normal impact energy for
each RF entry phase.  This module lifts every phase into an effective Child
potential and evaluates the ion-neutral linear Boltzmann expansion with fixed
Gaussian quadrature.  It is deliberately neither a particle Monte Carlo nor a
two-Gaussian fit:

* collision positions are deterministic Gauss--Legendre nodes in optical
  depth;
* elastic and resonant charge-exchange channels are separate;
* equal-mass two-body kinematics conserves energy at every collision;
* every stopped collision order is reported as unresolved probability rather
  than silently assigned to a fitted tail; and
* an exact JVP of the discrete operator with respect to neutral density is
  propagated alongside the probability weights.

The current provider is the first collisional-sheath rung.  It resolves the
wafer ion energy-angle distribution and transports the no-further-collision
branch of every born fast neutral to the wafer.  Subsequent neutral-neutral
collisions remain an explicit unresolved ledger; the effective
phase-conditioned Child profile is not yet a self-consistent moving sheath.
Those omissions keep ``supports_feature_depth`` false.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.optimize import brentq
from scipy.special import k0

from ..sheath import CollisionlessWaveformSheath, PeriodicSheathVoltage
from .transport import phelps_argon_momentum_transfer_cross_section_m2
from .wafer_sheath_transfer import PowerClosedRFSheathProjection


BOLTZMANN_J_K = 1.380649e-23
E_CHARGE_C = 1.602176634e-19
BOHR_RADIUS_M = 5.29177210903e-11


def _readonly_1d(value, *, nonnegative: bool = False) -> np.ndarray:
    array = np.asarray(value, dtype=float).copy()
    if (
        array.ndim != 1
        or np.any(~np.isfinite(array))
        or (nonnegative and np.any(array < 0.0))
    ):
        raise ValueError("expected a finite one-dimensional array")
    array.setflags(write=False)
    return array


def _readonly_velocity(value) -> np.ndarray:
    array = np.asarray(value, dtype=float).copy()
    if (
        array.ndim != 2
        or array.shape[1] != 3
        or np.any(~np.isfinite(array))
        or np.any(array[:, 2] < 0.0)
    ):
        raise ValueError("invalid incident velocity quadrature")
    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class ArgonBornMayerPhelpsCollisionModel:
    """Equal-mass Ar+--Ar elastic/CX collision model at energetic support.

    The integrated cross section is Phelps' published Ar+--Ar momentum-transfer
    law evaluated in the center-of-mass frame.  Khrabrov and Kaganovich's
    Born--Mayer construction makes the ion-atom collision region symmetric:
    elastic and identity-switch charge exchange each carry probability 1/2.

    Their arXiv source does not print its fitted ``V0`` and ``a`` constants.
    The defaults below are therefore a declared two-point source inference:
    ``epsilon=E_com/V0 ~= 0.06`` at 1 keV laboratory energy fixes ``V0``, and
    their ``theta_lab=0.05 degree at b=5.4 a0`` statement, inserted into their
    equations (8)--(10), fixes ``a``.  The independent published
    ``R_cx ~= 8.09 a0`` statement is retained as a test, not used in the fit.

    Below 400 eV the source explicitly declines to claim Born--Mayer angular
    validity.  The model retains resonant identity switching there but applies
    zero resolved elastic angle, exposing the count in solution provenance.
    """

    born_mayer_potential_eV: float = 500.0 / 0.06
    born_mayer_range_m: float = 0.5116381650132525 * BOHR_RADIUS_M
    born_mayer_minimum_lab_energy_eV: float = 400.0
    born_mayer_maximum_lab_energy_eV: float = 10_000.0
    elastic_probability: float = 0.5
    charge_exchange_probability: float = 0.5
    source: str = (
        "Phelps 1994 DOI 10.1063/1.357820 momentum-transfer law; "
        "Khrabrov & Kaganovich arXiv:2604.04214v2 equations 8--10"
    )
    evidence_kind: str = "published_model_inference"
    provenance: Mapping[str, object] = field(default_factory=lambda: {
        "born_mayer_parameter_inference": (
            "V0 from E_com/V0~=0.06 at 1 keV lab; range from "
            "theta_lab=0.05 deg at b=5.4 a0"
        ),
        "independent_check_not_fit": "R_cx~=8.09 a0 at 1 keV lab",
        "channel_rule": (
            "equal-mass ion-atom collision; elastic and identity-switch "
            "charge exchange each probability 1/2"
        ),
        "low_energy_rule": (
            "below 400 eV: identity-switch energy channel retained, "
            "Born-Mayer angular deflection unresolved and set to zero"
        ),
        "coefficient_selection_target": None,
    })

    def __post_init__(self):
        values = np.asarray([
            self.born_mayer_potential_eV,
            self.born_mayer_range_m,
            self.born_mayer_minimum_lab_energy_eV,
            self.born_mayer_maximum_lab_energy_eV,
            self.elastic_probability,
            self.charge_exchange_probability,
        ], dtype=float)
        if (
            np.any(~np.isfinite(values))
            or np.any(values <= 0.0)
            or not math.isclose(
                self.elastic_probability + self.charge_exchange_probability,
                1.0,
                rel_tol=0.0,
                abs_tol=1.0e-14,
            )
            or self.born_mayer_minimum_lab_energy_eV
            >= self.born_mayer_maximum_lab_energy_eV
            or not str(self.source).strip()
            or self.evidence_kind != "published_model_inference"
        ):
            raise ValueError("invalid Ar Born-Mayer collision model")
        object.__setattr__(
            self, "provenance", MappingProxyType(dict(self.provenance)))

    def total_cross_section_m2(self, laboratory_energy_eV) -> np.ndarray | float:
        energy = np.asarray(laboratory_energy_eV, dtype=float)
        if np.any(~np.isfinite(energy)) or np.any(energy <= 0.0):
            raise ValueError("laboratory collision energy must be positive")
        result = phelps_argon_momentum_transfer_cross_section_m2(0.5 * energy)
        return float(result) if energy.ndim == 0 else np.asarray(result)

    def collision_radius_m(self, laboratory_energy_eV: float) -> float:
        return float(np.sqrt(
            self.total_cross_section_m2(laboratory_energy_eV) / np.pi))

    def center_of_mass_scattering_angle_rad(
        self,
        laboratory_energy_eV: float,
        impact_parameter_m: float,
    ) -> float:
        energy = float(laboratory_energy_eV)
        impact = float(impact_parameter_m)
        if (
            not math.isfinite(energy)
            or energy <= 0.0
            or not math.isfinite(impact)
            or impact < 0.0
        ):
            raise ValueError("invalid Born-Mayer collision state")
        if energy < self.born_mayer_minimum_lab_energy_eV:
            return 0.0
        if energy > self.born_mayer_maximum_lab_energy_eV:
            raise ValueError("collision energy exceeds Born-Mayer evidence support")
        epsilon = 0.5 * energy / self.born_mayer_potential_eV
        if not 0.0 < epsilon < 1.0:
            raise ValueError("Born-Mayer small-epsilon condition is violated")
        head_on = -math.log(epsilon)
        beta = impact / self.born_mayer_range_m
        if beta == 0.0:
            return math.pi

        def scaled_impact(rho: float) -> float:
            radicand = max(0.0, 1.0 - math.exp(-rho) / epsilon)
            return rho * math.sqrt(radicand)

        lower = np.nextafter(head_on, np.inf)
        upper = max(head_on + 2.0, beta + head_on + 2.0)
        while scaled_impact(upper) < beta:
            upper *= 2.0
        rho = brentq(
            lambda value: scaled_impact(value) - beta,
            lower,
            upper,
            xtol=2.0e-13,
            rtol=2.0e-13,
        )
        vartheta = (
            rho ** 3
            * float(k0(rho))
            / (epsilon * beta ** 2)
        )
        return float(2.0 * np.arcsin(1.0 / (1.0 + 2.0 / vartheta)))

    def _impact_quadrature(
        self,
        laboratory_energy_eV: float,
        order: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        node, weight = np.polynomial.legendre.leggauss(int(order))
        area_fraction = 0.5 * (node + 1.0)
        weight = 0.5 * weight
        radius = self.collision_radius_m(laboratory_energy_eV)
        angles = np.asarray([
            self.center_of_mass_scattering_angle_rad(
                laboratory_energy_eV,
                radius * math.sqrt(float(fraction)),
            )
            for fraction in area_fraction
        ])
        return angles, np.asarray(weight)

    def impact_quadrature(
        self,
        laboratory_energy_eV: float,
        order: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        if int(order) < 1:
            raise ValueError("impact quadrature order must be positive")
        return self._impact_quadrature(
            float(laboratory_energy_eV), int(order))


@dataclass(frozen=True)
class CollisionalIonEnergyAngleDistribution:
    """Resolved wafer-arriving ion measure, conditional on arrival."""

    velocity_sqrt_eV: np.ndarray
    weight: np.ndarray
    entry_phase_rad: np.ndarray

    def __post_init__(self):
        velocity = _readonly_velocity(self.velocity_sqrt_eV)
        weight = _readonly_1d(self.weight, nonnegative=True)
        phase = _readonly_1d(self.entry_phase_rad)
        if (
            velocity.shape[0] == 0
            or weight.shape != (velocity.shape[0],)
            or phase.shape != weight.shape
            or not math.isclose(
                float(np.sum(weight)), 1.0, rel_tol=0.0, abs_tol=2.0e-12)
        ):
            raise ValueError("invalid collisional ion distribution")
        object.__setattr__(self, "velocity_sqrt_eV", velocity)
        object.__setattr__(self, "weight", weight)
        object.__setattr__(self, "entry_phase_rad", phase)

    @property
    def energy_eV(self) -> np.ndarray:
        return np.einsum(
            "ij,ij->i", self.velocity_sqrt_eV, self.velocity_sqrt_eV)

    @property
    def polar_angle_deg(self) -> np.ndarray:
        transverse = np.linalg.norm(self.velocity_sqrt_eV[:, :2], axis=1)
        return np.rad2deg(np.arctan2(
            transverse, self.velocity_sqrt_eV[:, 2]))

    @property
    def mean_energy_eV(self) -> float:
        return float(np.dot(self.weight, self.energy_eV))

    @property
    def mean_squared_polar_angle_rad2(self) -> float:
        angle = np.deg2rad(self.polar_angle_deg)
        return float(np.dot(self.weight, angle * angle))


@dataclass(frozen=True)
class DeterministicCollisionalSheathSolution:
    distribution: CollisionalIonEnergyAngleDistribution
    resolved_fast_neutral_distribution: (
        CollisionalIonEnergyAngleDistribution | None)
    source_ion_flux_m2_s: float
    arriving_ion_flux_m2_s: float
    resolved_fast_neutral_flux_m2_s: float
    ion_arrival_probability: float
    unresolved_probability: float
    escaped_probability: float
    uncollided_arrival_probability: float
    expected_collision_count_lower_bound: float
    expected_charge_exchange_count_lower_bound: float
    expected_fast_neutral_birth_count_lower_bound: float
    expected_fast_neutral_birth_energy_lower_bound_eV_per_source_ion: float
    resolved_fast_neutral_arrivals_per_source_ion: float
    unresolved_fast_neutral_collisions_per_source_ion: float
    escaped_fast_neutrals_per_source_ion: float
    fast_neutral_lineage_ledger_relative_residual: float
    maximum_resolved_energy_ledger_relative_residual: float
    probability_ledger_relative_residual: float
    collisionless_reference_mean_normal_energy_eV: float
    mean_total_optical_depth: float
    maximum_total_optical_depth: float
    below_born_mayer_support_collision_probability_lower_bound: float
    model_source: str
    provenance: Mapping[str, object]
    supports_density_jvp: bool = True
    supports_fast_neutral_wafer_flux: bool = False
    supports_feature_depth: bool = False

    def __post_init__(self):
        scalars = np.asarray([
            self.source_ion_flux_m2_s,
            self.arriving_ion_flux_m2_s,
            self.resolved_fast_neutral_flux_m2_s,
            self.ion_arrival_probability,
            self.unresolved_probability,
            self.escaped_probability,
            self.uncollided_arrival_probability,
            self.expected_collision_count_lower_bound,
            self.expected_charge_exchange_count_lower_bound,
            self.expected_fast_neutral_birth_count_lower_bound,
            self.expected_fast_neutral_birth_energy_lower_bound_eV_per_source_ion,
            self.resolved_fast_neutral_arrivals_per_source_ion,
            self.unresolved_fast_neutral_collisions_per_source_ion,
            self.escaped_fast_neutrals_per_source_ion,
            self.fast_neutral_lineage_ledger_relative_residual,
            self.maximum_resolved_energy_ledger_relative_residual,
            self.probability_ledger_relative_residual,
            self.collisionless_reference_mean_normal_energy_eV,
            self.mean_total_optical_depth,
            self.maximum_total_optical_depth,
            self.below_born_mayer_support_collision_probability_lower_bound,
        ], dtype=float)
        if (
            not isinstance(
                self.distribution, CollisionalIonEnergyAngleDistribution)
            or np.any(~np.isfinite(scalars))
            or np.any(scalars < 0.0)
            or self.source_ion_flux_m2_s <= 0.0
            or (
                self.resolved_fast_neutral_distribution is None
                and self.resolved_fast_neutral_flux_m2_s != 0.0
            )
            or (
                self.resolved_fast_neutral_distribution is not None
                and not isinstance(
                    self.resolved_fast_neutral_distribution,
                    CollisionalIonEnergyAngleDistribution,
                )
            )
            or not math.isclose(
                self.resolved_fast_neutral_flux_m2_s,
                self.source_ion_flux_m2_s
                * self.resolved_fast_neutral_arrivals_per_source_ion,
                rel_tol=2.0e-12,
                abs_tol=0.0,
            )
            or not math.isclose(
                self.arriving_ion_flux_m2_s,
                self.source_ion_flux_m2_s * self.ion_arrival_probability,
                rel_tol=2.0e-12,
                abs_tol=0.0,
            )
            or self.probability_ledger_relative_residual > 2.0e-11
            or self.maximum_resolved_energy_ledger_relative_residual > 2.0e-10
            or self.fast_neutral_lineage_ledger_relative_residual > 2.0e-11
            or not str(self.model_source).strip()
            or not self.supports_density_jvp
            or self.supports_fast_neutral_wafer_flux
            or self.supports_feature_depth
        ):
            raise ValueError("collisional sheath certification failed")
        object.__setattr__(
            self, "provenance", MappingProxyType(dict(self.provenance)))

    def to_boundary_state(
        self,
        *,
        ion_name: str,
        ion_mass_amu: float,
        reference_plane_m: float = 0.0,
    ):
        """Return the common forward boundary contract for resolved ions."""
        from ..boundary_state import PlasmaBoundaryState, SpeciesBoundaryState

        ion = SpeciesBoundaryState(
            name=str(ion_name),
            charge_number=1,
            mass_amu=float(ion_mass_amu),
            flux_m2_s=self.arriving_ion_flux_m2_s,
            velocity_sqrt_eV=self.distribution.velocity_sqrt_eV,
            weight=self.distribution.weight,
            phase_rad=self.distribution.entry_phase_rad,
            provenance={
                "model": "deterministic_collision_order_rf_sheath",
                "unresolved_probability": self.unresolved_probability,
                "fast_neutral_wafer_flux_closed": False,
            },
        )
        species = [ion]
        if self.resolved_fast_neutral_distribution is not None:
            neutral = self.resolved_fast_neutral_distribution
            species.append(SpeciesBoundaryState(
                name=f"{str(ion_name).removesuffix('+')}_fast_neutral",
                charge_number=0,
                mass_amu=float(ion_mass_amu),
                flux_m2_s=self.resolved_fast_neutral_flux_m2_s,
                velocity_sqrt_eV=neutral.velocity_sqrt_eV,
                weight=neutral.weight,
                phase_rad=neutral.entry_phase_rad,
                provenance={
                    "model": "unscattered_fast_neutral_lower_bound",
                    "unresolved_neutral_collisions_per_source_ion": (
                        self.unresolved_fast_neutral_collisions_per_source_ion),
                },
            ))
        return PlasmaBoundaryState(
            species=tuple(species),
            reference_plane_m=float(reference_plane_m),
            provenance={
                "source": self.model_source,
                "supports_feature_depth": False,
                "fast_neutral_boundary_is_lower_bound": True,
            },
        )


@dataclass(frozen=True)
class CollisionalSheathDensityTangent:
    gas_number_density_tangent_m3: float
    distribution_weight_tangent: np.ndarray
    ion_arrival_probability_tangent: float
    unresolved_probability_tangent: float
    escaped_probability_tangent: float
    uncollided_arrival_probability_tangent: float
    expected_collision_count_tangent: float
    expected_charge_exchange_count_tangent: float
    expected_fast_neutral_birth_count_tangent: float
    expected_fast_neutral_birth_energy_tangent_eV: float
    resolved_fast_neutral_arrivals_tangent: float
    unresolved_fast_neutral_collisions_tangent: float
    escaped_fast_neutrals_tangent: float
    fast_neutral_lineage_ledger_tangent_residual: float
    mean_impact_energy_tangent_eV: float
    probability_ledger_tangent_residual: float

    def __post_init__(self):
        weight = _readonly_1d(self.distribution_weight_tangent)
        scalars = np.asarray([
            self.gas_number_density_tangent_m3,
            self.ion_arrival_probability_tangent,
            self.unresolved_probability_tangent,
            self.escaped_probability_tangent,
            self.uncollided_arrival_probability_tangent,
            self.expected_collision_count_tangent,
            self.expected_charge_exchange_count_tangent,
            self.expected_fast_neutral_birth_count_tangent,
            self.expected_fast_neutral_birth_energy_tangent_eV,
            self.resolved_fast_neutral_arrivals_tangent,
            self.unresolved_fast_neutral_collisions_tangent,
            self.escaped_fast_neutrals_tangent,
            self.fast_neutral_lineage_ledger_tangent_residual,
            self.mean_impact_energy_tangent_eV,
            self.probability_ledger_tangent_residual,
        ], dtype=float)
        if (
            np.any(~np.isfinite(scalars))
            or abs(self.probability_ledger_tangent_residual) > 2.0e-10
            or abs(
                self.fast_neutral_lineage_ledger_tangent_residual
            ) > 2.0e-10
        ):
            raise ValueError("invalid collisional-sheath density tangent")
        object.__setattr__(self, "distribution_weight_tangent", weight)


@dataclass(frozen=True)
class PowerClosedArgonCollisionalSheathProjection:
    """Opt-in Ar wafer projection preserving the upstream power ledger.

    Ion-neutral collisions only redistribute the energy gained from the
    electric field between the ion and fast-neutral lineages.  They do not
    change the upstream ``A e Gamma_i <Delta E>`` closure.  This contract
    therefore carries that closure verbatim while exposing the collisional
    truncation separately.
    """

    collisionless: PowerClosedRFSheathProjection
    collisional: DeterministicCollisionalSheathSolution
    pressure_Pa: float
    gas_temperature_K: float
    neutral_density_m3: float
    collisionless_reference_relative_residual: float
    provenance: Mapping[str, object]
    supports_feature_depth: bool = False

    def __post_init__(self):
        scalars = np.asarray([
            self.pressure_Pa,
            self.gas_temperature_K,
            self.neutral_density_m3,
            self.collisionless_reference_relative_residual,
        ], dtype=float)
        if (
            not isinstance(self.collisionless, PowerClosedRFSheathProjection)
            or not isinstance(
                self.collisional, DeterministicCollisionalSheathSolution)
            or set(self.collisionless.distributions) != {"Ar+"}
            or np.any(~np.isfinite(scalars))
            or np.any(scalars[:3] <= 0.0)
            or abs(self.collisionless_reference_relative_residual) > 5.0e-4
            or self.supports_feature_depth
        ):
            raise ValueError("invalid power-closed Ar collisional projection")
        object.__setattr__(
            self, "provenance", MappingProxyType(dict(self.provenance)))

    def to_boundary_state(self, *, reference_plane_m: float = 0.0):
        return self.collisional.to_boundary_state(
            ion_name="Ar+",
            ion_mass_amu=self.collisionless.distributions["Ar+"].ion_mass_amu,
            reference_plane_m=reference_plane_m,
        )


@dataclass(frozen=True)
class DeterministicArgonCollisionalSheathTransfer:
    """Lift a power-closed pure-Ar sheath projection through collisions.

    The deliberately narrow ``Ar+``/Ar contract prevents the common but
    unjustified step of applying an argon charge-exchange law to molecular
    ions in fluorocarbon or chlorine plasmas.
    """

    collision_model: ArgonBornMayerPhelpsCollisionModel = field(
        default_factory=ArgonBornMayerPhelpsCollisionModel)
    initial_thermal_radial_order: int = 2
    initial_thermal_azimuth_order: int = 4
    position_quadrature_order: int = 5
    hazard_quadrature_order: int = 6
    impact_quadrature_order: int = 3
    collision_azimuth_order: int = 4
    maximum_collision_order: int = 2
    solver_kind: str = "implicit_discrete_ordinates"
    potential_node_count: int = 9
    total_energy_node_count: int = 9
    transverse_fraction_node_count: int = 13
    steps_per_period: int = 256
    steps_per_transit: int = 256

    def __post_init__(self):
        if not isinstance(
            self.collision_model, ArgonBornMayerPhelpsCollisionModel
        ):
            raise TypeError("an Ar+--Ar collision model is required")
        orders = (
            self.initial_thermal_radial_order,
            self.initial_thermal_azimuth_order,
            self.position_quadrature_order,
            self.hazard_quadrature_order,
            self.impact_quadrature_order,
            self.collision_azimuth_order,
            self.maximum_collision_order,
            self.potential_node_count,
            self.total_energy_node_count,
            self.transverse_fraction_node_count,
            self.steps_per_period,
            self.steps_per_transit,
        )
        if (
            any(int(value) < 1 for value in orders)
            or self.solver_kind not in {
                "implicit_discrete_ordinates", "collision_order_reference"
            }
        ):
            raise ValueError("collisional transfer orders must be positive")

    def project(
        self,
        collisionless: PowerClosedRFSheathProjection,
        *,
        pressure_Pa: float,
        gas_temperature_K: float,
    ) -> PowerClosedArgonCollisionalSheathProjection:
        if (
            not isinstance(collisionless, PowerClosedRFSheathProjection)
            or set(collisionless.distributions) != {"Ar+"}
        ):
            raise ValueError(
                "the collisional Ar transfer requires one Ar+ projection")
        pressure = float(pressure_Pa)
        gas_temperature = float(gas_temperature_K)
        if (
            not math.isfinite(pressure)
            or pressure <= 0.0
            or not math.isfinite(gas_temperature)
            or gas_temperature <= 0.0
        ):
            raise ValueError("positive pressure and gas temperature required")
        source = collisionless.distributions["Ar+"]
        waveform = PeriodicSheathVoltage.sinusoidal(
            dc_v=collisionless.sheath_dc_v,
            amplitude_v=collisionless.sheath_rf_amplitude_v,
            frequency_hz=collisionless.frequency_hz,
            source=(
                f"{collisionless.source}; deterministic Ar collision lift"),
            evidence_kind="assumed",
        )
        sheath = CollisionlessWaveformSheath(
            waveform=waveform,
            Te_eV=collisionless.electron_temperature_eV,
            ion_mass_amu=source.ion_mass_amu,
            thickness_m=(
                collisionless.sheath_thickness_m_by_species["Ar+"]),
        )
        density = pressure / (BOLTZMANN_J_K * gas_temperature)
        common = {
            "sheath": sheath,
            "collision_model": self.collision_model,
            "gas_number_density_m3": density,
            "neutral_gas_temperature_K": gas_temperature,
            "source_ion_flux_m2_s": source.flux_m2_s,
            "phase_count": source.energy_eV.size,
            "initial_thermal_radial_order": int(
                self.initial_thermal_radial_order),
            "initial_thermal_azimuth_order": int(
                self.initial_thermal_azimuth_order),
            "position_quadrature_order": int(self.position_quadrature_order),
            "hazard_quadrature_order": int(self.hazard_quadrature_order),
            "impact_quadrature_order": int(self.impact_quadrature_order),
            "collision_azimuth_order": int(self.collision_azimuth_order),
            "steps_per_period": int(self.steps_per_period),
            "steps_per_transit": int(self.steps_per_transit),
        }
        if self.solver_kind == "implicit_discrete_ordinates":
            # Local import avoids a module cycle: the bounded operator reuses
            # the source-audited collision law and solution contract here.
            from .collisional_sheath_discrete_ordinates import (
                DeterministicDiscreteOrdinatesRFSheath,
            )
            collisional = DeterministicDiscreteOrdinatesRFSheath(
                **common,
                potential_node_count=int(self.potential_node_count),
                total_energy_node_count=int(
                    self.total_energy_node_count),
                transverse_fraction_node_count=int(
                    self.transverse_fraction_node_count),
            ).solve()
        else:
            collisional = DeterministicCollisionalRFSheath(
                **common,
                maximum_collision_order=int(self.maximum_collision_order),
            ).solve()
        reference = collisional.collisionless_reference_mean_normal_energy_eV
        residual = (
            reference - source.mean_energy_eV
        ) / max(abs(reference), abs(source.mean_energy_eV), 1.0)
        return PowerClosedArgonCollisionalSheathProjection(
            collisionless=collisionless,
            collisional=collisional,
            pressure_Pa=pressure,
            gas_temperature_K=gas_temperature,
            neutral_density_m3=density,
            collisionless_reference_relative_residual=residual,
            provenance={
                "operator": (
                    "deterministic_Ar_implicit_discrete_ordinates_lift"
                    if self.solver_kind == "implicit_discrete_ordinates"
                    else "deterministic_Ar_collision_order_reference_lift"
                ),
                "solver_kind": self.solver_kind,
                "upstream_power_closure_relative_residual": (
                    collisionless.power_closure_relative_residual),
                "upstream_delivered_bias_power_W": (
                    collisionless.delivered_bias_power_W),
                "feature_depth_used": False,
                "fluorocarbon_transfer_authorized": False,
            },
            supports_feature_depth=False,
        )


def _orthogonal_frame(direction: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    direction = np.asarray(direction, dtype=float)
    if abs(direction[2]) < 0.9:
        reference = np.array([0.0, 0.0, 1.0])
    else:
        reference = np.array([1.0, 0.0, 0.0])
    first = np.cross(direction, reference)
    first /= np.linalg.norm(first)
    second = np.cross(direction, first)
    return first, second


def _equal_mass_collision_velocities(
    incoming_velocity_sqrt_eV: np.ndarray,
    center_of_mass_angle_rad: float,
    azimuth_rad: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return projectile and initially stationary target after elastic impact."""
    incoming = np.asarray(incoming_velocity_sqrt_eV, dtype=float)
    speed = float(np.linalg.norm(incoming))
    if speed <= 0.0:
        raise ValueError("collision requires a moving projectile")
    direction = incoming / speed
    first, second = _orthogonal_frame(direction)
    transverse = (
        math.cos(azimuth_rad) * first + math.sin(azimuth_rad) * second)
    half = 0.5 * float(center_of_mass_angle_rad)
    projectile = speed * math.cos(half) * (
        math.cos(half) * direction + math.sin(half) * transverse)
    target = speed * math.sin(half) * (
        math.sin(half) * direction - math.cos(half) * transverse)
    return projectile, target


def _canonical_axisymmetric_velocity(velocity: np.ndarray) -> np.ndarray:
    value = np.asarray(velocity, dtype=float)
    return np.array([
        float(np.linalg.norm(value[:2])),
        0.0,
        float(value[2]),
    ])


@dataclass(frozen=True)
class DeterministicCollisionalRFSheath:
    """Collision-order quadrature around a finite-transit RF sheath."""

    sheath: CollisionlessWaveformSheath
    collision_model: ArgonBornMayerPhelpsCollisionModel
    gas_number_density_m3: float
    neutral_gas_temperature_K: float
    source_ion_flux_m2_s: float
    phase_count: int = 24
    initial_thermal_radial_order: int = 2
    initial_thermal_azimuth_order: int = 4
    position_quadrature_order: int = 5
    hazard_quadrature_order: int = 6
    impact_quadrature_order: int = 3
    collision_azimuth_order: int = 4
    maximum_collision_order: int = 2
    steps_per_period: int = 256
    steps_per_transit: int = 256

    def __post_init__(self):
        if not isinstance(self.sheath, CollisionlessWaveformSheath):
            raise TypeError("a CollisionlessWaveformSheath is required")
        if not isinstance(
            self.collision_model, ArgonBornMayerPhelpsCollisionModel
        ):
            raise TypeError("an Ar Born-Mayer collision model is required")
        values = np.asarray([
            self.gas_number_density_m3,
            self.neutral_gas_temperature_K,
            self.source_ion_flux_m2_s,
        ], dtype=float)
        orders = (
            self.phase_count,
            self.initial_thermal_radial_order,
            self.initial_thermal_azimuth_order,
            self.position_quadrature_order,
            self.hazard_quadrature_order,
            self.impact_quadrature_order,
            self.collision_azimuth_order,
            self.maximum_collision_order,
            self.steps_per_period,
            self.steps_per_transit,
        )
        if (
            np.any(~np.isfinite(values))
            or self.gas_number_density_m3 < 0.0
            or self.neutral_gas_temperature_K <= 0.0
            or self.source_ion_flux_m2_s <= 0.0
            or any(int(value) < 1 for value in orders)
            or int(self.phase_count) < 4
            or int(self.initial_thermal_azimuth_order) < 2
            or int(self.collision_azimuth_order) % 2 != 0
        ):
            raise ValueError("invalid deterministic collisional sheath")

    @property
    def tangential_temperature_eV(self) -> float:
        return (
            BOLTZMANN_J_K * self.neutral_gas_temperature_K / E_CHARGE_C)

    @staticmethod
    def _potential_shape(position_fraction: float) -> float:
        return float(max(position_fraction, 0.0) ** (4.0 / 3.0))

    def _propagate_velocity(
        self,
        velocity: np.ndarray,
        start_fraction: float,
        end_fraction: float,
        phase_gain_eV: float,
    ) -> np.ndarray | None:
        result = np.asarray(velocity, dtype=float).copy()
        normal_energy = (
            result[2] ** 2
            + phase_gain_eV
            * (
                self._potential_shape(end_fraction)
                - self._potential_shape(start_fraction)
            )
        )
        if normal_energy < -1.0e-11:
            return None
        result[2] = math.sqrt(max(normal_energy, 0.0))
        return result

    def _hazard_rate_per_fraction(
        self,
        velocity: np.ndarray,
        density: float,
    ) -> float:
        normal = float(velocity[2])
        if normal <= 0.0:
            return math.inf
        speed = float(np.linalg.norm(velocity))
        energy = speed * speed
        return float(
            density
            * self.sheath.thickness
            * self.collision_model.total_cross_section_m2(energy)
            * speed
            / normal
        )

    def _integrated_hazard(
        self,
        velocity: np.ndarray,
        start_fraction: float,
        end_fraction: float,
        phase_gain_eV: float,
        density: float,
    ) -> float:
        if end_fraction <= start_fraction or density == 0.0:
            return 0.0
        node, weight = np.polynomial.legendre.leggauss(
            int(self.hazard_quadrature_order))
        half = 0.5 * (end_fraction - start_fraction)
        center = 0.5 * (end_fraction + start_fraction)
        result = 0.0
        for local, quadrature_weight in zip(node, weight):
            position = center + half * float(local)
            propagated = self._propagate_velocity(
                velocity, start_fraction, position, phase_gain_eV)
            if propagated is None:
                return math.inf
            result += float(quadrature_weight) * (
                self._hazard_rate_per_fraction(propagated, density))
        return float(half * result)

    def _initial_velocity_quadrature(self):
        temperature = self.tangential_temperature_eV
        radial_node, radial_weight = np.polynomial.laguerre.laggauss(
            int(self.initial_thermal_radial_order))
        bohm_energy = 0.5 * float(self.sheath.Te_eV)
        for node, weight in zip(radial_node, radial_weight):
            transverse = math.sqrt(temperature * float(node))
            # The internal collision operator is axisymmetric.  A single
            # canonical transverse direction avoids repeating rotated copies
            # of the same tree; wafer arrivals are lifted back to a uniform
            # 3-D azimuthal quadrature below.
            yield (
                np.array([transverse, 0.0, math.sqrt(bohm_energy)]),
                float(weight),
            )

    def _solve(
        self,
        *,
        gas_number_density_tangent_m3: float | None,
    ) -> tuple[
        DeterministicCollisionalSheathSolution,
        CollisionalSheathDensityTangent | None,
    ]:
        density = float(self.gas_number_density_m3)
        tangent_enabled = gas_number_density_tangent_m3 is not None
        density_tangent = (
            0.0 if gas_number_density_tangent_m3 is None
            else float(gas_number_density_tangent_m3)
        )
        if not math.isfinite(density_tangent):
            raise ValueError("neutral-density tangent must be finite")
        phases = (
            2.0 * np.pi
            * (np.arange(int(self.phase_count), dtype=float) + 0.5)
            / int(self.phase_count)
        )
        collisionless_energy = self.sheath.ion_impact_energies(
            phases,
            steps_per_period=int(self.steps_per_period),
            steps_per_transit=int(self.steps_per_transit),
        )
        bohm_energy = 0.5 * float(self.sheath.Te_eV)
        gain = np.maximum(collisionless_energy - bohm_energy, 0.0)
        position_node, position_weight = np.polynomial.legendre.leggauss(
            int(self.position_quadrature_order))

        arriving_velocity: list[np.ndarray] = []
        arriving_weight: list[float] = []
        arriving_weight_tangent: list[float] = []
        arriving_phase: list[float] = []
        arrived = 0.0
        arrived_tangent = 0.0
        unresolved = 0.0
        unresolved_tangent = 0.0
        escaped = 0.0
        escaped_tangent = 0.0
        uncollided = 0.0
        uncollided_tangent = 0.0
        collision_count = 0.0
        collision_count_tangent = 0.0
        charge_exchange_count = 0.0
        charge_exchange_count_tangent = 0.0
        neutral_birth_count = 0.0
        neutral_birth_count_tangent = 0.0
        neutral_birth_energy = 0.0
        neutral_birth_energy_tangent = 0.0
        neutral_arrival = 0.0
        neutral_arrival_tangent = 0.0
        neutral_unresolved = 0.0
        neutral_unresolved_tangent = 0.0
        neutral_escaped = 0.0
        neutral_escaped_tangent = 0.0
        neutral_velocity_arriving: list[np.ndarray] = []
        neutral_weight_arriving: list[float] = []
        neutral_phase_arriving: list[float] = []
        below_support = 0.0
        below_support_tangent = 0.0
        maximum_energy_residual = 0.0
        optical_depth: list[float] = []
        impact_cache: dict[float, tuple[np.ndarray, np.ndarray]] = {}

        output_azimuth = (
            2.0 * np.pi
            * (np.arange(int(self.initial_thermal_azimuth_order)) + 0.5)
            / int(self.initial_thermal_azimuth_order)
        )

        def append_axisymmetric_arrival(
            velocity: np.ndarray,
            weight: float,
            tangent: float,
            phase: float,
            *,
            neutral: bool,
        ) -> None:
            transverse = float(np.linalg.norm(velocity[:2]))
            normal = float(velocity[2])
            for angle in output_azimuth:
                lifted = np.array([
                    transverse * math.cos(float(angle)),
                    transverse * math.sin(float(angle)),
                    normal,
                ])
                if neutral:
                    neutral_velocity_arriving.append(lifted)
                    neutral_weight_arriving.append(
                        weight / len(output_azimuth))
                    neutral_phase_arriving.append(phase)
                else:
                    arriving_velocity.append(lifted)
                    arriving_weight.append(weight / len(output_azimuth))
                    arriving_weight_tangent.append(
                        tangent / len(output_azimuth))
                    arriving_phase.append(phase)

        def transport_born_fast_neutral(
            velocity: np.ndarray,
            position_fraction: float,
            phase: float,
            weight: float,
            weight_tangent: float,
        ) -> None:
            """Transport the no-further-collision neutral branch to wafer.

            Khrabrov--Kaganovich state that the Ar--Ar differential law is
            the ion--Ar law without the charge-exchange channel and that its
            total cross section is twice the CX cross section--equal to the
            full ion-event cross section used here.  Subsequent neutral
            collisions are therefore exposed as unresolved mass at this rung.
            """
            nonlocal neutral_arrival, neutral_arrival_tangent
            nonlocal neutral_unresolved, neutral_unresolved_tangent
            nonlocal neutral_escaped, neutral_escaped_tangent
            speed = float(np.linalg.norm(velocity))
            normal = float(velocity[2])
            if speed <= 0.0 or normal <= 0.0:
                neutral_escaped += weight
                neutral_escaped_tangent += weight_tangent
                return
            path_factor = (
                self.sheath.thickness
                * (1.0 - position_fraction)
                * speed
                / normal
            )
            unit_hazard = (
                path_factor
                * self.collision_model.total_cross_section_m2(speed * speed)
            )
            hazard = density * unit_hazard
            hazard_tangent = density_tangent * unit_hazard
            survival = math.exp(-hazard)
            survival_tangent = -survival * hazard_tangent
            arrival_weight = weight * survival
            arrival_tangent = (
                weight_tangent * survival + weight * survival_tangent)
            collision_weight = weight * (1.0 - survival)
            collision_tangent = (
                weight_tangent * (1.0 - survival)
                - weight * survival_tangent)
            neutral_arrival += arrival_weight
            neutral_arrival_tangent += arrival_tangent
            neutral_unresolved += collision_weight
            neutral_unresolved_tangent += collision_tangent
            if arrival_weight > np.finfo(float).tiny:
                append_axisymmetric_arrival(
                    velocity,
                    arrival_weight,
                    arrival_tangent,
                    phase,
                    neutral=True,
                )

        def add_probability(
            name: str,
            value: float,
            tangent: float,
        ) -> None:
            nonlocal arrived, arrived_tangent, unresolved
            nonlocal unresolved_tangent, escaped, escaped_tangent
            if name == "arrived":
                arrived += value
                arrived_tangent += tangent
            elif name == "unresolved":
                unresolved += value
                unresolved_tangent += tangent
            elif name == "escaped":
                escaped += value
                escaped_tangent += tangent
            else:  # pragma: no cover - internal closed set
                raise RuntimeError("unknown probability ledger")

        def recurse(
            velocity: np.ndarray,
            start_fraction: float,
            phase_gain: float,
            phase: float,
            weight: float,
            weight_tangent: float,
            collision_order: int,
            neutral_energy_accumulated: float,
            total_energy_target: float,
        ) -> None:
            nonlocal uncollided, uncollided_tangent
            nonlocal collision_count, collision_count_tangent
            nonlocal charge_exchange_count, charge_exchange_count_tangent
            nonlocal neutral_birth_energy, neutral_birth_energy_tangent
            nonlocal below_support, below_support_tangent
            nonlocal maximum_energy_residual

            unit_hazard_end = self._integrated_hazard(
                velocity,
                start_fraction,
                1.0,
                phase_gain,
                1.0,
            )
            hazard_end = density * unit_hazard_end
            if not math.isfinite(hazard_end):
                add_probability("escaped", weight, weight_tangent)
                return
            hazard_tangent = density_tangent * unit_hazard_end
            survival = math.exp(-hazard_end)
            survival_tangent = -survival * hazard_tangent
            resolved_weight = weight * survival
            resolved_tangent = (
                weight_tangent * survival + weight * survival_tangent)
            final_velocity = self._propagate_velocity(
                velocity, start_fraction, 1.0, phase_gain)
            if final_velocity is None or final_velocity[2] < 0.0:
                add_probability(
                    "escaped", resolved_weight, resolved_tangent)
            else:
                add_probability(
                    "arrived", resolved_weight, resolved_tangent)
                append_axisymmetric_arrival(
                    final_velocity,
                    resolved_weight,
                    resolved_tangent,
                    phase,
                    neutral=False,
                )
                if collision_order == 0:
                    uncollided += resolved_weight
                    uncollided_tangent += resolved_tangent
                accounted = (
                    float(np.dot(final_velocity, final_velocity))
                    + neutral_energy_accumulated
                )
                maximum_energy_residual = max(
                    maximum_energy_residual,
                    abs(accounted - total_energy_target)
                    / max(abs(total_energy_target), 1.0),
                )

            collision_probability = -math.expm1(-hazard_end)
            collision_probability_tangent = -survival_tangent
            colliding_weight = weight * collision_probability
            colliding_tangent = (
                weight_tangent * collision_probability
                + weight * collision_probability_tangent
            )
            if (
                collision_probability <= np.finfo(float).tiny
                and abs(colliding_tangent) <= np.finfo(float).tiny
            ):
                return
            if collision_order >= int(self.maximum_collision_order):
                add_probability(
                    "unresolved", colliding_weight, colliding_tangent)
                # This is a known additional collision even though its output
                # state lies beyond the retained collision-order expansion.
                collision_count += colliding_weight
                collision_count_tangent += colliding_tangent
                return

            half = 0.5 * (1.0 - start_fraction)
            center = 0.5 * (1.0 + start_fraction)
            event_records = []
            raw_total = 0.0
            raw_total_tangent = 0.0
            for local, quadrature_weight in zip(
                position_node, position_weight
            ):
                position = center + half * float(local)
                event_velocity = self._propagate_velocity(
                    velocity, start_fraction, position, phase_gain)
                if event_velocity is None:
                    continue
                unit_hazard = self._integrated_hazard(
                    velocity,
                    start_fraction,
                    position,
                    phase_gain,
                    1.0,
                )
                hazard = density * unit_hazard
                hazard_direction = density_tangent * unit_hazard
                unit_rate = self._hazard_rate_per_fraction(
                    event_velocity, 1.0)
                # The common density factor cancels from the first-event
                # conditional measure.  Removing it analytically keeps the
                # JVP well-defined at the collisionless boundary n_g=0.
                raw = (
                    half
                    * float(quadrature_weight)
                    * unit_rate
                    * math.exp(-hazard)
                )
                raw_tangent = -raw * hazard_direction
                event_records.append((
                    position,
                    event_velocity,
                    raw,
                    raw_tangent,
                ))
                raw_total += raw
                raw_total_tangent += raw_tangent
            if raw_total <= 0.0:
                # Positive integrated hazard with a zero event quadrature is a
                # numerical failure, never a collisionless reinterpretation.
                raise RuntimeError("first-collision quadrature lost its mass")

            for position, event_velocity, raw, raw_tangent in event_records:
                conditional = raw / raw_total
                conditional_tangent = (
                    raw_tangent * raw_total
                    - raw * raw_total_tangent
                ) / raw_total ** 2
                event_weight = colliding_weight * conditional
                event_tangent = (
                    colliding_tangent * conditional
                    + colliding_weight * conditional_tangent
                )
                event_energy = float(np.dot(event_velocity, event_velocity))
                cache_key = round(event_energy, 8)
                if cache_key not in impact_cache:
                    impact_cache[cache_key] = (
                        self.collision_model.impact_quadrature(
                            event_energy, int(self.impact_quadrature_order)))
                angles, impact_weight = impact_cache[cache_key]
                # The incoming global azimuth is redundant, but collision
                # azimuths are not: the wafer normal distinguishes scattering
                # toward from scattering away from the electrode.  Retain the
                # full azimuth quadrature before canonicalizing each output.
                azimuths = (
                    2.0 * np.pi
                    * (np.arange(int(self.collision_azimuth_order)) + 0.5)
                    / int(self.collision_azimuth_order)
                )

                def continue_branch(
                    ion_velocity: np.ndarray,
                    neutral_velocity: np.ndarray,
                    branch_fraction: float,
                    *,
                    charge_exchange: bool,
                ) -> None:
                    nonlocal collision_count, collision_count_tangent
                    nonlocal charge_exchange_count
                    nonlocal charge_exchange_count_tangent
                    nonlocal neutral_birth_count
                    nonlocal neutral_birth_count_tangent
                    nonlocal neutral_birth_energy
                    nonlocal neutral_birth_energy_tangent
                    nonlocal below_support, below_support_tangent
                    branch_weight = event_weight * branch_fraction
                    branch_tangent = event_tangent * branch_fraction
                    collision_count += branch_weight
                    collision_count_tangent += branch_tangent
                    neutral_birth_count += branch_weight
                    neutral_birth_count_tangent += branch_tangent
                    if charge_exchange:
                        charge_exchange_count += branch_weight
                        charge_exchange_count_tangent += branch_tangent
                    born_energy = float(np.dot(
                        neutral_velocity, neutral_velocity))
                    neutral_birth_energy += branch_weight * born_energy
                    neutral_birth_energy_tangent += (
                        branch_tangent * born_energy)
                    transport_born_fast_neutral(
                        neutral_velocity,
                        position,
                        phase,
                        branch_weight,
                        branch_tangent,
                    )
                    if (
                        event_energy
                        < self.collision_model
                        .born_mayer_minimum_lab_energy_eV
                    ):
                        below_support += branch_weight
                        below_support_tangent += branch_tangent
                    if ion_velocity[2] < 0.0:
                        add_probability(
                            "escaped", branch_weight, branch_tangent)
                        return
                    recurse(
                        ion_velocity,
                        position,
                        phase_gain,
                        phase,
                        branch_weight,
                        branch_tangent,
                        collision_order + 1,
                        neutral_energy_accumulated + born_energy,
                        total_energy_target,
                    )

                # Elastic scattering retains the Born--Mayer differential
                # angle.  The projectile remains the ion and the target atom
                # carries recoil energy.
                for angle, angular_weight in zip(angles, impact_weight):
                    for azimuth in azimuths:
                        projectile, target = _equal_mass_collision_velocities(
                            event_velocity, float(angle), float(azimuth))
                        continue_branch(
                            _canonical_axisymmetric_velocity(projectile),
                            _canonical_axisymmetric_velocity(target),
                            self.collision_model.elastic_probability
                            * float(angular_weight)
                            / len(azimuths),
                            charge_exchange=False,
                        )

                        # Charge exchange uses the same Born--Mayer encounter
                        # followed by a charge-label swap.  Equivalently the
                        # ion COM angle is pi-theta; the outgoing fast neutral
                        # has the elastic projectile's laboratory angle
                        # theta/2.  This is Khrabrov--Kaganovich's declared
                        # finite-angle replacement for a delta-function
                        # identity switch.
                        continue_branch(
                            _canonical_axisymmetric_velocity(target),
                            _canonical_axisymmetric_velocity(projectile),
                            self.collision_model.charge_exchange_probability
                            * float(angular_weight)
                            / len(azimuths),
                            charge_exchange=True,
                        )

        initial_nodes = tuple(self._initial_velocity_quadrature())
        for phase, phase_gain in zip(phases, gain):
            phase_optical_depth = []
            for initial_velocity, thermal_weight in initial_nodes:
                initial_weight = (
                    float(thermal_weight) / int(self.phase_count))
                target_energy = (
                    float(np.dot(initial_velocity, initial_velocity))
                    + float(phase_gain)
                )
                phase_optical_depth.append(density * self._integrated_hazard(
                    initial_velocity,
                    0.0,
                    1.0,
                    float(phase_gain),
                    1.0,
                ))
                recurse(
                    initial_velocity,
                    0.0,
                    float(phase_gain),
                    float(phase),
                    initial_weight,
                    0.0,
                    0,
                    0.0,
                    target_energy,
                )
            optical_depth.append(float(np.average(
                phase_optical_depth,
                weights=[item[1] for item in initial_nodes],
            )))

        ledger = arrived + unresolved + escaped
        ledger_tangent = (
            arrived_tangent + unresolved_tangent + escaped_tangent)
        probability_residual = abs(ledger - 1.0)
        neutral_lineage_residual = abs(
            neutral_arrival + neutral_unresolved + neutral_escaped
            - neutral_birth_count
        )
        neutral_lineage_tangent_residual = (
            neutral_arrival_tangent
            + neutral_unresolved_tangent
            + neutral_escaped_tangent
            - neutral_birth_count_tangent
        )
        if arrived <= 0.0:
            raise RuntimeError("collisional sheath delivered no resolved wafer ions")
        raw_weight = np.asarray(arriving_weight)
        raw_weight_tangent = np.asarray(arriving_weight_tangent)
        normalized = raw_weight / arrived
        normalized_tangent = (
            raw_weight_tangent * arrived
            - raw_weight * arrived_tangent
        ) / arrived ** 2
        distribution = CollisionalIonEnergyAngleDistribution(
            velocity_sqrt_eV=np.asarray(arriving_velocity),
            weight=normalized,
            entry_phase_rad=np.asarray(arriving_phase),
        )
        neutral_distribution = None
        if neutral_arrival > 0.0:
            neutral_distribution = CollisionalIonEnergyAngleDistribution(
                velocity_sqrt_eV=np.asarray(neutral_velocity_arriving),
                weight=(
                    np.asarray(neutral_weight_arriving) / neutral_arrival),
                entry_phase_rad=np.asarray(neutral_phase_arriving),
            )
        solution = DeterministicCollisionalSheathSolution(
            distribution=distribution,
            resolved_fast_neutral_distribution=neutral_distribution,
            source_ion_flux_m2_s=float(self.source_ion_flux_m2_s),
            arriving_ion_flux_m2_s=(
                float(self.source_ion_flux_m2_s) * arrived),
            resolved_fast_neutral_flux_m2_s=(
                float(self.source_ion_flux_m2_s) * neutral_arrival),
            ion_arrival_probability=arrived,
            unresolved_probability=unresolved,
            escaped_probability=escaped,
            uncollided_arrival_probability=uncollided,
            expected_collision_count_lower_bound=collision_count,
            expected_charge_exchange_count_lower_bound=charge_exchange_count,
            expected_fast_neutral_birth_count_lower_bound=(
                neutral_birth_count),
            expected_fast_neutral_birth_energy_lower_bound_eV_per_source_ion=(
                neutral_birth_energy),
            resolved_fast_neutral_arrivals_per_source_ion=neutral_arrival,
            unresolved_fast_neutral_collisions_per_source_ion=(
                neutral_unresolved),
            escaped_fast_neutrals_per_source_ion=neutral_escaped,
            fast_neutral_lineage_ledger_relative_residual=(
                neutral_lineage_residual),
            maximum_resolved_energy_ledger_relative_residual=(
                maximum_energy_residual),
            probability_ledger_relative_residual=probability_residual,
            collisionless_reference_mean_normal_energy_eV=float(
                np.mean(collisionless_energy)),
            mean_total_optical_depth=float(np.mean(optical_depth)),
            maximum_total_optical_depth=float(np.max(optical_depth)),
            below_born_mayer_support_collision_probability_lower_bound=(
                below_support),
            model_source=self.collision_model.source,
            provenance={
                **dict(self.collision_model.provenance),
                "solver": "deterministic_collision_order_gauss_quadrature",
                "phase_count": int(self.phase_count),
                "position_quadrature_order": int(
                    self.position_quadrature_order),
                "hazard_quadrature_order": int(
                    self.hazard_quadrature_order),
                "impact_quadrature_order": int(
                    self.impact_quadrature_order),
                "collision_azimuth_order": int(
                    self.collision_azimuth_order),
                "maximum_collision_order": int(
                    self.maximum_collision_order),
                "neutral_gas_temperature_K": float(
                    self.neutral_gas_temperature_K),
                "gas_number_density_m3": density,
                "fast_neutral_transport_closed": False,
                "fast_neutral_no_further_collision_branch_resolved": True,
                "moving_sheath_self_consistency_closed": False,
                "feature_depth_used": False,
            },
        )
        tangent = None
        if tangent_enabled:
            energy = distribution.energy_eV
            mean_energy_tangent = float(
                np.dot(normalized_tangent, energy))
            tangent = CollisionalSheathDensityTangent(
                gas_number_density_tangent_m3=density_tangent,
                distribution_weight_tangent=normalized_tangent,
                ion_arrival_probability_tangent=arrived_tangent,
                unresolved_probability_tangent=unresolved_tangent,
                escaped_probability_tangent=escaped_tangent,
                uncollided_arrival_probability_tangent=uncollided_tangent,
                expected_collision_count_tangent=collision_count_tangent,
                expected_charge_exchange_count_tangent=(
                    charge_exchange_count_tangent),
                expected_fast_neutral_birth_count_tangent=(
                    neutral_birth_count_tangent),
                expected_fast_neutral_birth_energy_tangent_eV=(
                    neutral_birth_energy_tangent),
                resolved_fast_neutral_arrivals_tangent=(
                    neutral_arrival_tangent),
                unresolved_fast_neutral_collisions_tangent=(
                    neutral_unresolved_tangent),
                escaped_fast_neutrals_tangent=neutral_escaped_tangent,
                fast_neutral_lineage_ledger_tangent_residual=(
                    neutral_lineage_tangent_residual),
                mean_impact_energy_tangent_eV=mean_energy_tangent,
                probability_ledger_tangent_residual=ledger_tangent,
            )
        return solution, tangent

    def solve(self) -> DeterministicCollisionalSheathSolution:
        solution, _ = self._solve(gas_number_density_tangent_m3=None)
        return solution

    def density_jvp(
        self,
        gas_number_density_tangent_m3: float,
    ) -> tuple[
        DeterministicCollisionalSheathSolution,
        CollisionalSheathDensityTangent,
    ]:
        """Exact JVP of the retained discrete collision-order operator."""
        solution, tangent = self._solve(
            gas_number_density_tangent_m3=(
                gas_number_density_tangent_m3))
        assert tangent is not None
        return solution, tangent
