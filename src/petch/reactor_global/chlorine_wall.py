"""Local, atom-conserving chlorine wall-recombination physics.

This module deliberately stops at the plasma-facing wall.  It converts a
local Cl number density and gas temperature into the isotropic Maxwellian
impingement flux and applies a condition-scoped recombination probability:

    2 Cl(wall) -> Cl2(gas)

It does not turn the local wall flux into a volume-averaged loss frequency.
That step requires neutral diffusion/Knudsen transport and the distribution
of conditioned surface states.  Omitting it here prevents a ballistic
well-mixed assumption from silently becoming reactor closure.
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
    "published_range_member",
    "assumed",
    "sensitivity",
})
_PREDICTIVE_EVIDENCE_KINDS = frozenset({"measured", "validated_model"})


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
    ) -> None:
        values = {
            "Cl/Cl2 ratio": (
                float(cl_to_cl2_ratio), self.valid_cl_to_cl2_ratio),
            "pressure": (float(pressure_Pa), self.valid_pressure_Pa),
            "ICP power": (float(icp_power_W), self.valid_icp_power_W),
        }
        for name, (value, domain) in values.items():
            if not np.isfinite(value) or not domain[0] <= value <= domain[1]:
                raise ValueError(
                    f"{name} is outside the wall-boundary evidence domain")

    def evaluate(
        self,
        *,
        chlorine_atom_density_m3: float,
        gas_temperature_K: float,
        cl_to_cl2_ratio: float,
        pressure_Pa: float,
        icp_power_W: float,
    ) -> "ChlorineWallFlux":
        self.require_applicable(
            cl_to_cl2_ratio=cl_to_cl2_ratio,
            pressure_Pa=pressure_Pa,
            icp_power_W=icp_power_W,
        )
        density = float(chlorine_atom_density_m3)
        if not np.isfinite(density) or density < 0.0:
            raise ValueError(
                "chlorine atom density must be finite and nonnegative")
        mean_speed = chlorine_atom_mean_thermal_speed_m_s(
            gas_temperature_K)
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
        )


@dataclass(frozen=True)
class ChlorineWallFlux:
    """Species fluxes at one plasma-facing wall location."""

    incident_cl_atom_flux_m2_s: float
    recombined_cl_atom_flux_m2_s: float
    returned_cl2_molecule_flux_m2_s: float
    mean_cl_speed_m_s: float
    source: str
    evidence_kind: str

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
