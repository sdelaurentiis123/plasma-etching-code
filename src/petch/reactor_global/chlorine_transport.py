"""Condition-scoped chlorine neutral transport for cylindrical reactors."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .chlorine_wall import (
    ChlorineWallRecombinationBoundary,
    chlorine_atom_mean_thermal_speed_m_s,
)
from .geometry import CylindricalReactor
from .neutral_transport import (
    CylindricalNeutralWallLoss,
    NeutralDiffusivityState,
    ReducedNeutralDiffusivity,
    solve_cylindrical_neutral_wall_loss,
)

ECONOMOU_CHLORINE_REDUCED_DIFFUSIVITY_M_INV_S = 6.21e20


def lymberopoulos_economou_1995_chlorine_diffusivity(
) -> ReducedNeutralDiffusivity:
    """Return the 500 K published-model Cl diffusivity.

    The source reports ``N D_Cl = 6.21e18 cm^-1 s^-1`` at its 500 K base
    condition.  The value is converted exactly to SI.  No temperature scaling
    is supplied because the Lennard--Jones parameters and collision integral
    used in the Chapman--Enskog calculation are not tabulated in the paper.
    """
    return ReducedNeutralDiffusivity(
        reduced_diffusivity_m_inv_s=(
            ECONOMOU_CHLORINE_REDUCED_DIFFUSIVITY_M_INV_S),
        reference_temperature_K=500.0,
        valid_temperature_K=(500.0, 500.0),
        source="lymberopoulos-economou-1995-chlorine-transport",
        evidence_kind="published_model",
        relative_uncertainty=None,
        provenance={
            "doi": "10.1109/27.467977",
            "source_value_cm-1_s-1": 6.21e18,
            "source_method": (
                "Chapman-Enskog using unlisted Lennard-Jones parameters"),
            "coefficient_selection_target": None,
            "temperature_conflict": (
                "Ramamurthi-Economou 2002 reused the same coefficient "
                "at 300 K"),
        },
    )


def ramamurthi_economou_2002_chlorine_diffusivity(
) -> ReducedNeutralDiffusivity:
    """Return the 300 K reproduction value, quarantined from prediction.

    The 2002 paper prints the same reduced diffusivity as the 500 K 1995
    source and cites that paper.  It is retained to reproduce the published
    model, not as evidence that the physical diffusivity is temperature
    independent.
    """
    return ReducedNeutralDiffusivity(
        reduced_diffusivity_m_inv_s=(
            ECONOMOU_CHLORINE_REDUCED_DIFFUSIVITY_M_INV_S),
        reference_temperature_K=300.0,
        valid_temperature_K=(300.0, 300.0),
        source="ramamurthi-economou-2002-chlorine-transport",
        evidence_kind="published_model",
        relative_uncertainty=None,
        provenance={
            "doi": "10.1116/1.1450581",
            "source_value_cm-1_s-1": 6.21e18,
            "source_reference": (
                "Lymberopoulos-Economou 1995, DOI 10.1109/27.467977"),
            "coefficient_selection_target": None,
            "temperature_conflict": (
                "same coefficient appears at 500 K in the cited source"),
        },
    )


@dataclass(frozen=True)
class ChlorineNeutralWallTransport:
    """Bulk-to-wall chlorine transport with explicit evidence composition."""

    diffusivity: NeutralDiffusivityState
    wall_loss: CylindricalNeutralWallLoss
    wall_boundary: ChlorineWallRecombinationBoundary

    def __post_init__(self):
        if (
            not isinstance(self.diffusivity, NeutralDiffusivityState)
            or not isinstance(self.wall_loss, CylindricalNeutralWallLoss)
            or not isinstance(
                self.wall_boundary, ChlorineWallRecombinationBoundary)
        ):
            raise TypeError("invalid chlorine neutral wall-transport state")
        if not np.isclose(
            self.diffusivity.diffusivity_m2_s,
            self.wall_loss.diffusivity_m2_s,
            rtol=1.0e-14,
            atol=0.0,
        ):
            raise ValueError("wall loss does not use the supplied diffusivity")
        if not np.isclose(
            self.wall_boundary.recombination_probability,
            self.wall_loss.wall_reaction_probability,
            rtol=0.0,
            atol=0.0,
        ):
            raise ValueError(
                "wall loss does not use the supplied recombination boundary")

    @property
    def supports_prediction(self) -> bool:
        return (
            self.diffusivity.supports_prediction
            and self.wall_boundary.supports_local_prediction
            and self.wall_loss.numerical_closure_passes
        )

    def evaluate_volume_rates(
        self,
        chlorine_atom_density_m3: float,
    ) -> "ChlorineVolumeWallRates":
        density = float(chlorine_atom_density_m3)
        if not np.isfinite(density) or density < 0.0:
            raise ValueError(
                "chlorine atom density must be finite and nonnegative")
        atom_loss = (
            self.wall_loss.exact_loss_frequency_s_inv * density)
        return ChlorineVolumeWallRates(
            chlorine_atom_loss_m3_s=atom_loss,
            chlorine_molecule_return_m3_s=0.5 * atom_loss,
        )


@dataclass(frozen=True)
class ChlorineVolumeWallRates:
    """Atom-conserving volume rates for ``2 Cl -> Cl2`` at the wall."""

    chlorine_atom_loss_m3_s: float
    chlorine_molecule_return_m3_s: float

    def __post_init__(self):
        values = np.asarray([
            self.chlorine_atom_loss_m3_s,
            self.chlorine_molecule_return_m3_s,
        ], dtype=float)
        if np.any(~np.isfinite(values)) or np.any(values < 0.0):
            raise ValueError("invalid chlorine volume wall rates")
        object.__setattr__(
            self, "chlorine_atom_loss_m3_s", float(values[0]))
        object.__setattr__(
            self, "chlorine_molecule_return_m3_s", float(values[1]))
        if not np.isclose(
            self.chlorine_atom_loss_m3_s,
            2.0 * self.chlorine_molecule_return_m3_s,
            rtol=1.0e-14,
            atol=0.0,
        ):
            raise ValueError("chlorine volume wall rates do not conserve atoms")

    @property
    def chlorine_atom_inventory_residual_m3_s(self) -> float:
        return float(
            -self.chlorine_atom_loss_m3_s
            + 2.0 * self.chlorine_molecule_return_m3_s)


def solve_chlorine_neutral_wall_transport(
    *,
    geometry: CylindricalReactor,
    wall_boundary: ChlorineWallRecombinationBoundary,
    diffusivity_model: ReducedNeutralDiffusivity,
    total_neutral_density_m3: float,
    gas_temperature_K: float,
    cl_to_cl2_ratio: float,
    pressure_Pa: float,
    icp_power_W: float,
) -> ChlorineNeutralWallTransport:
    """Compose source-scoped diffusivity, wall state, and exact geometry."""
    if not isinstance(wall_boundary, ChlorineWallRecombinationBoundary):
        raise TypeError("chlorine wall boundary is required")
    if not isinstance(diffusivity_model, ReducedNeutralDiffusivity):
        raise TypeError("reduced neutral diffusivity is required")
    wall_boundary.require_applicable(
        cl_to_cl2_ratio=cl_to_cl2_ratio,
        pressure_Pa=pressure_Pa,
        icp_power_W=icp_power_W,
        gas_temperature_K=gas_temperature_K,
    )
    diffusivity = diffusivity_model.evaluate(
        total_neutral_density_m3=total_neutral_density_m3,
        gas_temperature_K=gas_temperature_K,
    )
    wall_loss = solve_cylindrical_neutral_wall_loss(
        geometry=geometry,
        diffusivity_m2_s=diffusivity.diffusivity_m2_s,
        mean_thermal_speed_m_s=chlorine_atom_mean_thermal_speed_m_s(
            gas_temperature_K),
        wall_reaction_probability=(
            wall_boundary.recombination_probability),
    )
    return ChlorineNeutralWallTransport(
        diffusivity=diffusivity,
        wall_loss=wall_loss,
        wall_boundary=wall_boundary,
    )
