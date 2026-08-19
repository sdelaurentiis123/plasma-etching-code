"""Deterministic depth/mask integration for a square-pillar transport board.

The 3-D transport solver supplies depth-indexed, dimensionless floor and mask
dose factors.  These are normalized per projected horizontal area and can be
slightly above one when a resolved edge/bevel presents additional collecting
surface; they are not probabilities.  This module performs only the
atom-counted time integration.
It deliberately keeps the effective TiO2 removal yield and TiO2:Cr selectivity
as caller-supplied evidence axes; neither is inferred from a held-out profile.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .tio2_ion_dose import tio2_formula_unit_density_m3


def _positive_finite(name: str, value: float) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return result


def _transmission_curve(depth_nodes_nm, values, film_thickness_nm):
    depth = np.asarray(depth_nodes_nm, dtype=float)
    transmission = np.asarray(values, dtype=float)
    film = _positive_finite("film thickness", film_thickness_nm)
    if (
        depth.ndim != 1
        or transmission.shape != depth.shape
        or depth.size < 2
        or np.any(~np.isfinite(depth))
        or np.any(~np.isfinite(transmission))
        or depth[0] != 0.0
        or np.any(np.diff(depth) <= 0.0)
        or depth[-1] >= film
        or np.any(transmission <= 0.0)
        or np.any(transmission > 2.0)
    ):
        raise ValueError("invalid depth-indexed transmission curve")
    # The last transport snapshot intentionally retains a finite TiO2 floor.
    # Its limiting value is extended over the final unresolved sliver instead
    # of evaluating a post-clear fused-silica surface as though it were TiO2.
    extended_depth = np.append(depth, film)
    extended_transmission = np.append(transmission, transmission[-1])
    return extended_depth, extended_transmission


@dataclass(frozen=True)
class SquarePillarDepthResult:
    """One conditional 20-minute depth/mask trajectory."""

    blanket_rate_nm_s: float
    mask_pinned_depth_nm: float
    controlled_depth_nm: float
    mask_survives_duration: bool
    mask_exhaustion_time_s: float | None
    depth_at_mask_exhaustion_nm: float | None
    residual_mask_nm: float
    integration_step_s: float


def integrate_square_pillar_depth(
    *,
    depth_nodes_nm,
    floor_transmission,
    mask_transmission,
    film_thickness_nm: float,
    mask_thickness_nm: float,
    positive_ion_flux_m2_s: float,
    duration_s: float,
    mass_density_kg_m3: float,
    formula_units_per_incident_ion: float,
    tio2_to_cr_selectivity: float,
    integration_step_s: float = 1.0,
) -> SquarePillarDepthResult:
    """Integrate ideal-floor relief and Cr consumption in physical time.

    ``mask_pinned_depth_nm`` answers the transport/surface-dose problem with an
    indefinitely surviving mask.  ``controlled_depth_nm`` stops when the Cr
    thickness is exhausted.  No continuation after mask loss is manufactured:
    the post-exhaustion geometry is not identified by the supplied evidence.
    """
    flux = _positive_finite("positive-ion flux", positive_ion_flux_m2_s)
    density = _positive_finite("mass density", mass_density_kg_m3)
    removal_yield = _positive_finite(
        "formula-unit removal yield", formula_units_per_incident_ion)
    number_density = tio2_formula_unit_density_m3(density)
    blanket_rate = removal_yield * flux / number_density * 1.0e9
    return integrate_square_pillar_depth_from_blanket_rate(
        depth_nodes_nm=depth_nodes_nm,
        floor_transmission=floor_transmission,
        mask_transmission=mask_transmission,
        film_thickness_nm=film_thickness_nm,
        mask_thickness_nm=mask_thickness_nm,
        blanket_tio2_rate_nm_s=blanket_rate,
        duration_s=duration_s,
        tio2_to_cr_selectivity=tio2_to_cr_selectivity,
        integration_step_s=integration_step_s,
    )


def integrate_square_pillar_depth_from_blanket_rate(
    *,
    depth_nodes_nm,
    floor_transmission,
    mask_transmission,
    film_thickness_nm: float,
    mask_thickness_nm: float,
    blanket_tio2_rate_nm_s: float,
    duration_s: float,
    tio2_to_cr_selectivity: float,
    integration_step_s: float = 1.0,
) -> SquarePillarDepthResult:
    """Integrate from a directly supplied blanket TiO2 rate.

    This entry point is for a measured or explicitly cross-machine blanket-rate
    analog.  It prevents an observed rate from being disguised as a microscopic
    yield/wafer-flux pair.
    """
    film = _positive_finite("film thickness", film_thickness_nm)
    mask = _positive_finite("mask thickness", mask_thickness_nm)
    blanket_rate = _positive_finite("blanket TiO2 rate", blanket_tio2_rate_nm_s)
    duration = _positive_finite("duration", duration_s)
    selectivity = _positive_finite("TiO2:Cr selectivity", tio2_to_cr_selectivity)
    step = _positive_finite("integration step", integration_step_s)
    if step > duration:
        raise ValueError("integration step cannot exceed duration")
    count = int(round(duration / step))
    if not math.isclose(count * step, duration, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("duration must be an integer number of integration steps")

    depth_axis, floor_axis = _transmission_curve(
        depth_nodes_nm, floor_transmission, film)
    mask_depth, mask_axis = _transmission_curve(
        depth_nodes_nm, mask_transmission, film)
    if not np.array_equal(depth_axis, mask_depth):
        raise ValueError("floor and mask transmission nodes differ")

    def advance_depth(value, dt):
        if value >= film:
            return film
        # Midpoint update is deterministic, second order, and needs no
        # stochastic trajectory samples.
        first = float(np.interp(value, depth_axis, floor_axis))
        midpoint = min(film, value + 0.5 * dt * blanket_rate * first)
        second = float(np.interp(midpoint, depth_axis, floor_axis))
        return min(film, value + dt * blanket_rate * second)

    pinned_depth = 0.0
    controlled_depth = 0.0
    mask_loss = 0.0
    exhaustion_time = None
    exhaustion_depth = None
    for index in range(count):
        pinned_depth = advance_depth(pinned_depth, step)
        if exhaustion_time is not None:
            continue

        local_mask = float(np.interp(controlled_depth, depth_axis, mask_axis))
        next_depth = advance_depth(controlled_depth, step)
        next_mask = float(np.interp(next_depth, depth_axis, mask_axis))
        loss_increment = (
            blanket_rate * step * 0.5 * (local_mask + next_mask) / selectivity
        )
        if mask_loss + loss_increment >= mask:
            fraction = (mask - mask_loss) / loss_increment
            exhaustion_time = (index + fraction) * step
            exhaustion_depth = controlled_depth + fraction * (
                next_depth - controlled_depth)
            controlled_depth = exhaustion_depth
            mask_loss = mask
        else:
            mask_loss += loss_increment
            controlled_depth = next_depth

    survives = exhaustion_time is None
    return SquarePillarDepthResult(
        blanket_rate_nm_s=blanket_rate,
        mask_pinned_depth_nm=pinned_depth,
        controlled_depth_nm=controlled_depth,
        mask_survives_duration=survives,
        mask_exhaustion_time_s=exhaustion_time,
        depth_at_mask_exhaustion_nm=exhaustion_depth,
        residual_mask_nm=max(0.0, mask - mask_loss),
        integration_step_s=step,
    )
