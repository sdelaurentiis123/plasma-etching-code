"""Atom-counted TiO2 depth arithmetic for reactor-to-feature audits.

This module intentionally stops before choosing a TiO2 surface yield.  It
converts between positive-ion dose, formula-unit removal yield, material
density, and depth.  That separation is useful for blind validation: a
reactor solution can be asked what surface yield it *requires* without using
the withheld profile to select the yield.

``feature_transmission`` is the run-averaged fraction of the wafer-plane
positive-ion dose delivered to the evolving feature floor.  It therefore
collects geometric shadowing, angular transport, and charging, but not surface
chemistry.  Keeping it explicit prevents a blanket etch yield from silently
becoming a high-aspect-ratio feature yield.
"""
from __future__ import annotations

from dataclasses import dataclass
import math


AVOGADRO_PER_MOL = 6.02214076e23  # exact in the SI
TIO2_MOLAR_MASS_KG_MOL = 79.866e-3


def _positive_finite(name: str, value: float) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return result


def tio2_formula_unit_density_m3(
    mass_density_kg_m3: float,
    *,
    molar_mass_kg_mol: float = TIO2_MOLAR_MASS_KG_MOL,
) -> float:
    """Return TiO2 formula units per cubic metre."""
    density = _positive_finite("mass density", mass_density_kg_m3)
    molar_mass = _positive_finite("molar mass", molar_mass_kg_mol)
    return density / molar_mass * AVOGADRO_PER_MOL


def required_formula_units_per_incident_ion(
    depth_nm: float,
    mass_density_kg_m3: float,
    positive_ion_flux_m2_s: float,
    duration_s: float,
    *,
    feature_transmission: float = 1.0,
) -> float:
    """Yield required to remove ``depth_nm`` at the supplied incident dose.

    The result is formula units removed per positive ion that reaches the
    surface reaction operator.  A feature transmission below one increases
    the required microscopic yield relative to a blanket surface.
    """
    depth_m = _positive_finite("depth", depth_nm) * 1.0e-9
    flux = _positive_finite("positive-ion flux", positive_ion_flux_m2_s)
    duration = _positive_finite("duration", duration_s)
    transmission = _positive_finite(
        "feature transmission", feature_transmission)
    if transmission > 1.0:
        raise ValueError("feature transmission cannot exceed one")
    number_density = tio2_formula_unit_density_m3(mass_density_kg_m3)
    delivered_ion_dose = flux * duration * transmission
    return depth_m * number_density / delivered_ion_dose


def depth_nm_from_positive_ion_dose(
    formula_units_per_incident_ion: float,
    mass_density_kg_m3: float,
    positive_ion_flux_m2_s: float,
    duration_s: float,
    *,
    feature_transmission: float = 1.0,
    film_thickness_nm: float | None = None,
) -> float:
    """Convert a positive-ion dose and effective removal yield to depth."""
    yield_value = _positive_finite(
        "formula-unit removal yield", formula_units_per_incident_ion)
    flux = _positive_finite("positive-ion flux", positive_ion_flux_m2_s)
    duration = _positive_finite("duration", duration_s)
    transmission = _positive_finite(
        "feature transmission", feature_transmission)
    if transmission > 1.0:
        raise ValueError("feature transmission cannot exceed one")
    density = tio2_formula_unit_density_m3(mass_density_kg_m3)
    depth_nm = yield_value * flux * duration * transmission / density * 1.0e9
    if film_thickness_nm is not None:
        cap = _positive_finite("film thickness", film_thickness_nm)
        depth_nm = min(depth_nm, cap)
    return depth_nm


def minimum_feature_transmission_for_depth(
    depth_nm: float,
    mass_density_kg_m3: float,
    positive_ion_flux_m2_s: float,
    duration_s: float,
    formula_units_per_incident_ion: float,
) -> float:
    """Run-averaged floor transmission needed for a candidate surface yield."""
    yield_value = _positive_finite(
        "formula-unit removal yield", formula_units_per_incident_ion)
    blanket_requirement = required_formula_units_per_incident_ion(
        depth_nm,
        mass_density_kg_m3,
        positive_ion_flux_m2_s,
        duration_s,
    )
    return blanket_requirement / yield_value


@dataclass(frozen=True)
class Tio2ClearanceGate:
    """Dimensional ledger for one reactor-dose sensitivity state."""

    depth_nm: float
    mass_density_kg_m3: float
    positive_ion_flux_m2_s: float
    duration_s: float
    feature_transmission: float
    formula_unit_density_m3: float
    incident_positive_ion_dose_m2: float
    delivered_positive_ion_dose_m2: float
    required_formula_units_per_incident_ion: float


def build_clearance_gate(
    depth_nm: float,
    mass_density_kg_m3: float,
    positive_ion_flux_m2_s: float,
    duration_s: float,
    *,
    feature_transmission: float = 1.0,
) -> Tio2ClearanceGate:
    """Build an auditable TiO2 atom/dose ledger."""
    flux = _positive_finite("positive-ion flux", positive_ion_flux_m2_s)
    duration = _positive_finite("duration", duration_s)
    transmission = _positive_finite(
        "feature transmission", feature_transmission)
    if transmission > 1.0:
        raise ValueError("feature transmission cannot exceed one")
    return Tio2ClearanceGate(
        depth_nm=_positive_finite("depth", depth_nm),
        mass_density_kg_m3=_positive_finite(
            "mass density", mass_density_kg_m3),
        positive_ion_flux_m2_s=flux,
        duration_s=duration,
        feature_transmission=transmission,
        formula_unit_density_m3=tio2_formula_unit_density_m3(
            mass_density_kg_m3),
        incident_positive_ion_dose_m2=flux * duration,
        delivered_positive_ion_dose_m2=flux * duration * transmission,
        required_formula_units_per_incident_ion=(
            required_formula_units_per_incident_ion(
                depth_nm,
                mass_density_kg_m3,
                flux,
                duration,
                feature_transmission=transmission,
            )
        ),
    )
