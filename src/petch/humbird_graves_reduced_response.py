"""Reduced transient response fitted to Humbird--Graves Si/CF2/F/Ar+ MD.

This is the kinetic companion to :mod:`petch.stratified_fluorocarbon_si`.
It predicts the observables exposed by the 2004 primary-author seminar:
retained C and F, instantaneous Si yield, and cumulative Si removal versus
CF2 fluence.  Its coefficients were regressed only on the protocol's
calibration panels (20 and 200 eV CF2/Ar+=9/1, plus 200 eV
CF2/F/Ar+=8/1/1).  The 7/2/1 case is calibration-excluded.

The model form is collision based:

* carbon incorporation is thermal sticking plus an NRT-displacement term;
* C loss is written with CF2 and atomic-F collision cross sections;
* retained F obeys an exact first-order site-capacity balance;
* the carbon-rich layer suppresses Si removal with two fitted percolation
  response functions; and
* direct atomic-F response is capped by one ion-renewed site,
  ``g(R)=2R/(1+R)``, where ``R=F/Ar+``.

All fitted numbers are surface-response coefficients.  No feature depth,
reactor flux normalization, Krueger endpoint, or optimizer-selected reactor
yield enters.  The response is a reduced representation of one classical-MD
potential and is not transferable to SiO2 or arbitrary reactor ions without a
new chemistry card and direct validation.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .ion_energy_deposition import nuclear_energy_in_layer_eV
from .stratified_fluorocarbon_si import SILICON_CARBIDE


_TRAPEZOID = getattr(np, "trapezoid", None) or np.trapz


@dataclass(frozen=True)
class HumbirdGravesResponseParameters:
    """Directly derived constants plus calibration-regressed coefficients."""

    monolayer_areal_density_m2: float = 6.78e18
    source_fluence_unit_m2: float = 1.0e19
    carbon_capacity_fraction_of_transport_layer: float = 0.5
    displacement_energy_eV: float = 25.0

    thermal_cf2_carbon_sticking: float = 5.97221050e-3
    carbon_incorporation_per_nrt_opportunity: float = 6.08583540e-2
    cf2_carbon_turnover_cross_section_m2: float = 7.50947163e-23
    atomic_f_carbon_turnover_cross_section_m2: float = 1.13069296e-21

    atomic_f_apparent_front_sticking: float = 1.0
    cf2_f_turnover_cross_section_m2: float = 9.53753699e-23
    ion_f_turnover_cross_section_m2: float = 1.40572631e-21
    atomic_f_f_turnover_cross_section_m2: float = 1.87036901e-20

    baseline_si_yield: float = 1.6858490e-1
    baseline_carbon_half_response_ml: float = 1.123515466e1
    baseline_carbon_percolation_exponent: float = 4.35429311
    atomic_f_si_yield_increment: float = 1.3325600e-1
    atomic_f_carbon_half_response_ml: float = 9.1815723
    atomic_f_carbon_percolation_exponent: float = 5.82307943

    evidence: Mapping[str, str] | None = None

    def __post_init__(self):
        positive = (
            self.monolayer_areal_density_m2,
            self.source_fluence_unit_m2,
            self.carbon_capacity_fraction_of_transport_layer,
            self.displacement_energy_eV,
            self.thermal_cf2_carbon_sticking,
            self.carbon_incorporation_per_nrt_opportunity,
            self.cf2_carbon_turnover_cross_section_m2,
            self.atomic_f_carbon_turnover_cross_section_m2,
            self.atomic_f_apparent_front_sticking,
            self.cf2_f_turnover_cross_section_m2,
            self.ion_f_turnover_cross_section_m2,
            self.atomic_f_f_turnover_cross_section_m2,
            self.baseline_si_yield,
            self.baseline_carbon_half_response_ml,
            self.baseline_carbon_percolation_exponent,
            self.atomic_f_si_yield_increment,
            self.atomic_f_carbon_half_response_ml,
            self.atomic_f_carbon_percolation_exponent,
        )
        if any(not np.isfinite(value) or value <= 0.0 for value in positive):
            raise ValueError("invalid Humbird--Graves response parameters")
        if not 0.0 < self.carbon_capacity_fraction_of_transport_layer <= 1.0:
            raise ValueError("carbon capacity fraction must lie in (0, 1]")
        if not 0.0 < self.atomic_f_apparent_front_sticking <= 1.0:
            raise ValueError("apparent F sticking must lie in (0, 1]")
        evidence = dict(self.evidence or {
            "monolayer_areal_density_m2": (
                "Si(100) surface atomic density used by the source ML axis"),
            "source_fluence_unit_m2": (
                "axis unit 10^15 cm^-2 converted to 10^19 m^-2"),
            "carbon_capacity_fraction_of_transport_layer": (
                "one C site per Si site in a dense Si-C layer"),
            "displacement_energy_eV": (
                "25 eV declared central value of the 10--80 eV "
                "polymer/near-surface displacement band"),
            "carbon_kinetics": (
                "regressed on calibration C-inventory panels only; "
                "20% F panel excluded"),
            "fluorine_kinetics": (
                "regressed on calibration F-inventory panels only; "
                "20% F panel excluded"),
            "silicon_yield": (
                "regressed on calibration instantaneous-yield and cumulative-"
                "Si panels only; no feature depth"),
            "atomic_f_renewal_response": (
                "one ion-renewed site cap: 2R/(1+R), fixed before evaluating "
                "the calibration-excluded R=2 case"),
        })
        required = {
            "monolayer_areal_density_m2",
            "source_fluence_unit_m2",
            "carbon_capacity_fraction_of_transport_layer",
            "displacement_energy_eV",
            "carbon_kinetics",
            "fluorine_kinetics",
            "silicon_yield",
            "atomic_f_renewal_response",
        }
        if set(evidence) != required or any(
                not isinstance(value, str) or not value
                for value in evidence.values()):
            raise ValueError("incomplete reduced-response provenance")
        object.__setattr__(self, "evidence", MappingProxyType(evidence))


class HumbirdGravesReducedResponse:
    """Evaluate the source-domain transient response without feature physics."""

    _DEPTH_BY_ENERGY_NM = {20.0: 1.5, 200.0: 3.0}

    def __init__(self, parameters=None):
        self.parameters = parameters or HumbirdGravesResponseParameters()
        self._nrt_cache = {}

    @staticmethod
    def _validate_boundary(
            energy_eV, cf2_per_ion, atomic_f_per_ion, fluence):
        energy = float(energy_eV)
        if energy not in HumbirdGravesReducedResponse._DEPTH_BY_ENERGY_NM:
            raise ValueError(
                "strict seminar response is declared only at 20 or 200 eV")
        if (not np.isfinite(cf2_per_ion) or cf2_per_ion <= 0.0
                or not np.isfinite(atomic_f_per_ion)
                or atomic_f_per_ion < 0.0):
            raise ValueError("invalid CF2/F/ion boundary ratio")
        fluence = np.asarray(fluence, dtype=float)
        if np.any(~np.isfinite(fluence)) or np.any(fluence < 0.0):
            raise ValueError("fluence must be finite and nonnegative")
        return energy, float(cf2_per_ion), float(atomic_f_per_ion), fluence

    def nrt_displacements_per_ion(self, energy_eV):
        energy = float(energy_eV)
        if energy not in self._DEPTH_BY_ENERGY_NM:
            raise ValueError("NRT source depth is declared only at 20 or 200 eV")
        cached = self._nrt_cache.get(energy)
        if cached is None:
            deposited = nuclear_energy_in_layer_eV(
                energy,
                1.0,
                self._DEPTH_BY_ENERGY_NM[energy],
                18,
                39.948,
                SILICON_CARBIDE,
            )
            cached = deposited / (
                2.0 * self.parameters.displacement_energy_eV)
            self._nrt_cache[energy] = cached
        return cached

    def carbon_capacity_ml(self, energy_eV):
        energy = float(energy_eV)
        if energy not in self._DEPTH_BY_ENERGY_NM:
            raise ValueError("source capacity is declared only at 20 or 200 eV")
        depth_m = self._DEPTH_BY_ENERGY_NM[energy] * 1.0e-9
        return (
            self.parameters.carbon_capacity_fraction_of_transport_layer
            * SILICON_CARBIDE.atom_density_m3
            * depth_m
            / self.parameters.monolayer_areal_density_m2
        )

    def _carbon_coefficients(
            self, energy_eV, cf2_per_ion, atomic_f_per_ion):
        par = self.parameters
        conversion = (
            par.source_fluence_unit_m2
            / par.monolayer_areal_density_m2
        )
        displacements = self.nrt_displacements_per_ion(energy_eV)
        source_ml_per_axis_unit = conversion * (
            par.thermal_cf2_carbon_sticking
            + par.carbon_incorporation_per_nrt_opportunity
            * displacements
            / cf2_per_ion
        )
        turnover_per_axis_unit = (
            par.cf2_carbon_turnover_cross_section_m2
            * par.source_fluence_unit_m2
            + par.atomic_f_carbon_turnover_cross_section_m2
            * par.source_fluence_unit_m2
            * atomic_f_per_ion
            / cf2_per_ion
        )
        capacity = self.carbon_capacity_ml(energy_eV)
        relaxation = source_ml_per_axis_unit / capacity + turnover_per_axis_unit
        return source_ml_per_axis_unit, relaxation, capacity

    def carbon_inventory_ml(
            self, cf2_fluence_1e15_cm2, *, energy_eV,
            cf2_per_ion, atomic_f_per_ion):
        energy, ratio, f_ratio, fluence = self._validate_boundary(
            energy_eV, cf2_per_ion, atomic_f_per_ion,
            cf2_fluence_1e15_cm2)
        source, relaxation, _ = self._carbon_coefficients(
            energy, ratio, f_ratio)
        return (
            source / relaxation
            * (1.0 - np.exp(-relaxation * fluence))
        )

    def fluorine_inventory_ml(
            self, cf2_fluence_1e15_cm2, *, energy_eV,
            cf2_per_ion, atomic_f_per_ion):
        energy, ratio, f_ratio, fluence = self._validate_boundary(
            energy_eV, cf2_per_ion, atomic_f_per_ion,
            cf2_fluence_1e15_cm2)
        par = self.parameters
        carbon_source, carbon_relaxation, capacity = (
            self._carbon_coefficients(energy, ratio, f_ratio))
        carbon_asymptote = carbon_source / carbon_relaxation
        conversion = (
            par.source_fluence_unit_m2
            / par.monolayer_areal_density_m2
        )
        direct_f_source = (
            par.atomic_f_apparent_front_sticking
            * conversion
            * f_ratio
            / ratio
        )
        # Two F sites per carbon capacity plus four active F bonds per one
        # source monolayer of interfacial Si.
        f_capacity_ml = 2.0 * capacity + 4.0
        f_relaxation = (
            par.cf2_f_turnover_cross_section_m2
            * par.source_fluence_unit_m2
            + par.ion_f_turnover_cross_section_m2
            * par.source_fluence_unit_m2
            * self.nrt_displacements_per_ion(energy)
            / ratio
            + par.atomic_f_f_turnover_cross_section_m2
            * par.source_fluence_unit_m2
            * f_ratio
            / ratio
            + direct_f_source / f_capacity_ml
        )
        constant_source = (
            2.0 * carbon_source
            * (1.0 - carbon_asymptote / capacity)
            + direct_f_source
        )
        exponential_source = (
            2.0 * carbon_source * carbon_asymptote / capacity)
        base = (
            constant_source / f_relaxation
            * (1.0 - np.exp(-f_relaxation * fluence))
        )
        difference = f_relaxation - carbon_relaxation
        if abs(difference) < 1.0e-12:
            transient = (
                exponential_source
                * fluence
                * np.exp(-f_relaxation * fluence)
            )
        else:
            transient = (
                exponential_source
                * (
                    np.exp(-carbon_relaxation * fluence)
                    - np.exp(-f_relaxation * fluence)
                )
                / difference
            )
        return base + transient

    def si_yield_per_ion(
            self, cf2_fluence_1e15_cm2, *, energy_eV,
            cf2_per_ion, atomic_f_per_ion):
        energy, ratio, f_ratio, fluence = self._validate_boundary(
            energy_eV, cf2_per_ion, atomic_f_per_ion,
            cf2_fluence_1e15_cm2)
        par = self.parameters
        carbon = self.carbon_inventory_ml(
            fluence,
            energy_eV=energy,
            cf2_per_ion=ratio,
            atomic_f_per_ion=f_ratio,
        )
        if energy <= 20.0:
            energy_factor = 0.0
        else:
            energy_factor = (
                (np.sqrt(energy) - np.sqrt(20.0))
                / (np.sqrt(200.0) - np.sqrt(20.0))
            )
        renewal = (
            0.0
            if f_ratio == 0.0
            else 2.0 * f_ratio / (1.0 + f_ratio)
        )
        baseline = (
            par.baseline_si_yield
            / (
                1.0
                + (
                    carbon / par.baseline_carbon_half_response_ml
                ) ** par.baseline_carbon_percolation_exponent
            )
        )
        atomic_f = (
            par.atomic_f_si_yield_increment
            * renewal
            / (
                1.0
                + (
                    carbon / par.atomic_f_carbon_half_response_ml
                ) ** par.atomic_f_carbon_percolation_exponent
            )
        )
        return energy_factor * (baseline + atomic_f)

    def cumulative_si_etch_ml(
            self, cf2_fluence_1e15_cm2, *, energy_eV,
            cf2_per_ion, atomic_f_per_ion, quadrature_points=4097):
        energy, ratio, f_ratio, fluence = self._validate_boundary(
            energy_eV, cf2_per_ion, atomic_f_per_ion,
            cf2_fluence_1e15_cm2)
        if int(quadrature_points) != quadrature_points or quadrature_points < 65:
            raise ValueError("quadrature_points must be an integer >= 65")
        flat = fluence.ravel()
        result = np.zeros(flat.shape, dtype=float)
        conversion = (
            self.parameters.source_fluence_unit_m2
            / ratio
            / self.parameters.monolayer_areal_density_m2
        )
        for index, endpoint in enumerate(flat):
            if endpoint == 0.0:
                continue
            grid = np.linspace(0.0, endpoint, int(quadrature_points))
            yield_values = self.si_yield_per_ion(
                grid,
                energy_eV=energy,
                cf2_per_ion=ratio,
                atomic_f_per_ion=f_ratio,
            )
            result[index] = _TRAPEZOID(yield_values, grid) * conversion
        result = result.reshape(fluence.shape)
        return float(result) if result.ndim == 0 else result
