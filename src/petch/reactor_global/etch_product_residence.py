"""Zero-dimensional etch-product residence transfer to a wafer boundary.

Lee, Graves, and Lieberman (1996) identify long-lived SiClx etch products as
the mechanism coupling flow, reactor walls, and poly-Si etch rate.  They use
reactive- and reflective-wall limiting cases.  This module implements that
same identifiability boundary analytically: an independently supplied center
SiCl2/ion ratio fixes the unknown product-generation normalization, while
gross Si removal and the 0-D residence operator transfer it across recipes.

No feature depth or etch-rate observation is accepted by the class.  The
caller supplies a *predicted* gross Si source rate, so the enclosing pipeline
can solve the reactor/surface feedback as a deterministic scalar fixed point.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

from .wafer_power_transfer import isotropic_thermal_particle_flux_m2_s


SICL2_MASS_AMU = 98.991


@dataclass(frozen=True)
class EtchProductResidenceProjection:
    sicl2_flux_m2_s: float
    residence_time_s: float
    exhaust_loss_frequency_s_inv: float
    wall_loss_frequency_s_inv: float
    gross_si_source_scale: float
    residence_time_scale: float
    reference_sicl2_to_total_ion_ratio: float
    evidence_supports_prediction: bool = False

    def __post_init__(self):
        values = (
            self.sicl2_flux_m2_s,
            self.residence_time_s,
            self.exhaust_loss_frequency_s_inv,
            self.wall_loss_frequency_s_inv,
            self.gross_si_source_scale,
            self.residence_time_scale,
            self.reference_sicl2_to_total_ion_ratio,
        )
        if (
            any(not math.isfinite(value) or value < 0.0 for value in values)
            or self.residence_time_s <= 0.0
            or self.residence_time_scale <= 0.0
            or self.evidence_supports_prediction
        ):
            raise ValueError("invalid etch-product residence projection")


@dataclass(frozen=True)
class DiagnosticConditionedEtchProductResidenceTransfer:
    """Transfer a center product-flux estimate through a 0-D loss operator."""

    reference_total_ion_flux_m2_s: float
    reference_sicl2_to_total_ion_ratio: float
    reference_gross_si_source_rate_m2_s: float
    reference_exhaust_loss_frequency_s_inv: float
    reactor_volume_m3: float
    reactor_physical_area_m2: float
    gas_temperature_K: float
    wall_reactivity: float
    sicl2_sticking_coefficient: float = 0.3
    source: str = (
        "Mahorowala Chapter 4 source-plane SiCl2/ion estimate with "
        "Lee-Graves-Lieberman reactive/reflective wall limits"
    )

    def __post_init__(self):
        values = (
            self.reference_total_ion_flux_m2_s,
            self.reference_sicl2_to_total_ion_ratio,
            self.reference_gross_si_source_rate_m2_s,
            self.reference_exhaust_loss_frequency_s_inv,
            self.reactor_volume_m3,
            self.reactor_physical_area_m2,
            self.gas_temperature_K,
            self.wall_reactivity,
            self.sicl2_sticking_coefficient,
        )
        if (
            any(not math.isfinite(value) for value in values)
            or self.reference_total_ion_flux_m2_s <= 0.0
            or self.reference_sicl2_to_total_ion_ratio < 0.0
            or self.reference_gross_si_source_rate_m2_s <= 0.0
            or self.reference_exhaust_loss_frequency_s_inv <= 0.0
            or self.reactor_volume_m3 <= 0.0
            or self.reactor_physical_area_m2 <= 0.0
            or self.gas_temperature_K <= 0.0
            or not 0.0 <= self.wall_reactivity <= 1.0
            or not 0.0 <= self.sicl2_sticking_coefficient <= 1.0
            or not str(self.source).strip()
        ):
            raise ValueError("invalid etch-product residence transfer")

    @property
    def wall_loss_frequency_s_inv(self) -> float:
        unit_density_flux = isotropic_thermal_particle_flux_m2_s(
            1.0,
            self.gas_temperature_K,
            mass_amu=SICL2_MASS_AMU,
        )
        return float(
            self.wall_reactivity
            * self.sicl2_sticking_coefficient
            * unit_density_flux
            * self.reactor_physical_area_m2
            / self.reactor_volume_m3
        )

    def residence_time_s(self, exhaust_loss_frequency_s_inv: float) -> float:
        exhaust = float(exhaust_loss_frequency_s_inv)
        if not math.isfinite(exhaust) or exhaust <= 0.0:
            raise ValueError("exhaust loss frequency must be positive")
        return float(1.0 / (exhaust + self.wall_loss_frequency_s_inv))

    @property
    def reference_residence_time_s(self) -> float:
        return self.residence_time_s(
            self.reference_exhaust_loss_frequency_s_inv
        )

    def predict(
        self,
        *,
        gross_si_source_rate_m2_s: float,
        exhaust_loss_frequency_s_inv: float,
    ) -> EtchProductResidenceProjection:
        source_rate = float(gross_si_source_rate_m2_s)
        exhaust = float(exhaust_loss_frequency_s_inv)
        if not math.isfinite(source_rate) or source_rate < 0.0:
            raise ValueError("gross Si source rate must be nonnegative")
        residence = self.residence_time_s(exhaust)
        source_scale = source_rate / self.reference_gross_si_source_rate_m2_s
        residence_scale = residence / self.reference_residence_time_s
        reference_flux = (
            self.reference_sicl2_to_total_ion_ratio
            * self.reference_total_ion_flux_m2_s
        )
        return EtchProductResidenceProjection(
            sicl2_flux_m2_s=(
                reference_flux * source_scale * residence_scale
            ),
            residence_time_s=residence,
            exhaust_loss_frequency_s_inv=exhaust,
            wall_loss_frequency_s_inv=self.wall_loss_frequency_s_inv,
            gross_si_source_scale=source_scale,
            residence_time_scale=residence_scale,
            reference_sicl2_to_total_ion_ratio=(
                self.reference_sicl2_to_total_ion_ratio
            ),
        )
