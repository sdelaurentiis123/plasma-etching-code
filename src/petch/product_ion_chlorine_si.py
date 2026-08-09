"""Auditable SiClx+ surface limits for chlorine/poly-Si depth closure.

Lee, Graves, and Lieberman include every ``SiClx+`` ion in the total ion flux
that drives their wafer etch law, but report no species-resolved product-ion
yield.  They separately bound redeposition with reflective and reactive wall
limits.  This module makes that unresolved branch explicit: it transfers the
measured Chang Cl+ full-chlorination yield card to product ions, integrates
each deterministic IEAD, and subtracts only the source-declared reactive-wall
deposition of Si+ and SiCl+.

The result is a source-model sensitivity, never a beam-validated prediction.
Its purpose is to quantify the exact experiment needed to decide whether the
missing product-ion channel can explain a depth residual.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .chang_sawin_chlorine_si import ChangSawinClIonSiParameters
from .reactor_global.wafer_sheath_transfer import (
    SpeciesResolvedIonEnergyDistribution,
)


PRODUCT_ION_MASS_AMU = MappingProxyType({
    "Si+": 28.0855,
    "SiCl+": 28.0855 + 35.453,
    "SiCl2+": 28.0855 + 2.0 * 35.453,
    "SiCl3+": 28.0855 + 3.0 * 35.453,
    "SiCl4+": 28.0855 + 4.0 * 35.453,
})


@dataclass(frozen=True)
class ProductIonSiSurfaceLimit:
    wall_limit: str
    gross_removal_rate_si_m2_s: Mapping[str, float]
    deposition_rate_si_m2_s: Mapping[str, float]
    net_removal_rate_si_m2_s: float
    source: str = (
        "Lee-Graves-Lieberman Eq. 3 all-ion etch participation plus Section "
        "2.3.1 redeposition limits; Chang Table 5.2 Cl+ yield transferred"
    )
    supports_prediction: bool = False

    def __post_init__(self):
        gross = {str(k): float(v) for k, v in self.gross_removal_rate_si_m2_s.items()}
        deposition = {str(k): float(v) for k, v in self.deposition_rate_si_m2_s.items()}
        if (
            self.wall_limit not in {"reflective", "reactive"}
            or set(gross) != set(PRODUCT_ION_MASS_AMU)
            or set(deposition) != set(PRODUCT_ION_MASS_AMU)
            or any(not math.isfinite(v) or v < 0.0
                   for v in (*gross.values(), *deposition.values()))
            or not math.isfinite(self.net_removal_rate_si_m2_s)
            or not str(self.source).strip()
            or self.supports_prediction
        ):
            raise ValueError("invalid product-ion Si surface limit")
        object.__setattr__(self, "gross_removal_rate_si_m2_s", MappingProxyType(gross))
        object.__setattr__(self, "deposition_rate_si_m2_s", MappingProxyType(deposition))

    @property
    def total_gross_removal_rate_si_m2_s(self) -> float:
        return float(sum(self.gross_removal_rate_si_m2_s.values()))

    @property
    def total_deposition_rate_si_m2_s(self) -> float:
        return float(sum(self.deposition_rate_si_m2_s.values()))


class LeeChangProductIonSiSurfaceSensitivity:
    """Integrate the printed all-ion limit using the measured Cl+ cards."""

    def __init__(self, parameters: ChangSawinClIonSiParameters | None = None):
        self.parameters = (
            ChangSawinClIonSiParameters.chang_thesis_table_5_2()
            if parameters is None else parameters
        )
        if not isinstance(self.parameters, ChangSawinClIonSiParameters):
            raise TypeError("parameters must be ChangSawinClIonSiParameters")

    def evaluate(
        self,
        distributions: Mapping[str, SpeciesResolvedIonEnergyDistribution],
        *,
        chlorination_fraction: float,
        wall_limit: str,
        allow_energy_extrapolation: bool = True,
    ) -> ProductIonSiSurfaceLimit:
        distributions = dict(distributions)
        theta = float(chlorination_fraction)
        if (
            set(distributions) != set(PRODUCT_ION_MASS_AMU)
            or any(
                not isinstance(item, SpeciesResolvedIonEnergyDistribution)
                or item.species != name
                or not math.isclose(
                    item.ion_mass_amu,
                    PRODUCT_ION_MASS_AMU[name],
                    rel_tol=2.0e-6,
                    abs_tol=0.0,
                )
                for name, item in distributions.items()
            )
            or not math.isfinite(theta)
            or not 0.0 <= theta <= 1.0
            or wall_limit not in {"reflective", "reactive"}
        ):
            raise ValueError("invalid product-ion surface input")
        gross = {}
        deposition = {}
        for name, distribution in distributions.items():
            _, beta = self.parameters.coefficients(
                distribution.energy_eV,
                allow_extrapolation=allow_energy_extrapolation,
            )
            gross[name] = float(
                distribution.flux_m2_s
                * theta
                * np.sum(distribution.weight * beta)
            )
            sticking = (
                0.5
                if wall_limit == "reactive" and name in {"Si+", "SiCl+"}
                else 0.0
            )
            deposition[name] = float(sticking * distribution.flux_m2_s)
        net = float(sum(gross.values()) - sum(deposition.values()))
        return ProductIonSiSurfaceLimit(
            wall_limit=wall_limit,
            gross_removal_rate_si_m2_s=gross,
            deposition_rate_si_m2_s=deposition,
            net_removal_rate_si_m2_s=net,
        )
