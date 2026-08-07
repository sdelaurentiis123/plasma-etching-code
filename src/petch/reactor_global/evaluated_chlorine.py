"""Evaluated electron-collision data for predictive chlorine reactors.

This module is intentionally separate from the Lee--Lieberman reproduction
deck.  The values below are measured/evaluated cross sections, not parameters
tuned to a reactor state or an etched profile.
"""
from __future__ import annotations

from .network import ElectronMaxwellianCrossSectionRateCoefficient

ATOMIC_CHLORINE_IONIZATION_THRESHOLD_EV = 12.967633

_HAYES_ATOMIC_CHLORINE_IONIZATION_ENERGY_EV = (
    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
    22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 38.0, 40.0,
    45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0,
    95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0,
    135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 170.0, 180.0,
    190.0, 200.0,
)

_HAYES_ATOMIC_CHLORINE_IONIZATION_CROSS_SECTION_1E20_M2 = (
    0.00, 0.01, 0.02, 0.24, 0.52, 0.74, 1.01, 1.27, 1.50, 1.65,
    1.99, 2.34, 2.59, 2.80, 2.96, 3.16, 3.20, 3.27, 3.35, 3.35,
    3.43, 3.44, 3.47, 3.49, 3.49, 3.47, 3.44, 3.43, 3.43, 3.37,
    3.34, 3.31, 3.23, 3.20, 3.21, 3.15, 3.13, 3.07, 3.05, 3.01,
    2.97, 2.96, 2.91, 2.85, 2.81, 2.72, 2.68, 2.63,
)


def nist_hayes_atomic_chlorine_ionization_rate(
        *, maximum_kernel_tail_fraction: float = 1.0e-6
) -> ElectronMaxwellianCrossSectionRateCoefficient:
    """Return the NIST-recommended measured Cl-I ionization rate provider.

    Christophorou and Olthoff reproduce the selected Hayes et al. measurements
    in their Table 25 and quote a +/-14% absolute uncertainty.  NIST ASD's
    ground-state ionization limit supplies the physical threshold.  The raw
    table's tiny sub-threshold values are retained in the source transcription
    but excluded from the physical Maxwellian integral.
    """
    return ElectronMaxwellianCrossSectionRateCoefficient(
        electron_energy_eV=(
            _HAYES_ATOMIC_CHLORINE_IONIZATION_ENERGY_EV),
        cross_section_m2=tuple(
            value * 1.0e-20
            for value in (
                _HAYES_ATOMIC_CHLORINE_IONIZATION_CROSS_SECTION_1E20_M2)
        ),
        threshold_eV=ATOMIC_CHLORINE_IONIZATION_THRESHOLD_EV,
        relative_uncertainty=0.14,
        source=(
            "christophorou-olthoff-1999-cl2 Table 25; Hayes et al. "
            "measurements; NIST ASD Cl I threshold"
        ),
        evidence_kind="measured",
        maximum_kernel_tail_fraction=maximum_kernel_tail_fraction,
    )
