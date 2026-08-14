"""Source-backed transport closures for the Lee--Lieberman argon model.

The ion closure retains the energy dependence of the Phelps Ar+--Ar momentum
transfer cross section and averages ``sigma * v`` over the relative-energy
distribution.  The metastable closure combines evaluated dilute-gas
self-diffusion with the Lee--Lieberman Knudsen limit.  None of these
coefficients is selected from a plasma-density, flux, or etch-depth target.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.integrate import quad

from .argon import ARGON_MASS_AMU
from .model import (
    ArgonGlobalCondition,
    ArgonTransportState,
    BOLTZMANN_J_K,
    PASCAL_PER_MTORR,
)
from .network import E_CHARGE_C

ATOMIC_MASS_UNIT_KG = 1.66053906660e-27
STANDARD_PRESSURE_PA = 101_325.0
ARGON_MASS_KG = ARGON_MASS_AMU * ATOMIC_MASS_UNIT_KG

_PHELPS_1994 = "phelps-1994-ar-ion-scattering"
_LEE_1994 = "lee-lieberman-1994-global"
_NIST_DIFFUSION_2024 = "nist-tn-2279-gas-diffusion"


def phelps_argon_momentum_transfer_cross_section_m2(
        laboratory_energy_eV: float | np.ndarray,
        ) -> float | np.ndarray:
    """Return Phelps' Ar+--Ar momentum-transfer cross section.

    Energy is the projectile laboratory energy for a stationary Ar target.
    For equal masses, ``E_lab = 2 E_cm``.  This convention matches the Phelps
    LXCat Ar+/Ar records, whose analytic component expressions evaluate
    ``Qm`` at ``2 E_cm``.

    The analytic approximation is Phelps' 1994 suggested momentum-transfer
    law,

    ``Qm = 1.15e-18 E^-0.1 (1 + 0.015/E)^0.6 m2``.
    """
    energy = np.asarray(laboratory_energy_eV, dtype=float)
    if np.any(~np.isfinite(energy)) or np.any(energy <= 0.0):
        raise ValueError(
            "laboratory collision energy must be positive and finite")
    cross_section = (
        1.15e-18
        * energy ** -0.1
        * (1.0 + 0.015 / energy) ** 0.6
    )
    if np.ndim(laboratory_energy_eV) == 0:
        return float(cross_section)
    return cross_section


def argon_relative_temperature_eV(
        ion_temperature_eV: float, gas_temperature_K: float) -> float:
    """Equivalent Maxwellian relative temperature for equal-mass Ar+ and Ar."""
    values = np.asarray(
        [ion_temperature_eV, gas_temperature_K], dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("ion and gas temperatures must be positive and finite")
    gas_temperature_eV = BOLTZMANN_J_K * gas_temperature_K / E_CHARGE_C
    return float(0.5 * (ion_temperature_eV + gas_temperature_eV))


@lru_cache(maxsize=256)
def phelps_argon_momentum_transfer_rate_m3_s(
        ion_temperature_eV: float, gas_temperature_K: float) -> float:
    """Maxwellian-average ``Qm(2 E_rel) * v_rel`` for Ar+ in Ar.

    The relative-energy distribution is integrated in dimensionless energy
    ``x = E_rel / T_rel``.  For equal ion and neutral masses the reduced mass
    is one half of the argon mass, and the equivalent stationary-target
    projectile energy used by the Phelps law is ``2 E_rel``.
    """
    relative_temperature = argon_relative_temperature_eV(
        ion_temperature_eV, gas_temperature_K)
    reduced_mass_kg = 0.5 * ARGON_MASS_KG
    speed_scale = np.sqrt(
        2.0 * E_CHARGE_C * relative_temperature / reduced_mass_kg)

    def integrand(x: float) -> float:
        if x == 0.0:
            return 0.0
        return float(
            x
            * np.exp(-x)
            * phelps_argon_momentum_transfer_cross_section_m2(
                2.0 * relative_temperature * x)
        )

    integral, error = quad(
        integrand,
        0.0,
        np.inf,
        epsabs=0.0,
        epsrel=2.0e-11,
        limit=200,
    )
    rate = 2.0 / np.sqrt(np.pi) * speed_scale * integral
    if (
        not np.isfinite(rate)
        or rate <= 0.0
        or not np.isfinite(error)
        or error > max(abs(integral), np.finfo(float).tiny) * 1.0e-8
    ):
        raise RuntimeError("argon momentum-transfer quadrature failed")
    return float(rate)


def lee_lieberman_argon_ion_temperature_eV(
        pressure_Pa: float, gas_temperature_K: float) -> float:
    """Return the pressure-dependent ion temperature used by Lee--Lieberman.

    The source holds ``Ti = 0.5 eV`` through 1 mTorr.  Above 1 mTorr,
    ``Ti - Tg`` decreases in proportion to ``1/p`` and asymptotes to the
    neutral temperature.
    """
    values = np.asarray([pressure_Pa, gas_temperature_K], dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("pressure and gas temperature must be positive")
    gas_temperature_eV = BOLTZMANN_J_K * gas_temperature_K / E_CHARGE_C
    if pressure_Pa <= PASCAL_PER_MTORR:
        return 0.5
    return float(
        gas_temperature_eV
        + (0.5 - gas_temperature_eV)
        * PASCAL_PER_MTORR
        / pressure_Pa
    )


def nist_argon_self_diffusion_m2_s(
        gas_temperature_K: float, pressure_Pa: float) -> float:
    """Evaluate NIST TN 2279's Ar-in-Ar dilute-gas correlation.

    NIST reports ``ln(D / (cm2 s-1)) = A + B/T + C ln(T)`` at standard
    atmospheric pressure, with ``A=-11.097``, ``B=-45.486 K``, and
    ``C=1.676``.  Dilute-gas diffusivity is scaled inversely with pressure.

    The NIST fit is tabulated over 235--418 K.  Lee--Lieberman's 600 K
    condition is therefore an explicit temperature extrapolation and is
    recorded as such in transport provenance.
    """
    values = np.asarray([gas_temperature_K, pressure_Pa], dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("gas temperature and pressure must be positive")
    diffusion_cm2_s_at_standard_pressure = np.exp(
        -11.097 - 45.486 / gas_temperature_K
        + 1.676 * np.log(gas_temperature_K)
    )
    return float(
        diffusion_cm2_s_at_standard_pressure
        * 1.0e-4
        * STANDARD_PRESSURE_PA
        / pressure_Pa
    )


@dataclass(frozen=True)
class LeeLiebermanArgonTransportProvider:
    """No-fit transport provider for published-model reproduction.

    This provider is deliberately classified ``published_model`` rather than
    ``validated_model``.  In particular, the 600 K Ar self-diffusion value is
    extrapolated beyond the evaluated NIST correlation's stated range.
    """

    name: str = "lee_lieberman_phelps_nist_argon_transport"
    version: str = "1"

    def predict(
            self, condition: ArgonGlobalCondition,
            electron_temperature_eV: float) -> ArgonTransportState:
        if not isinstance(condition, ArgonGlobalCondition):
            raise TypeError("argon condition is required")
        if (
            not np.isfinite(electron_temperature_eV)
            or electron_temperature_eV <= 0.0
        ):
            raise ValueError("electron temperature must be positive")

        ion_temperature_eV = lee_lieberman_argon_ion_temperature_eV(
            condition.pressure_Pa, condition.gas_temperature_K)
        momentum_rate_m3_s = phelps_argon_momentum_transfer_rate_m3_s(
            ion_temperature_eV, condition.gas_temperature_K)
        collision_frequency_s = (
            condition.neutral_ground_density_m3 * momentum_rate_m3_s)
        ion_thermal_speed_m_s = np.sqrt(
            E_CHARGE_C * ion_temperature_eV / ARGON_MASS_KG)
        ion_mean_free_path_m = (
            ion_thermal_speed_m_s / collision_frequency_s)
        ambipolar_diffusion_m2_s = (
            E_CHARGE_C * electron_temperature_eV
            / (ARGON_MASS_KG * collision_frequency_s)
        )

        bulk_diffusion_m2_s = nist_argon_self_diffusion_m2_s(
            condition.gas_temperature_K, condition.pressure_Pa)
        neutral_thermal_speed_m_s = np.sqrt(
            BOLTZMANN_J_K
            * condition.gas_temperature_K
            / ARGON_MASS_KG
        )
        knudsen_diffusion_m2_s = (
            neutral_thermal_speed_m_s
            * condition.geometry.diffusion_length_m
            / 3.0
        )
        effective_diffusion_m2_s = 1.0 / (
            1.0 / bulk_diffusion_m2_s
            + 1.0 / knudsen_diffusion_m2_s
        )
        nist_extrapolated = not (
            235.0 <= condition.gas_temperature_K <= 418.0)

        return ArgonTransportState(
            ion_mean_free_path_m=ion_mean_free_path_m,
            ambipolar_diffusion_m2_s=ambipolar_diffusion_m2_s,
            metastable_effective_diffusion_m2_s=effective_diffusion_m2_s,
            source=(
                f"{_PHELPS_1994}; {_LEE_1994}; "
                f"{_NIST_DIFFUSION_2024}"
            ),
            evidence_kind="published_model",
            provenance={
                "ion_temperature_eV": ion_temperature_eV,
                "relative_temperature_eV": argon_relative_temperature_eV(
                    ion_temperature_eV, condition.gas_temperature_K),
                "momentum_transfer_rate_m3_s": momentum_rate_m3_s,
                "collision_frequency_s-1": collision_frequency_s,
                "bulk_metastable_diffusion_m2_s": bulk_diffusion_m2_s,
                "knudsen_metastable_diffusion_m2_s":
                    knudsen_diffusion_m2_s,
                "nist_temperature_range_K": [235.0, 418.0],
                "nist_temperature_extrapolated": nist_extrapolated,
                "coefficient_selection_target": None,
            },
        )
