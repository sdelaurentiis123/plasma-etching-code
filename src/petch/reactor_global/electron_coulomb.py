"""Deterministic differentiable Coulomb kinetics in electron-energy space."""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from .electron_kinetics import (
    ELECTRON_SPEED_PER_SQRT_EV_M_S,
    ElectronEnergyDistribution,
    ElectronEnergyGrid,
)
from .network import E_CHARGE_C


VACUUM_PERMITTIVITY_F_M = 8.8541878128e-12


def _scalar_or_array(value):
    array = np.asarray(value, dtype=float)
    return float(array) if array.ndim == 0 else array


@dataclass(frozen=True)
class ElectronElectronCoulombCoefficients:
    """Nonlinear isotropic Coulomb drift/diffusion on energy faces."""

    drift_eV_m3_s: np.ndarray
    diffusion_eV2_m3_s: np.ndarray
    coulomb_logarithm: float | np.ndarray
    kinetic_temperature_eV: float | np.ndarray
    lower_population_integral: np.ndarray
    lower_energy_integral_eV: np.ndarray
    upper_eepf_integral_eV_minus_1_over_2: np.ndarray
    model_id: str = "hagelaar_pitchford_2005_isotropic_ee_classical_log"
    supports_reactor_state_prediction: bool = False
    supports_wafer_flux: bool = False
    supports_feature_depth: bool = False


@dataclass(frozen=True)
class IsotropicElectronElectronCoulombKernel:
    """Hagelaar--Pitchford (2005) equations 34--41.

    The kernel evaluates exact piecewise-constant finite-volume moments
    ``A1``, ``A2``, and ``A3`` at every energy face and returns the nonlinear
    electron--electron contributions to reduced drift and diffusion. The
    classical Debye Coulomb logarithm is used exactly as printed in equation
    37. Electron--ion and anisotropic electron--electron momentum terms are
    deliberately outside this provider.
    """

    grid: ElectronEnergyGrid

    def __post_init__(self):
        if not isinstance(self.grid, ElectronEnergyGrid):
            raise TypeError("an ElectronEnergyGrid is required")

    @staticmethod
    def _validated_plasma_state(
        electron_to_neutral_density_ratio,
        gas_number_density_m3,
        batch_shape: tuple[int, ...],
    ) -> tuple[np.ndarray, np.ndarray]:
        ratio = np.asarray(electron_to_neutral_density_ratio, dtype=float)
        gas_density = np.asarray(gas_number_density_m3, dtype=float)
        try:
            ratio = np.broadcast_to(ratio, batch_shape)
            gas_density = np.broadcast_to(gas_density, batch_shape)
        except ValueError as exc:
            raise ValueError(
                "Coulomb plasma state cannot broadcast to EEPF batch"
            ) from exc
        if (
            np.any(~np.isfinite(ratio))
            or np.any(ratio <= 0.0)
            or np.any(~np.isfinite(gas_density))
            or np.any(gas_density <= 0.0)
        ):
            raise ValueError("Coulomb plasma densities must be positive")
        return np.asarray(ratio), np.asarray(gas_density)

    def _moments(self, values: np.ndarray) -> tuple[
        np.ndarray, np.ndarray, np.ndarray, np.ndarray
    ]:
        normalization_cells = values * self.grid.normalization_weights
        energy_cells = values * self.grid.mean_energy_weights_eV
        plain_cells = values * self.grid.cell_widths_eV
        lower_population = np.cumsum(
            normalization_cells, axis=-1)[..., :-1]
        lower_energy = np.cumsum(energy_cells, axis=-1)[..., :-1]
        upper_plain = np.flip(
            np.cumsum(np.flip(plain_cells, axis=-1), axis=-1),
            axis=-1,
        )[..., 1:]
        mean_energy = np.sum(energy_cells, axis=-1)
        return lower_population, lower_energy, upper_plain, mean_energy

    @staticmethod
    def _coulomb_state(
        mean_energy_eV: np.ndarray,
        ratio: np.ndarray,
        gas_density_m3: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        kinetic_temperature_eV = (2.0 / 3.0) * mean_energy_eV
        electron_density_m3 = ratio * gas_density_m3
        debye_particle_number = (
            12.0 * math.pi
            * (VACUUM_PERMITTIVITY_F_M
               * E_CHARGE_C * kinetic_temperature_eV) ** 1.5
            / (E_CHARGE_C ** 3 * np.sqrt(electron_density_m3))
        )
        if (
            np.any(~np.isfinite(debye_particle_number))
            or np.any(debye_particle_number <= 1.0)
        ):
            raise ValueError(
                "classical Debye Coulomb logarithm is outside its domain")
        coulomb_logarithm = np.log(debye_particle_number)
        base_coefficient = (
            E_CHARGE_C ** 2
            * ELECTRON_SPEED_PER_SQRT_EV_M_S
            / (24.0 * math.pi * VACUUM_PERMITTIVITY_F_M ** 2)
        )
        reduced_strength = base_coefficient * ratio * coulomb_logarithm
        return kinetic_temperature_eV, coulomb_logarithm, reduced_strength

    def evaluate(
        self,
        distribution: ElectronEnergyDistribution,
        *,
        electron_to_neutral_density_ratio,
        gas_number_density_m3,
    ) -> ElectronElectronCoulombCoefficients:
        if not isinstance(distribution, ElectronEnergyDistribution):
            raise TypeError("an ElectronEnergyDistribution is required")
        if distribution.grid != self.grid:
            raise ValueError("Coulomb kernel and EEPF grids differ")
        ratio, gas_density = self._validated_plasma_state(
            electron_to_neutral_density_ratio,
            gas_number_density_m3,
            distribution.batch_shape,
        )
        values = distribution.eepf_eV_minus_3_over_2
        lower_population, lower_energy, upper_plain, mean_energy = (
            self._moments(values))
        temperature, logarithm, strength = self._coulomb_state(
            mean_energy, ratio, gas_density)
        faces = self.grid.boundaries[1:-1]
        diffusion_shape = lower_energy + faces ** 1.5 * upper_plain
        drift = -3.0 * strength[..., np.newaxis] * lower_population
        diffusion = 2.0 * strength[..., np.newaxis] * diffusion_shape
        if np.any(diffusion <= 0.0) or np.any(~np.isfinite(diffusion)):
            raise FloatingPointError("invalid electron-electron diffusion")
        return ElectronElectronCoulombCoefficients(
            drift_eV_m3_s=np.asarray(drift),
            diffusion_eV2_m3_s=np.asarray(diffusion),
            coulomb_logarithm=_scalar_or_array(logarithm),
            kinetic_temperature_eV=_scalar_or_array(temperature),
            lower_population_integral=np.asarray(lower_population),
            lower_energy_integral_eV=np.asarray(lower_energy),
            upper_eepf_integral_eV_minus_1_over_2=np.asarray(upper_plain),
        )

    def jvp(
        self,
        distribution: ElectronEnergyDistribution,
        eepf_tangent,
        *,
        electron_to_neutral_density_ratio,
        electron_to_neutral_density_ratio_tangent=0.0,
        gas_number_density_m3,
        gas_number_density_m3_tangent=0.0,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Exact JVP of drift and diffusion with respect to all inputs."""
        state = self.evaluate(
            distribution,
            electron_to_neutral_density_ratio=(
                electron_to_neutral_density_ratio),
            gas_number_density_m3=gas_number_density_m3,
        )
        values = distribution.eepf_eV_minus_3_over_2
        tangent = np.asarray(eepf_tangent, dtype=float)
        if tangent.shape != values.shape or np.any(~np.isfinite(tangent)):
            raise ValueError("Coulomb EEPF tangent shape differs")
        ratio, gas_density = self._validated_plasma_state(
            electron_to_neutral_density_ratio,
            gas_number_density_m3,
            distribution.batch_shape,
        )
        ratio_tangent = np.broadcast_to(
            np.asarray(
                electron_to_neutral_density_ratio_tangent, dtype=float),
            distribution.batch_shape,
        )
        gas_tangent = np.broadcast_to(
            np.asarray(gas_number_density_m3_tangent, dtype=float),
            distribution.batch_shape,
        )
        if np.any(~np.isfinite(ratio_tangent)) or np.any(
            ~np.isfinite(gas_tangent)
        ):
            raise ValueError("Coulomb plasma-state tangent must be finite")
        d_a1, d_a2, d_a3, d_mean = self._moments(tangent)
        temperature = np.asarray(state.kinetic_temperature_eV)
        logarithm = np.asarray(state.coulomb_logarithm)
        d_temperature = (2.0 / 3.0) * d_mean
        d_logarithm = (
            1.5 * d_temperature / temperature
            - 0.5 * ratio_tangent / ratio
            - 0.5 * gas_tangent / gas_density
        )
        base_coefficient = (
            E_CHARGE_C ** 2
            * ELECTRON_SPEED_PER_SQRT_EV_M_S
            / (24.0 * math.pi * VACUUM_PERMITTIVITY_F_M ** 2)
        )
        strength = base_coefficient * ratio * logarithm
        d_strength = base_coefficient * (
            ratio_tangent * logarithm + ratio * d_logarithm)
        faces = self.grid.boundaries[1:-1]
        a1 = state.lower_population_integral
        diffusion_shape = (
            state.lower_energy_integral_eV
            + faces ** 1.5
            * state.upper_eepf_integral_eV_minus_1_over_2)
        d_diffusion_shape = d_a2 + faces ** 1.5 * d_a3
        d_drift = -3.0 * (
            d_strength[..., np.newaxis] * a1
            + strength[..., np.newaxis] * d_a1)
        d_diffusion = 2.0 * (
            d_strength[..., np.newaxis] * diffusion_shape
            + strength[..., np.newaxis] * d_diffusion_shape)
        return np.asarray(d_drift), np.asarray(d_diffusion)

    def vjp(
        self,
        distribution: ElectronEnergyDistribution,
        drift_cotangent,
        diffusion_cotangent,
        *,
        electron_to_neutral_density_ratio,
        gas_number_density_m3,
    ) -> tuple[np.ndarray, float | np.ndarray, float | np.ndarray]:
        """Exact VJP returning EEPF, ionization-degree, and gas gradients."""
        state = self.evaluate(
            distribution,
            electron_to_neutral_density_ratio=(
                electron_to_neutral_density_ratio),
            gas_number_density_m3=gas_number_density_m3,
        )
        drift_bar = np.asarray(drift_cotangent, dtype=float)
        diffusion_bar = np.asarray(diffusion_cotangent, dtype=float)
        expected = distribution.batch_shape + (self.grid.cell_count - 1,)
        if (
            drift_bar.shape != expected
            or diffusion_bar.shape != expected
            or np.any(~np.isfinite(drift_bar))
            or np.any(~np.isfinite(diffusion_bar))
        ):
            raise ValueError("Coulomb coefficient cotangent shape differs")
        ratio, gas_density = self._validated_plasma_state(
            electron_to_neutral_density_ratio,
            gas_number_density_m3,
            distribution.batch_shape,
        )
        logarithm = np.asarray(state.coulomb_logarithm)
        temperature = np.asarray(state.kinetic_temperature_eV)
        base_coefficient = (
            E_CHARGE_C ** 2
            * ELECTRON_SPEED_PER_SQRT_EV_M_S
            / (24.0 * math.pi * VACUUM_PERMITTIVITY_F_M ** 2)
        )
        strength = base_coefficient * ratio * logarithm
        faces = self.grid.boundaries[1:-1]
        a1 = state.lower_population_integral
        diffusion_shape = (
            state.lower_energy_integral_eV
            + faces ** 1.5
            * state.upper_eepf_integral_eV_minus_1_over_2)
        a1_bar = -3.0 * strength[..., np.newaxis] * drift_bar
        diffusion_shape_bar = (
            2.0 * strength[..., np.newaxis] * diffusion_bar)
        strength_bar = np.sum(
            -3.0 * a1 * drift_bar
            + 2.0 * diffusion_shape * diffusion_bar,
            axis=-1,
        )
        reverse_a1 = np.flip(
            np.cumsum(np.flip(a1_bar, axis=-1), axis=-1), axis=-1)
        reverse_a1 = np.concatenate((
            reverse_a1,
            np.zeros(distribution.batch_shape + (1,)),
        ), axis=-1)
        reverse_a2 = np.flip(
            np.cumsum(
                np.flip(diffusion_shape_bar, axis=-1), axis=-1),
            axis=-1,
        )
        reverse_a2 = np.concatenate((
            reverse_a2,
            np.zeros(distribution.batch_shape + (1,)),
        ), axis=-1)
        a3_bar = faces ** 1.5 * diffusion_shape_bar
        forward_a3 = np.concatenate((
            np.zeros(distribution.batch_shape + (1,)),
            np.cumsum(a3_bar, axis=-1),
        ), axis=-1)
        eepf_bar = (
            self.grid.normalization_weights * reverse_a1
            + self.grid.mean_energy_weights_eV * reverse_a2
            + self.grid.cell_widths_eV * forward_a3
        )
        temperature_bar = (
            strength_bar * base_coefficient * ratio * 1.5 / temperature)
        eepf_bar += (
            (2.0 / 3.0) * temperature_bar[..., np.newaxis]
            * self.grid.mean_energy_weights_eV)
        ratio_bar = strength_bar * base_coefficient * (
            logarithm - 0.5)
        gas_density_bar = (
            -0.5 * strength_bar * base_coefficient * ratio / gas_density)
        return (
            np.asarray(eepf_bar),
            _scalar_or_array(ratio_bar),
            _scalar_or_array(gas_density_bar),
        )
