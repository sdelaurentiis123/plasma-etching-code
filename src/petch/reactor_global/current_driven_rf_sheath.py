"""Deterministic current-driven RF sheath with a moving electron front.

This module implements the Turner--Chabert arbitrary-waveform sheath model
(Applied Physics Letters 104, 164102, 2014; arXiv:1212.2612).  Unlike the
phase-conditioned effective Child lift used by the first collisional-transport
rung, the electric field here is explicitly a function of position and time.

The input is the *sheath RF current-density waveform*, not generator power.
That distinction is physical: a matching network and the plasma/electrode
impedances are required to infer this waveform from a generator setpoint.
Consequently this provider closes current -> moving sheath -> collisionless
IED, but deliberately does not claim generator-to-wafer closure.

The implementation is deterministic and differentiable.  Fourier current
waveforms give an analytic charge trajectory.  The nonlinear sheath scale is
then obtained algebraically from the same Bohm-current, Poisson, and
time-averaged Child equations used in the source.  Ion trajectories use a
fixed, vectorized velocity-Verlet integration through the moving field.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from types import MappingProxyType
from typing import Mapping

import numpy as np

from ..sheath import AMU, ECHARGE, EPS0, PeriodicSheathVoltage, bohm_speed


_EVIDENCE_KINDS = {
    "measured_sheath_current",
    "validated_circuit_model",
    "assumed",
}


def _readonly_1d(value, *, dtype=float) -> np.ndarray:
    array = np.asarray(value, dtype=dtype).copy()
    if array.ndim != 1 or np.any(~np.isfinite(array)):
        raise ValueError("expected a finite one-dimensional array")
    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class PeriodicCurrentDensity:
    """Zero-mean Fourier representation of sheath RF current density.

    The waveform is

    ``J(t) = sum[J_s,n sin(n omega t) + J_c,n cos(n omega t)]``.

    DC ion and electron conduction currents are not included in this object;
    their cycle averages cancel in the periodic sheath state.  The current
    must be the current at the sheath/electrode reference plane, after any
    matching-network and stray-current de-embedding.
    """

    fundamental_frequency_hz: float
    harmonic_number: np.ndarray
    sine_A_m2: np.ndarray
    cosine_A_m2: np.ndarray
    source: str
    evidence_kind: str = "assumed"

    def __post_init__(self):
        harmonic = _readonly_1d(self.harmonic_number, dtype=int)
        sine = _readonly_1d(self.sine_A_m2)
        cosine = _readonly_1d(self.cosine_A_m2)
        if (
            not math.isfinite(float(self.fundamental_frequency_hz))
            or self.fundamental_frequency_hz <= 0.0
            or harmonic.size == 0
            or sine.shape != harmonic.shape
            or cosine.shape != harmonic.shape
            or np.any(harmonic <= 0)
            or len(np.unique(harmonic)) != harmonic.size
            or float(np.max(np.hypot(sine, cosine))) <= 0.0
            or not str(self.source).strip()
            or self.evidence_kind not in _EVIDENCE_KINDS
        ):
            raise ValueError("invalid periodic sheath-current waveform")
        order = np.argsort(harmonic)
        object.__setattr__(self, "harmonic_number", harmonic[order])
        object.__setattr__(self, "sine_A_m2", sine[order])
        object.__setattr__(self, "cosine_A_m2", cosine[order])

    @property
    def period_s(self) -> float:
        return 1.0 / float(self.fundamental_frequency_hz)

    @property
    def maximum_frequency_hz(self) -> float:
        return float(
            self.fundamental_frequency_hz * self.harmonic_number[-1])

    @property
    def supports_predictive_boundary(self) -> bool:
        return self.evidence_kind in {
            "measured_sheath_current", "validated_circuit_model"}

    def current_density_A_m2(self, time_s) -> np.ndarray | float:
        time = np.asarray(time_s, dtype=float)
        angle = (
            2.0
            * np.pi
            * self.fundamental_frequency_hz
            * time[..., None]
            * self.harmonic_number
        )
        result = np.sum(
            self.sine_A_m2 * np.sin(angle)
            + self.cosine_A_m2 * np.cos(angle),
            axis=-1,
        )
        return float(result) if time.ndim == 0 else result

    def charge_primitive_C_m2(self, time_s) -> np.ndarray | float:
        """Periodic antiderivative whose time derivative is ``J(t)``."""
        time = np.asarray(time_s, dtype=float)
        omega = 2.0 * np.pi * self.fundamental_frequency_hz
        angle = omega * time[..., None] * self.harmonic_number
        divisor = omega * self.harmonic_number
        result = np.sum(
            -self.sine_A_m2 / divisor * np.cos(angle)
            + self.cosine_A_m2 / divisor * np.sin(angle),
            axis=-1,
        )
        return float(result) if time.ndim == 0 else result

    def scaled(self, factor: float, *, source: str | None = None):
        value = float(factor)
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError("current-density scale must be positive")
        return PeriodicCurrentDensity(
            fundamental_frequency_hz=self.fundamental_frequency_hz,
            harmonic_number=self.harmonic_number,
            sine_A_m2=value * self.sine_A_m2,
            cosine_A_m2=value * self.cosine_A_m2,
            source=self.source if source is None else str(source),
            evidence_kind=self.evidence_kind,
        )


@dataclass(frozen=True)
class TurnerChabertCurrentSheathTangent:
    """Exact tangent for a uniform multiplicative RF-current perturbation."""

    logarithmic_current_scale_tangent: float
    maximum_width_tangent_m: float
    maximum_voltage_tangent_v: float
    mean_voltage_tangent_v: float
    charge_excursion_tangent_C_m2: float

    def __post_init__(self):
        values = np.asarray([
            self.logarithmic_current_scale_tangent,
            self.maximum_width_tangent_m,
            self.maximum_voltage_tangent_v,
            self.mean_voltage_tangent_v,
            self.charge_excursion_tangent_C_m2,
        ])
        if np.any(~np.isfinite(values)):
            raise ValueError("invalid current-driven sheath tangent")


@dataclass(frozen=True)
class TurnerChabertCurrentDrivenSheath:
    """Self-consistent arbitrary-current RF sheath and moving Poisson field."""

    current: PeriodicCurrentDensity
    electron_temperature_eV: float
    ion_mass_amu: float
    sheath_edge_density_m3: float
    phase_quadrature_count: int = 4096
    source: str = "turner-chabert-2014-rf-sheath"
    provenance: Mapping[str, object] | None = None
    _primitive_min_C_m2: float = field(init=False, repr=False, compare=False)
    _charge_excursion_C_m2: float = field(
        init=False, repr=False, compare=False)
    _xi: float = field(init=False, repr=False, compare=False)
    _maximum_width_m: float = field(init=False, repr=False, compare=False)
    _maximum_voltage_v: float = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        values = np.asarray([
            self.electron_temperature_eV,
            self.ion_mass_amu,
            self.sheath_edge_density_m3,
        ])
        if (
            not isinstance(self.current, PeriodicCurrentDensity)
            or np.any(~np.isfinite(values))
            or np.any(values <= 0.0)
            or int(self.phase_quadrature_count) < 256
            or not str(self.source).strip()
        ):
            raise ValueError("invalid Turner-Chabert sheath condition")
        object.__setattr__(
            self,
            "provenance",
            MappingProxyType({} if self.provenance is None
                             else dict(self.provenance)),
        )
        # Cache the nonlinear cycle closure once.  Moving-field kinetic
        # characteristics call ``electric_field_V_m`` thousands of times;
        # rebuilding the 4096-point normalization grid there would preserve
        # the equations but destroy the intended reduced-model performance.
        phase = (
            2.0
            * np.pi
            * np.arange(int(self.phase_quadrature_count), dtype=float)
            / int(self.phase_quadrature_count)
        )
        time = phase / (
            2.0 * np.pi * self.current.fundamental_frequency_hz)
        primitive = np.asarray(self.current.charge_primitive_C_m2(time))
        lower = float(np.min(primitive))
        excursion = float(np.max(primitive) - lower)
        if excursion <= 0.0:
            raise ValueError("RF current does not produce a valid sheath cycle")
        normalized = np.clip((primitive - lower) / excursion, 0.0, 1.0)
        xi = float(np.mean(self._voltage_shape(normalized)))
        if not 0.0 < xi < 1.0:
            raise ValueError("RF current does not produce a valid sheath cycle")
        mass = self.ion_mass_amu * AMU
        ion_current = (
            ECHARGE
            * self.sheath_edge_density_m3
            * bohm_speed(self.electron_temperature_eV, self.ion_mass_amu)
        )
        child_coefficient = 4.0 / (9.0 * xi)
        coefficient = (
            child_coefficient
            * EPS0
            * math.sqrt(2.0 * ECHARGE / mass)
            * (xi * 3.0 * excursion / (4.0 * EPS0)) ** 1.5
        )
        maximum_width = float((coefficient / ion_current) ** 2)
        maximum_voltage = float(
            3.0 * maximum_width * excursion / (4.0 * EPS0))
        object.__setattr__(self, "_primitive_min_C_m2", lower)
        object.__setattr__(self, "_charge_excursion_C_m2", excursion)
        object.__setattr__(self, "_xi", xi)
        object.__setattr__(self, "_maximum_width_m", maximum_width)
        object.__setattr__(self, "_maximum_voltage_v", maximum_voltage)

    @property
    def Te_eV(self) -> float:
        """Compatibility alias for the common sheath transport interface."""
        return float(self.electron_temperature_eV)

    @property
    def density_m3(self) -> float:
        return float(self.sheath_edge_density_m3)

    @property
    def phase_grid_rad(self) -> np.ndarray:
        return (
            2.0
            * np.pi
            * np.arange(int(self.phase_quadrature_count), dtype=float)
            / int(self.phase_quadrature_count)
        )

    def _cycle_primitive(self) -> np.ndarray:
        return np.asarray(self.current.charge_primitive_C_m2(
            self.phase_grid_rad
            / (2.0 * np.pi * self.current.fundamental_frequency_hz)
        ))

    @property
    def charge_excursion_C_m2(self) -> float:
        return self._charge_excursion_C_m2

    def normalized_charge(self, time_s) -> np.ndarray | float:
        primitive = np.asarray(self.current.charge_primitive_C_m2(time_s))
        normalized = np.clip(
            (primitive - self._primitive_min_C_m2)
            / self._charge_excursion_C_m2,
            0.0,
            1.0,
        )
        return float(normalized) if primitive.ndim == 0 else normalized

    @staticmethod
    def _voltage_shape(normalized_charge) -> np.ndarray | float:
        y = np.asarray(normalized_charge, dtype=float)
        value = 1.0 - (4.0 / 3.0) * y + (1.0 / 3.0) * y ** 4
        return float(value) if y.ndim == 0 else value

    @property
    def xi(self) -> float:
        return self._xi

    @property
    def ion_current_density_A_m2(self) -> float:
        return float(
            ECHARGE
            * self.sheath_edge_density_m3
            * bohm_speed(self.electron_temperature_eV, self.ion_mass_amu)
        )

    @property
    def maximum_width_m(self) -> float:
        """Maximum sheath width from source equations (2)--(5), (13)."""
        return self._maximum_width_m

    @property
    def thickness(self) -> float:
        """Compatibility alias: the transport domain is the maximum width."""
        return self.maximum_width_m

    @property
    def maximum_voltage_v(self) -> float:
        return self._maximum_voltage_v

    @property
    def mean_voltage_v(self) -> float:
        return float(self.xi * self.maximum_voltage_v)

    @property
    def maximum_frequency_hz(self) -> float:
        return self.current.maximum_frequency_hz

    @property
    def period_s(self) -> float:
        return self.current.period_s

    @property
    def supports_predictive_boundary(self) -> bool:
        return self.current.supports_predictive_boundary

    def voltage(self, time_s) -> np.ndarray | float:
        return self.maximum_voltage_v * self._voltage_shape(
            self.normalized_charge(time_s))

    def electron_front_fraction(self, time_s) -> np.ndarray | float:
        """Return ``s(t)/s_max`` from source equation (14)."""
        y = np.asarray(self.normalized_charge(time_s))
        result = y ** 3
        return float(result) if y.ndim == 0 else result

    def potential_drop_v(self, position_m, time_s) -> np.ndarray | float:
        """Potential drop from the instantaneous electron front to ``x``."""
        position, time = np.broadcast_arrays(
            np.asarray(position_m, dtype=float), np.asarray(time_s, dtype=float))
        if np.any(~np.isfinite(position)) or np.any(~np.isfinite(time)):
            raise ValueError("non-finite sheath field coordinate")
        r = np.clip(position / self.maximum_width_m, 0.0, 1.0)
        y = np.asarray(self.normalized_charge(time))
        depleted = r >= y ** 3
        value = self.maximum_voltage_v * (
            r ** (4.0 / 3.0)
            - (4.0 / 3.0) * y * r
            + (1.0 / 3.0) * y ** 4
        )
        result = np.where(depleted, np.maximum(value, 0.0), 0.0)
        return float(result) if result.ndim == 0 else result

    def electric_field_V_m(self, position_m, time_s) -> np.ndarray | float:
        """Positive field magnitude accelerating ions toward the electrode."""
        position, time = np.broadcast_arrays(
            np.asarray(position_m, dtype=float), np.asarray(time_s, dtype=float))
        if np.any(~np.isfinite(position)) or np.any(~np.isfinite(time)):
            raise ValueError("non-finite sheath field coordinate")
        r = np.clip(position / self.maximum_width_m, 0.0, 1.0)
        y = np.asarray(self.normalized_charge(time))
        depleted = r >= y ** 3
        field = (
            (4.0 / 3.0)
            * self.maximum_voltage_v
            / self.maximum_width_m
            * (np.cbrt(r) - y)
        )
        result = np.where(depleted, np.maximum(field, 0.0), 0.0)
        return float(result) if result.ndim == 0 else result

    @property
    def child_current_relative_residual(self) -> float:
        mass = self.ion_mass_amu * AMU
        reconstructed = (
            (4.0 / (9.0 * self.xi))
            * EPS0
            / self.maximum_width_m ** 2
            * math.sqrt(2.0 * ECHARGE / mass)
            * self.mean_voltage_v ** 1.5
        )
        return float(abs(reconstructed - self.ion_current_density_A_m2)
                     / self.ion_current_density_A_m2)

    @property
    def charge_voltage_relative_residual(self) -> float:
        reconstructed = (
            4.0
            * EPS0
            * self.maximum_voltage_v
            / (3.0 * self.maximum_width_m)
        )
        return float(abs(reconstructed - self.charge_excursion_C_m2)
                     / self.charge_excursion_C_m2)

    def current_scale_jvp(
        self, logarithmic_current_scale_tangent: float
    ) -> TurnerChabertCurrentSheathTangent:
        """Exact JVP for ``J -> exp(alpha) J`` at ``alpha=0``.

        The normalized charge waveform and ``xi`` are invariant to a common
        current scale.  Source equations then give ``s_max ~ J^3`` and
        ``V_max ~ J^4`` exactly.
        """
        direction = float(logarithmic_current_scale_tangent)
        if not math.isfinite(direction):
            raise ValueError("current-scale tangent must be finite")
        return TurnerChabertCurrentSheathTangent(
            logarithmic_current_scale_tangent=direction,
            maximum_width_tangent_m=(
                3.0 * self.maximum_width_m * direction),
            maximum_voltage_tangent_v=(
                4.0 * self.maximum_voltage_v * direction),
            mean_voltage_tangent_v=(
                4.0 * self.mean_voltage_v * direction),
            charge_excursion_tangent_C_m2=(
                self.charge_excursion_C_m2 * direction),
        )

    def ion_impact_energies(
        self,
        phases,
        steps_per_period: int = 512,
        steps_per_transit: int = 512,
        max_periods: float = 100.0,
    ) -> np.ndarray:
        """Integrate Bohm ions through the moving time-dependent field."""
        phases = np.asarray(phases, dtype=float)
        shape = phases.shape
        phase = phases.ravel()
        if phase.size == 0 or np.any(~np.isfinite(phase)):
            raise ValueError("entry phases must be a nonempty finite array")
        if (
            int(steps_per_period) < 16
            or int(steps_per_transit) < 16
            or not math.isfinite(float(max_periods))
            or max_periods <= 0.0
        ):
            raise ValueError("invalid moving-sheath trajectory resolution")
        mass = self.ion_mass_amu * AMU
        v0 = bohm_speed(self.electron_temperature_eV, self.ion_mass_amu)
        vmax = math.sqrt(
            v0 ** 2 + 2.0 * ECHARGE * self.maximum_voltage_v / mass)
        transit_estimate = (
            2.0 * self.maximum_width_m / max(v0 + vmax, 1.0e-30))
        dt = min(
            1.0 / self.maximum_frequency_hz / int(steps_per_period),
            transit_estimate / int(steps_per_transit),
        )
        max_steps = int(math.ceil(max_periods * self.period_s / dt))
        entry_time = (
            np.mod(phase, 2.0 * np.pi)
            * self.period_s
            / (2.0 * np.pi)
        )
        position = np.zeros(phase.size)
        velocity = np.full(phase.size, v0)
        time = entry_time.copy()
        active = np.ones(phase.size, dtype=bool)
        energy = np.full(phase.size, np.nan)

        def acceleration(x, t):
            return ECHARGE * self.electric_field_V_m(x, t) / mass

        for _ in range(max_steps):
            if not active.any():
                break
            index = np.where(active)[0]
            x0 = position[index]
            v_start = velocity[index]
            t0 = time[index]
            a0 = acceleration(x0, t0)
            v_half = v_start + 0.5 * a0 * dt
            x1 = x0 + v_half * dt
            t1 = t0 + dt
            a1 = acceleration(x1, t1)
            v1 = v_half + 0.5 * a1 * dt
            crossed = x1 >= self.maximum_width_m
            if crossed.any():
                fraction = np.clip(
                    (self.maximum_width_m - x0[crossed])
                    / np.maximum(x1[crossed] - x0[crossed], 1.0e-30),
                    0.0,
                    1.0,
                )
                impact_velocity = (
                    v_start[crossed]
                    + fraction * (v1[crossed] - v_start[crossed])
                )
                hit = index[crossed]
                energy[hit] = (
                    0.5 * mass * impact_velocity ** 2 / ECHARGE)
                active[hit] = False
            keep = ~crossed
            position[index[keep]] = x1[keep]
            velocity[index[keep]] = v1[keep]
            time[index[keep]] = t1[keep]
        if active.any():
            raise RuntimeError(
                "ion moving-sheath transit did not finish within max_periods")
        return energy.reshape(shape)

    def voltage_fourier_projection(
        self, harmonic_count: int = 32
    ) -> PeriodicSheathVoltage:
        """Project the nonlinear voltage onto the common Fourier contract."""
        count = int(harmonic_count)
        sample_count = int(self.phase_quadrature_count)
        if count < 1 or 2 * count >= sample_count:
            raise ValueError("invalid voltage Fourier projection order")
        phase = self.phase_grid_rad
        time = phase / (
            2.0 * np.pi * self.current.fundamental_frequency_hz)
        voltage = np.asarray(self.voltage(time))
        harmonic = np.arange(1, count + 1, dtype=int)
        sine = np.asarray([
            2.0 * np.mean(voltage * np.sin(number * phase))
            for number in harmonic
        ])
        cosine = np.asarray([
            2.0 * np.mean(voltage * np.cos(number * phase))
            for number in harmonic
        ])
        return PeriodicSheathVoltage(
            fundamental_frequency_hz=(
                self.current.fundamental_frequency_hz),
            dc_v=float(np.mean(voltage)),
            harmonic_number=harmonic,
            sine_v=sine,
            cosine_v=cosine,
            source=self.source,
            evidence_kind=(
                "validated_reactor_model"
                if self.current.supports_predictive_boundary
                else "assumed"
            ),
        )

    def certification_receipt(self) -> Mapping[str, object]:
        return MappingProxyType({
            "model": "Turner-Chabert arbitrary-current moving RF sheath",
            "model_source": (
                "Turner & Chabert, APL 104, 164102 (2014), "
                "DOI 10.1063/1.4872172, equations 1-19"
            ),
            "current_source": self.current.source,
            "current_evidence_kind": self.current.evidence_kind,
            "charge_excursion_C_m2": self.charge_excursion_C_m2,
            "xi": self.xi,
            "maximum_width_m": self.maximum_width_m,
            "maximum_voltage_v": self.maximum_voltage_v,
            "mean_voltage_v": self.mean_voltage_v,
            "child_current_relative_residual": (
                self.child_current_relative_residual),
            "charge_voltage_relative_residual": (
                self.charge_voltage_relative_residual),
            "moving_electron_front_resolved": True,
            "time_dependent_poisson_field_resolved": True,
            "supports_collisionless_ion_trajectories": True,
            "supports_collisional_ion_trajectories": False,
            "supports_generator_power_inversion": False,
            "supports_external_circuit_prediction": False,
            "supports_feature_depth": False,
            "feature_depth_used": False,
            **dict(self.provenance),
        })
