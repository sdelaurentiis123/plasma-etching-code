"""Species, rate-law, and reaction-network primitives for 0-D reactors.

All densities are in ``m^-3``, time is in seconds, electron temperature is in
eV, and reaction event rates are in ``m^-3 s^-1``.  Electrons are explicit in
reaction stoichiometry so charge conservation cannot be hidden by customary
plasma-chemistry shorthand.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping

import numpy as np

CM3_TO_M3 = 1.0e-6
E_CHARGE_C = 1.602176634e-19
# 2022 CODATA recommended electron mass in kg.
ELECTRON_MASS_KG = 9.1093837139e-31
_SPECIES_ROLES = {
    "electron",
    "neutral",
    "excited_neutral",
    "positive_ion",
    "negative_ion",
}
_EVIDENCE_KINDS = {
    "measured",
    "regressed",
    "semi_empirical",
    "estimated",
    "derived",
    "published_compilation",
}
INCIDENT_ELECTRON_KINETIC_ENERGY_MOMENT = (
    "incident_electron_kinetic_energy"
)
STATIONARY_TARGET_ELASTIC_ENERGY_MOMENT = (
    "stationary_target_elastic_momentum_transfer_energy"
)
KEMANECI_ELASTIC_ENERGY_APPROXIMATION = (
    "kemaneci_maxwellian_elastic_energy_approximation"
)


def _immutable_numeric_mapping(
        values: Mapping[str, float], *, field_name: str,
        positive: bool = False, integral: bool = False):
    output: dict[str, float | int] = {}
    for raw_name, raw_value in dict(values).items():
        name = str(raw_name)
        value = float(raw_value)
        if (
            not name.strip()
            or not np.isfinite(value)
            or (positive and value <= 0.0)
            or (integral and value != int(value))
        ):
            raise ValueError(f"invalid {field_name}")
        output[name] = int(value) if integral else value
    return MappingProxyType(output)


@dataclass(frozen=True)
class Species:
    """One gas-phase species with elemental and charge bookkeeping."""

    name: str
    mass_amu: float
    charge_number: int
    composition: Mapping[str, int]
    role: str
    source: str
    evidence_kind: str = "measured"

    def __post_init__(self):
        composition = _immutable_numeric_mapping(
            self.composition,
            field_name="species composition",
            positive=True,
            integral=True,
        )
        if (
            not str(self.name).strip()
            or not np.isfinite(self.mass_amu)
            or self.mass_amu <= 0.0
            or int(self.charge_number) != self.charge_number
            or self.role not in _SPECIES_ROLES
            or not str(self.source).strip()
            or self.evidence_kind not in _EVIDENCE_KINDS
        ):
            raise ValueError("invalid reactor species")
        if self.role == "electron":
            if self.charge_number != -1 or composition:
                raise ValueError("an electron must have charge -1 and no nuclei")
        elif not composition:
            raise ValueError("a heavy species must declare elemental composition")
        if self.role in {"neutral", "excited_neutral"} and self.charge_number != 0:
            raise ValueError("a neutral species must have zero charge")
        if self.role == "positive_ion" and self.charge_number <= 0:
            raise ValueError("a positive ion must have positive charge")
        if self.role == "negative_ion" and self.charge_number >= 0:
            raise ValueError("a negative ion must have negative charge")
        object.__setattr__(self, "mass_amu", float(self.mass_amu))
        object.__setattr__(self, "charge_number", int(self.charge_number))
        object.__setattr__(self, "composition", composition)


@dataclass(frozen=True)
class RateContext:
    """Thermodynamic state supplied to a rate coefficient."""

    electron_temperature_eV: float
    gas_temperature_K: float | None = None
    scalars: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self):
        if (
            not np.isfinite(self.electron_temperature_eV)
            or self.electron_temperature_eV <= 0.0
            or (
                self.gas_temperature_K is not None
                and (
                    not np.isfinite(self.gas_temperature_K)
                    or self.gas_temperature_K <= 0.0
                )
            )
        ):
            raise ValueError("invalid rate context")
        scalars = _immutable_numeric_mapping(
            self.scalars, field_name="rate-context scalar")
        object.__setattr__(
            self, "electron_temperature_eV",
            float(self.electron_temperature_eV))
        if self.gas_temperature_K is not None:
            object.__setattr__(
                self, "gas_temperature_K", float(self.gas_temperature_K))
        object.__setattr__(self, "scalars", scalars)


@dataclass(frozen=True)
class ConstantRateCoefficient:
    """A constant SI rate coefficient of declared density order."""

    value_si: float
    density_order: float
    source: str
    source_units: str
    evidence_kind: str

    def __post_init__(self):
        if (
            not np.isfinite(self.value_si)
            or self.value_si <= 0.0
            or not np.isfinite(self.density_order)
            or self.density_order <= 0.0
            or not str(self.source).strip()
            or not str(self.source_units).strip()
            or self.evidence_kind not in _EVIDENCE_KINDS
        ):
            raise ValueError("invalid constant rate coefficient")
        object.__setattr__(self, "value_si", float(self.value_si))
        object.__setattr__(self, "density_order", float(self.density_order))

    def coefficient_si(self, context: RateContext) -> float:
        if not isinstance(context, RateContext):
            raise TypeError("rate context is required")
        return self.value_si

    @classmethod
    def from_per_second(
            cls, value: float, *, source: str,
            evidence_kind: str = "measured"):
        return cls(
            value_si=float(value),
            density_order=1.0,
            source=source,
            source_units="s^-1",
            evidence_kind=evidence_kind,
        )

    @classmethod
    def from_cm3_per_s(
            cls, value: float, *, source: str,
            evidence_kind: str = "measured"):
        return cls(
            value_si=float(value) * CM3_TO_M3,
            density_order=2.0,
            source=source,
            source_units="cm^3 s^-1",
            evidence_kind=evidence_kind,
        )


@dataclass(frozen=True)
class ElectronArrheniusRateCoefficient:
    """``A Te^b exp(-E/Te)`` electron-impact coefficient in SI units."""

    prefactor_si: float
    activation_eV: float
    temperature_power: float
    density_order: float
    source: str
    source_units: str
    evidence_kind: str

    def __post_init__(self):
        if (
            not np.isfinite(self.prefactor_si)
            or self.prefactor_si <= 0.0
            or not np.isfinite(self.activation_eV)
            or self.activation_eV < 0.0
            or not np.isfinite(self.temperature_power)
            or not np.isfinite(self.density_order)
            or self.density_order <= 0.0
            or not str(self.source).strip()
            or not str(self.source_units).strip()
            or self.evidence_kind not in _EVIDENCE_KINDS
        ):
            raise ValueError("invalid electron Arrhenius rate coefficient")
        object.__setattr__(self, "prefactor_si", float(self.prefactor_si))
        object.__setattr__(self, "activation_eV", float(self.activation_eV))
        object.__setattr__(
            self, "temperature_power", float(self.temperature_power))
        object.__setattr__(self, "density_order", float(self.density_order))

    def coefficient_si(self, context: RateContext) -> float:
        if not isinstance(context, RateContext):
            raise TypeError("rate context is required")
        temperature = context.electron_temperature_eV
        value = (
            self.prefactor_si
            * temperature ** self.temperature_power
            * np.exp(-self.activation_eV / temperature)
        )
        if not np.isfinite(value) or value <= 0.0:
            raise FloatingPointError("nonpositive or nonfinite rate coefficient")
        return float(value)

    @classmethod
    def from_cm3_per_s(
            cls, prefactor: float, *, activation_eV: float,
            temperature_power: float = 0.0, source: str,
            evidence_kind: str = "regressed"):
        return cls(
            prefactor_si=float(prefactor) * CM3_TO_M3,
            activation_eV=activation_eV,
            temperature_power=temperature_power,
            density_order=2.0,
            source=source,
            source_units="cm^3 s^-1; Te in eV",
            evidence_kind=evidence_kind,
        )


@dataclass(frozen=True)
class ElectronInverseTemperaturePolynomialRateCoefficient:
    """``A exp(sum(c_n / Te**n))`` electron-impact coefficient.

    Several published molecular-plasma fits use inverse-temperature
    polynomials rather than a one-threshold Arrhenius form.  Coefficients are
    indexed from ``n=1`` and retain the source's eV convention for ``Te``.
    """

    prefactor_si: float
    inverse_temperature_coefficients: tuple[float, ...]
    density_order: float
    source: str
    source_units: str
    evidence_kind: str

    def __post_init__(self):
        coefficients = tuple(
            float(value) for value in self.inverse_temperature_coefficients)
        if (
            not np.isfinite(self.prefactor_si)
            or self.prefactor_si <= 0.0
            or not coefficients
            or np.any(~np.isfinite(np.asarray(coefficients)))
            or not np.isfinite(self.density_order)
            or self.density_order <= 0.0
            or not str(self.source).strip()
            or not str(self.source_units).strip()
            or self.evidence_kind not in _EVIDENCE_KINDS
        ):
            raise ValueError(
                "invalid inverse-temperature polynomial rate coefficient")
        object.__setattr__(self, "prefactor_si", float(self.prefactor_si))
        object.__setattr__(
            self, "inverse_temperature_coefficients", coefficients)
        object.__setattr__(self, "density_order", float(self.density_order))

    def coefficient_si(self, context: RateContext) -> float:
        if not isinstance(context, RateContext):
            raise TypeError("rate context is required")
        temperature = context.electron_temperature_eV
        exponent = sum(
            coefficient / temperature ** order
            for order, coefficient in enumerate(
                self.inverse_temperature_coefficients, start=1)
        )
        value = self.prefactor_si * np.exp(exponent)
        if not np.isfinite(value) or value <= 0.0:
            raise FloatingPointError("nonpositive or nonfinite rate coefficient")
        return float(value)

    @classmethod
    def from_cm3_per_s(
            cls, prefactor: float, *,
            inverse_temperature_coefficients: tuple[float, ...],
            source: str, evidence_kind: str = "regressed"):
        return cls(
            prefactor_si=float(prefactor) * CM3_TO_M3,
            inverse_temperature_coefficients=(
                inverse_temperature_coefficients),
            density_order=2.0,
            source=source,
            source_units="cm^3 s^-1; Te in eV",
            evidence_kind=evidence_kind,
        )


@dataclass(frozen=True)
class ElectronLogTemperatureInversePolynomialRateCoefficient:
    """``exp(A) Te**B exp(C/Te + D/Te**2 + ...)`` electron rate.

    This is the rate form used by Kokkoris et al. (J. Phys. D 42,
    055209, 2009), whose printed equation is

    ``ln(k) = A + B ln(Te) + C/Te + D/Te**2 + E/Te**3``.

    The leading coefficient is stored after exponentiating ``A`` so the
    object remains usable for other sources with the same analytic form.
    ``Te`` is in eV and the coefficient is in SI units.
    """

    prefactor_si: float
    temperature_power: float
    inverse_temperature_coefficients: tuple[float, ...]
    density_order: float
    source: str
    source_units: str
    evidence_kind: str

    def __post_init__(self):
        coefficients = tuple(
            float(value) for value in self.inverse_temperature_coefficients)
        if (
            not np.isfinite(self.prefactor_si)
            or self.prefactor_si <= 0.0
            or not np.isfinite(self.temperature_power)
            or not coefficients
            or np.any(~np.isfinite(np.asarray(coefficients)))
            or not np.isfinite(self.density_order)
            or self.density_order <= 0.0
            or not str(self.source).strip()
            or not str(self.source_units).strip()
            or self.evidence_kind not in _EVIDENCE_KINDS
        ):
            raise ValueError(
                "invalid log-temperature inverse-polynomial rate coefficient")
        object.__setattr__(self, "prefactor_si", float(self.prefactor_si))
        object.__setattr__(self, "temperature_power", float(self.temperature_power))
        object.__setattr__(
            self, "inverse_temperature_coefficients", coefficients)
        object.__setattr__(self, "density_order", float(self.density_order))

    def coefficient_si(self, context: RateContext) -> float:
        if not isinstance(context, RateContext):
            raise TypeError("rate context is required")
        temperature = context.electron_temperature_eV
        exponent = sum(
            coefficient / temperature ** order
            for order, coefficient in enumerate(
                self.inverse_temperature_coefficients, start=1)
        )
        value = (
            self.prefactor_si
            * temperature ** self.temperature_power
            * np.exp(exponent)
        )
        if not np.isfinite(value) or value <= 0.0:
            raise FloatingPointError("nonpositive or nonfinite rate coefficient")
        return float(value)

    @classmethod
    def from_log_si(
            cls, log_prefactor: float, *, temperature_power: float,
            inverse_temperature_coefficients: tuple[float, ...],
            source: str, evidence_kind: str = "regressed"):
        return cls(
            prefactor_si=float(np.exp(log_prefactor)),
            temperature_power=temperature_power,
            inverse_temperature_coefficients=(
                inverse_temperature_coefficients),
            density_order=2.0,
            source=source,
            source_units="m^3 s^-1; Te in eV",
            evidence_kind=evidence_kind,
        )


@dataclass(frozen=True)
class ElectronBase10LogPolynomialRateCoefficient:
    """Lennon Eq.-6 base-10-log electron-impact coefficient.

    The evaluated form is

    ``(Te/Tref)^b exp(-E/Te) sum(c_n log10(Te/Tref)^n)``.

    The polynomial coefficients carry the rate units.  Lennon et al.
    explicitly define Eq. 6 with ``log10`` and bound it to
    ``0.1 < Te/Tref < 10``.  Both details are enforced here rather than
    inferred from Lee--Lieberman's abbreviated Table-2 transcription.
    """

    polynomial_coefficients_si: tuple[float, ...]
    reference_temperature_eV: float
    activation_eV: float
    temperature_power: float
    density_order: float
    source: str
    source_units: str
    evidence_kind: str

    def __post_init__(self):
        coefficients = tuple(
            float(value) for value in self.polynomial_coefficients_si)
        values = np.asarray([
            self.reference_temperature_eV,
            self.activation_eV,
            self.temperature_power,
            self.density_order,
        ], dtype=float)
        if (
            not coefficients
            or np.any(~np.isfinite(np.asarray(coefficients)))
            or np.any(~np.isfinite(values))
            or self.reference_temperature_eV <= 0.0
            or self.activation_eV < 0.0
            or self.density_order <= 0.0
            or not str(self.source).strip()
            or not str(self.source_units).strip()
            or self.evidence_kind not in _EVIDENCE_KINDS
        ):
            raise ValueError("invalid electron log-polynomial rate coefficient")
        object.__setattr__(
            self, "polynomial_coefficients_si", coefficients)
        object.__setattr__(
            self, "reference_temperature_eV",
            float(self.reference_temperature_eV))
        object.__setattr__(self, "activation_eV", float(self.activation_eV))
        object.__setattr__(
            self, "temperature_power", float(self.temperature_power))
        object.__setattr__(self, "density_order", float(self.density_order))

    def coefficient_si(self, context: RateContext) -> float:
        if not isinstance(context, RateContext):
            raise TypeError("rate context is required")
        temperature = context.electron_temperature_eV
        ratio = temperature / self.reference_temperature_eV
        if not 0.1 < ratio < 10.0:
            raise ValueError(
                "temperature is outside the Lennon Eq.-6 validity domain")
        log_ratio = np.log10(ratio)
        polynomial = sum(
            coefficient * log_ratio ** order
            for order, coefficient in enumerate(
                self.polynomial_coefficients_si)
        )
        value = (
            ratio ** self.temperature_power
            * np.exp(-self.activation_eV / temperature)
            * polynomial
        )
        if not np.isfinite(value) or value <= 0.0:
            raise FloatingPointError("nonpositive or nonfinite rate coefficient")
        return float(value)

    @classmethod
    def from_cm3_per_s(
            cls, polynomial_coefficients: tuple[float, ...], *,
            reference_temperature_eV: float, activation_eV: float,
            temperature_power: float, source: str,
            evidence_kind: str = "published_compilation"):
        return cls(
            polynomial_coefficients_si=tuple(
                float(value) * CM3_TO_M3
                for value in polynomial_coefficients),
            reference_temperature_eV=reference_temperature_eV,
            activation_eV=activation_eV,
            temperature_power=temperature_power,
            density_order=2.0,
            source=source,
            source_units=(
                "cm^3 s^-1; Te in eV; base-10 log; "
                "0.1 < Te/Tref < 10"
            ),
            evidence_kind=evidence_kind,
        )


@dataclass(frozen=True)
class ElectronMaxwellianCrossSectionRateCoefficient:
    """Maxwellian rate coefficient from a tabulated electron cross section.

    The tabulated cross section is linearly interpolated in energy and
    integrated analytically over each segment:

    ``<sigma v> = sqrt(8 e / (pi m_e)) Te^(-3/2)
                  integral sigma(E) E exp(-E/Te) dE``.

    Both ``E`` and ``Te`` are in eV.  A separately sourced physical threshold
    forces the cross section to zero below threshold even when finite
    experimental energy resolution leaves small sub-threshold table entries.
    The coefficient fails closed when the unmeasured high-energy support
    contains more than ``maximum_kernel_tail_fraction`` of the rate kernel.
    Its separate incident-energy moment applies the stricter
    ``maximum_energy_kernel_tail_fraction`` to the ``E^2 exp(-E/Te)`` kernel.
    """

    electron_energy_eV: tuple[float, ...]
    cross_section_m2: tuple[float, ...]
    threshold_eV: float
    relative_uncertainty: float | None
    source: str
    evidence_kind: str
    maximum_kernel_tail_fraction: float = 1.0e-6
    maximum_energy_kernel_tail_fraction: float = 1.0e-6
    density_order: float = field(init=False, default=2.0)
    source_units: str = field(
        init=False,
        default="tabulated E in eV and cross section in m^2; Maxwellian EEDF",
    )

    def __post_init__(self):
        energies = tuple(float(value) for value in self.electron_energy_eV)
        cross_sections = tuple(float(value) for value in self.cross_section_m2)
        uncertainty = self.relative_uncertainty
        if uncertainty is not None:
            uncertainty = float(uncertainty)
        if (
            len(energies) < 2
            or len(energies) != len(cross_sections)
            or np.any(~np.isfinite(np.asarray(energies)))
            or np.any(~np.isfinite(np.asarray(cross_sections)))
            or np.any(np.diff(np.asarray(energies)) <= 0.0)
            or energies[0] < 0.0
            or np.any(np.asarray(cross_sections) < 0.0)
            or not any(value > 0.0 for value in cross_sections)
            or not np.isfinite(self.threshold_eV)
            or not 0.0 <= self.threshold_eV < energies[-1]
            or (
                uncertainty is not None
                and (
                    not np.isfinite(uncertainty)
                    or not 0.0 <= uncertainty < 1.0
                )
            )
            or not np.isfinite(self.maximum_kernel_tail_fraction)
            or not 0.0 < self.maximum_kernel_tail_fraction < 1.0
            or not np.isfinite(self.maximum_energy_kernel_tail_fraction)
            or not 0.0 < self.maximum_energy_kernel_tail_fraction < 1.0
            or not str(self.source).strip()
            or self.evidence_kind not in _EVIDENCE_KINDS
        ):
            raise ValueError(
                "invalid tabulated Maxwellian cross-section coefficient")
        object.__setattr__(self, "electron_energy_eV", energies)
        object.__setattr__(self, "cross_section_m2", cross_sections)
        object.__setattr__(self, "threshold_eV", float(self.threshold_eV))
        object.__setattr__(self, "relative_uncertainty", uncertainty)
        object.__setattr__(
            self,
            "maximum_kernel_tail_fraction",
            float(self.maximum_kernel_tail_fraction),
        )
        object.__setattr__(
            self,
            "maximum_energy_kernel_tail_fraction",
            float(self.maximum_energy_kernel_tail_fraction),
        )

    def maxwellian_kernel_tail_fraction(self, temperature_eV: float) -> float:
        """Return the constant-cross-section kernel above table support."""
        temperature = float(temperature_eV)
        if not np.isfinite(temperature) or temperature <= 0.0:
            raise ValueError("electron temperature must be finite and positive")
        ratio = self.electron_energy_eV[-1] / temperature
        return float((ratio + 1.0) * np.exp(-ratio))

    def incident_energy_kernel_tail_fraction(
            self, temperature_eV: float) -> float:
        """Return the constant-cross-section energy kernel above support."""
        temperature = float(temperature_eV)
        if not np.isfinite(temperature) or temperature <= 0.0:
            raise ValueError("electron temperature must be finite and positive")
        ratio = self.electron_energy_eV[-1] / temperature
        return float(
            0.5
            * (ratio ** 2 + 2.0 * ratio + 2.0)
            * np.exp(-ratio)
        )

    @staticmethod
    def _first_moment(
            lower_eV: np.ndarray, upper_eV: np.ndarray,
            temperature_eV: float) -> np.ndarray:
        lower = lower_eV / temperature_eV
        upper = upper_eV / temperature_eV
        return temperature_eV ** 2 * (
            (lower + 1.0) * np.exp(-lower)
            - (upper + 1.0) * np.exp(-upper)
        )

    @staticmethod
    def _second_moment(
            lower_eV: np.ndarray, upper_eV: np.ndarray,
            temperature_eV: float) -> np.ndarray:
        lower = lower_eV / temperature_eV
        upper = upper_eV / temperature_eV
        return temperature_eV ** 3 * (
            (lower ** 2 + 2.0 * lower + 2.0) * np.exp(-lower)
            - (upper ** 2 + 2.0 * upper + 2.0) * np.exp(-upper)
        )

    @staticmethod
    def _third_moment(
            lower_eV: np.ndarray, upper_eV: np.ndarray,
            temperature_eV: float) -> np.ndarray:
        lower = lower_eV / temperature_eV
        upper = upper_eV / temperature_eV
        return temperature_eV ** 4 * (
            (
                lower ** 3 + 3.0 * lower ** 2 + 6.0 * lower + 6.0
            ) * np.exp(-lower)
            - (
                upper ** 3 + 3.0 * upper ** 2 + 6.0 * upper + 6.0
            ) * np.exp(-upper)
        )

    def _thresholded_support(self) -> tuple[np.ndarray, np.ndarray]:
        """Return source nodes at/above threshold plus an unsampled zero node.

        If the source tabulates the threshold itself, its printed value is
        retained: the physical rule is zero *below* threshold.  If threshold
        lies between samples or below the first sample, a zero-cross-section
        threshold node is inserted and the unknown interval is represented by
        the explicit linear rise to the first measured point.
        """
        source_energies = np.asarray(self.electron_energy_eV)
        source_cross_sections = np.asarray(self.cross_section_m2)
        first = int(np.searchsorted(
            source_energies, self.threshold_eV, side="left"))
        if (
            first < source_energies.size
            and source_energies[first] == self.threshold_eV
        ):
            return (
                source_energies[first:],
                source_cross_sections[first:],
            )
        return (
            np.concatenate((
                np.array([self.threshold_eV]),
                source_energies[first:],
            )),
            np.concatenate((
                np.array([0.0]),
                source_cross_sections[first:],
            )),
        )

    def _piecewise_moments(
            self, context: RateContext) -> tuple[float, float]:
        if not isinstance(context, RateContext):
            raise TypeError("rate context is required")
        temperature = context.electron_temperature_eV
        energies, cross_sections = self._thresholded_support()
        lower = energies[:-1]
        upper = energies[1:]
        lower_cross_section = cross_sections[:-1]
        slopes = np.diff(cross_sections) / np.diff(energies)
        first_moment = self._first_moment(lower, upper, temperature)
        second_moment = self._second_moment(lower, upper, temperature)
        third_moment = self._third_moment(lower, upper, temperature)
        rate_integral = np.sum(
            lower_cross_section * first_moment
            + slopes * (second_moment - lower * first_moment)
        )
        energy_integral = np.sum(
            lower_cross_section * second_moment
            + slopes * (third_moment - lower * second_moment)
        )
        prefactor = np.sqrt(
            8.0 * E_CHARGE_C / (np.pi * ELECTRON_MASS_KG)
        ) / temperature ** 1.5
        rate = prefactor * rate_integral
        energy = prefactor * energy_integral
        if (
            not np.isfinite(rate)
            or rate <= 0.0
            or not np.isfinite(energy)
            or energy <= 0.0
        ):
            raise FloatingPointError(
                "nonpositive or nonfinite cross-section moment")
        return float(rate), float(energy)

    def coefficient_si(self, context: RateContext) -> float:
        if not isinstance(context, RateContext):
            raise TypeError("rate context is required")
        temperature = context.electron_temperature_eV
        tail_fraction = self.maxwellian_kernel_tail_fraction(temperature)
        if tail_fraction > self.maximum_kernel_tail_fraction:
            raise ValueError(
                "unmeasured cross-section tail exceeds declared Maxwellian "
                "kernel tolerance"
            )

        return self._piecewise_moments(context)[0]

    def incident_energy_moment_eV_m3_s(
            self, context: RateContext) -> float:
        """Return ``<sigma v E>`` for the same evaluated collision table."""
        if not isinstance(context, RateContext):
            raise TypeError("rate context is required")
        tail_fraction = self.incident_energy_kernel_tail_fraction(
            context.electron_temperature_eV)
        if tail_fraction > self.maximum_energy_kernel_tail_fraction:
            raise ValueError(
                "unmeasured cross-section tail exceeds declared incident-"
                "energy Maxwellian kernel tolerance"
            )
        return self._piecewise_moments(context)[1]


@dataclass(frozen=True)
class ElectronTabulatedCrossSectionSupport:
    """Exact Maxwellian moments over only the tabulated energy support.

    Some evaluated collision tables do not span either a physical threshold
    or the full high-energy tail needed by a reactor EEDF.  Extending those
    tables to zero or infinity would add an unmeasured closure.  This class
    instead integrates the piecewise-linear source table only between its
    first and last samples and reports the missing constant-cross-section
    kernel fractions on both sides.

    ``tabulated_rate_coefficient_si`` returns the supported contribution to
    ``<sigma v>``.  ``tabulated_incident_energy_moment_eV_m3_s`` returns the
    supported contribution to ``<sigma v E>``.  The latter is the required
    electron-energy moment for a particle-removing event such as
    dissociative attachment; it must not be replaced by an Arrhenius fit
    exponent.
    """

    electron_energy_eV: tuple[float, ...]
    cross_section_m2: tuple[float, ...]
    relative_uncertainty: float | None
    source: str
    evidence_kind: str
    density_order: float = field(init=False, default=2.0)
    source_units: str = field(
        init=False,
        default=(
            "tabulated E in eV and cross section in m^2; support-only "
            "Maxwellian moments"
        ),
    )

    def __post_init__(self):
        energies = tuple(float(value) for value in self.electron_energy_eV)
        cross_sections = tuple(float(value) for value in self.cross_section_m2)
        uncertainty = self.relative_uncertainty
        if uncertainty is not None:
            uncertainty = float(uncertainty)
        if (
            len(energies) < 2
            or len(energies) != len(cross_sections)
            or np.any(~np.isfinite(np.asarray(energies)))
            or np.any(~np.isfinite(np.asarray(cross_sections)))
            or np.any(np.diff(np.asarray(energies)) <= 0.0)
            or energies[0] < 0.0
            or np.any(np.asarray(cross_sections) < 0.0)
            or not any(value > 0.0 for value in cross_sections)
            or (
                uncertainty is not None
                and (
                    not np.isfinite(uncertainty)
                    or not 0.0 <= uncertainty < 1.0
                )
            )
            or not str(self.source).strip()
            or self.evidence_kind not in _EVIDENCE_KINDS
        ):
            raise ValueError("invalid tabulated cross-section support")
        object.__setattr__(self, "electron_energy_eV", energies)
        object.__setattr__(self, "cross_section_m2", cross_sections)
        object.__setattr__(self, "relative_uncertainty", uncertainty)

    @staticmethod
    def _first_moment(
            lower_eV: np.ndarray, upper_eV: np.ndarray,
            temperature_eV: float) -> np.ndarray:
        return ElectronMaxwellianCrossSectionRateCoefficient._first_moment(
            lower_eV, upper_eV, temperature_eV)

    @staticmethod
    def _second_moment(
            lower_eV: np.ndarray, upper_eV: np.ndarray,
            temperature_eV: float) -> np.ndarray:
        return ElectronMaxwellianCrossSectionRateCoefficient._second_moment(
            lower_eV, upper_eV, temperature_eV)

    @staticmethod
    def _third_moment(
            lower_eV: np.ndarray, upper_eV: np.ndarray,
            temperature_eV: float) -> np.ndarray:
        lower = lower_eV / temperature_eV
        upper = upper_eV / temperature_eV
        return temperature_eV ** 4 * (
            (
                lower ** 3 + 3.0 * lower ** 2 + 6.0 * lower + 6.0
            ) * np.exp(-lower)
            - (
                upper ** 3 + 3.0 * upper ** 2 + 6.0 * upper + 6.0
            ) * np.exp(-upper)
        )

    @staticmethod
    def _prefactor(temperature_eV: float) -> float:
        return float(
            np.sqrt(8.0 * E_CHARGE_C / (np.pi * ELECTRON_MASS_KG))
            / temperature_eV ** 1.5
        )

    def _piecewise_moments(
            self, context: RateContext) -> tuple[float, float]:
        if not isinstance(context, RateContext):
            raise TypeError("rate context is required")
        temperature = context.electron_temperature_eV
        energies = np.asarray(self.electron_energy_eV)
        cross_sections = np.asarray(self.cross_section_m2)
        lower = energies[:-1]
        upper = energies[1:]
        lower_cross_section = cross_sections[:-1]
        slopes = np.diff(cross_sections) / np.diff(energies)
        first = self._first_moment(lower, upper, temperature)
        second = self._second_moment(lower, upper, temperature)
        third = self._third_moment(lower, upper, temperature)
        rate_integral = np.sum(
            lower_cross_section * first
            + slopes * (second - lower * first)
        )
        energy_integral = np.sum(
            lower_cross_section * second
            + slopes * (third - lower * second)
        )
        prefactor = self._prefactor(temperature)
        rate = prefactor * rate_integral
        energy = prefactor * energy_integral
        if (
            not np.isfinite(rate)
            or rate <= 0.0
            or not np.isfinite(energy)
            or energy <= 0.0
        ):
            raise FloatingPointError("nonpositive or nonfinite support moment")
        return float(rate), float(energy)

    def tabulated_rate_coefficient_si(self, context: RateContext) -> float:
        """Return the ``<sigma v>`` contribution on measured support."""
        return self._piecewise_moments(context)[0]

    def tabulated_incident_energy_moment_eV_m3_s(
            self, context: RateContext) -> float:
        """Return the ``<sigma v E>`` contribution on measured support."""
        return self._piecewise_moments(context)[1]

    def rate_kernel_missing_fractions(
            self, temperature_eV: float) -> tuple[float, float]:
        """Return missing fractions of the ``E exp(-E/Te)`` kernel."""
        temperature = float(temperature_eV)
        if not np.isfinite(temperature) or temperature <= 0.0:
            raise ValueError(
                "electron temperature must be finite and positive")
        lower = self.electron_energy_eV[0] / temperature
        upper = self.electron_energy_eV[-1] / temperature
        lower_fraction = 1.0 - (lower + 1.0) * np.exp(-lower)
        upper_fraction = (upper + 1.0) * np.exp(-upper)
        return float(lower_fraction), float(upper_fraction)

    def incident_energy_kernel_missing_fractions(
            self, temperature_eV: float) -> tuple[float, float]:
        """Return missing fractions of the ``E^2 exp(-E/Te)`` kernel."""
        temperature = float(temperature_eV)
        if not np.isfinite(temperature) or temperature <= 0.0:
            raise ValueError(
                "electron temperature must be finite and positive")
        lower = self.electron_energy_eV[0] / temperature
        upper = self.electron_energy_eV[-1] / temperature

        def survival(ratio: float) -> float:
            return float(
                0.5
                * (ratio ** 2 + 2.0 * ratio + 2.0)
                * np.exp(-ratio)
            )

        return float(1.0 - survival(lower)), survival(upper)


@dataclass(frozen=True)
class ElectronTemperatureTabulatedRateCoefficient:
    """Fast bounded electron-temperature rate table.

    Positive coefficients are interpolated linearly in ``log(k)`` against
    ``1 / Te``.  Threshold-dominated electron-impact rates are approximately
    Arrhenius in that coordinate, which gives a compact, positive numerical
    representation without fitting a reactor observable.  Extrapolation is
    forbidden.
    """

    electron_temperature_eV: tuple[float, ...]
    coefficient_m3_s: tuple[float, ...]
    source: str
    evidence_kind: str
    relative_uncertainty: float | None = None
    density_order: float = field(init=False, default=2.0)
    _inverse_temperature: np.ndarray = field(
        init=False, repr=False, compare=False)
    _log_coefficient: np.ndarray = field(
        init=False, repr=False, compare=False)
    source_units: str = field(
        init=False,
        default=(
            "tabulated Te in eV and coefficient in m^3 s^-1; "
            "log(k) linear in 1/Te"
        ),
    )

    def __post_init__(self):
        temperatures = tuple(
            float(value) for value in self.electron_temperature_eV)
        coefficients = tuple(float(value) for value in self.coefficient_m3_s)
        uncertainty = self.relative_uncertainty
        if uncertainty is not None:
            uncertainty = float(uncertainty)
        if (
            len(temperatures) < 2
            or len(temperatures) != len(coefficients)
            or np.any(~np.isfinite(np.asarray(temperatures)))
            or np.any(~np.isfinite(np.asarray(coefficients)))
            or np.any(np.asarray(temperatures) <= 0.0)
            or np.any(np.diff(np.asarray(temperatures)) <= 0.0)
            or np.any(np.asarray(coefficients) <= 0.0)
            or (
                uncertainty is not None
                and (
                    not np.isfinite(uncertainty)
                    or not 0.0 <= uncertainty < 1.0
                )
            )
            or not str(self.source).strip()
            or self.evidence_kind not in _EVIDENCE_KINDS
        ):
            raise ValueError(
                "invalid tabulated electron-temperature coefficient")
        object.__setattr__(self, "electron_temperature_eV", temperatures)
        object.__setattr__(self, "coefficient_m3_s", coefficients)
        object.__setattr__(self, "relative_uncertainty", uncertainty)
        inverse_temperature = 1.0 / np.asarray(temperatures)
        log_coefficient = np.log(np.asarray(coefficients))
        inverse_temperature.setflags(write=False)
        log_coefficient.setflags(write=False)
        object.__setattr__(
            self, "_inverse_temperature", inverse_temperature)
        object.__setattr__(self, "_log_coefficient", log_coefficient)

    def coefficient_si(self, context: RateContext) -> float:
        if not isinstance(context, RateContext):
            raise TypeError("rate context is required")
        temperature = context.electron_temperature_eV
        lower = self.electron_temperature_eV[0]
        upper = self.electron_temperature_eV[-1]
        if not lower <= temperature <= upper:
            raise ValueError(
                "electron temperature is outside the tabulated rate domain")
        value = np.exp(np.interp(
            1.0 / temperature,
            self._inverse_temperature[::-1],
            self._log_coefficient[::-1],
        ))
        if not np.isfinite(value) or value <= 0.0:
            raise FloatingPointError("nonpositive or nonfinite rate coefficient")
        return float(value)


@dataclass(frozen=True)
class ElectronAnalyticRateTerm:
    """One term in a bounded published electron-temperature fit.

    The represented form is

    ``A Te^b exp(sum(c_n / Te^n) + c_s / (Te + d)
                   - (ln(Te) + mu)^2 / (2 width^2))``.

    Zero-valued exponent groups are omitted.  This is deliberately a numeric
    expression object rather than an evaluated string, so a literature table
    can be reproduced without ``eval`` or silently changing its algebra.
    """

    prefactor_si: float
    temperature_power: float = 0.0
    inverse_temperature_coefficients: tuple[float, ...] = ()
    shifted_inverse_coefficient_eV: float = 0.0
    shifted_inverse_offset_eV: float | None = None
    log_temperature_shift: float | None = None
    log_temperature_width: float | None = None

    def __post_init__(self):
        coefficients = tuple(
            float(value) for value in self.inverse_temperature_coefficients)
        values = np.asarray([
            self.prefactor_si,
            self.temperature_power,
            self.shifted_inverse_coefficient_eV,
        ], dtype=float)
        paired_shift = (
            self.shifted_inverse_offset_eV is not None
            or self.shifted_inverse_coefficient_eV != 0.0
        )
        paired_log = (
            self.log_temperature_shift is not None
            or self.log_temperature_width is not None
        )
        if (
            np.any(~np.isfinite(values))
            or self.prefactor_si <= 0.0
            or np.any(~np.isfinite(np.asarray(coefficients)))
            or (
                paired_shift
                and (
                    self.shifted_inverse_offset_eV is None
                    or not np.isfinite(self.shifted_inverse_offset_eV)
                    or self.shifted_inverse_offset_eV <= 0.0
                )
            )
            or (
                paired_log
                and (
                    self.log_temperature_shift is None
                    or self.log_temperature_width is None
                    or not np.isfinite(self.log_temperature_shift)
                    or not np.isfinite(self.log_temperature_width)
                    or self.log_temperature_width <= 0.0
                )
            )
        ):
            raise ValueError("invalid analytic electron-rate term")
        object.__setattr__(self, "prefactor_si", float(self.prefactor_si))
        object.__setattr__(
            self, "temperature_power", float(self.temperature_power))
        object.__setattr__(
            self, "inverse_temperature_coefficients", coefficients)
        object.__setattr__(
            self,
            "shifted_inverse_coefficient_eV",
            float(self.shifted_inverse_coefficient_eV),
        )
        if self.shifted_inverse_offset_eV is not None:
            object.__setattr__(
                self,
                "shifted_inverse_offset_eV",
                float(self.shifted_inverse_offset_eV),
            )
        if self.log_temperature_shift is not None:
            object.__setattr__(
                self,
                "log_temperature_shift",
                float(self.log_temperature_shift),
            )
            object.__setattr__(
                self,
                "log_temperature_width",
                float(self.log_temperature_width),
            )

    def value_si(self, temperature_eV: float) -> float:
        temperature = float(temperature_eV)
        if not np.isfinite(temperature) or temperature <= 0.0:
            raise ValueError("electron temperature must be finite and positive")
        exponent = sum(
            coefficient / temperature ** order
            for order, coefficient in enumerate(
                self.inverse_temperature_coefficients, start=1)
        )
        if self.shifted_inverse_offset_eV is not None:
            exponent += self.shifted_inverse_coefficient_eV / (
                temperature + self.shifted_inverse_offset_eV)
        if self.log_temperature_shift is not None:
            exponent -= (
                (np.log(temperature) + self.log_temperature_shift) ** 2
                / (2.0 * self.log_temperature_width ** 2)
            )
        value = (
            self.prefactor_si
            * temperature ** self.temperature_power
            * np.exp(exponent)
        )
        if not np.isfinite(value) or value <= 0.0:
            raise FloatingPointError("nonpositive analytic electron-rate term")
        return float(value)


@dataclass(frozen=True)
class ElectronCompositeRateCoefficient:
    """A sum of analytic terms on an enforced published ``Te`` domain."""

    terms: tuple[ElectronAnalyticRateTerm, ...]
    minimum_temperature_eV: float
    maximum_temperature_eV: float
    density_order: float
    source: str
    source_units: str
    evidence_kind: str

    def __post_init__(self):
        terms = tuple(self.terms)
        values = np.asarray([
            self.minimum_temperature_eV,
            self.maximum_temperature_eV,
            self.density_order,
        ], dtype=float)
        if (
            not terms
            or any(not isinstance(term, ElectronAnalyticRateTerm)
                   for term in terms)
            or np.any(~np.isfinite(values))
            or self.minimum_temperature_eV <= 0.0
            or self.maximum_temperature_eV <= self.minimum_temperature_eV
            or self.density_order <= 0.0
            or not str(self.source).strip()
            or not str(self.source_units).strip()
            or self.evidence_kind not in _EVIDENCE_KINDS
        ):
            raise ValueError("invalid composite electron-rate coefficient")
        object.__setattr__(self, "terms", terms)
        object.__setattr__(
            self,
            "minimum_temperature_eV",
            float(self.minimum_temperature_eV),
        )
        object.__setattr__(
            self,
            "maximum_temperature_eV",
            float(self.maximum_temperature_eV),
        )
        object.__setattr__(self, "density_order", float(self.density_order))

    def coefficient_si(self, context: RateContext) -> float:
        if not isinstance(context, RateContext):
            raise TypeError("rate context is required")
        temperature = context.electron_temperature_eV
        if not (
            self.minimum_temperature_eV
            <= temperature
            <= self.maximum_temperature_eV
        ):
            raise ValueError(
                "electron temperature is outside the published fit domain")
        value = sum(term.value_si(temperature) for term in self.terms)
        if not np.isfinite(value) or value <= 0.0:
            raise FloatingPointError(
                "nonpositive composite electron-rate coefficient")
        return float(value)


_DETAILED_BALANCE_FORWARD_TYPES = (
    ConstantRateCoefficient,
    ElectronArrheniusRateCoefficient,
    ElectronInverseTemperaturePolynomialRateCoefficient,
    ElectronBase10LogPolynomialRateCoefficient,
    ElectronMaxwellianCrossSectionRateCoefficient,
    ElectronTemperatureTabulatedRateCoefficient,
    ElectronCompositeRateCoefficient,
)


@dataclass(frozen=True)
class ElectronDetailedBalanceRateCoefficient:
    """Maxwellian electron de-excitation coefficient from detailed balance.

    For lower and upper target levels ``l`` and ``u`` separated by
    ``delta_E``, equilibrium requires

    ``k_ul = k_lu (g_l / g_u) exp(delta_E / Te)``.

    The statistical weights and level gap are explicit because omitting either
    can still reproduce a legacy implementation while violating the physical
    Boltzmann population ratio.  The wrapped forward coefficient retains its
    own temperature-domain checks and density order.
    """

    forward_rate_coefficient: object
    excitation_energy_eV: float
    lower_statistical_weight: float
    upper_statistical_weight: float
    source: str
    evidence_kind: str = "derived"

    def __post_init__(self):
        values = np.asarray([
            self.excitation_energy_eV,
            self.lower_statistical_weight,
            self.upper_statistical_weight,
        ], dtype=float)
        if (
            not isinstance(
                self.forward_rate_coefficient,
                _DETAILED_BALANCE_FORWARD_TYPES,
            )
            or np.any(~np.isfinite(values))
            or self.excitation_energy_eV <= 0.0
            or self.lower_statistical_weight <= 0.0
            or self.upper_statistical_weight <= 0.0
            or not str(self.source).strip()
            or self.evidence_kind not in _EVIDENCE_KINDS
        ):
            raise ValueError("invalid detailed-balance rate coefficient")
        for name in (
            "excitation_energy_eV",
            "lower_statistical_weight",
            "upper_statistical_weight",
        ):
            object.__setattr__(self, name, float(getattr(self, name)))

    @property
    def density_order(self) -> float:
        return float(self.forward_rate_coefficient.density_order)

    @property
    def source_units(self) -> str:
        return str(self.forward_rate_coefficient.source_units)

    def coefficient_si(self, context: RateContext) -> float:
        if not isinstance(context, RateContext):
            raise TypeError("rate context is required")
        forward = self.forward_rate_coefficient.coefficient_si(context)
        factor = (
            self.lower_statistical_weight
            / self.upper_statistical_weight
            * np.exp(
                self.excitation_energy_eV
                / context.electron_temperature_eV
            )
        )
        value = forward * factor
        if not np.isfinite(value) or value <= 0.0:
            raise FloatingPointError(
                "nonpositive or nonfinite detailed-balance coefficient")
        return float(value)


@dataclass(frozen=True)
class GasTemperatureArrheniusRateCoefficient:
    """``A (Tg/Tref)^b exp(C/Tg)`` heavy-particle coefficient."""

    prefactor_si: float
    temperature_power: float
    activation_temperature_K: float
    reference_temperature_K: float
    density_order: float
    source: str
    source_units: str
    evidence_kind: str

    def __post_init__(self):
        values = np.asarray([
            self.prefactor_si,
            self.temperature_power,
            self.activation_temperature_K,
            self.reference_temperature_K,
            self.density_order,
        ], dtype=float)
        if (
            np.any(~np.isfinite(values))
            or self.prefactor_si <= 0.0
            or self.reference_temperature_K <= 0.0
            or self.density_order <= 0.0
            or not str(self.source).strip()
            or not str(self.source_units).strip()
            or self.evidence_kind not in _EVIDENCE_KINDS
        ):
            raise ValueError("invalid gas-temperature Arrhenius coefficient")
        for name in (
            "prefactor_si",
            "temperature_power",
            "activation_temperature_K",
            "reference_temperature_K",
            "density_order",
        ):
            object.__setattr__(self, name, float(getattr(self, name)))

    def coefficient_si(self, context: RateContext) -> float:
        if not isinstance(context, RateContext):
            raise TypeError("rate context is required")
        if context.gas_temperature_K is None:
            raise ValueError("gas temperature is required for this rate")
        temperature = context.gas_temperature_K
        value = (
            self.prefactor_si
            * (temperature / self.reference_temperature_K) ** (
                self.temperature_power)
            * np.exp(self.activation_temperature_K / temperature)
        )
        if not np.isfinite(value) or value <= 0.0:
            raise FloatingPointError(
                "nonpositive gas-temperature rate coefficient")
        return float(value)


_RATE_COEFFICIENT_TYPES = (
    ConstantRateCoefficient,
    ElectronArrheniusRateCoefficient,
    ElectronInverseTemperaturePolynomialRateCoefficient,
    ElectronLogTemperatureInversePolynomialRateCoefficient,
    ElectronBase10LogPolynomialRateCoefficient,
    ElectronMaxwellianCrossSectionRateCoefficient,
    ElectronTemperatureTabulatedRateCoefficient,
    ElectronCompositeRateCoefficient,
    ElectronDetailedBalanceRateCoefficient,
    GasTemperatureArrheniusRateCoefficient,
)


@dataclass(frozen=True)
class Reaction:
    """A closed, atom- and charge-conserving reaction event.

    ``kinetic_orders`` are deliberately independent of stoichiometric
    coefficients. This supports reduced wall-return events without pretending
    they are elementary volume collisions.
    """

    name: str
    reactants: Mapping[str, float]
    products: Mapping[str, float]
    kinetic_orders: Mapping[str, float]
    rate_coefficient: (
        ConstantRateCoefficient
        | ElectronArrheniusRateCoefficient
        | ElectronInverseTemperaturePolynomialRateCoefficient
        | ElectronLogTemperatureInversePolynomialRateCoefficient
        | ElectronBase10LogPolynomialRateCoefficient
        | ElectronMaxwellianCrossSectionRateCoefficient
        | ElectronTemperatureTabulatedRateCoefficient
        | ElectronCompositeRateCoefficient
        | ElectronDetailedBalanceRateCoefficient
        | GasTemperatureArrheniusRateCoefficient
    )
    electron_energy_loss_eV: float | None
    source: str
    electron_energy_loss_moment: str | None = None
    domain: str = "volume"
    elastic_target_mass_kg: float | None = None

    def __post_init__(self):
        reactants = _immutable_numeric_mapping(
            self.reactants, field_name="reaction reactant", positive=True)
        products = _immutable_numeric_mapping(
            self.products, field_name="reaction product", positive=True)
        orders = _immutable_numeric_mapping(
            self.kinetic_orders, field_name="reaction kinetic order", positive=True)
        if (
            not str(self.name).strip()
            or not reactants
            or not products
            or not orders
            or self.domain not in {"volume", "closed_wall_return"}
            or not isinstance(self.rate_coefficient, _RATE_COEFFICIENT_TYPES)
            or (
                self.electron_energy_loss_eV is not None
                and not np.isfinite(self.electron_energy_loss_eV)
            )
            or self.electron_energy_loss_moment not in {
                None,
                INCIDENT_ELECTRON_KINETIC_ENERGY_MOMENT,
                STATIONARY_TARGET_ELASTIC_ENERGY_MOMENT,
                KEMANECI_ELASTIC_ENERGY_APPROXIMATION,
            }
            or (
                self.electron_energy_loss_eV is not None
                and self.electron_energy_loss_moment is not None
            )
            or (
                self.elastic_target_mass_kg is not None
                and (
                    not np.isfinite(self.elastic_target_mass_kg)
                    or self.elastic_target_mass_kg <= ELECTRON_MASS_KG
                )
            )
            or not str(self.source).strip()
        ):
            raise ValueError("invalid reactor reaction")
        if not np.isclose(
                sum(orders.values()),
                self.rate_coefficient.density_order,
                rtol=0.0,
                atol=1.0e-14):
            raise ValueError("rate coefficient units disagree with kinetic order")
        object.__setattr__(self, "reactants", reactants)
        object.__setattr__(self, "products", products)
        object.__setattr__(self, "kinetic_orders", orders)
        if self.electron_energy_loss_eV is not None:
            object.__setattr__(
                self, "electron_energy_loss_eV",
                float(self.electron_energy_loss_eV))
        if (
            self.electron_energy_loss_moment
            == INCIDENT_ELECTRON_KINETIC_ENERGY_MOMENT
        ):
            if (
                not isinstance(
                    self.rate_coefficient,
                    ElectronMaxwellianCrossSectionRateCoefficient,
                )
                or reactants.get("e") != 1.0
                or products.get("e", 0.0) != 0.0
                or orders.get("e") != 1.0
                or self.domain != "volume"
            ):
                raise ValueError(
                    "incident-electron energy moment requires a one-electron "
                    "removal reaction driven by the same Maxwellian cross "
                    "section"
                )
            if self.elastic_target_mass_kg is not None:
                raise ValueError(
                    "an incident-electron removal moment has no elastic target"
                )
        elif self.electron_energy_loss_moment in {
            STATIONARY_TARGET_ELASTIC_ENERGY_MOMENT,
            KEMANECI_ELASTIC_ENERGY_APPROXIMATION,
        }:
            heavy_reactants = {
                name: amount for name, amount in reactants.items()
                if name != "e"
            }
            heavy_products = {
                name: amount for name, amount in products.items()
                if name != "e"
            }
            if (
                not isinstance(
                    self.rate_coefficient,
                    ElectronMaxwellianCrossSectionRateCoefficient,
                )
                or reactants.get("e") != 1.0
                or products.get("e") != 1.0
                or orders.get("e") != 1.0
                or heavy_reactants != heavy_products
                or len(heavy_reactants) != 1
                or orders != reactants
                or self.elastic_target_mass_kg is None
                or self.domain != "volume"
            ):
                raise ValueError(
                    "elastic electron-energy transfer requires one unchanged "
                    "heavy target, its mass, and the same Maxwellian "
                    "momentum-transfer cross section"
                )
            object.__setattr__(
                self,
                "elastic_target_mass_kg",
                float(self.elastic_target_mass_kg),
            )
        elif self.elastic_target_mass_kg is not None:
            raise ValueError(
                "an elastic target mass requires an elastic energy moment")

    def event_rate_m3_s(
            self, densities_m3: Mapping[str, float],
            context: RateContext) -> float:
        rate = self.rate_coefficient.coefficient_si(context)
        for name, order in self.kinetic_orders.items():
            if name not in densities_m3:
                raise KeyError(f"missing density for {name}")
            density = float(densities_m3[name])
            if not np.isfinite(density) or density < 0.0:
                raise ValueError("densities must be finite and nonnegative")
            rate *= density ** order
        if not np.isfinite(rate) or rate < 0.0:
            raise FloatingPointError("nonfinite reaction event rate")
        return float(rate)

    def electron_energy_loss_rate_eV_m3_s(
            self, densities_m3: Mapping[str, float],
            context: RateContext) -> float:
        """Return this reaction's signed free-electron energy loss rate."""
        if self.electron_energy_loss_eV is not None:
            return float(
                self.electron_energy_loss_eV
                * self.event_rate_m3_s(densities_m3, context)
            )
        if (
            self.electron_energy_loss_moment
            == INCIDENT_ELECTRON_KINETIC_ENERGY_MOMENT
        ):
            coefficient = self.rate_coefficient
            moment = coefficient.incident_energy_moment_eV_m3_s(context)
            for name, order in self.kinetic_orders.items():
                if name not in densities_m3:
                    raise KeyError(f"missing density for {name}")
                density = float(densities_m3[name])
                if not np.isfinite(density) or density < 0.0:
                    raise ValueError(
                        "densities must be finite and nonnegative")
                moment *= density ** order
            if not np.isfinite(moment) or moment < 0.0:
                raise FloatingPointError("nonfinite electron-energy loss rate")
            return float(moment)
        if (
            self.electron_energy_loss_moment
            == STATIONARY_TARGET_ELASTIC_ENERGY_MOMENT
        ):
            coefficient = self.rate_coefficient
            target_mass = self.elastic_target_mass_kg
            transfer_fraction = (
                2.0 * ELECTRON_MASS_KG * target_mass
                / (ELECTRON_MASS_KG + target_mass) ** 2
            )
            moment = (
                transfer_fraction
                * coefficient.incident_energy_moment_eV_m3_s(context)
            )
            for name, order in self.kinetic_orders.items():
                if name not in densities_m3:
                    raise KeyError(f"missing density for {name}")
                density = float(densities_m3[name])
                if not np.isfinite(density) or density < 0.0:
                    raise ValueError(
                        "densities must be finite and nonnegative")
                moment *= density ** order
            if not np.isfinite(moment) or moment < 0.0:
                raise FloatingPointError("nonfinite electron-energy loss rate")
            return float(moment)
        if (
            self.electron_energy_loss_moment
            == KEMANECI_ELASTIC_ENERGY_APPROXIMATION
        ):
            energy_per_event = (
                3.0
                * context.electron_temperature_eV
                * ELECTRON_MASS_KG
                / self.elastic_target_mass_kg
            )
            return float(
                energy_per_event
                * self.event_rate_m3_s(densities_m3, context)
            )
        raise ValueError(
            f"electron-energy ledger is incomplete for reaction {self.name}")


@dataclass(frozen=True)
class ReactionNetwork:
    """Immutable closed reaction network with exact conservation audits."""

    species: tuple[Species, ...]
    reactions: tuple[Reaction, ...]

    def __post_init__(self):
        species = tuple(self.species)
        reactions = tuple(self.reactions)
        if (
            not species
            or not reactions
            or any(not isinstance(item, Species) for item in species)
            or any(not isinstance(item, Reaction) for item in reactions)
            or len({item.name for item in species}) != len(species)
            or len({item.name for item in reactions}) != len(reactions)
        ):
            raise ValueError("invalid reaction network")
        names = {item.name for item in species}
        for reaction in reactions:
            participants = (
                set(reaction.reactants)
                | set(reaction.products)
                | set(reaction.kinetic_orders)
            )
            missing = participants - names
            if missing:
                raise ValueError(
                    f"reaction {reaction.name} uses unknown species {sorted(missing)}")
        object.__setattr__(self, "species", species)
        object.__setattr__(self, "reactions", reactions)
        self.assert_closed_conservation()

    @property
    def species_names(self) -> tuple[str, ...]:
        return tuple(item.name for item in self.species)

    @property
    def elements(self) -> tuple[str, ...]:
        return tuple(sorted({
            element
            for species in self.species
            for element in species.composition
        }))

    @property
    def stoichiometric_matrix(self) -> np.ndarray:
        index = {name: position for position, name in enumerate(self.species_names)}
        matrix = np.zeros((len(self.species), len(self.reactions)), dtype=float)
        for column, reaction in enumerate(self.reactions):
            for name, coefficient in reaction.products.items():
                matrix[index[name], column] += coefficient
            for name, coefficient in reaction.reactants.items():
                matrix[index[name], column] -= coefficient
        matrix.setflags(write=False)
        return matrix

    @property
    def elemental_matrix(self) -> np.ndarray:
        matrix = np.array([
            [species.composition.get(element, 0) for species in self.species]
            for element in self.elements
        ], dtype=float)
        matrix.setflags(write=False)
        return matrix

    @property
    def charge_vector(self) -> np.ndarray:
        vector = np.array(
            [species.charge_number for species in self.species], dtype=float)
        vector.setflags(write=False)
        return vector

    def reaction_conservation_residuals(self) -> dict[str, dict[str, object]]:
        elemental = self.elemental_matrix @ self.stoichiometric_matrix
        charge = self.charge_vector @ self.stoichiometric_matrix
        return {
            reaction.name: {
                "elements": {
                    element: float(elemental[row, column])
                    for row, element in enumerate(self.elements)
                },
                "charge_number": float(charge[column]),
            }
            for column, reaction in enumerate(self.reactions)
        }

    def assert_closed_conservation(self, *, tolerance: float = 1.0e-14) -> None:
        for name, residual in self.reaction_conservation_residuals().items():
            values = tuple(residual["elements"].values()) + (
                residual["charge_number"],)
            if any(abs(value) > tolerance for value in values):
                raise ValueError(
                    f"reaction {name} does not conserve atoms and charge: {residual}")

    def event_rates_m3_s(
            self, densities_m3: Mapping[str, float],
            context: RateContext) -> np.ndarray:
        extra = set(densities_m3) - set(self.species_names)
        if extra:
            raise KeyError(f"unknown density species {sorted(extra)}")
        rates = np.array([
            reaction.event_rate_m3_s(densities_m3, context)
            for reaction in self.reactions
        ])
        rates.setflags(write=False)
        return rates

    @property
    def has_complete_electron_energy_ledger(self) -> bool:
        """Whether every reaction declares its electron-energy exchange."""
        return all(
            (
                reaction.electron_energy_loss_eV is not None
                or reaction.electron_energy_loss_moment is not None
            )
            for reaction in self.reactions
        )

    def source_vector_m3_s(
            self, densities_m3: Mapping[str, float],
            context: RateContext) -> np.ndarray:
        source = self.stoichiometric_matrix @ self.event_rates_m3_s(
            densities_m3, context)
        source.setflags(write=False)
        return source

    def electron_power_loss_density_W_m3(
            self, densities_m3: Mapping[str, float],
            context: RateContext) -> float:
        extra = set(densities_m3) - set(self.species_names)
        if extra:
            raise KeyError(f"unknown density species {sorted(extra)}")
        if not self.has_complete_electron_energy_ledger:
            missing = [
                reaction.name
                for reaction in self.reactions
                if (
                    reaction.electron_energy_loss_eV is None
                    and reaction.electron_energy_loss_moment is None
                )
            ]
            raise ValueError(
                "electron-energy ledger is incomplete for reactions "
                f"{missing}")
        losses_eV_m3_s = sum(
            reaction.electron_energy_loss_rate_eV_m3_s(
                densities_m3, context)
            for reaction in self.reactions
        )
        return float(E_CHARGE_C * losses_eV_m3_s)

    def source_conservation_report(
            self, densities_m3: Mapping[str, float],
            context: RateContext) -> dict[str, object]:
        source = self.source_vector_m3_s(densities_m3, context)
        element_source = self.elemental_matrix @ source
        charge_source = float(np.dot(self.charge_vector, source))
        scale = max(float(np.max(np.abs(source))), 1.0)
        return {
            "element_source_m3_s": {
                element: float(element_source[index])
                for index, element in enumerate(self.elements)
            },
            "charge_number_source_m3_s": charge_source,
            "normalized_maximum_residual": float(max(
                np.max(np.abs(element_source), initial=0.0),
                abs(charge_source),
            ) / scale),
        }
