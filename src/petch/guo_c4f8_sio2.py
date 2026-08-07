"""Guo--Kwon translating-mixed-layer closure for FC/Ar etching of SiO2.

This module implements the equations printed in Guo's 2009 MIT thesis,
Chapter 4, using the twenty-reaction Table 4.1 deck.  It deliberately remains
separate from :mod:`petch.mixed_layer`: the latter is a Krueger/Huang reduced
projection with a different state topology and different fitted constants.

Evidence boundary
-----------------
The reaction coefficients were regressed by Guo against Yin's C4F8/Ar oxide
yield corpus.  This is therefore an L1, experiment-regressed surface closure,
not an interatomic potential and not independent validation on that corpus.
Its value is that no feature depth was used.  A feature-profile or different
reactor can be a transfer test only when its incident species, flux, energy and
angle boundary is supplied independently.

The source leaves some implementation details implicit.  This implementation
uses Kwon's printed translating-layer elemental balance: incoming and outgoing
atoms are compensated by a stoichiometric SiO2 movement flux so the real-atom
fractions sum to one.  Vacancy is the source's massless dangling-bond state and
is advanced separately.  Adsorption consumes the vacancy appearing in the
printed adsorption rate.  Both choices are exposed in provenance and tested by
elemental ledgers.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import re
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.integrate import solve_ivp


ELEMENTS = ("Si", "O", "C", "F")
VALENCE = MappingProxyType({"Si": 4.0, "O": 2.0, "C": 4.0, "F": 1.0, "V": 1.0})
_FORMULA_TOKEN = re.compile(r"(Si|C|O|F)([0-9]*)")


class GuoSourceLawUnderspecified(ValueError):
    """The requested condition crosses a source ambiguity that needs data."""


def formula_atoms(formula: str) -> Mapping[str, float]:
    """Return elemental stoichiometry for one neutral or ion formula."""
    formula = str(formula).strip()
    if formula in {"Ar", "I", ""}:
        return MappingProxyType({})
    if formula.endswith("+"):
        formula = formula[:-1]
    counts = {element: 0.0 for element in ELEMENTS}
    position = 0
    for match in _FORMULA_TOKEN.finditer(formula):
        if match.start() != position:
            raise ValueError(f"unsupported chemical formula: {formula!r}")
        element, count = match.groups()
        counts[element] += float(count or "1")
        position = match.end()
    if position != len(formula):
        raise ValueError(f"unsupported chemical formula: {formula!r}")
    return MappingProxyType({
        element: count for element, count in counts.items() if count})


@dataclass(frozen=True)
class GuoIonQuadrature:
    """Flux-normalized incident-ion energy/cosine quadrature."""

    energy_eV: np.ndarray
    cosine_incidence: np.ndarray
    weight: np.ndarray

    def __post_init__(self):
        energy = np.asarray(self.energy_eV, dtype=float).reshape(-1).copy()
        cosine = np.asarray(
            self.cosine_incidence, dtype=float).reshape(-1).copy()
        weight = np.asarray(self.weight, dtype=float).reshape(-1).copy()
        if (
            not energy.size
            or energy.shape != cosine.shape
            or energy.shape != weight.shape
            or np.any(~np.isfinite(energy))
            or np.any(energy < 0.0)
            or np.any(~np.isfinite(cosine))
            or np.any((cosine < 0.0) | (cosine > 1.0))
            or np.any(~np.isfinite(weight))
            or np.any(weight < 0.0)
            or not np.isfinite(weight.sum())
            or weight.sum() <= 0.0
        ):
            raise ValueError("invalid Guo ion quadrature")
        weight /= weight.sum()
        for array in (energy, cosine, weight):
            array.setflags(write=False)
        object.__setattr__(self, "energy_eV", energy)
        object.__setattr__(self, "cosine_incidence", cosine)
        object.__setattr__(self, "weight", weight)

    @classmethod
    def monoenergetic(cls, energy_eV: float, cosine_incidence: float = 1.0):
        return cls(
            np.asarray([energy_eV]),
            np.asarray([cosine_incidence]),
            np.asarray([1.0]),
        )


@dataclass(frozen=True)
class GuoIncidentComposition:
    """Incident chemistry normalized to total positive-ion flux.

    ``neutral_flux_ratio`` values are molecular flux / total positive-ion
    flux.  ``ion_fraction`` values are fractions of that positive-ion flux;
    their sum may be below one, with the remainder explicitly interpreted as
    non-incorporating inert ions (normally Ar+).
    """

    neutral_flux_ratio: Mapping[str, float]
    ion_fraction: Mapping[str, float]

    def __post_init__(self):
        neutrals = {
            str(name): float(value)
            for name, value in self.neutral_flux_ratio.items()
        }
        ions = {
            str(name).removesuffix("+"): float(value)
            for name, value in self.ion_fraction.items()
        }
        for values in (neutrals, ions):
            if any(
                not name
                or not np.isfinite(value)
                or value < 0.0
                for name, value in values.items()
            ):
                raise ValueError("invalid Guo incident composition")
            for name in values:
                formula_atoms(name)
        if sum(ions.values()) > 1.0 + 1.0e-12:
            raise ValueError("Guo ion fractions exceed total positive-ion flux")
        object.__setattr__(
            self, "neutral_flux_ratio", MappingProxyType(neutrals))
        object.__setattr__(self, "ion_fraction", MappingProxyType(ions))

    @property
    def inert_ion_fraction(self) -> float:
        return max(1.0 - sum(self.ion_fraction.values()), 0.0)


@dataclass(frozen=True)
class GuoTmlState:
    """Real-atom fractions plus massless vacancy density per real atom."""

    si: float
    o: float
    c: float
    f: float
    vacancy: float

    def __post_init__(self):
        values = np.asarray(
            [self.si, self.o, self.c, self.f, self.vacancy], dtype=float)
        if (
            np.any(~np.isfinite(values))
            or np.any(values < -2.0e-9)
            or self.vacancy < -2.0e-9
            or not np.isclose(values[:4].sum(), 1.0, atol=2.0e-8)
        ):
            raise ValueError("invalid Guo translating-layer state")

    @classmethod
    def oxide(cls):
        return cls(1.0 / 3.0, 2.0 / 3.0, 0.0, 0.0, 0.0)

    def as_array(self) -> np.ndarray:
        return np.asarray([self.si, self.o, self.c, self.f, self.vacancy])

    def fractions(self) -> Mapping[str, float]:
        return MappingProxyType({
            "Si": self.si,
            "O": self.o,
            "C": self.c,
            "F": self.f,
            "V": self.vacancy,
        })


@dataclass(frozen=True)
class GuoTmlEvaluation:
    state: GuoTmlState
    sio2_yield_per_ion: float
    movement_atoms_per_ion: float
    reaction_yields: Mapping[str, float]
    incoming_atoms_per_ion: Mapping[str, float]
    removed_atoms_per_ion: Mapping[str, float]
    elemental_derivative_residual: float
    steady_state_residual: float
    integration_coordinate: float
    source_extrapolation: Mapping[str, object]


def nearest_neighbor_probability(
    state: GuoTmlState, first: str, second: str
) -> float:
    """Guo Eq. (1), including its identical-species factor of two."""
    if first not in VALENCE or second not in VALENCE:
        raise ValueError("unknown Guo bond species")
    fractions = state.fractions()
    denominator = sum(VALENCE[name] * fractions[name] for name in VALENCE)
    if denominator <= 0.0:
        return 0.0
    identical = 1.0 if first == second else 0.0
    return (
        VALENCE[first] * fractions[first]
        * VALENCE[second] * fractions[second]
        / ((1.0 + identical) * denominator)
    )


def _theta_deg(cosine):
    return np.rad2deg(np.arccos(np.clip(cosine, 0.0, 1.0)))


def ion_incorporation_angular(cosine):
    return np.clip(np.asarray(cosine, dtype=float), 0.0, 1.0)


def physical_sputtering_angular(cosine):
    """Guo Table 4.2 polynomial with its evident typesetting error repaired.

    The table prints both the -421.98 and +95.31 terms as ``cos^2``.  That
    literal reading is negative from about 25--90 degrees.  A sixth-degree
    polynomial containing the otherwise complete descending degree sequence,
    together with the source's plotted positive sputter-like curve, implies
    that the +95.31 term is linear in cosine.  This is a declared inference,
    not a silently corrected transcription.  We retain every printed
    coefficient, make that single exponent repair, and clip the tiny negative
    grazing tail to zero.
    """
    cosine = np.clip(np.asarray(cosine, dtype=float), 0.0, 1.0)
    value = (
        -141.29 * cosine ** 6
        + 641.11 * cosine ** 5
        - 1111.3 * cosine ** 4
        + 944.63 * cosine ** 3
        - 421.98 * cosine ** 2
        + 95.31 * cosine
        - 5.46
    )
    return np.maximum(value, 0.0)


def physical_sputtering_angular_literal(cosine):
    """Return the impossible literal Table 4.2 polynomial for source auditing.

    This function is intentionally not used by the mechanism.  It preserves
    the duplicated ``cos^2`` exponent so tests and reports can demonstrate why
    the printed expression cannot be a physical sputtering multiplier.
    """
    cosine = np.clip(np.asarray(cosine, dtype=float), 0.0, 1.0)
    return (
        -141.29 * cosine ** 6
        + 641.11 * cosine ** 5
        - 1111.3 * cosine ** 4
        + 944.63 * cosine ** 3
        - 421.98 * cosine ** 2
        + 95.31 * cosine ** 2
        - 5.46
    )


def ion_enhanced_angular(cosine):
    """Guo Table 4.2 oxide ion-enhanced/vacancy angular law."""
    theta = _theta_deg(np.asarray(cosine, dtype=float))
    high = (
        (110.0 - theta) / 85.0
        - ((theta - 25.0) * (theta - 90.0)) / 5000.0
    )
    return np.maximum(np.where(theta < 25.0, 1.0, high), 0.0)


class GuoC4F8ArSiO2Mechanism:
    """Twenty-reaction Guo translating-layer mechanism at one surface point."""

    # Yin identifies DC 350 V as approximately 370 eV at the sample.  That is
    # the highest energy in the yield corpus used by Guo's Figure 4.2 fit.
    source_energy_max_eV = 370.0
    maximum_printed_threshold_eV = 72.0

    def __init__(
        self,
        composition: GuoIncidentComposition,
        ions: GuoIonQuadrature,
    ):
        if np.min(ions.energy_eV) < self.maximum_printed_threshold_eV:
            raise GuoSourceLawUnderspecified(
                "Table 4.1 does not assign an incorporation threshold to the "
                "generic ion I+; conditions below the largest printed "
                "elemental sputter threshold (72 eV) need a species-resolved "
                "threshold closure"
            )
        self.composition = composition
        self.ions = ions
        self._moments = MappingProxyType(self._build_ion_moments())

    def _average(self, value) -> float:
        return float(np.dot(self.ions.weight, np.asarray(value, dtype=float)))

    def _sqrt_moment(self, coefficient, threshold, angular) -> float:
        energetic = np.maximum(
            np.sqrt(self.ions.energy_eV) - math.sqrt(threshold), 0.0)
        return self._average(coefficient * energetic * angular)

    def _build_ion_moments(self) -> dict[str, float]:
        physical = physical_sputtering_angular(
            self.ions.cosine_incidence)
        enhanced = ion_enhanced_angular(self.ions.cosine_incidence)
        moments = {
            "incorporation": self._average(
                ion_incorporation_angular(self.ions.cosine_incidence)),
            "annihilation": self._average(enhanced),
        }
        for name, coefficient, threshold in (
            ("sputter_Si", 0.042, 44.0),
            ("sputter_O", 0.018, 60.0),
            ("sputter_C", 0.009, 72.0),
            ("sputter_F", 0.023, 55.0),
        ):
            moments[name] = self._sqrt_moment(
                coefficient, threshold, physical)
        for name, coefficient, threshold in (
            ("SiF2", 6.75, 22.0),
            ("F2", 0.0, 20.0),
            ("O2", 0.22, 10.0),
            ("SiO", 0.007, 20.0),
            ("CO", 0.24, 35.0),
            ("CO2", 0.95, 14.0),
            ("CF2", 2.0, 0.0),
            ("vacancy_creation", 0.14, 40.0),
            ("densification", 1.66, 20.0),
        ):
            moments[name] = self._sqrt_moment(
                coefficient, threshold, enhanced)
        return moments

    def _rates_and_derivative(self, values):
        real = np.maximum(np.asarray(values[:4], dtype=float), 0.0)
        total = real.sum()
        if total <= 0.0:
            real = np.asarray([1.0 / 3.0, 2.0 / 3.0, 0.0, 0.0])
        else:
            real /= total
        state = GuoTmlState(
            float(real[0]), float(real[1]), float(real[2]), float(real[3]),
            max(float(values[4]), 0.0),
        )
        j = lambda a, b: nearest_neighbor_probability(state, a, b)
        rates = {
            "r5_sputter_Si": self._moments["sputter_Si"] * state.si,
            "r6_sputter_O": self._moments["sputter_O"] * state.o,
            "r7_sputter_C": self._moments["sputter_C"] * state.c,
            "r8_sputter_F": self._moments["sputter_F"] * state.f,
            "r9_SiF2": self._moments["SiF2"] * j("Si", "F") ** 2,
            "r10_F2": self._moments["F2"] * j("F", "F"),
            "r11_O2": self._moments["O2"] * j("O", "O"),
            "r12_SiO": self._moments["SiO"] * j("Si", "O"),
            "r13_CO": self._moments["CO"] * j("C", "O"),
            "r14_CO2": self._moments["CO2"] * j("C", "O") ** 2,
            "r15_CF2": self._moments["CF2"] * j("C", "F") ** 2,
            "r16_vacancy_creation": self._moments["vacancy_creation"],
            "r17_densification": (
                self._moments["densification"] * state.vacancy),
            "r18_annihilation": (
                10000.0 * self._moments["annihilation"]
                * j("C", "V") ** 2
            ),
        }

        added = {element: 0.0 for element in ELEMENTS}
        removed = {element: 0.0 for element in ELEMENTS}
        ion_scale = self._moments["incorporation"]
        for species, fraction in self.composition.ion_fraction.items():
            for element, count in formula_atoms(species).items():
                added[element] += ion_scale * fraction * count

        atomic_f_ratio = self.composition.neutral_flux_ratio.get("F", 0.0)
        rates["r2_F_adsorption"] = 20.0 * j("Si", "V") * atomic_f_ratio
        added["F"] += rates["r2_F_adsorption"]

        generic_adsorption = 0.0
        for species, flux_ratio in self.composition.neutral_flux_ratio.items():
            if species == "F":
                continue
            rate = (
                3.5 * j("C", "V") + 1.8 * j("O", "V")
            ) * flux_ratio
            rates[f"r3_4_adsorption_{species}"] = rate
            generic_adsorption += rate
            for element, count in formula_atoms(species).items():
                added[element] += rate * count

        rates["r19_CF4_recombination"] = (
            0.0
            * state.f
            * self.composition.neutral_flux_ratio.get("CF3", 0.0)
        )
        rates["r20_spontaneous_SiF4"] = (
            2.99e-5 * state.si * atomic_f_ratio)

        removed["Si"] += (
            rates["r5_sputter_Si"]
            + rates["r9_SiF2"]
            + rates["r12_SiO"]
            + rates["r20_spontaneous_SiF4"]
        )
        removed["O"] += (
            rates["r6_sputter_O"]
            + 2.0 * rates["r11_O2"]
            + rates["r12_SiO"]
            + rates["r13_CO"]
            + 2.0 * rates["r14_CO2"]
        )
        removed["C"] += (
            rates["r7_sputter_C"]
            + rates["r13_CO"]
            + rates["r14_CO2"]
            + rates["r15_CF2"]
            + rates["r19_CF4_recombination"]
        )
        removed["F"] += (
            rates["r8_sputter_F"]
            + 2.0 * rates["r9_SiF2"]
            + 2.0 * rates["r10_F2"]
            + 2.0 * rates["r15_CF2"]
            + rates["r19_CF4_recombination"]
            + 4.0 * rates["r20_spontaneous_SiF4"]
        )

        added_total = sum(added.values())
        removed_total = sum(removed.values())
        movement_atoms = removed_total - added_total
        movement_sio2 = movement_atoms / 3.0
        derivative = np.asarray([
            added["Si"] - removed["Si"] + movement_sio2,
            added["O"] - removed["O"] + 2.0 * movement_sio2,
            added["C"] - removed["C"],
            added["F"] - removed["F"],
            (
                rates["r16_vacancy_creation"]
                - rates["r17_densification"]
                - 2.0 * rates["r18_annihilation"]
                - rates["r2_F_adsorption"]
                - generic_adsorption
            ),
        ])
        return state, rates, added, removed, movement_atoms, derivative

    def derivative(self, coordinate, values):
        del coordinate
        return self._rates_and_derivative(values)[-1]

    def solve_steady_state(
        self,
        initial_state: GuoTmlState | None = None,
        *,
        maximum_coordinate: float = 1.0e5,
        relative_tolerance: float = 2.0e-9,
        absolute_tolerance: float = 2.0e-11,
        residual_tolerance: float = 2.0e-8,
    ) -> GuoTmlEvaluation:
        initial_state = initial_state or GuoTmlState.oxide()
        solution = solve_ivp(
            self.derivative,
            (0.0, float(maximum_coordinate)),
            initial_state.as_array(),
            method="BDF",
            rtol=float(relative_tolerance),
            atol=float(absolute_tolerance),
        )
        if not solution.success:
            raise RuntimeError(
                f"Guo translating-layer integration failed: {solution.message}")
        state, rates, added, removed, movement_atoms, derivative = (
            self._rates_and_derivative(solution.y[:, -1]))
        residual = float(np.max(np.abs(derivative)))
        if residual > residual_tolerance:
            raise RuntimeError(
                "Guo translating-layer steady state did not converge: "
                f"max residual {residual:.6g}")
        return GuoTmlEvaluation(
            state=state,
            sio2_yield_per_ion=movement_atoms / 3.0,
            movement_atoms_per_ion=movement_atoms,
            reaction_yields=MappingProxyType(dict(rates)),
            incoming_atoms_per_ion=MappingProxyType(dict(added)),
            removed_atoms_per_ion=MappingProxyType(dict(removed)),
            elemental_derivative_residual=float(
                sum(derivative[:4])),
            steady_state_residual=residual,
            integration_coordinate=float(solution.t[-1]),
            source_extrapolation=MappingProxyType({
                "source_fit_energy_max_eV": self.source_energy_max_eV,
                "quadrature_max_energy_eV": float(
                    np.max(self.ions.energy_eV)),
                "beyond_source_fit_energy": bool(
                    np.max(self.ions.energy_eV)
                    > self.source_energy_max_eV
                ),
                "physical_angular_polynomial_interpretation": (
                    "declared_degree_sequence_repair_95p31_cosine_power_1"
                ),
                "surface_evidence_ceiling": "L1_yield_regressed",
            }),
        )
