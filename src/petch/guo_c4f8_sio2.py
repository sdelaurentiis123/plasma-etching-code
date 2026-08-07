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
from scipy.optimize import least_squares


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


@dataclass(frozen=True)
class GuoTmlTransientEvaluation:
    """Finite-fluence translating-layer evolution.

    Guo Eq. (17) advances atom fractions with reaction yields, whose units are
    atoms per incident ion.  Its dimensionless integration coordinate is
    therefore incident ions per real atom in the translating layer.  A
    feature adapter must supply an independently declared areal layer
    capacity before converting physical time to this coordinate.
    """

    state: GuoTmlState
    incident_ions_per_tml_atom: float
    integrated_removed_movement_atoms_per_tml_atom: float
    integrated_deposited_movement_atoms_per_tml_atom: float
    average_sio2_removal_yield_per_ion: float
    average_deposition_atoms_per_ion: float
    average_net_movement_atoms_per_ion: float
    final_net_movement_atoms_per_ion: float
    final_state_derivative_residual: float
    maximum_atom_ledger_residual_atoms_per_ion: float
    solver_step_count: int
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

    def _rates_and_derivative(self, values, *, movement_regime=None):
        if movement_regime not in {None, "etch", "deposition"}:
            raise ValueError("unknown Guo movement regime")
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
        etch_regime = (
            movement_atoms >= 0.0
            if movement_regime is None else movement_regime == "etch"
        )
        if etch_regime:
            # Net etch: the translating layer moves into SiO2 and receives
            # substrate atoms in the source stoichiometry.  Kwon's worked
            # movement-flux balance uses this same substrate-feed branch.
            movement = {
                "Si": movement_atoms / 3.0,
                "O": 2.0 * movement_atoms / 3.0,
                "C": 0.0,
                "F": 0.0,
            }
        else:
            # Net deposition: Guo Eq. (16) explicitly switches the movement
            # denominator from substrate to deposited film.  Atoms convected
            # out of a perfectly mixed translating layer therefore carry its
            # instantaneous composition.  Continuing to drain SiO2 here
            # makes a fluorocarbon-rich deposition state mathematically
            # impossible and contradicts the Film_or_Sub branch in Eq. (16).
            movement = {
                "Si": movement_atoms * state.si,
                "O": movement_atoms * state.o,
                "C": movement_atoms * state.c,
                "F": movement_atoms * state.f,
            }
        derivative = np.asarray([
            added["Si"] - removed["Si"] + movement["Si"],
            added["O"] - removed["O"] + movement["O"],
            added["C"] - removed["C"] + movement["C"],
            added["F"] - removed["F"] + movement["F"],
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

    def advance_fluence(
        self,
        initial_state: GuoTmlState,
        incident_ions_per_tml_atom: float,
        *,
        relative_tolerance: float = 2.0e-9,
        absolute_tolerance: float = 2.0e-11,
    ) -> GuoTmlTransientEvaluation:
        """Advance Guo Eq. (17) over a finite, dimensionless ion fluence.

        Two accumulator states integrate the gross etch-branch and
        deposition-branch movement separately.  This preserves their
        complementarity even if a transient crosses zero net movement; it
        does not infer a film composition or add a kinetic coefficient.
        """
        if not isinstance(initial_state, GuoTmlState):
            raise TypeError("initial_state must be a GuoTmlState")
        coordinate = float(incident_ions_per_tml_atom)
        if not np.isfinite(coordinate) or coordinate < 0.0:
            raise ValueError(
                "incident_ions_per_tml_atom must be finite and nonnegative")

        initial_values = initial_state.as_array()

        def augmented_derivative(_coordinate, augmented):
            (
                _state,
                _rates,
                added,
                removed,
                movement_atoms,
                state_derivative,
            ) = self._rates_and_derivative(augmented[:5])
            # The pointwise atom ledger is exact by construction:
            # movement = removed - incoming.  Keep gross positive and
            # negative movement in separate accumulators so a branch crossing
            # cannot cancel physical recession and deposition.
            ledger_residual = (
                movement_atoms - sum(removed.values()) + sum(added.values()))
            if abs(ledger_residual) > 5.0e-13:
                raise RuntimeError(
                    "Guo transient atom ledger lost closure: "
                    f"{ledger_residual:.6g} atoms/ion")
            return np.concatenate((
                state_derivative,
                [max(movement_atoms, 0.0), max(-movement_atoms, 0.0)],
            ))

        if coordinate == 0.0:
            (
                _state,
                _rates,
                added,
                removed,
                movement_atoms,
                state_derivative,
            ) = self._rates_and_derivative(initial_values)
            positive = max(movement_atoms, 0.0)
            negative = max(-movement_atoms, 0.0)
            maximum_ledger_residual = abs(
                movement_atoms - sum(removed.values()) + sum(added.values()))
            final_state = initial_state
            solver_step_count = 0
            integrated_positive = 0.0
            integrated_negative = 0.0
            average_positive = positive
            average_negative = negative
        else:
            solution = solve_ivp(
                augmented_derivative,
                (0.0, coordinate),
                np.concatenate((initial_values, [0.0, 0.0])),
                method="BDF",
                rtol=float(relative_tolerance),
                atol=float(absolute_tolerance),
            )
            if not solution.success:
                raise RuntimeError(
                    "Guo translating-layer transient integration failed: "
                    f"{solution.message}")
            final_values = np.asarray(solution.y[:5, -1], dtype=float)
            if np.any(final_values < -2.0e-9):
                raise RuntimeError(
                    "Guo transient left the nonnegative state domain")
            # BDF can return boundary-active species a few ulps below zero.
            # Canonicalize only values already inside GuoTmlState's numerical
            # feasibility tolerance, then restore the exact real-atom
            # normalization required by the feature-state contract.
            final_values = np.maximum(final_values, 0.0)
            real_total = float(np.sum(final_values[:4]))
            if real_total <= 0.0:
                raise RuntimeError("Guo transient erased all real atoms")
            final_values[:4] /= real_total
            final_state = GuoTmlState(*final_values)
            solver_step_count = max(int(solution.t.size) - 1, 0)
            integrated_positive = float(solution.y[5, -1])
            integrated_negative = float(solution.y[6, -1])
            if integrated_positive < -1.0e-10 or integrated_negative < -1.0e-10:
                raise RuntimeError(
                    "Guo transient produced a negative gross movement integral")
            integrated_positive = max(integrated_positive, 0.0)
            integrated_negative = max(integrated_negative, 0.0)
            average_positive = integrated_positive / coordinate
            average_negative = integrated_negative / coordinate
            maximum_ledger_residual = 0.0
            for values in solution.y[:5].T:
                (
                    _state,
                    _rates,
                    added,
                    removed,
                    movement_atoms,
                    _state_derivative,
                ) = self._rates_and_derivative(values)
                maximum_ledger_residual = max(
                    maximum_ledger_residual,
                    abs(
                        movement_atoms
                        - sum(removed.values())
                        + sum(added.values())
                    ),
                )

        (
            _state,
            _rates,
            _added,
            _removed,
            final_movement,
            final_derivative,
        ) = self._rates_and_derivative(final_state.as_array())
        return GuoTmlTransientEvaluation(
            state=final_state,
            incident_ions_per_tml_atom=coordinate,
            integrated_removed_movement_atoms_per_tml_atom=(
                integrated_positive),
            integrated_deposited_movement_atoms_per_tml_atom=(
                integrated_negative),
            average_sio2_removal_yield_per_ion=average_positive / 3.0,
            average_deposition_atoms_per_ion=average_negative,
            average_net_movement_atoms_per_ion=(
                average_positive - average_negative),
            final_net_movement_atoms_per_ion=float(final_movement),
            final_state_derivative_residual=float(
                np.max(np.abs(final_derivative))),
            maximum_atom_ledger_residual_atoms_per_ion=float(
                maximum_ledger_residual),
            solver_step_count=solver_step_count,
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
                "transient_coordinate": (
                    "incident_ions_per_translating_layer_real_atom"
                ),
                "transient_solver": "BDF_finite_fluence",
            }),
        )

    def _steady_evaluation(
        self,
        values,
        *,
        integration_coordinate: float,
        residual_tolerance: float,
        solver_name: str,
    ) -> GuoTmlEvaluation:
        state, rates, added, removed, movement_atoms, derivative = (
            self._rates_and_derivative(values))
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
            elemental_derivative_residual=float(sum(derivative[:4])),
            steady_state_residual=residual,
            integration_coordinate=float(integration_coordinate),
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
                "steady_state_solver": solver_name,
            }),
        )

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
        return self._steady_evaluation(
            solution.y[:, -1],
            integration_coordinate=float(solution.t[-1]),
            residual_tolerance=float(residual_tolerance),
            solver_name="BDF_fluence_integration",
        )

    def solve_steady_state_complementarity_bdf(
        self,
        initial_state: GuoTmlState | None = None,
        *,
        maximum_coordinate: float = 1.0e5,
        relative_tolerance: float = 2.0e-9,
        absolute_tolerance: float = 2.0e-11,
        residual_tolerance: float = 2.0e-8,
    ) -> GuoTmlEvaluation:
        """Integrate each smooth movement branch and enforce its sign.

        Guo Eq. (16) changes the translating-layer movement denominator
        between substrate and deposited film.  Switching that branch inside a
        stiff ODE step can make BDF chatter at zero net movement.  Each branch
        is smooth on its own, so integrate the branch continuous with the
        incoming state first, accept it only if its converged movement sign is
        complementary, and otherwise integrate the other branch.  This is the
        differential-equation analogue of :meth:`solve_steady_state_algebraic`
        and adds no smoothing or fitted transition.
        """
        initial_state = initial_state or GuoTmlState.oxide()
        initial_values = initial_state.as_array()
        initial_movement = self._rates_and_derivative(initial_values)[-2]
        preferred_regime = (
            "etch" if initial_movement >= 0.0 else "deposition")
        regime_order = (
            preferred_regime,
            "deposition" if preferred_regime == "etch" else "etch",
        )
        failures = []
        for regime in regime_order:
            solution = solve_ivp(
                lambda coordinate, values: self._rates_and_derivative(
                    values, movement_regime=regime)[-1],
                (0.0, float(maximum_coordinate)),
                initial_values,
                method="BDF",
                rtol=float(relative_tolerance),
                atol=float(absolute_tolerance),
            )
            if not solution.success:
                failures.append(f"{regime}: {solution.message}")
                continue
            final_values = solution.y[:, -1]
            movement = self._rates_and_derivative(final_values)[-2]
            sign_consistent = (
                movement >= -residual_tolerance
                if regime == "etch"
                else movement <= residual_tolerance
            )
            if not sign_consistent:
                failures.append(
                    f"{regime}: converged movement {movement:.6g} "
                    "violates complementarity")
                continue
            try:
                return self._steady_evaluation(
                    final_values,
                    integration_coordinate=float(solution.t[-1]),
                    residual_tolerance=float(residual_tolerance),
                    solver_name=f"BDF_{regime}_complementarity_integration",
                )
            except RuntimeError as error:
                failures.append(f"{regime}: {error}")
        raise RuntimeError(
            "Guo translating-layer complementarity integration found no "
            "sign-consistent steady branch; " + "; ".join(failures))

    def solve_steady_state_algebraic(
        self,
        initial_state: GuoTmlState | None = None,
        *,
        residual_tolerance: float = 2.0e-8,
        maximum_function_evaluations: int = 3000,
    ) -> GuoTmlEvaluation:
        """Solve the same steady balances as an explicitly constrained root.

        The first four elemental derivatives sum to zero by construction.
        Therefore three elemental balances, the vacancy balance, and the
        real-atom normalization form five independent residuals for the five
        state variables.  Nonnegative bounded least squares is much faster
        than integrating to a very large fluence at every feature face.  It
        introduces no new physics or fitted quantity; tests require agreement
        with the independent BDF fluence integration.
        """
        initial_state = initial_state or GuoTmlState.oxide()
        incident_elements = set()
        incident_atom_weight = {element: 0.0 for element in ELEMENTS}
        for species, value in (
            list(self.composition.neutral_flux_ratio.items())
            + list(self.composition.ion_fraction.items())
        ):
            if value > 0.0:
                atoms = formula_atoms(species)
                incident_elements.update(atoms)
                for element, count in atoms.items():
                    incident_atom_weight[element] += value * count

        candidates = []
        attempted_residuals = []
        initial_values = initial_state.as_array()
        initial_movement = self._rates_and_derivative(initial_values)[-2]
        preferred_regime = (
            "etch" if initial_movement >= 0.0 else "deposition")
        regime_order = (
            preferred_regime,
            "deposition" if preferred_regime == "etch" else "etch",
        )
        for regime in regime_order:
            active_elements = set(incident_elements)
            if regime == "etch":
                active_elements.update({"Si", "O"})
            if not active_elements:
                continue
            active_real = [
                index for index, element in enumerate(ELEMENTS)
                if element in active_elements
            ]

            def unpack(transformed):
                logits = np.concatenate((
                    np.asarray(transformed[:-1], dtype=float),
                    np.zeros(1),
                ))
                logits = logits - np.max(logits)
                active_fraction = np.exp(logits)
                active_fraction /= active_fraction.sum()
                values = np.zeros(5)
                values[active_real] = active_fraction
                values[4] = np.exp(transformed[-1])
                return values

            def independent_residual(transformed):
                values = unpack(transformed)
                derivative = self._rates_and_derivative(
                    values, movement_regime=regime)[-1]
                return np.asarray([
                    *[
                        derivative[index]
                        for index in active_real[:-1]
                    ],
                    derivative[4],
                ])

            def transformed_seed(real, vacancy):
                real = np.maximum(
                    np.asarray(real, dtype=float), 1.0e-12)
                real /= real.sum()
                logits = np.clip(
                    np.log(real[:-1] / real[-1]), -29.0, 29.0)
                return np.concatenate((
                    logits,
                    [np.clip(
                        np.log(max(float(vacancy), 1.0e-12)),
                        -29.0,
                        4.0,
                    )],
                ))

            seed_real = initial_values[active_real]
            seeds = [transformed_seed(seed_real, initial_values[4])]
            incident_seed = np.asarray([
                incident_atom_weight[ELEMENTS[index]]
                for index in active_real
            ])
            if incident_seed.sum() > 0.0:
                seeds.append(transformed_seed(incident_seed, 0.01))
            seeds.append(transformed_seed(
                np.ones(len(active_real)), 0.01))
            if len(active_real) > 1:
                for dominant_index in range(len(active_real)):
                    dominant = np.full(
                        len(active_real), 0.01 / (len(active_real) - 1))
                    dominant[dominant_index] = 0.99
                    seeds.append(transformed_seed(dominant, 0.01))

            for seed_index, transformed_initial in enumerate(seeds):
                solution = least_squares(
                    independent_residual,
                    transformed_initial,
                    bounds=(
                        np.asarray(
                            [-30.0] * (len(active_real) - 1) + [-30.0]),
                        np.asarray(
                            [30.0] * (len(active_real) - 1) + [5.0]),
                    ),
                    # The independently enforced full-balance acceptance gate
                    # is 2e-8 by default.  Driving the optimizer's internal
                    # step/cost/gradient criteria to 1e-14 spent most feature
                    # runtime resolving digits that are discarded by that
                    # physical residual gate.  These tolerances remain two
                    # orders tighter than the accepted balance residual.
                    xtol=min(1.0e-10, residual_tolerance * 0.01),
                    ftol=min(1.0e-10, residual_tolerance * 0.01),
                    gtol=min(1.0e-10, residual_tolerance * 0.01),
                    # Log-ratio composition coordinates enforce exact
                    # nonnegativity and normalization while resolving roots
                    # on the Si=0 or C=F=0 active boundaries.  Log vacancy
                    # removes its two-order scale mismatch.
                    x_scale="jac",
                    max_nfev=int(maximum_function_evaluations),
                )
                values = unpack(solution.x)
                auto = self._rates_and_derivative(values)
                movement_atoms = auto[-2]
                full_residual = float(np.max(np.abs(auto[-1])))
                attempted_residuals.append(full_residual)
                sign_consistent = (
                    movement_atoms >= -residual_tolerance
                    if regime == "etch"
                    else movement_atoms <= residual_tolerance
                )
                if (
                    solution.success
                    and sign_consistent
                    and full_residual <= residual_tolerance
                ):
                    if regime == preferred_regime and seed_index == 0:
                        # The source integrates from the current substrate or
                        # film composition.  A converged root on that same
                        # movement branch is the local continuation and cannot
                        # be displaced by a second, disconnected steady root.
                        return self._steady_evaluation(
                            values,
                            integration_coordinate=0.0,
                            residual_tolerance=float(residual_tolerance),
                            solver_name=(
                                f"bounded_{regime}_"
                                "complementarity_root"
                            ),
                        )
                    continuation_distance = float(
                        np.sum(np.abs(values[:4] - initial_values[:4]))
                        + 0.1 * abs(values[4] - initial_values[4])
                    )
                    candidates.append((
                        continuation_distance,
                        full_residual,
                        regime,
                        values,
                    ))

        if not candidates:
            best = min(attempted_residuals, default=float("inf"))
            raise RuntimeError(
                "Guo translating-layer complementarity solve found no "
                f"sign-consistent branch; best full residual {best:.6g}")
        (
            _continuation_distance,
            residual,
            regime,
            values,
        ) = min(candidates, key=lambda item: (item[0], item[1]))
        if residual > residual_tolerance:
            raise RuntimeError(
                "Guo translating-layer complementarity solve did not "
                f"converge: best full residual {residual:.6g}")
        return self._steady_evaluation(
            values,
            integration_coordinate=0.0,
            residual_tolerance=float(residual_tolerance),
            solver_name=f"bounded_{regime}_complementarity_root",
        )
