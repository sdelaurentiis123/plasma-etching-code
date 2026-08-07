"""Feature-engine adapter for the Guo/Kwon C4F8/Ar--SiO2 closure.

The adapter preserves the source's distinction between a fitted adsorption
coefficient and a collision probability.  Guo's adsorption yields are

``20 J(Si,V) * Gamma_F / Gamma_i``

and

``[3.5 J(C,V) + 1.8 J(O,V)] * Gamma_N / Gamma_i``.

Consequently the bracketed, state-dependent factors are the per-collision
loss probabilities required by neutral radiosity.  The raw coefficients
20, 3.5, and 1.8 are never passed to transport as probabilities.

At each exposed face the exact energetic event measure and local incident
neutral fluxes define one Guo translating-layer steady state.  The bounded
algebraic root is mathematically identical to the independently tested BDF
fluence integration and avoids an arbitrary mixing-layer thickness in a
quasi-steady feature calculation.

This is a transfer-audit mechanism, not a declaration that Krueger's missing
ion composition or C4F6 boundary has been identified.  Its validity record
remains false outside Guo's measured/regressed board even when explicitly
allowed to execute for sensitivity scoring.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .guo_c4f8_sio2 import (
    GuoC4F8ArSiO2Mechanism,
    GuoIncidentComposition,
    GuoIonQuadrature,
    GuoSourceLawUnderspecified,
    GuoTmlState,
    formula_atoms,
)
from .surface_exchange import (
    SurfaceMaterialExchange,
    unresolved_surface_exchange,
)
from .surface_kinetics import (
    EnergeticFlux,
    FaceResolvedEnergeticFlux,
    MechanismValidity,
    SurfaceFluxes,
)


_REAL_FIELDS = ("si", "o", "c", "f")
_STATE_FIELDS = _REAL_FIELDS + ("vacancy", "removed_sio2_formula_units_m2")
_SOURCE_NEUTRALS = frozenset({"C", "F", "O", "CF", "CF2", "CF3"})


@dataclass(frozen=True)
class GuoC4F8ArSiO2FeatureState:
    """Local intensive TML composition plus cumulative substrate removal."""

    si: np.ndarray | float
    o: np.ndarray | float
    c: np.ndarray | float
    f: np.ndarray | float
    vacancy: np.ndarray | float
    removed_sio2_formula_units_m2: np.ndarray | float = 0.0

    def __post_init__(self):
        supplied = np.broadcast_arrays(*[
            np.asarray(getattr(self, name), dtype=float)
            for name in _STATE_FIELDS
        ])
        values = [np.array(item, copy=True) for item in supplied]
        if any(np.any(~np.isfinite(item)) or np.any(item < 0.0)
               for item in values):
            raise ValueError("invalid Guo feature surface state")
        real_sum = sum(values[:4])
        if not np.allclose(real_sum, 1.0, rtol=0.0, atol=2.0e-8):
            raise ValueError("Guo feature real-atom fractions must sum to one")
        for name, value in zip(_STATE_FIELDS, values):
            value.setflags(write=False)
            object.__setattr__(self, name, value)

    @classmethod
    def oxide(cls, shape=()):
        one = np.ones(shape)
        zero = np.zeros(shape)
        return cls(one / 3.0, 2.0 * one / 3.0, zero, zero, zero, zero)

    @classmethod
    def bare(cls, shape=()):
        return cls.oxide(shape)

    def conservative_surface_fields(self):
        return {name: getattr(self, name) for name in _STATE_FIELDS}

    def conservative_surface_upper_bounds(self):
        return {
            "si": 1.0,
            "o": 1.0,
            "c": 1.0,
            "f": 1.0,
            "vacancy": None,
            "removed_sio2_formula_units_m2": None,
        }

    def surface_field_remap_modes(self):
        return {
            "si": "intensive",
            "o": "intensive",
            "c": "intensive",
            "f": "intensive",
            "vacancy": "intensive",
            "removed_sio2_formula_units_m2": "conservative",
        }

    def with_conservative_surface_fields(self, fields):
        fields = dict(fields)
        if set(fields) != set(_STATE_FIELDS):
            raise ValueError("Guo feature remap fields do not match")
        real = [np.asarray(fields[name], dtype=float) for name in _REAL_FIELDS]
        total = sum(real)
        if np.any(total <= 0.0):
            raise ValueError("Guo feature remap erased all real atoms")
        normalized = [value / total for value in real]
        return type(self)(
            *normalized,
            fields["vacancy"],
            fields["removed_sio2_formula_units_m2"],
        )


@dataclass(frozen=True)
class GuoC4F8ArSiO2FeatureStepResult:
    state: GuoC4F8ArSiO2FeatureState
    etch_velocity_m_s: np.ndarray
    normal_growth_velocity_m_s: np.ndarray
    removed_sio2_formula_units_m2: np.ndarray
    deposited_mixed_layer_atoms_m2: np.ndarray
    sio2_yield_per_ion: np.ndarray
    net_movement_atoms_per_ion: np.ndarray
    steady_state_residual: np.ndarray
    atom_ledger_residual_atoms_per_ion: np.ndarray
    bdf_fallback_face_count: int
    maximum_neutral_reaction_probability: float
    material_exchange: SurfaceMaterialExchange
    validity: MechanismValidity
    product_populations: tuple = ()


class GuoC4F8ArSiO2FeatureMechanism:
    """Quasi-steady local Guo surface law behind the common 3-D engine."""

    quasi_steady_surface_state = True

    def __init__(
        self,
        *,
        neutral_species=(
            "C3F4", "C2F3", "CF", "CF2", "CF3", "O",
        ),
        ion_species_mapping: Mapping[str, str | None] | None = None,
        bulk_sio2_formula_density_m3: float = 2.2e28,
        deposited_film_atom_density_m3: float = 7.5e28,
        allow_out_of_board_transfer_audit: bool = False,
    ):
        neutral = tuple(str(name) for name in neutral_species)
        if not neutral or len(set(neutral)) != len(neutral):
            raise ValueError("Guo feature neutral species must be unique")
        for name in neutral:
            formula_atoms(name)
        mapping = dict(
            {"ions": None}
            if ion_species_mapping is None else ion_species_mapping
        )
        if not mapping:
            raise ValueError("Guo feature requires an energetic species mapping")
        normalized_mapping = {}
        for population, formula in mapping.items():
            population = str(population)
            if not population:
                raise ValueError("empty Guo energetic population name")
            if formula is not None:
                formula = str(formula).removesuffix("+")
                formula_atoms(formula)
            normalized_mapping[population] = formula
        density = float(bulk_sio2_formula_density_m3)
        film_density = float(deposited_film_atom_density_m3)
        if (not np.isfinite(density) or density <= 0.0
                or not np.isfinite(film_density) or film_density <= 0.0):
            raise ValueError("invalid SiO2 or deposited-film density")

        self.neutral_species = neutral
        self.ion_species_mapping = MappingProxyType(normalized_mapping)
        self.bulk_sio2_formula_density_m3 = density
        self.deposited_film_atom_density_m3 = film_density
        self.allow_out_of_board_transfer_audit = bool(
            allow_out_of_board_transfer_audit)
        self.provenance = MappingProxyType({
            "model": "guo-kwon-translating-mixed-layer-feature-adapter-v1",
            "source_surface_model": (
                "Guo MIT PhD thesis 2009 Table 4.1/4.2 and Kwon MIT ScD "
                "thesis 2004 translating-layer balances"
            ),
            "neutral_collision_probability": {
                "F": "20*J(Si,V)",
                "generic_N": "3.5*J(C,V)+1.8*J(O,V)",
                "policy": (
                    "state-dependent products are checked in [0,1]; raw "
                    "coefficients are never treated as probabilities"
                ),
            },
            "steady_state_solver": (
                "bounded algebraic root, independently pinned against BDF "
                "fluence integration"
            ),
            "neutral_species": list(neutral),
            "ion_species_mapping": {
                name: formula for name, formula in normalized_mapping.items()
            },
            "bulk_sio2_formula_density_m3": density,
            "deposited_film_atom_density_m3": film_density,
            "allow_out_of_board_transfer_audit":
                self.allow_out_of_board_transfer_audit,
            "calibration": {
                "feature_depth_used": False,
                "adjustable_surface_parameters": [],
            },
            "claim_status": (
                "source-law transfer audit; not an identified C4F6 reactor "
                "boundary and not atomic-level first-principles accuracy"
            ),
        })

    @classmethod
    def krueger_2024_transfer_audit(cls):
        return cls(allow_out_of_board_transfer_audit=True)

    def initial_state(self, shape=()):
        return GuoC4F8ArSiO2FeatureState.oxide(shape)

    @staticmethod
    def _bond_probability_arrays(state, first):
        valence = {"Si": 4.0, "O": 2.0, "C": 4.0}
        fraction = {
            "Si": np.asarray(state.si, dtype=float),
            "O": np.asarray(state.o, dtype=float),
            "C": np.asarray(state.c, dtype=float),
        }[first]
        denominator = (
            4.0 * state.si + 2.0 * state.o + 4.0 * state.c
            + state.f + state.vacancy
        )
        return (
            valence[first] * fraction * state.vacancy
            / np.maximum(denominator, np.finfo(float).tiny)
        )

    def neutral_reaction_probability(
        self, state: GuoC4F8ArSiO2FeatureState
    ):
        if not isinstance(state, GuoC4F8ArSiO2FeatureState):
            raise TypeError("Guo feature probabilities require Guo feature state")
        atomic_f = 20.0 * self._bond_probability_arrays(state, "Si")
        generic = (
            3.5 * self._bond_probability_arrays(state, "C")
            + 1.8 * self._bond_probability_arrays(state, "O")
        )
        output = {
            name: np.asarray(atomic_f if name == "F" else generic, dtype=float)
            for name in self.neutral_species
        }
        for name, value in output.items():
            if (np.any(~np.isfinite(value))
                    or np.any((value < 0.0) | (value > 1.0))):
                raise GuoSourceLawUnderspecified(
                    f"Guo state makes {name} adsorption yield leave [0,1]; "
                    "a collision-level saturation law is not published")
        return MappingProxyType(output)

    @staticmethod
    def _positive_energy_angle(population):
        if isinstance(population, FaceResolvedEnergeticFlux):
            selected = population.event_flux_m2_s > 0.0
            return (
                population.event_energy_eV[selected],
                population.event_cosine_incidence[selected],
            )
        if isinstance(population, EnergeticFlux):
            if not np.any(np.asarray(population.flux_m2_s) > 0.0):
                return np.empty(0), np.empty(0)
            selected = population.weight > 0.0
            return (
                population.energy_eV[selected],
                population.cosine_incidence[selected],
            )
        raise TypeError(type(population).__name__)  # pragma: no cover

    def validity(self, state, fluxes):
        if not isinstance(state, GuoC4F8ArSiO2FeatureState):
            raise TypeError("Guo feature validity requires Guo feature state")
        unsupported_neutral = tuple(sorted(
            name for name, value in fluxes.neutral_flux_m2_s.items()
            if name not in self.neutral_species
            and np.any(np.asarray(value, dtype=float) > 0.0)
        ))
        unsupported_energetic = tuple(sorted({
            population.name for population in fluxes.energetic_fluxes
            if population.name not in self.ion_species_mapping
            and np.any(np.asarray(population.flux_m2_s, dtype=float) > 0.0)
        }))
        outside_source_neutrals = tuple(sorted(
            name for name in self.neutral_species
            if name not in _SOURCE_NEUTRALS
            and np.any(np.asarray(
                fluxes.neutral_flux_m2_s.get(name, 0.0),
                dtype=float,
            ) > 0.0)
        ))
        leaves_energy = False
        below_defined_threshold = False
        off_normal_uses_inferred_polynomial = False
        aggregate_ion_composition = False
        for population in fluxes.energetic_fluxes:
            if population.name not in self.ion_species_mapping:
                continue
            energy, cosine = self._positive_energy_angle(population)
            leaves_energy |= bool(np.any(
                energy > GuoC4F8ArSiO2Mechanism.source_energy_max_eV))
            below_defined_threshold |= bool(np.any(
                energy
                < GuoC4F8ArSiO2Mechanism.maximum_printed_threshold_eV))
            off_normal_uses_inferred_polynomial |= bool(
                np.any(cosine < 1.0 - 1.0e-12))
            aggregate_ion_composition |= bool(
                energy.size
                and self.ion_species_mapping[population.name] is None)
        reasons = []
        if unsupported_neutral or unsupported_energetic:
            reasons.append("positive incident flux has no declared Guo channel")
        if outside_source_neutrals:
            reasons.append(
                "neutral species leave Guo's printed source list: "
                + ", ".join(outside_source_neutrals))
        if leaves_energy:
            reasons.append(
                "ion energy leaves the <=370 eV Guo/Yin regression board")
        if below_defined_threshold:
            reasons.append(
                "generic ion incorporation threshold is undefined below 72 eV")
        if off_normal_uses_inferred_polynomial:
            reasons.append(
                "off-normal physical sputtering uses the declared Table-4.2 "
                "typesetting repair, not an independently traced original fit")
        if aggregate_ion_composition:
            reasons.append(
                "aggregate ion population has no species-resolved composition")
        return MechanismValidity(
            within_declared_scope=not reasons,
            reasons=tuple(reasons),
            unsupported_neutral_species=unsupported_neutral,
            known_model_form_omissions=(
                "C4F6 parent flux and parent/ion co-incidence are unpublished",
                "C2F3/C3F4 use the source's generic-neutral topology outside "
                "its printed species list",
                "the Guo/Yin deck is beam-regressed, not an elementary "
                "reaction network or interatomic potential",
                "surface-product launch identities and energy-angle laws are "
                "unresolved, so removed oxide cannot be return-transported",
                "ion-free spontaneous etching is not advanced on faces with "
                "zero energetic flux",
            ),
            parameter_evidence_supports_prediction=False,
            nonpredictive_parameters=(
                "aggregate positive-ion composition",
                "C2F3/C3F4 generic-neutral transfer",
                "reaction coefficients above 370 eV",
                "Table-4.2 physical angular exponent repair",
            ),
        )

    @staticmethod
    def _broadcast_flux(value, shape):
        return np.broadcast_to(np.asarray(value, dtype=float), shape)

    def _face_ion_measure(self, fluxes, shape):
        size = int(np.prod(shape, dtype=int)) if shape else 1
        energies = [[] for _ in range(size)]
        cosines = [[] for _ in range(size)]
        measures = [[] for _ in range(size)]
        formula_flux = [dict() for _ in range(size)]
        total_flux = np.zeros(size)

        for population in fluxes.energetic_fluxes:
            if population.name not in self.ion_species_mapping:
                continue
            formula = self.ion_species_mapping[population.name]
            if isinstance(population, FaceResolvedEnergeticFlux):
                if shape != (population.face_count,):
                    raise ValueError(
                        "face-resolved Guo ion measure requires one state per face")
                for event in range(population.event_face.size):
                    contribution = float(population.event_flux_m2_s[event])
                    if contribution <= 0.0:
                        continue
                    face = int(population.event_face[event])
                    energies[face].append(
                        float(population.event_energy_eV[event]))
                    cosines[face].append(
                        float(population.event_cosine_incidence[event]))
                    measures[face].append(contribution)
                    total_flux[face] += contribution
                    if formula is not None:
                        formula_flux[face][formula] = (
                            formula_flux[face].get(formula, 0.0)
                            + contribution
                        )
            elif isinstance(population, EnergeticFlux):
                face_flux = self._broadcast_flux(
                    population.flux_m2_s, shape).reshape(-1)
                for face, supplied in enumerate(face_flux):
                    if supplied <= 0.0:
                        continue
                    contribution = supplied * population.weight
                    positive = contribution > 0.0
                    energies[face].extend(
                        population.energy_eV[positive].tolist())
                    cosines[face].extend(
                        population.cosine_incidence[positive].tolist())
                    measures[face].extend(contribution[positive].tolist())
                    total_flux[face] += float(supplied)
                    if formula is not None:
                        formula_flux[face][formula] = (
                            formula_flux[face].get(formula, 0.0)
                            + float(supplied)
                        )
            else:  # pragma: no cover - SurfaceFluxes validates.
                raise TypeError(type(population).__name__)
        return total_flux, energies, cosines, measures, formula_flux

    def advance(
        self,
        state: GuoC4F8ArSiO2FeatureState,
        fluxes: SurfaceFluxes,
        duration_s: float,
        *,
        strict: bool = True,
    ):
        if not isinstance(state, GuoC4F8ArSiO2FeatureState):
            raise TypeError("Guo feature mechanism requires Guo feature state")
        duration = float(duration_s)
        if not np.isfinite(duration) or duration < 0.0:
            raise ValueError("duration_s must be finite and nonnegative")
        validity = self.validity(state, fluxes)
        if (strict and not validity.within_declared_scope
                and not self.allow_out_of_board_transfer_audit):
            raise ValueError(
                "surface mechanism outside declared scope: "
                + "; ".join(validity.reasons))

        shape = np.asarray(state.si).shape
        (
            ion_flux,
            face_energy,
            face_cosine,
            face_measure,
            face_formula_flux,
        ) = self._face_ion_measure(fluxes, shape)
        neutral_flux = {
            name: self._broadcast_flux(
                fluxes.neutral_flux_m2_s.get(name, 0.0), shape
            ).reshape(-1)
            for name in self.neutral_species
        }
        output = {
            name: np.asarray(getattr(state, name), dtype=float).reshape(-1).copy()
            for name in _REAL_FIELDS + ("vacancy",)
        }
        yield_per_ion = np.zeros_like(ion_flux)
        residual = np.zeros_like(ion_flux)
        atom_residual = np.zeros_like(ion_flux)
        movement_atoms_per_ion = np.zeros_like(ion_flux)
        bdf_fallback_face_count = 0

        for face, supplied_ion_flux in enumerate(ion_flux):
            if supplied_ion_flux <= 0.0:
                continue
            measure = np.asarray(face_measure[face], dtype=float)
            quadrature = GuoIonQuadrature(
                np.asarray(face_energy[face], dtype=float),
                np.asarray(face_cosine[face], dtype=float),
                measure,
            )
            composition = GuoIncidentComposition(
                {
                    name: values[face] / supplied_ion_flux
                    for name, values in neutral_flux.items()
                    if values[face] > 0.0
                },
                {
                    formula: value / supplied_ion_flux
                    for formula, value in face_formula_flux[face].items()
                },
            )
            local = GuoC4F8ArSiO2Mechanism(composition, quadrature)
            initial = GuoTmlState(
                output["si"][face],
                output["o"][face],
                output["c"][face],
                output["f"][face],
                output["vacancy"][face],
            )
            try:
                solved = local.solve_steady_state_algebraic(initial)
            except RuntimeError as algebraic_error:
                try:
                    # Guo integrated the differential balances from the
                    # substrate composition.  The BDF path is therefore the
                    # source-faithful fallback when the faster constrained
                    # root misses a bound-active deposition solution.
                    solved = local.solve_steady_state(initial)
                except RuntimeError as bdf_error:
                    ratios = {
                        name: values[face] / supplied_ion_flux
                        for name, values in neutral_flux.items()
                        if values[face] > 0.0
                    }
                    raise RuntimeError(
                        "Guo local feature steady state failed at "
                        f"face={face}, "
                        f"ion_flux_m2_s={supplied_ion_flux:.9g}, "
                        f"energy_eV=[{min(face_energy[face]):.9g},"
                        f"{max(face_energy[face]):.9g}], "
                        f"neutral_to_ion={ratios}; algebraic="
                        f"{algebraic_error}; BDF={bdf_error}"
                    ) from bdf_error
                bdf_fallback_face_count += 1
            for name, value in zip(
                _REAL_FIELDS + ("vacancy",),
                solved.state.as_array(),
            ):
                output[name][face] = value
            yield_per_ion[face] = solved.sio2_yield_per_ion
            movement_atoms_per_ion[face] = solved.movement_atoms_per_ion
            residual[face] = solved.steady_state_residual
            atom_residual[face] = (
                solved.movement_atoms_per_ion
                - sum(solved.removed_atoms_per_ion.values())
                + sum(solved.incoming_atoms_per_ion.values())
            )

        movement_atom_rate = movement_atoms_per_ion * ion_flux
        removal_rate = np.maximum(movement_atom_rate, 0.0) / 3.0
        deposition_rate = np.maximum(-movement_atom_rate, 0.0)
        removed = removal_rate * duration
        deposited = deposition_rate * duration
        next_state = GuoC4F8ArSiO2FeatureState(
            *[
                output[name].reshape(shape)
                for name in _REAL_FIELDS + ("vacancy",)
            ],
            (
                np.asarray(state.removed_sio2_formula_units_m2)
                + removed.reshape(shape)
            ),
        )
        probability = self.neutral_reaction_probability(next_state)
        maximum_probability = max(
            (float(np.max(value)) for value in probability.values()),
            default=0.0,
        )
        exchange = unresolved_surface_exchange(
            removed_units_m2={
                "SiO2_formula_unit": removed.reshape(shape),
            },
            deposited_units_m2={
                "unresolved_mixed_layer_atom":
                    deposited.reshape(shape),
            },
            limitations=(
                "Guo product identities are atom-counted internally but "
                "product launch energy-angle laws are unresolved",
                "net translating-layer deposition composition is unresolved",
            ),
        )
        return GuoC4F8ArSiO2FeatureStepResult(
            state=next_state,
            etch_velocity_m_s=(
                removal_rate / self.bulk_sio2_formula_density_m3
            ).reshape(shape),
            normal_growth_velocity_m_s=(
                deposition_rate / self.deposited_film_atom_density_m3
            ).reshape(shape),
            removed_sio2_formula_units_m2=removed.reshape(shape),
            deposited_mixed_layer_atoms_m2=deposited.reshape(shape),
            sio2_yield_per_ion=yield_per_ion.reshape(shape),
            net_movement_atoms_per_ion=
                movement_atoms_per_ion.reshape(shape),
            steady_state_residual=residual.reshape(shape),
            atom_ledger_residual_atoms_per_ion=atom_residual.reshape(shape),
            bdf_fallback_face_count=bdf_fallback_face_count,
            maximum_neutral_reaction_probability=maximum_probability,
            material_exchange=exchange,
            validity=validity,
        )


__all__ = [
    "GuoC4F8ArSiO2FeatureMechanism",
    "GuoC4F8ArSiO2FeatureState",
    "GuoC4F8ArSiO2FeatureStepResult",
]
