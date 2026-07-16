"""Material-ID routing for one common transport/charging/profile engine."""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .surface_exchange import (
    SurfaceMaterialExchange, SurfaceProductPopulation, validate_surface_product_routing,
)
from .surface_kinetics import (
    EnergeticFlux, FaceResolvedEnergeticFlux, MechanismValidity, SurfaceFluxes,
)


def _state_prefix(material_id, field_name):
    if (not isinstance(field_name, str) or not field_name
            or any(character not in "abcdefghijklmnopqrstuvwxyz0123456789_"
                   for character in field_name)):
        raise ValueError("material mechanism state fields must be safe lowercase identifiers")
    return f"m{int(material_id)}__{field_name}"


def _provenance_value(value, path="provenance"):
    """Freeze the router evidence to a JSON-compatible value at construction time."""
    if value is None or isinstance(value, (str, bool, int, float)):
        if isinstance(value, float) and not np.isfinite(value):
            raise ValueError(f"{path} contains a non-finite value")
        return value
    if isinstance(value, np.generic):
        return _provenance_value(value.item(), path)
    if isinstance(value, Mapping):
        output = {}
        for key, item in value.items():
            if not isinstance(key, (str, int)):
                raise ValueError(f"{path} contains an invalid mapping key")
            output[str(key)] = _provenance_value(item, f"{path}.{key}")
        return output
    if isinstance(value, (tuple, list)):
        return [_provenance_value(item, f"{path}[{index}]")
                for index, item in enumerate(value)]
    raise ValueError(f"{path} must be machine-readable, not {type(value).__name__}")


@dataclass(frozen=True)
class MaterialSurfaceState3D:
    """Full active-face fields namespaced by material for conservative generic remap."""

    fields: Mapping[str, np.ndarray]
    upper_bounds: Mapping[str, float | None]
    remap_modes: Mapping[str, str] | None = None

    def __post_init__(self):
        fields = {}
        for name, supplied in dict(self.fields).items():
            value = np.asarray(supplied, dtype=float).copy()
            if (not isinstance(name, str) or not name or not name.startswith("m")
                    or "__" not in name or np.any(~np.isfinite(value))
                    or np.any(value < 0.0)):
                raise ValueError("invalid material surface-state field")
            value.setflags(write=False)
            fields[name] = value
        upper = dict(self.upper_bounds)
        if upper and set(upper) != set(fields):
            raise ValueError("material surface-state upper bounds must match its fields")
        for name, value in upper.items():
            if value is not None and (not np.isfinite(value) or value <= 0.0):
                raise ValueError(f"invalid upper bound for {name}")
        modes = ({name: "conservative" for name in fields}
                 if self.remap_modes is None else dict(self.remap_modes))
        if (set(modes) != set(fields)
                or any(value not in {"conservative", "intensive"}
                       for value in modes.values())):
            raise ValueError("material surface-state remap modes must match its fields")
        object.__setattr__(self, "fields", MappingProxyType(fields))
        object.__setattr__(self, "upper_bounds", MappingProxyType(upper))
        object.__setattr__(self, "remap_modes", MappingProxyType(modes))

    @classmethod
    def bare(cls, shape=()):
        # Safe-checkpoint restoration supplies the explicitly serialized names immediately after.
        return cls({}, {}, {})

    def conservative_surface_fields(self):
        return dict(self.fields)

    def conservative_surface_upper_bounds(self):
        return ({name: None for name in self.fields}
                if not self.upper_bounds else dict(self.upper_bounds))

    def surface_field_remap_modes(self):
        return dict(self.remap_modes)

    def with_conservative_surface_fields(self, fields):
        fields = dict(fields)
        if self.fields and set(fields) != set(self.fields):
            raise ValueError("material remap fields do not match its state contract")
        upper = (dict(self.upper_bounds) if self.upper_bounds
                 else {name: None for name in fields})
        modes = (dict(self.remap_modes) if self.remap_modes
                 else {name: "conservative" for name in fields})
        return type(self)(fields, upper, modes)


@dataclass(frozen=True)
class MaterialSurfaceStepResult3D:
    state: MaterialSurfaceState3D
    etch_velocity_m_s: np.ndarray
    normal_growth_velocity_m_s: np.ndarray
    material_exchange: SurfaceMaterialExchange
    product_populations: tuple[SurfaceProductPopulation, ...]
    validity: MechanismValidity
    material_results: Mapping[int, object]

    def __post_init__(self):
        velocity = np.asarray(self.etch_velocity_m_s, dtype=float).copy()
        growth = np.asarray(self.normal_growth_velocity_m_s, dtype=float).copy()
        if (not isinstance(self.state, MaterialSurfaceState3D)
                or velocity.ndim != 1 or np.any(~np.isfinite(velocity))
                or np.any(velocity < 0.0)
                or growth.shape != velocity.shape or np.any(~np.isfinite(growth))
                or np.any(growth < 0.0)
                or not isinstance(self.material_exchange, SurfaceMaterialExchange)
                or any(not isinstance(item, SurfaceProductPopulation)
                       for item in self.product_populations)
                or not isinstance(self.validity, MechanismValidity)):
            raise ValueError("invalid material-routed surface result")
        velocity.setflags(write=False)
        growth.setflags(write=False)
        object.__setattr__(self, "etch_velocity_m_s", velocity)
        object.__setattr__(self, "normal_growth_velocity_m_s", growth)
        object.__setattr__(
            self, "product_populations",
            validate_surface_product_routing(
                self.material_exchange, tuple(self.product_populations)))
        object.__setattr__(
            self, "material_results", MappingProxyType(dict(self.material_results)))


def _subset_fluxes(fluxes, selected, face_count):
    selected = np.asarray(selected, dtype=int)
    old_to_new = np.full(int(face_count), -1, dtype=int)
    old_to_new[selected] = np.arange(len(selected))
    neutral = {
        name: np.asarray(value)[selected]
        for name, value in fluxes.neutral_flux_m2_s.items()}
    energetic = []
    for population in fluxes.energetic_fluxes:
        if isinstance(population, FaceResolvedEnergeticFlux):
            mapped = old_to_new[population.event_face]
            retained = mapped >= 0
            energetic.append(FaceResolvedEnergeticFlux(
                population.name, len(selected), mapped[retained],
                population.event_flux_m2_s[retained], population.event_energy_eV[retained],
                population.event_cosine_incidence[retained],
                event_position=(None if population.event_position is None
                                else population.event_position[retained]),
                event_incident_direction=(
                    None if population.event_incident_direction is None
                    else population.event_incident_direction[retained])))
        elif isinstance(population, EnergeticFlux):
            flux = np.asarray(population.flux_m2_s)
            energetic.append(EnergeticFlux(
                population.name, flux if flux.ndim == 0 else flux[selected],
                population.energy_eV, population.cosine_incidence, population.weight))
        else:  # pragma: no cover - SurfaceFluxes validates this contract.
            raise TypeError(type(population).__name__)
    return SurfaceFluxes(neutral, tuple(energetic))


def _expanded_inventory(target, local, selected, face_count):
    for name, value in local.items():
        supplied = np.asarray(value, dtype=float)
        try:
            supplied = np.broadcast_to(supplied, (len(selected),))
        except ValueError as error:
            raise ValueError(f"material exchange inventory {name!r} has the wrong shape") from error
        if name not in target:
            target[name] = np.zeros(face_count)
        target[name][selected] += supplied


class MaterialMechanismRouter3D:
    """Dispatch exposed faces to independent material laws without branching the engine."""

    def __init__(self, mechanisms: Mapping[int, object], *, provenance: Mapping[int, object]):
        mechanisms = {int(key): value for key, value in dict(mechanisms).items()}
        evidence = {int(key): value for key, value in dict(provenance).items()}
        if (not mechanisms or any(key <= 0 for key in mechanisms)
                or set(evidence) != set(mechanisms)
                or any(value is None for value in mechanisms.values())
                or any(value in (None, "", {}) for value in evidence.values())):
            raise ValueError("every routed material requires a mechanism and provenance")
        self.mechanisms = MappingProxyType(mechanisms)
        self.provenance = MappingProxyType(dict(
            model="material-mechanism-router-3d-v1",
            materials={str(key): dict(
                mechanism=type(mechanisms[key]).__name__,
                evidence=_provenance_value(
                    evidence[key], f"material_router.materials.{key}.evidence"))
                for key in sorted(mechanisms)}))

    def _validate_materials(self, face_material_id):
        material = np.asarray(face_material_id, dtype=int)
        if (material.ndim != 1 or np.any(material <= 0)
                or not set(np.unique(material)).issubset(self.mechanisms)):
            missing = set(np.unique(material)) - set(self.mechanisms)
            raise ValueError(f"material router has no mechanism for ids {sorted(missing)}")
        return material

    def initial_state_by_material(self, face_material_id):
        material = self._validate_materials(face_material_id)
        fields = {}
        upper = {}
        modes = {}
        for material_id in sorted(set(material)):
            selected = np.where(material == material_id)[0]
            mechanism = self.mechanisms[int(material_id)]
            if not hasattr(mechanism, "initial_state"):
                raise TypeError(f"material {material_id} mechanism has no initial_state")
            state = mechanism.initial_state((len(selected),))
            if (not hasattr(state, "conservative_surface_fields")
                    or not hasattr(state, "conservative_surface_upper_bounds")
                    or not hasattr(state, "with_conservative_surface_fields")):
                raise TypeError(f"material {material_id} state is not conservatively remappable")
            local_fields = dict(state.conservative_surface_fields())
            local_upper = dict(state.conservative_surface_upper_bounds())
            local_modes = (
                dict(state.surface_field_remap_modes())
                if hasattr(state, "surface_field_remap_modes")
                else {name: "conservative" for name in local_fields})
            if (not local_fields or set(local_fields) != set(local_upper)
                    or set(local_fields) != set(local_modes)):
                raise ValueError(f"material {material_id} state contract is incomplete")
            for name, value in local_fields.items():
                key = _state_prefix(material_id, name)
                fields[key] = np.zeros(len(material))
                fields[key][selected] = np.asarray(value, dtype=float)
                upper[key] = local_upper[name]
                modes[key] = local_modes[name]
        return MaterialSurfaceState3D(fields, upper, modes)

    def _local_state(self, state, material_id, selected):
        mechanism = self.mechanisms[int(material_id)]
        template = mechanism.initial_state((len(selected),))
        expected = dict(template.conservative_surface_fields())
        local = {}
        upper = {}
        for name in expected:
            key = _state_prefix(material_id, name)
            if key not in state.fields:
                raise ValueError(f"material state is missing {key}")
            local[name] = state.fields[key][selected]
            upper[key] = template.conservative_surface_upper_bounds()[name]
        return template.with_conservative_surface_fields(local), upper

    def neutral_reaction_probability_by_material(
            self, state: MaterialSurfaceState3D, face_material_id):
        material = self._validate_materials(face_material_id)
        if not isinstance(state, MaterialSurfaceState3D):
            raise TypeError("material router requires MaterialSurfaceState3D")
        output = {}
        for material_id in sorted(set(material)):
            selected = np.where(material == material_id)[0]
            mechanism = self.mechanisms[int(material_id)]
            if not hasattr(mechanism, "neutral_reaction_probability"):
                continue
            local_state, _ = self._local_state(state, material_id, selected)
            for name, value in mechanism.neutral_reaction_probability(local_state).items():
                if name not in output:
                    output[name] = np.zeros(len(material))
                output[name][selected] = np.broadcast_to(value, (len(selected),))
        return output

    def advance_by_material(
            self, state: MaterialSurfaceState3D, fluxes: SurfaceFluxes,
            duration_s: float, face_material_id):
        material = self._validate_materials(face_material_id)
        if not isinstance(state, MaterialSurfaceState3D):
            raise TypeError("material router requires MaterialSurfaceState3D")
        face_count = len(material)
        expected_state = self.initial_state_by_material(material)
        if set(state.fields) != set(expected_state.fields) or any(
                np.asarray(value).shape != (face_count,) for value in state.fields.values()):
            raise ValueError("material surface state does not match the active material mesh")
        output_fields = {name: np.asarray(value).copy() for name, value in state.fields.items()}
        upper = dict(expected_state.upper_bounds)
        modes = dict(expected_state.remap_modes)
        velocity = np.zeros(face_count)
        growth = np.zeros(face_count)
        removed = {}; outgoing = {}; unresolved = {}; deposited = {}
        products = []; results = {}; reasons = []; unsupported = []; omissions = []
        nonpredictive = []; evidence_supports = True; exchange_limitations = []
        for material_id in sorted(set(material)):
            selected = np.where(material == material_id)[0]
            mechanism = self.mechanisms[int(material_id)]
            local_state, _ = self._local_state(state, material_id, selected)
            result = mechanism.advance(
                local_state, _subset_fluxes(fluxes, selected, face_count), float(duration_s))
            results[int(material_id)] = result
            local_velocity = np.asarray(result.etch_velocity_m_s, dtype=float)
            velocity[selected] = np.broadcast_to(local_velocity, (len(selected),))
            local_growth = np.asarray(
                getattr(result, "normal_growth_velocity_m_s", 0.0), dtype=float)
            growth[selected] = np.broadcast_to(local_growth, (len(selected),))
            for name, value in result.state.conservative_surface_fields().items():
                output_fields[_state_prefix(material_id, name)][selected] = value
            exchange = getattr(result, "material_exchange", None)
            if not isinstance(exchange, SurfaceMaterialExchange):
                raise TypeError(f"material {material_id} result lacks a material-exchange ledger")
            _expanded_inventory(removed, exchange.removed_units_m2, selected, face_count)
            _expanded_inventory(outgoing, exchange.outgoing_units_m2, selected, face_count)
            _expanded_inventory(unresolved, exchange.unresolved_units_m2, selected, face_count)
            _expanded_inventory(deposited, exchange.deposited_units_m2, selected, face_count)
            exchange_limitations.extend(
                f"material {material_id}: {item}" for item in exchange.known_limitations)
            for population in tuple(getattr(result, "product_populations", ())):
                count = np.zeros(face_count)
                count[selected] = np.broadcast_to(
                    population.integrated_particle_count_m2, (len(selected),))
                products.append(SurfaceProductPopulation(
                    name=f"material_{material_id}:{population.name}",
                    source_inventory=population.source_inventory,
                    integrated_particle_count_m2=count,
                    material_units_per_particle=population.material_units_per_particle,
                    mass_amu=population.mass_amu, angular_model=population.angular_model,
                    energy_model=population.energy_model,
                    energy_parameters=population.energy_parameters,
                    provenance=dict(population.provenance, material_id=int(material_id)),
                    relative_standard_uncertainty=population.relative_standard_uncertainty))
            validity = result.validity
            reasons.extend(f"material {material_id}: {item}" for item in validity.reasons)
            unsupported.extend(validity.unsupported_neutral_species)
            omissions.extend(
                f"material {material_id}: {item}"
                for item in validity.known_model_form_omissions)
            evidence_supports &= validity.parameter_evidence_supports_prediction
            nonpredictive.extend(
                f"material_{material_id}.{item}" for item in validity.nonpredictive_parameters)
        exchange = SurfaceMaterialExchange(
            removed, outgoing, unresolved, deposited, tuple(exchange_limitations))
        validity = MechanismValidity(
            within_declared_scope=not reasons, reasons=tuple(reasons),
            unsupported_neutral_species=tuple(sorted(set(unsupported))),
            known_model_form_omissions=tuple(omissions),
            parameter_evidence_supports_prediction=evidence_supports,
            nonpredictive_parameters=tuple(nonpredictive))
        return MaterialSurfaceStepResult3D(
            state=MaterialSurfaceState3D(output_fields, upper, modes),
            etch_velocity_m_s=velocity, normal_growth_velocity_m_s=growth,
            material_exchange=exchange, product_populations=tuple(products),
            validity=validity, material_results=results)
