"""Embedded surface-state / diffuse-radiosity coupling on one fixed 3-D mesh.

This module advances the existing material mechanism while repeatedly solving the
existing diffuse-neutral radiosity equation.  Direct transport and geometric form
factors are immutable cached inputs; only state-dependent reaction probabilities and
the resulting radiosity are recomputed.  It deliberately contains no geometry move,
benchmark parameters, transport tracing, or profile-driver policy.

The v1 contract is intentionally narrow:

* every surface face must be routed to a material mechanism;
* emitted product populations are refused rather than silently discarded;
* one-full-step versus two-half-step error estimates cover every conservative state
  increment, every exchange inventory, and per-face recession/growth;
* rejected trials do not mutate state or contribute to the accepted ledger; and
* a caller-supplied exact operator identity is checked on every use of the cache.

Radiosity and material ledgers are each certified exactly by their authoritative
operators.  A generic neutral-species-to-material stoichiometric cross-ledger is not
claimed because the material-mechanism API does not expose that mapping.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
import json
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .material_mechanism_3d import (
    MaterialSurfaceState3D,
    MaterialSurfaceStepResult3D,
)
from .neutral_radiosity_3d import (
    DiffuseFormFactors3D,
    solve_diffuse_neutral_radiosity_3d,
)
from .surface_exchange import SurfaceMaterialExchange
from .surface_kinetics import (
    EnergeticFlux,
    FaceResolvedEnergeticFlux,
    SurfaceFluxes,
)


class SurfaceRadiosityCouplingRefusal(RuntimeError):
    """The declared fixed-mesh coupling contract could not be satisfied."""


class SurfaceRadiosityCacheIdentityError(SurfaceRadiosityCouplingRefusal):
    """A cached geometric/operator epoch was used with a different identity."""


def _identity_jsonable(value, path="identity"):
    """Return a deterministic JSON value or refuse an ambiguous identity input."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, (float, np.floating)):
        value = float(value)
        if not np.isfinite(value):
            raise ValueError(f"{path} contains a non-finite float")
        return value
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, Mapping):
        output = {}
        for key, item in value.items():
            if not isinstance(key, str) or not key:
                raise ValueError(f"{path} keys must be nonempty strings")
            output[key] = _identity_jsonable(item, f"{path}.{key}")
        return output
    if isinstance(value, (tuple, list)):
        return [
            _identity_jsonable(item, f"{path}[{index}]")
            for index, item in enumerate(value)
        ]
    raise ValueError(
        f"{path} must be JSON-compatible evidence, not {type(value).__name__}")


def _canonical_identity(payload):
    if not isinstance(payload, Mapping) or not payload:
        raise ValueError("operator identity payload must be a nonempty mapping")
    value = _identity_jsonable(payload)
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False).encode("ascii")
    return value, encoded


def _freeze_json(value):
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze_json(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze_json(item) for item in value)
    return value


def _digest_array(digest, name, supplied):
    array = np.ascontiguousarray(np.asarray(supplied))
    digest.update(name.encode("utf-8") + b"\0")
    digest.update(array.dtype.str.encode("ascii") + b"\0")
    digest.update(json.dumps(array.shape, separators=(",", ":")).encode("ascii") + b"\0")
    digest.update(array.tobytes())


def _copy_surface_fluxes(fluxes):
    if not isinstance(fluxes, SurfaceFluxes):
        raise TypeError("direct surface flux must be SurfaceFluxes")
    energetic = []
    for population in fluxes.energetic_fluxes:
        if isinstance(population, EnergeticFlux):
            energetic.append(EnergeticFlux(
                population.name, population.flux_m2_s, population.energy_eV,
                population.cosine_incidence, population.weight))
        elif isinstance(population, FaceResolvedEnergeticFlux):
            energetic.append(FaceResolvedEnergeticFlux(
                population.name, population.face_count, population.event_face,
                population.event_flux_m2_s, population.event_energy_eV,
                population.event_cosine_incidence,
                event_position=population.event_position,
                event_incident_direction=population.event_incident_direction))
        else:  # pragma: no cover - SurfaceFluxes validates this already.
            raise TypeError(type(population).__name__)
    return SurfaceFluxes(dict(fluxes.neutral_flux_m2_s), tuple(energetic))


def _surface_flux_digest(digest, fluxes):
    for name, value in sorted(fluxes.neutral_flux_m2_s.items()):
        _digest_array(digest, f"neutral:{name}", value)
    for index, population in enumerate(fluxes.energetic_fluxes):
        prefix = f"energetic:{index}:{type(population).__name__}:{population.name}"
        if isinstance(population, EnergeticFlux):
            for name in ("flux_m2_s", "energy_eV", "cosine_incidence", "weight"):
                _digest_array(digest, f"{prefix}:{name}", getattr(population, name))
        else:
            digest.update(f"{prefix}:face_count={population.face_count}\n".encode("utf-8"))
            for name in (
                    "event_face", "event_flux_m2_s", "event_energy_eV",
                    "event_cosine_incidence"):
                _digest_array(digest, f"{prefix}:{name}", getattr(population, name))
            for name in ("event_position", "event_incident_direction"):
                value = getattr(population, name)
                digest.update(f"{prefix}:{name}:none={value is None}\n".encode("utf-8"))
                if value is not None:
                    _digest_array(digest, f"{prefix}:{name}", value)


@dataclass(frozen=True)
class RadiositySpeciesDiagnostic3D:
    source_rate_s: float
    reacted_rate_s: float
    escaped_rate_s: float
    relative_balance_error: float
    relative_linear_residual: float
    solver_method: str
    iteration_count: int
    inactive_face_count: int
    repeated_incident_flux_elided: bool


@dataclass(frozen=True)
class CachedRadiosityEvaluation3D:
    surface_fluxes: SurfaceFluxes
    reaction_probability: Mapping[str, np.ndarray]
    species_diagnostics: Mapping[str, RadiositySpeciesDiagnostic3D]
    maximum_relative_balance_error: float

    def __post_init__(self):
        if not isinstance(self.surface_fluxes, SurfaceFluxes):
            raise TypeError("cached radiosity evaluation requires SurfaceFluxes")
        probability = {}
        for name, supplied in dict(self.reaction_probability).items():
            value = np.asarray(supplied, dtype=float).copy()
            if (not name or value.ndim != 1 or np.any(~np.isfinite(value))
                    or np.any((value < 0.0) | (value > 1.0))):
                raise ValueError("invalid cached reaction probability")
            value.setflags(write=False)
            probability[name] = value
        diagnostics = dict(self.species_diagnostics)
        if (set(diagnostics) != set(probability)
                or any(not isinstance(value, RadiositySpeciesDiagnostic3D)
                       for value in diagnostics.values())):
            raise ValueError("radiosity diagnostics must cover every neutral species")
        maximum = float(self.maximum_relative_balance_error)
        if not np.isfinite(maximum) or maximum < 0.0:
            raise ValueError("invalid maximum radiosity balance")
        object.__setattr__(self, "reaction_probability", MappingProxyType(probability))
        object.__setattr__(self, "species_diagnostics", MappingProxyType(diagnostics))
        object.__setattr__(self, "maximum_relative_balance_error", maximum)


@dataclass(frozen=True)
class SurfaceRadiosityOperatorCache3D:
    """Immutable direct-transport and geometric-radiosity operator for one epoch."""

    direct_surface_fluxes: SurfaceFluxes
    form_factors: DiffuseFormFactors3D
    face_area_m2: np.ndarray
    active_face_index: np.ndarray
    face_material_id: np.ndarray
    species_role: Mapping[str, str]
    operator_identity_payload: Mapping[str, object]
    relative_tolerance: float = 1.0e-10
    maximum_iterations: int = 500
    operator_identity_sha256: str = field(init=False)
    cache_sha256: str = field(init=False)
    _identity_encoding: bytes = field(init=False, repr=False, compare=False)

    def __post_init__(self):
        direct = _copy_surface_fluxes(self.direct_surface_fluxes)
        if not isinstance(self.form_factors, DiffuseFormFactors3D):
            raise TypeError("form_factors must be DiffuseFormFactors3D")
        factors = DiffuseFormFactors3D(
            self.form_factors.face_count, self.form_factors.source_face,
            self.form_factors.target_face, self.form_factors.transfer_fraction,
            self.form_factors.escape_fraction, self.form_factors.rays_per_face)
        area = np.asarray(self.face_area_m2, dtype=float).copy()
        active = np.asarray(self.active_face_index, dtype=int).copy()
        material = np.asarray(self.face_material_id, dtype=int).copy()
        n_face = factors.face_count
        if (area.shape != (n_face,) or np.any(~np.isfinite(area))
                or np.any(area <= 0.0)):
            raise ValueError("face areas must be finite, positive, and match form factors")
        if (active.shape != (n_face,)
                or not np.array_equal(active, np.arange(n_face))):
            raise SurfaceRadiosityCouplingRefusal(
                "surface-radiosity coupling v1 refuses incomplete active-face routing")
        if (material.shape != (n_face,) or np.any(material <= 0)):
            raise ValueError("every cached face requires a positive material id")
        for name, value in direct.neutral_flux_m2_s.items():
            if np.asarray(value).shape != (n_face,):
                raise ValueError(f"direct neutral flux {name!r} must be face-resolved")
        species_names = set(direct.neutral_flux_m2_s)
        for population in direct.energetic_fluxes:
            species_names.add(population.name)
            if (isinstance(population, EnergeticFlux)
                    and np.asarray(population.flux_m2_s).ndim
                    and np.asarray(population.flux_m2_s).shape != (n_face,)):
                raise ValueError(
                    f"energetic flux {population.name!r} must be scalar or face-resolved")
            if (isinstance(population, FaceResolvedEnergeticFlux)
                    and population.face_count != n_face):
                raise ValueError(
                    f"energetic events {population.name!r} use a different face mesh")
        role = dict(self.species_role)
        if set(role) != species_names:
            raise SurfaceRadiosityCouplingRefusal(
                "species routing must cover every cached direct population exactly")
        for name in direct.neutral_flux_m2_s:
            if role[name] != "neutral_reactant":
                raise SurfaceRadiosityCouplingRefusal(
                    f"neutral population {name!r} lacks neutral_reactant routing")
        for population in direct.energetic_fluxes:
            if role[population.name] != "energetic_bombardment":
                raise SurfaceRadiosityCouplingRefusal(
                    f"energetic population {population.name!r} lacks energetic routing")
        tolerance = float(self.relative_tolerance)
        iterations = int(self.maximum_iterations)
        if not np.isfinite(tolerance) or tolerance <= 0.0 or iterations <= 0:
            raise ValueError("invalid cached radiosity solver controls")
        identity, identity_encoding = _canonical_identity(
            self.operator_identity_payload)
        identity_digest = sha256(identity_encoding).hexdigest()

        digest = sha256()
        digest.update(b"petch.surface-radiosity-cache-3d.v1\0")
        digest.update(identity_encoding)
        _surface_flux_digest(digest, direct)
        for name in ("source_face", "target_face", "transfer_fraction", "escape_fraction"):
            _digest_array(digest, f"form_factors:{name}", getattr(factors, name))
        digest.update(f"rays_per_face={factors.rays_per_face}\n".encode("ascii"))
        _digest_array(digest, "face_area_m2", area)
        _digest_array(digest, "active_face_index", active)
        _digest_array(digest, "face_material_id", material)
        digest.update(json.dumps(role, sort_keys=True, separators=(",", ":")).encode("utf-8"))
        digest.update(f"relative_tolerance={tolerance.hex()}\n".encode("ascii"))
        digest.update(f"maximum_iterations={iterations}\n".encode("ascii"))

        for value in (area, active, material):
            value.setflags(write=False)
        object.__setattr__(self, "direct_surface_fluxes", direct)
        object.__setattr__(self, "form_factors", factors)
        object.__setattr__(self, "face_area_m2", area)
        object.__setattr__(self, "active_face_index", active)
        object.__setattr__(self, "face_material_id", material)
        object.__setattr__(self, "species_role", MappingProxyType(role))
        object.__setattr__(self, "operator_identity_payload", _freeze_json(identity))
        object.__setattr__(self, "relative_tolerance", tolerance)
        object.__setattr__(self, "maximum_iterations", iterations)
        object.__setattr__(self, "operator_identity_sha256", identity_digest)
        object.__setattr__(self, "cache_sha256", digest.hexdigest())
        object.__setattr__(self, "_identity_encoding", identity_encoding)

    def assert_identity(self, operator_identity_payload):
        _value, encoded = _canonical_identity(operator_identity_payload)
        if encoded != self._identity_encoding:
            raise SurfaceRadiosityCacheIdentityError(
                "cached surface-radiosity operator identity does not match the caller")

    @property
    def face_count(self):
        return int(self.form_factors.face_count)

    def _reaction_probability(self, state, mechanism):
        if not isinstance(state, MaterialSurfaceState3D):
            raise TypeError("surface-radiosity coupling requires MaterialSurfaceState3D")
        if not hasattr(mechanism, "neutral_reaction_probability_by_material"):
            raise TypeError(
                "surface-radiosity coupling requires material-routed reaction probabilities")
        supplied = mechanism.neutral_reaction_probability_by_material(
            state, self.face_material_id)
        expected = set(self.direct_surface_fluxes.neutral_flux_m2_s)
        if not expected.issubset(supplied):
            missing = sorted(expected - set(supplied))
            raise SurfaceRadiosityCouplingRefusal(
                "reaction-probability routing must cover every cached neutral species: "
                f"missing={missing}")
        # A reusable material mechanism may advertise supported channels that are absent from this
        # boundary deck. They carry no incident measure and are intentionally projected out here;
        # missing boundary species still refuse above, including inert species unless the mechanism
        # declares them explicitly with zero probability.
        output = {}
        for name in self.direct_surface_fluxes.neutral_flux_m2_s:
            value = np.asarray(supplied[name], dtype=float)
            try:
                value = np.broadcast_to(value, (self.face_count,)).copy()
            except ValueError as error:
                raise ValueError(
                    f"reaction probability {name!r} does not match cached faces") from error
            if (np.any(~np.isfinite(value))
                    or np.any((value < 0.0) | (value > 1.0))):
                raise ValueError(f"reaction probability {name!r} lies outside [0,1]")
            value.setflags(write=False)
            output[name] = value
        return output

    def solve(self, state, mechanism, *, operator_identity_payload):
        """Resolve state-dependent neutral incidence with the production equation."""
        self.assert_identity(operator_identity_payload)
        probability = self._reaction_probability(state, mechanism)
        neutral = {}
        diagnostics = {}
        for name, supplied in self.direct_surface_fluxes.neutral_flux_m2_s.items():
            direct = np.asarray(supplied, dtype=float)
            reaction = probability[name]
            if not np.any(reaction > 0.0):
                source_rate = float(np.dot(self.face_area_m2, direct))
                neutral[name] = direct
                diagnostics[name] = RadiositySpeciesDiagnostic3D(
                    source_rate, 0.0, source_rate, 0.0, 0.0,
                    "analytic_zero_reaction_elision", 0, self.face_count, True)
                continue
            solution = solve_diffuse_neutral_radiosity_3d(
                direct, self.face_area_m2, self.form_factors.source_face,
                self.form_factors.target_face,
                self.form_factors.transfer_fraction,
                self.form_factors.escape_fraction, reaction,
                relative_tolerance=self.relative_tolerance,
                maximum_iterations=self.maximum_iterations)
            neutral[name] = solution.incident_flux_m2_s
            diagnostics[name] = RadiositySpeciesDiagnostic3D(
                solution.source_rate_s, solution.reacted_rate_s,
                solution.escaped_rate_s, solution.relative_balance_error,
                solution.relative_linear_residual, solution.solver_method,
                solution.iteration_count, solution.inactive_face_count, False)
        maximum_balance = max(
            (value.relative_balance_error for value in diagnostics.values()),
            default=0.0)
        return CachedRadiosityEvaluation3D(
            SurfaceFluxes(neutral, self.direct_surface_fluxes.energetic_fluxes),
            probability, diagnostics, maximum_balance)


@dataclass(frozen=True)
class ArrayErrorNorm3D:
    maximum_absolute_error: float
    relative_linf_error: float
    area_weighted_relative_l1_error: float

    @property
    def maximum_relative_error(self):
        return max(self.relative_linf_error, self.area_weighted_relative_l1_error)


@dataclass(frozen=True)
class EmbeddedSurfaceError3D:
    state_increment: Mapping[str, ArrayErrorNorm3D]
    exchange_inventory: Mapping[str, Mapping[str, ArrayErrorNorm3D]]
    integrated_recession: ArrayErrorNorm3D
    integrated_growth: ArrayErrorNorm3D
    maximum_relative_error: float

    def __post_init__(self):
        state = dict(self.state_increment)
        exchange = {
            kind: MappingProxyType(dict(values))
            for kind, values in dict(self.exchange_inventory).items()
        }
        required = {"removed", "outgoing", "unresolved", "deposited"}
        if set(exchange) != required:
            raise ValueError("embedded error must cover all four exchange inventories")
        if (any(not isinstance(value, ArrayErrorNorm3D) for value in state.values())
                or any(
                    not isinstance(value, ArrayErrorNorm3D)
                    for values in exchange.values() for value in values.values())):
            raise TypeError("embedded error mappings require ArrayErrorNorm3D values")
        object.__setattr__(self, "state_increment", MappingProxyType(state))
        object.__setattr__(self, "exchange_inventory", MappingProxyType(exchange))


@dataclass(frozen=True)
class ProbabilityChange3D:
    by_species: Mapping[str, float]
    maximum: float

    def __post_init__(self):
        values = {str(name): float(value) for name, value in self.by_species.items()}
        if any(not np.isfinite(value) or value < 0.0 for value in values.values()):
            raise ValueError("invalid probability-change diagnostics")
        object.__setattr__(self, "by_species", MappingProxyType(values))


@dataclass(frozen=True)
class SurfaceCouplingTrialDiagnostic3D:
    start_time_s: float
    duration_s: float
    accepted: bool
    rejection_reason: str | None
    embedded_error: EmbeddedSurfaceError3D
    start_to_midpoint_probability_change: ProbabilityChange3D
    midpoint_to_final_probability_change: ProbabilityChange3D
    start_to_final_probability_change: ProbabilityChange3D
    maximum_path_probability_change: float
    maximum_radiosity_relative_balance_error: float


@dataclass(frozen=True)
class SurfaceRadiosityCouplingDiagnostics3D:
    cache_sha256: str
    operator_identity_sha256: str
    accepted_trial_count: int
    rejected_trial_count: int
    radiosity_solve_count: int
    minimum_accepted_chemistry_substep_s: float
    maximum_accepted_chemistry_substep_s: float
    maximum_accepted_embedded_relative_error: float
    maximum_accepted_path_probability_change: float
    maximum_radiosity_relative_balance_error: float
    rejected_trial_exchange_contribution_is_zero: bool
    trials: tuple[SurfaceCouplingTrialDiagnostic3D, ...]


@dataclass(frozen=True)
class SurfaceRadiosityCouplingResult3D:
    state: MaterialSurfaceState3D
    material_exchange: SurfaceMaterialExchange
    integrated_recession_m: np.ndarray
    integrated_growth_m: np.ndarray
    diagnostics: SurfaceRadiosityCouplingDiagnostics3D

    def __post_init__(self):
        if (not isinstance(self.state, MaterialSurfaceState3D)
                or not isinstance(self.material_exchange, SurfaceMaterialExchange)
                or not isinstance(self.diagnostics, SurfaceRadiosityCouplingDiagnostics3D)):
            raise TypeError("invalid surface-radiosity coupling result")
        recession = np.asarray(self.integrated_recession_m, dtype=float).copy()
        growth = np.asarray(self.integrated_growth_m, dtype=float).copy()
        if (recession.ndim != 1 or growth.shape != recession.shape
                or np.any(~np.isfinite(recession)) or np.any(recession < 0.0)
                or np.any(~np.isfinite(growth)) or np.any(growth < 0.0)):
            raise ValueError("invalid integrated surface displacement")
        recession.setflags(write=False)
        growth.setflags(write=False)
        object.__setattr__(self, "integrated_recession_m", recession)
        object.__setattr__(self, "integrated_growth_m", growth)


def _probability_change(left, right):
    by_species = {}
    maximum = 0.0
    for name in sorted(set(left) | set(right)):
        first, second = np.broadcast_arrays(
            np.asarray(left.get(name, 0.0), dtype=float),
            np.asarray(right.get(name, 0.0), dtype=float))
        value = float(np.max(np.abs(second - first))) if first.size else 0.0
        by_species[name] = value
        maximum = max(maximum, value)
    return ProbabilityChange3D(by_species, maximum)


def _array_error(left, right, area):
    left, right = np.broadcast_arrays(
        np.asarray(left, dtype=float), np.asarray(right, dtype=float))
    if left.shape != area.shape:
        raise ValueError("embedded comparison requires one value per cached face")
    difference = np.abs(left - right)
    maximum_absolute = float(np.max(difference)) if difference.size else 0.0
    linf_scale = max(
        float(np.max(np.abs(left))) if left.size else 0.0,
        float(np.max(np.abs(right))) if right.size else 0.0,
        np.finfo(float).tiny)
    weighted_scale = max(
        float(np.dot(area, np.abs(left))),
        float(np.dot(area, np.abs(right))),
        np.finfo(float).tiny)
    return ArrayErrorNorm3D(
        maximum_absolute,
        maximum_absolute / linf_scale,
        float(np.dot(area, difference)) / weighted_scale)


_EXCHANGE_KINDS = {
    "removed": "removed_units_m2",
    "outgoing": "outgoing_units_m2",
    "unresolved": "unresolved_units_m2",
    "deposited": "deposited_units_m2",
}


def _exchange_dict(exchange):
    return {
        kind: {name: np.asarray(value, dtype=float)
               for name, value in getattr(exchange, attribute).items()}
        for kind, attribute in _EXCHANGE_KINDS.items()
    }


def _sum_exchange_dicts(*exchanges):
    output = {kind: {} for kind in _EXCHANGE_KINDS}
    for exchange in exchanges:
        values = exchange if isinstance(exchange, dict) else _exchange_dict(exchange)
        for kind in _EXCHANGE_KINDS:
            for name, supplied in values[kind].items():
                value = np.asarray(supplied, dtype=float)
                if name not in output[kind]:
                    output[kind][name] = np.zeros_like(value, dtype=float)
                output[kind][name] += value
    return output


def _exchange_object(values, limitations=()):
    return SurfaceMaterialExchange(
        values["removed"], values["outgoing"], values["unresolved"],
        values["deposited"], tuple(limitations))


def _maximum_exact_ledger_residual(exchange):
    maximum = 0.0
    for name, removed in exchange.removed_units_m2.items():
        source, outgoing, unresolved = np.broadcast_arrays(
            np.asarray(removed, dtype=float),
            np.asarray(exchange.outgoing_units_m2.get(name, 0.0), dtype=float),
            np.asarray(exchange.unresolved_units_m2.get(name, 0.0), dtype=float))
        residual = source - outgoing - unresolved
        maximum = max(
            maximum,
            float(np.max(np.abs(residual))) if residual.size else 0.0)
    return maximum


def _embedded_error(
        initial_state, coarse, fine_final, coarse_exchange, fine_exchange,
        coarse_recession, fine_recession, coarse_growth, fine_growth, area):
    state_error = {}
    if (set(initial_state.fields) != set(coarse.state.fields)
            or set(initial_state.fields) != set(fine_final.state.fields)):
        raise SurfaceRadiosityCouplingRefusal(
            "material mechanism changed its conservative state-field contract")
    for name in sorted(initial_state.fields):
        initial = np.asarray(initial_state.fields[name], dtype=float)
        state_error[name] = _array_error(
            np.asarray(coarse.state.fields[name], dtype=float) - initial,
            np.asarray(fine_final.state.fields[name], dtype=float) - initial,
            area)
    exchange_error = {}
    for kind in _EXCHANGE_KINDS:
        exchange_error[kind] = {}
        names = sorted(set(coarse_exchange[kind]) | set(fine_exchange[kind]))
        for name in names:
            exchange_error[kind][name] = _array_error(
                coarse_exchange[kind].get(name, np.zeros_like(area)),
                fine_exchange[kind].get(name, np.zeros_like(area)), area)
    recession_error = _array_error(coarse_recession, fine_recession, area)
    growth_error = _array_error(coarse_growth, fine_growth, area)
    maximum = max(
        [value.maximum_relative_error for value in state_error.values()]
        + [value.maximum_relative_error
           for values in exchange_error.values() for value in values.values()]
        + [recession_error.maximum_relative_error, growth_error.maximum_relative_error],
        default=0.0)
    return EmbeddedSurfaceError3D(
        state_error, exchange_error, recession_error, growth_error, maximum)


def _advance_material(mechanism, state, fluxes, duration_s, material):
    result = mechanism.advance_by_material(state, fluxes, duration_s, material)
    if not isinstance(result, MaterialSurfaceStepResult3D):
        raise TypeError(
            "material-routed mechanism must return MaterialSurfaceStepResult3D")
    if result.product_populations:
        raise SurfaceRadiosityCouplingRefusal(
            "surface-radiosity coupling v1 refuses emitted product populations")
    if not result.validity.within_declared_scope:
        raise SurfaceRadiosityCouplingRefusal(
            "material mechanism left its declared scope: "
            + "; ".join(result.validity.reasons))
    if _maximum_exact_ledger_residual(result.material_exchange) != 0.0:
        raise SurfaceRadiosityCouplingRefusal(
            "material mechanism exchange ledger is not exactly closed")
    return result


def _input_state_snapshot(state):
    return {
        name: np.asarray(value).copy()
        for name, value in state.fields.items()
    }


def _assert_input_state_unchanged(state, snapshot):
    if (set(state.fields) != set(snapshot)
            or any(not np.array_equal(state.fields[name], value)
                   for name, value in snapshot.items())):
        raise RuntimeError("surface-radiosity coupling mutated its input state")


def integrate_surface_radiosity_chemistry_3d(
        cache: SurfaceRadiosityOperatorCache3D, mechanism,
        state: MaterialSurfaceState3D, duration_s: float, *,
        operator_identity_payload: Mapping[str, object],
        maximum_embedded_relative_error: float,
        maximum_local_reaction_probability_change: float,
        minimum_chemistry_substep_s: float,
        maximum_trial_count: int = 100000):
    """Advance chemistry on a fixed mesh with embedded radiosity feedback.

    Each trial compares one full chemistry step using start-state radiosity with
    two half steps whose second half uses midpoint radiosity.  Only the two-half
    path is accepted.  Candidate steps are halved deterministically until both
    the embedded error and the pathwise reaction-probability cap pass.
    """
    if not isinstance(cache, SurfaceRadiosityOperatorCache3D):
        raise TypeError("cache must be SurfaceRadiosityOperatorCache3D")
    if not isinstance(state, MaterialSurfaceState3D):
        raise TypeError("state must be MaterialSurfaceState3D")
    cache.assert_identity(operator_identity_payload)
    duration = float(duration_s)
    error_limit = float(maximum_embedded_relative_error)
    probability_limit = float(maximum_local_reaction_probability_change)
    minimum = float(minimum_chemistry_substep_s)
    trial_limit = int(maximum_trial_count)
    if (not np.isfinite(duration) or duration < 0.0
            or not np.isfinite(error_limit) or error_limit <= 0.0
            or not np.isfinite(probability_limit) or probability_limit <= 0.0
            or not np.isfinite(minimum) or minimum <= 0.0
            or trial_limit <= 0):
        raise ValueError("invalid surface-radiosity integration controls")
    if any(np.asarray(value).shape != (cache.face_count,) for value in state.fields.values()):
        raise ValueError("surface state does not match cached active faces")
    snapshot = _input_state_snapshot(state)
    flux_signature = cache.cache_sha256
    empty_exchange = SurfaceMaterialExchange({}, {}, {}, {})
    zero = np.zeros(cache.face_count)
    if duration == 0.0:
        diagnostics = SurfaceRadiosityCouplingDiagnostics3D(
            cache.cache_sha256, cache.operator_identity_sha256, 0, 0, 0,
            0.0, 0.0, 0.0, 0.0, 0.0, True, ())
        _assert_input_state_unchanged(state, snapshot)
        return SurfaceRadiosityCouplingResult3D(
            state, empty_exchange, zero, zero, diagnostics)
    if 0.5 * duration < minimum:
        raise SurfaceRadiosityCouplingRefusal(
            "duration cannot supply two half steps above the declared minimum")

    current = state
    elapsed = 0.0
    candidate = duration
    accepted = 0
    rejected = 0
    solve_count = 0
    minimum_accepted = np.inf
    maximum_accepted = 0.0
    maximum_accepted_error = 0.0
    maximum_accepted_dp = 0.0
    maximum_balance = 0.0
    trials = []
    cumulative_values = {kind: {} for kind in _EXCHANGE_KINDS}
    limitations = []
    recession_total = np.zeros(cache.face_count)
    growth_total = np.zeros(cache.face_count)
    roundoff = max(16.0 * np.finfo(float).eps * duration, 1.0e-18)

    while duration - elapsed > roundoff:
        if len(trials) >= trial_limit:
            raise SurfaceRadiosityCouplingRefusal(
                "surface-radiosity coupling exceeded maximum_trial_count")
        remaining = duration - elapsed
        candidate = min(candidate, remaining)
        if 0.5 * candidate < minimum:
            raise SurfaceRadiosityCouplingRefusal(
                "remaining interval fell below two declared minimum chemistry substeps")

        start_eval = cache.solve(
            current, mechanism, operator_identity_payload=operator_identity_payload)
        solve_count += 1
        coarse = _advance_material(
            mechanism, current, start_eval.surface_fluxes, candidate,
            cache.face_material_id)
        half = 0.5 * candidate
        fine_first = _advance_material(
            mechanism, current, start_eval.surface_fluxes, half,
            cache.face_material_id)
        midpoint_eval = cache.solve(
            fine_first.state, mechanism,
            operator_identity_payload=operator_identity_payload)
        solve_count += 1
        fine_second = _advance_material(
            mechanism, fine_first.state, midpoint_eval.surface_fluxes, half,
            cache.face_material_id)
        final_probability = cache._reaction_probability(fine_second.state, mechanism)

        start_midpoint = _probability_change(
            start_eval.reaction_probability, midpoint_eval.reaction_probability)
        midpoint_final = _probability_change(
            midpoint_eval.reaction_probability, final_probability)
        start_final = _probability_change(
            start_eval.reaction_probability, final_probability)
        path_dp = max(
            start_midpoint.maximum, midpoint_final.maximum, start_final.maximum)

        coarse_values = _exchange_dict(coarse.material_exchange)
        fine_values = _sum_exchange_dicts(
            fine_first.material_exchange, fine_second.material_exchange)
        fine_exchange = _exchange_object(
            fine_values,
            fine_first.material_exchange.known_limitations
            + fine_second.material_exchange.known_limitations)
        if _maximum_exact_ledger_residual(fine_exchange) != 0.0:
            raise SurfaceRadiosityCouplingRefusal(
                "combined fine-path exchange ledger is not exactly closed")
        coarse_recession = np.asarray(coarse.etch_velocity_m_s) * candidate
        coarse_growth = np.asarray(coarse.normal_growth_velocity_m_s) * candidate
        fine_recession = half * (
            np.asarray(fine_first.etch_velocity_m_s)
            + np.asarray(fine_second.etch_velocity_m_s))
        fine_growth = half * (
            np.asarray(fine_first.normal_growth_velocity_m_s)
            + np.asarray(fine_second.normal_growth_velocity_m_s))
        embedded = _embedded_error(
            current, coarse, fine_second, coarse_values, fine_values,
            coarse_recession, fine_recession, coarse_growth, fine_growth,
            cache.face_area_m2)
        balance = max(
            start_eval.maximum_relative_balance_error,
            midpoint_eval.maximum_relative_balance_error)
        maximum_balance = max(maximum_balance, balance)

        reason = None
        if path_dp > probability_limit:
            reason = "reaction_probability_safety_cap"
        elif embedded.maximum_relative_error > error_limit:
            reason = "embedded_step_error"
        accepted_trial = reason is None
        trials.append(SurfaceCouplingTrialDiagnostic3D(
            elapsed, candidate, accepted_trial, reason, embedded,
            start_midpoint, midpoint_final, start_final, path_dp, balance))
        if not accepted_trial:
            rejected += 1
            reduced = 0.5 * candidate
            if 0.5 * reduced < minimum:
                _assert_input_state_unchanged(state, snapshot)
                raise SurfaceRadiosityCouplingRefusal(
                    "coupling trial failed above the declared minimum chemistry substep: "
                    f"reason={reason}, candidate={candidate:.17g}, minimum={minimum:.17g}")
            candidate = reduced
            continue

        # This is the only mutation of accepted accumulators.  Coarse and rejected
        # fine trials are therefore incapable of leaking into the reported ledger.
        cumulative_values = _sum_exchange_dicts(cumulative_values, fine_values)
        limitations.extend(fine_exchange.known_limitations)
        cumulative_exchange = _exchange_object(
            cumulative_values, tuple(dict.fromkeys(limitations)))
        if _maximum_exact_ledger_residual(cumulative_exchange) != 0.0:
            raise SurfaceRadiosityCouplingRefusal(
                "accepted cumulative exchange ledger is not exactly closed")
        recession_total += fine_recession
        growth_total += fine_growth
        current = fine_second.state
        elapsed += candidate
        accepted += 1
        minimum_accepted = min(minimum_accepted, half)
        maximum_accepted = max(maximum_accepted, half)
        maximum_accepted_error = max(
            maximum_accepted_error, embedded.maximum_relative_error)
        maximum_accepted_dp = max(maximum_accepted_dp, path_dp)

    if abs(elapsed - duration) > roundoff:
        raise RuntimeError("surface-radiosity coupling did not cover the requested duration")
    _assert_input_state_unchanged(state, snapshot)
    if cache.cache_sha256 != flux_signature:
        raise RuntimeError("surface-radiosity coupling mutated its cached operator")
    cumulative_exchange = _exchange_object(
        cumulative_values, tuple(dict.fromkeys(limitations)))
    diagnostics = SurfaceRadiosityCouplingDiagnostics3D(
        cache.cache_sha256, cache.operator_identity_sha256,
        accepted, rejected, solve_count, float(minimum_accepted),
        maximum_accepted, maximum_accepted_error, maximum_accepted_dp,
        maximum_balance, True, tuple(trials))
    return SurfaceRadiosityCouplingResult3D(
        current, cumulative_exchange, recession_total, growth_total, diagnostics)
