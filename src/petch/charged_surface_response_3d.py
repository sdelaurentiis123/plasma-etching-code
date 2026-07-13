"""Conservative accounting contract for charged-particle surface response in 3-D.

This module does not choose reflection or emission physics.  It defines the rate measure emitted by
such a model and converts incident/outgoing charged events into signed deposited surface current.
Keeping this identity centralized prevents reflection, true-secondary emission, and future material
models from each inventing their own charge convention.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Protocol

import numpy as np
from scipy.stats import qmc

from .sheath import ECHARGE
from .surface_kinetics import FaceResolvedEnergeticFlux, ParameterEvidence


@dataclass(frozen=True)
class ChargedSurfaceContext3D:
    """Geometry/material fields presented unchanged to every charged response model."""

    face_area_m2: np.ndarray
    face_gas_normal: np.ndarray
    face_material_id: np.ndarray
    material_state: object | None = None

    def __post_init__(self):
        area = np.asarray(self.face_area_m2, dtype=float).copy()
        normal = np.asarray(self.face_gas_normal, dtype=float).copy()
        material = np.asarray(self.face_material_id).copy()
        if (area.ndim != 1 or area.size == 0 or normal.shape != (area.size, 3)
                or material.shape != area.shape or np.any(~np.isfinite(area))
                or np.any(area <= 0.0) or np.any(~np.isfinite(normal))
                or not np.allclose(np.linalg.norm(normal, axis=1), 1.0, rtol=0.0, atol=2e-6)):
            raise ValueError("invalid charged-surface geometry or material fields")
        for value in (area, normal, material):
            value.setflags(write=False)
        object.__setattr__(self, "face_area_m2", area)
        object.__setattr__(self, "face_gas_normal", normal)
        object.__setattr__(self, "face_material_id", material)


class ChargedSurfaceResponse3D(Protocol):
    """Material response contract shared by transient, PTC, and steady current audits."""

    def evaluate(
            self, incident_populations, charge_number_by_species,
            context: ChargedSurfaceContext3D) -> "ChargedSurfaceTransfer3D":
        ...


@dataclass(frozen=True)
class OutgoingChargedParticleEvents3D:
    """Sparse charged particles launched from surface impacts.

    ``event_rate_s`` is a particle rate, not a flux density. ``event_position`` uses the same mesh
    coordinate system as the incident triangle mesh. ``event_velocity_sqrt_eV`` is the engine's
    energy-coordinate velocity vector and must point from the surface into the gas when produced by
    a response model; that geometric condition is checked by the later launch operator, which owns
    the oriented face normals.
    """

    name: str
    charge_number: int
    face_count: int
    source_face: np.ndarray
    event_rate_s: np.ndarray
    event_position: np.ndarray
    event_velocity_sqrt_eV: np.ndarray

    def __post_init__(self):
        if (not self.name or int(self.charge_number) != self.charge_number
                or int(self.charge_number) == 0 or int(self.face_count) <= 0):
            raise ValueError("outgoing charged events require a name, charge, and surface")
        face = np.asarray(self.source_face, dtype=int).copy()
        rate = np.asarray(self.event_rate_s, dtype=float).copy()
        position = np.asarray(self.event_position, dtype=float).copy()
        velocity = np.asarray(self.event_velocity_sqrt_eV, dtype=float).copy()
        if (face.ndim != 1 or rate.shape != face.shape
                or position.shape != (face.size, 3) or velocity.shape != (face.size, 3)
                or np.any(face < 0) or np.any(face >= int(self.face_count))
                or np.any(~np.isfinite(rate)) or np.any(rate < 0.0)
                or np.any(~np.isfinite(position)) or np.any(~np.isfinite(velocity))
                or (velocity.size and np.any(np.linalg.norm(velocity, axis=1) <= 0.0))):
            raise ValueError("invalid outgoing charged-particle event measure")
        for value in (face, rate, position, velocity):
            value.setflags(write=False)
        object.__setattr__(self, "charge_number", int(self.charge_number))
        object.__setattr__(self, "face_count", int(self.face_count))
        object.__setattr__(self, "source_face", face)
        object.__setattr__(self, "event_rate_s", rate)
        object.__setattr__(self, "event_position", position)
        object.__setattr__(self, "event_velocity_sqrt_eV", velocity)

    @property
    def charge_rate_c_s(self):
        return float(ECHARGE * self.charge_number * np.sum(self.event_rate_s))


@dataclass(frozen=True)
class ChargedSurfaceTransfer3D:
    """Charge-conservative surface transfer before outgoing particles are re-transported."""

    positive_deposition_current_density_a_m2: np.ndarray
    negative_deposition_current_density_a_m2: np.ndarray
    face_current_density_a_m2: np.ndarray
    outgoing: tuple[OutgoingChargedParticleEvents3D, ...]
    incident_charge_rate_c_s: float
    outgoing_charge_rate_c_s: float
    deposited_charge_rate_c_s: float
    charge_balance_residual_c_s: float
    relative_charge_balance_error: float
    incident_kinetic_energy_rate_eV_s: float
    outgoing_kinetic_energy_rate_eV_s: float
    deposited_kinetic_energy_rate_eV_s: float
    kinetic_energy_balance_residual_eV_s: float
    relative_kinetic_energy_balance_error: float

    def __post_init__(self):
        positive = np.asarray(
            self.positive_deposition_current_density_a_m2, dtype=float).copy()
        negative = np.asarray(
            self.negative_deposition_current_density_a_m2, dtype=float).copy()
        signed = np.asarray(self.face_current_density_a_m2, dtype=float).copy()
        if (positive.ndim != 1 or negative.shape != positive.shape or signed.shape != positive.shape
                or np.any(~np.isfinite(positive)) or np.any(positive < 0.0)
                or np.any(~np.isfinite(negative)) or np.any(negative < 0.0)
                or np.any(~np.isfinite(signed))
                or not np.allclose(signed, positive - negative, rtol=1e-14, atol=0.0)):
            raise ValueError("invalid charged surface-transfer currents")
        outgoing = tuple(self.outgoing)
        if any(not isinstance(item, OutgoingChargedParticleEvents3D) for item in outgoing):
            raise TypeError("outgoing surface populations must contain charged events")
        rates = np.asarray([
            self.incident_charge_rate_c_s, self.outgoing_charge_rate_c_s,
            self.deposited_charge_rate_c_s, self.charge_balance_residual_c_s,
            self.relative_charge_balance_error], dtype=float)
        if (np.any(~np.isfinite(rates)) or self.relative_charge_balance_error < 0.0
                or any(item.face_count != positive.size for item in outgoing)):
            raise ValueError("invalid charged surface-transfer diagnostics")
        energy = np.asarray([
            self.incident_kinetic_energy_rate_eV_s,
            self.outgoing_kinetic_energy_rate_eV_s,
            self.deposited_kinetic_energy_rate_eV_s,
            self.kinetic_energy_balance_residual_eV_s,
            self.relative_kinetic_energy_balance_error], dtype=float)
        if (np.any(~np.isfinite(energy)) or self.incident_kinetic_energy_rate_eV_s < 0.0
                or self.outgoing_kinetic_energy_rate_eV_s < 0.0
                or self.relative_kinetic_energy_balance_error < 0.0):
            raise ValueError("invalid charged surface-transfer kinetic-energy ledger")
        for value in (positive, negative, signed):
            value.setflags(write=False)
        object.__setattr__(self, "positive_deposition_current_density_a_m2", positive)
        object.__setattr__(self, "negative_deposition_current_density_a_m2", negative)
        object.__setattr__(self, "face_current_density_a_m2", signed)
        object.__setattr__(self, "outgoing", outgoing)


def account_charged_surface_transfer_3d(
        incident_populations, charge_number_by_species, face_area_m2, *, outgoing=()):
    """Apply ``deposited = incident - outgoing`` to one surface event measure.

    Outgoing charge is accounted at its source face. Its later landing is a new incident population
    and must be passed through this same function again. The positive/negative arrays retain charge
    throughput separately, while ``face_current_density_a_m2`` is the physical net storage current.
    """
    incident = tuple(incident_populations)
    if (not incident or any(not isinstance(item, FaceResolvedEnergeticFlux) for item in incident)):
        raise ValueError("one or more face-resolved incident populations are required")
    face_count = incident[0].face_count
    if any(item.face_count != face_count for item in incident):
        raise ValueError("incident populations must share one surface mesh")
    names = [item.name for item in incident]
    if len(set(names)) != len(names):
        raise ValueError("incident charged species names must be unique")
    charge = dict(charge_number_by_species)
    if (set(charge) != set(names)
            or any(int(value) != value or int(value) == 0 for value in charge.values())):
        raise ValueError("charge_number_by_species must classify every incident population")
    area = np.asarray(face_area_m2, dtype=float)
    if (area.shape != (face_count,) or np.any(~np.isfinite(area)) or np.any(area <= 0.0)):
        raise ValueError("face_area_m2 must contain one positive physical area per face")
    outgoing = tuple(outgoing)
    if any(not isinstance(item, OutgoingChargedParticleEvents3D) for item in outgoing):
        raise TypeError("outgoing must contain OutgoingChargedParticleEvents3D")
    if any(item.face_count != face_count for item in outgoing):
        raise ValueError("outgoing populations must use the incident surface mesh")

    positive = np.zeros(face_count)
    negative = np.zeros(face_count)
    incident_charge_rate = 0.0
    incident_absolute_charge_rate = 0.0
    incident_energy_rate = 0.0
    for population in incident:
        event_current_density = (
            ECHARGE * int(charge[population.name]) * population.event_flux_m2_s)
        positive += np.bincount(
            population.event_face, weights=np.maximum(event_current_density, 0.0),
            minlength=face_count)
        negative += np.bincount(
            population.event_face, weights=np.maximum(-event_current_density, 0.0),
            minlength=face_count)
        incident_charge_rate += float(np.dot(event_current_density, area[population.event_face]))
        incident_absolute_charge_rate += float(np.dot(
            np.abs(event_current_density), area[population.event_face]))
        incident_energy_rate += float(np.dot(
            population.event_flux_m2_s * population.event_energy_eV,
            area[population.event_face]))

    outgoing_charge_rate = 0.0
    outgoing_absolute_charge_rate = 0.0
    outgoing_energy_rate = 0.0
    for population in outgoing:
        # Outgoing current leaves the surface, so its contribution to stored charge has the
        # opposite sign. Convert particle rate back to source-face current density exactly once.
        event_deposition_density = (
            -ECHARGE * population.charge_number * population.event_rate_s
            / area[population.source_face])
        positive += np.bincount(
            population.source_face, weights=np.maximum(event_deposition_density, 0.0),
            minlength=face_count)
        negative += np.bincount(
            population.source_face, weights=np.maximum(-event_deposition_density, 0.0),
            minlength=face_count)
        outgoing_charge_rate += population.charge_rate_c_s
        outgoing_absolute_charge_rate += float(
            ECHARGE * abs(population.charge_number) * np.sum(population.event_rate_s))
        outgoing_energy_rate += float(np.dot(
            population.event_rate_s,
            np.einsum("rc,rc->r", population.event_velocity_sqrt_eV,
                      population.event_velocity_sqrt_eV)))

    signed = positive - negative
    deposited_charge_rate = float(np.dot(signed, area))
    residual = incident_charge_rate - outgoing_charge_rate - deposited_charge_rate
    scale = max(
        incident_absolute_charge_rate, outgoing_absolute_charge_rate,
        abs(deposited_charge_rate), np.finfo(float).tiny)
    deposited_energy_rate = incident_energy_rate - outgoing_energy_rate
    energy_residual = (
        incident_energy_rate - outgoing_energy_rate - deposited_energy_rate)
    energy_scale = max(
        incident_energy_rate, outgoing_energy_rate, abs(deposited_energy_rate),
        np.finfo(float).tiny)
    return ChargedSurfaceTransfer3D(
        positive, negative, signed, outgoing,
        incident_charge_rate_c_s=incident_charge_rate,
        outgoing_charge_rate_c_s=outgoing_charge_rate,
        deposited_charge_rate_c_s=deposited_charge_rate,
        charge_balance_residual_c_s=float(residual),
        relative_charge_balance_error=float(abs(residual) / scale),
        incident_kinetic_energy_rate_eV_s=incident_energy_rate,
        outgoing_kinetic_energy_rate_eV_s=outgoing_energy_rate,
        deposited_kinetic_energy_rate_eV_s=deposited_energy_rate,
        kinetic_energy_balance_residual_eV_s=float(energy_residual),
        relative_kinetic_energy_balance_error=float(abs(energy_residual) / energy_scale))


def perfect_absorber_surface_transfer_3d(
        incident_populations, charge_number_by_species, face_area_m2):
    """Exact material-response identity for the engine's historical absorbing current law.

    This deliberately accumulates face flux before multiplying by charge, matching the original
    charging operator's floating-point order.  The general event accountant instead works event by
    event because an emitting response must preserve sparse outgoing-particle provenance.
    """
    incident = tuple(incident_populations)
    if (not incident or any(not isinstance(item, FaceResolvedEnergeticFlux) for item in incident)):
        raise ValueError("one or more face-resolved incident populations are required")
    face_count = incident[0].face_count
    if any(item.face_count != face_count for item in incident):
        raise ValueError("incident populations must share one surface mesh")
    names = [item.name for item in incident]
    if len(set(names)) != len(names):
        raise ValueError("incident charged species names must be unique")
    charge = dict(charge_number_by_species)
    if (set(charge) != set(names)
            or any(int(value) != value or int(value) == 0 for value in charge.values())):
        raise ValueError("charge_number_by_species must classify every incident population")
    area = np.asarray(face_area_m2, dtype=float)
    if (area.shape != (face_count,) or np.any(~np.isfinite(area)) or np.any(area <= 0.0)):
        raise ValueError("face_area_m2 must contain one positive physical area per face")

    positive = np.zeros(face_count)
    negative = np.zeros(face_count)
    for population in incident:
        current = ECHARGE * abs(float(charge[population.name])) * population.flux_m2_s
        if charge[population.name] > 0:
            positive += current
        else:
            negative += current
    signed = positive - negative
    deposited_charge_rate = float(np.dot(signed, area))
    incident_energy_rate = float(sum(
        np.dot(population.event_flux_m2_s * population.event_energy_eV,
               area[population.event_face])
        for population in incident))
    return ChargedSurfaceTransfer3D(
        positive, negative, signed, (),
        incident_charge_rate_c_s=deposited_charge_rate,
        outgoing_charge_rate_c_s=0.0,
        deposited_charge_rate_c_s=deposited_charge_rate,
        charge_balance_residual_c_s=0.0,
        relative_charge_balance_error=0.0,
        incident_kinetic_energy_rate_eV_s=incident_energy_rate,
        outgoing_kinetic_energy_rate_eV_s=0.0,
        deposited_kinetic_energy_rate_eV_s=incident_energy_rate,
        kinetic_energy_balance_residual_eV_s=0.0,
        relative_kinetic_energy_balance_error=0.0)


@dataclass(frozen=True)
class PerfectAbsorberChargedSurfaceResponse3D:
    """Material-independent identity response used by the historical charging operator."""

    def evaluate(
            self, incident_populations, charge_number_by_species,
            context: ChargedSurfaceContext3D):
        if not isinstance(context, ChargedSurfaceContext3D):
            raise TypeError("context must be ChargedSurfaceContext3D")
        return perfect_absorber_surface_transfer_3d(
            incident_populations, charge_number_by_species, context.face_area_m2)


@dataclass(frozen=True)
class GrazingSpecularIonReflection3D:
    """Declared v1 grazing-ion reflection sensitivity with exact weighted lineage.

    ``P_reflect = grazing_reflection_probability * (1 - cos(theta)**angular_exponent)``.
    The reflected event retains ``energy_retention_fraction`` of incident kinetic energy and uses
    the specular direction about the local gas normal.  The remaining weighted particles and energy
    are deposited.  This is a phenomenological, material-tagged sensitivity law: it is not promoted
    as a universal Si/SiO2 scattering distribution, and all three parameters carry provenance and
    declared bounds.
    """

    material_id: object
    ion_species_name: str
    grazing_reflection_probability: float
    angular_exponent: float
    energy_retention_fraction: float
    parameter_evidence: Mapping[str, ParameterEvidence]
    parameter_bounds: Mapping[str, tuple[float, float]]

    def __post_init__(self):
        parameters = {
            "grazing_reflection_probability": self.grazing_reflection_probability,
            "angular_exponent": self.angular_exponent,
            "energy_retention_fraction": self.energy_retention_fraction}
        evidence = dict(self.parameter_evidence)
        bounds = {name: tuple(value) for name, value in self.parameter_bounds.items()}
        if (not self.ion_species_name or set(evidence) != set(parameters)
                or set(bounds) != set(parameters)
                or any(not isinstance(value, ParameterEvidence) for value in evidence.values())
                or any(len(value) != 2 or not np.all(np.isfinite(value))
                       or value[0] > value[1] for value in bounds.values())
                or not 0.0 <= self.grazing_reflection_probability <= 1.0
                or not np.isfinite(self.angular_exponent) or self.angular_exponent <= 0.0
                or not 0.0 < self.energy_retention_fraction <= 1.0
                or any(not bounds[name][0] <= value <= bounds[name][1]
                       for name, value in parameters.items())):
            raise ValueError("invalid grazing specular-ion reflection model")
        object.__setattr__(self, "parameter_evidence", MappingProxyType(evidence))
        object.__setattr__(self, "parameter_bounds", MappingProxyType(bounds))

    @classmethod
    def literature_bounded_sensitivity(cls, material_id, ion_species_name="Ar+"):
        """Central v1 point inside published grazing-reflection observations."""
        sources = (
            "Helmer and Graves, JVST A 16, 3502 (1998), DOI 10.1116/1.580993; "
            "Hoekstra et al., JVST B 16, 2102 (1998); "
            "Du et al., JVST A 40, 053007 (2022), DOI 10.1116/6.0002008")
        evidence = {
            name: ParameterEvidence(
                sources, "literature_bounded_phenomenological_sensitivity",
                note=(
                    "Bounds span reported >90% grazing reflection, up-to-99% retained energy, "
                    "and modern MD angle/energy/material dependence; not a calibrated material law"))
            for name in (
                "grazing_reflection_probability", "angular_exponent",
                "energy_retention_fraction")}
        return cls(
            material_id=material_id, ion_species_name=ion_species_name,
            grazing_reflection_probability=0.95, angular_exponent=3.0,
            energy_retention_fraction=0.90, parameter_evidence=evidence,
            parameter_bounds={
                "grazing_reflection_probability": (0.80, 1.0),
                "angular_exponent": (2.0, 8.0),
                "energy_retention_fraction": (0.50, 0.99)})

    @property
    def provenance(self):
        return MappingProxyType(dict(
            model="weighted angle-dependent specular reflection with constant energy retention",
            material_id=self.material_id,
            ion_species_name=self.ion_species_name,
            parameters={
                "grazing_reflection_probability": self.grazing_reflection_probability,
                "angular_exponent": self.angular_exponent,
                "energy_retention_fraction": self.energy_retention_fraction},
            bounds=dict(self.parameter_bounds),
            evidence={name: value.source for name, value in self.parameter_evidence.items()},
            claim="literature-bounded sensitivity; not material-calibrated prediction"))

    def reflection_probability(self, cosine_incidence):
        cosine = np.asarray(cosine_incidence, dtype=float)
        if (np.any(~np.isfinite(cosine))
                or np.any((cosine < 0.0) | (cosine > 1.0))):
            raise ValueError("incidence cosine must lie in [0, 1]")
        return self.grazing_reflection_probability * (
            1.0 - cosine ** self.angular_exponent)

    def evaluate(
            self, incident_populations, charge_number_by_species,
            context: ChargedSurfaceContext3D):
        incident = tuple(incident_populations)
        charge = dict(charge_number_by_species)
        if not isinstance(context, ChargedSurfaceContext3D):
            raise TypeError("context must be ChargedSurfaceContext3D")
        selected_population = next(
            (item for item in incident if item.name == self.ion_species_name), None)
        outgoing = ()
        if selected_population is not None:
            ion_charge = charge.get(self.ion_species_name)
            if ion_charge is None or ion_charge <= 0:
                raise ValueError("ion reflection requires a positively charged incident species")
            if (selected_population.event_position is None
                    or selected_population.event_incident_direction is None):
                raise ValueError("ion reflection requires impact position and direction lineage")
            face = selected_population.event_face
            normal = context.face_gas_normal[face]
            geometric_cosine = -np.einsum(
                "rc,rc->r", selected_population.event_incident_direction, normal)
            cosine_difference = np.abs(
                geometric_cosine - selected_population.event_cosine_incidence)
            invalid_cosine = (geometric_cosine < -2e-6) | (cosine_difference > 2e-5)
            if np.any(invalid_cosine):
                event = int(np.flatnonzero(invalid_cosine)[
                    np.argmax(cosine_difference[invalid_cosine])])
                raise ValueError(
                    "ion impact cosine is inconsistent with direction and gas normal: "
                    f"event={event}, face={int(face[event])}, "
                    f"stored={selected_population.event_cosine_incidence[event]:.9g}, "
                    f"geometric={geometric_cosine[event]:.9g}, "
                    f"difference={cosine_difference[event]:.9g}, "
                    f"direction={selected_population.event_incident_direction[event].tolist()}, "
                    f"normal={normal[event].tolist()}")
            selected = context.face_material_id[face] == self.material_id
            probability = self.reflection_probability(
                selected_population.event_cosine_incidence)
            selected &= probability > 0.0
            if np.any(selected):
                normal = normal[selected]
                incident_direction = selected_population.event_incident_direction[selected]
                specular = incident_direction - 2.0 * np.einsum(
                    "rc,rc->r", incident_direction, normal)[:, None] * normal
                specular /= np.linalg.norm(specular, axis=1)[:, None]
                if np.any(np.einsum("rc,rc->r", specular, normal) <= 0.0):
                    raise ValueError("specular ion response did not point into the gas")
                event_rate = (
                    selected_population.event_flux_m2_s[selected]
                    * context.face_area_m2[face[selected]] * probability[selected])
                outgoing = (OutgoingChargedParticleEvents3D(
                    self.ion_species_name, int(ion_charge), len(context.face_area_m2),
                    face[selected], event_rate,
                    selected_population.event_position[selected],
                    np.sqrt(
                        self.energy_retention_fraction
                        * selected_population.event_energy_eV[selected])[:, None] * specular),)
        if not outgoing:
            return perfect_absorber_surface_transfer_3d(
                incident, charge, context.face_area_m2)
        return account_charged_surface_transfer_3d(
            incident, charge, context.face_area_m2, outgoing=outgoing)


def sobolewski_2021_ar_kinetic_see_yield(impact_energy_eV):
    """Recommended kinetic Ar+-induced electron yield for plasma-exposed SiO2.

    This is Eq. (8) of Sobolewski, Plasma Sources Sci. Technol. 30 025004 (2021),
    DOI 10.1088/1361-6595/abd61f.  It excludes the paper's separately recommended photon,
    metastable, and ion-potential emission terms.
    """
    energy = np.asarray(impact_energy_eV, dtype=float)
    if np.any(~np.isfinite(energy)) or np.any(energy < 0.0):
        raise ValueError("impact energy must be finite and nonnegative")
    return 0.030 * energy ** 2 / (200.0 + energy) ** 1.5


@dataclass(frozen=True)
class Sobolewski2021ArKineticSEE3D:
    """Material-tagged Ar+ kinetic SEE with deterministic Lambertian quadrature.

    The yield is sourced, but the cited yield paper does not determine a full emitted-electron
    spectrum.  ``emission_energy_eV`` and its evidence are therefore mandatory declared inputs.
    Huang and Kushner (JVST A 44, 023013, 2026) support a Lambertian angular law and an average
    energy of a few eV for ion-induced secondaries; callers must refine the declared energy and
    angular quadrature instead of treating either as a fitted charging knob.
    """

    material_id: object
    emission_energy_eV: float
    emission_energy_evidence: str
    ion_species_name: str = "Ar+"
    emitted_species_name: str = "secondary_electron"
    angular_log2_samples: int = 3
    angular_seed: int = 0
    minimum_supported_impact_energy_eV: float = 10.0
    maximum_supported_impact_energy_eV: float = 1.0e4

    def __post_init__(self):
        if (not self.ion_species_name or not self.emitted_species_name
                or self.ion_species_name == self.emitted_species_name
                or not np.isfinite(self.emission_energy_eV) or self.emission_energy_eV <= 0.0
                or not self.emission_energy_evidence
                or int(self.angular_log2_samples) != self.angular_log2_samples
                or self.angular_log2_samples < 0
                or int(self.angular_seed) != self.angular_seed
                or not np.isfinite(self.minimum_supported_impact_energy_eV)
                or not np.isfinite(self.maximum_supported_impact_energy_eV)
                or self.minimum_supported_impact_energy_eV < 0.0
                or self.maximum_supported_impact_energy_eV
                <= self.minimum_supported_impact_energy_eV):
            raise ValueError("invalid Sobolewski-2021 Ar+ kinetic-SEE response")
        object.__setattr__(self, "angular_log2_samples", int(self.angular_log2_samples))
        object.__setattr__(self, "angular_seed", int(self.angular_seed))

    @property
    def provenance(self):
        return MappingProxyType(dict(
            yield_source=(
                "Sobolewski, Plasma Sources Sci. Technol. 30 025004 (2021), Eq. 8; "
                "DOI 10.1088/1361-6595/abd61f"),
            yield_scope=(
                "kinetic Ar+-induced electron emission from plasma-exposed SiO2 only; "
                "photon, metastable, fast-neutral, and ion-potential channels excluded"),
            yield_fit_range="10 eV to 10 keV; values below 50 eV derive from literature, not the in-situ measurement",
            angular_source=(
                "Huang and Kushner, J. Vac. Sci. Technol. A 44, 023013 (2026), "
                "DOI 10.1116/6.0005187; Lambertian ion-induced secondary emission"),
            emission_energy_source=self.emission_energy_evidence,
        ))

    def evaluate(
            self, incident_populations, charge_number_by_species,
            context: ChargedSurfaceContext3D):
        incident = tuple(incident_populations)
        charge = dict(charge_number_by_species)
        if not isinstance(context, ChargedSurfaceContext3D):
            raise TypeError("context must be ChargedSurfaceContext3D")
        selected_population = next(
            (item for item in incident if item.name == self.ion_species_name), None)
        outgoing = ()
        if selected_population is not None:
            if charge.get(self.ion_species_name) != 1:
                raise ValueError("Sobolewski Ar+ response requires a singly charged positive ion")
            if (selected_population.event_position is None
                    or selected_population.event_incident_direction is None):
                raise ValueError("Ar+ SEE requires preserved impact position and direction")
            material_selected = (
                context.face_material_id[selected_population.event_face] == self.material_id)
            emitted_event = np.flatnonzero(material_selected)
            if emitted_event.size:
                impact_energy = selected_population.event_energy_eV[emitted_event]
                if (np.any(impact_energy < self.minimum_supported_impact_energy_eV)
                        or np.any(impact_energy > self.maximum_supported_impact_energy_eV)):
                    raise ValueError("Ar+ impact lies outside the declared Sobolewski fit range")
                emission_yield = sobolewski_2021_ar_kinetic_see_yield(impact_energy)
                if np.any(emission_yield * self.emission_energy_eV > impact_energy):
                    raise ValueError("declared emitted-electron energy exceeds incident energy budget")
                direction_count = 2 ** self.angular_log2_samples
                angular = qmc.Sobol(
                    2, scramble=True, seed=self.angular_seed).random_base2(
                        self.angular_log2_samples)
                cosine = np.sqrt(angular[:, 0])
                sine = np.sqrt(1.0 - angular[:, 0])
                base_azimuth = 2.0 * np.pi * angular[:, 1]
                face = selected_population.event_face[emitted_event]
                normal = context.face_gas_normal[face]
                reference = np.zeros_like(normal)
                use_z = np.abs(normal[:, 2]) < 0.9
                reference[use_z, 2] = 1.0
                reference[~use_z, 0] = 1.0
                tangent_a = np.cross(reference, normal)
                tangent_a /= np.linalg.norm(tangent_a, axis=1)[:, None]
                tangent_b = np.cross(normal, tangent_a)
                azimuth_rotation = 2.0 * np.pi * np.mod(
                    emitted_event * 0.6180339887498949, 1.0)
                azimuth = base_azimuth[None, :] + azimuth_rotation[:, None]
                direction = (
                    cosine[None, :, None] * normal[:, None, :]
                    + sine[None, :, None] * (
                        np.cos(azimuth)[:, :, None] * tangent_a[:, None, :]
                        + np.sin(azimuth)[:, :, None] * tangent_b[:, None, :]))
                incident_rate = (
                    selected_population.event_flux_m2_s[emitted_event]
                    * context.face_area_m2[face])
                event_rate = np.repeat(
                    incident_rate * emission_yield / direction_count, direction_count)
                outgoing = (OutgoingChargedParticleEvents3D(
                    self.emitted_species_name, -1, len(context.face_area_m2),
                    np.repeat(face, direction_count), event_rate,
                    np.repeat(
                        selected_population.event_position[emitted_event],
                        direction_count, axis=0),
                    np.sqrt(self.emission_energy_eV) * direction.reshape(-1, 3)),)
        if not outgoing:
            return perfect_absorber_surface_transfer_3d(
                incident, charge, context.face_area_m2)
        return account_charged_surface_transfer_3d(
            incident, charge, context.face_area_m2, outgoing=outgoing)
