"""Neutralized SiO2 forward scattering from Hwang--Giapis (1997).

This is separate from charged-particle reflection on purpose.  The incident Cl+
deposits its charge at SiO2, while the specularly scattered projectile is a fast
neutral Cl atom and therefore follows a straight, field-free second flight.
Treating that projectile as a charged reflection would change both the charging
ODE and the notch mechanism.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType

import numpy as np

from .boundary_transport_3d import (
    BoundaryTransport3DResult,
    _first_segment_triangle_hit_float64_3d,
    _inset_surface_launch_positions_3d,
)
from .charged_surface_response_3d import ChargedSurfaceContext3D
from .surface_kinetics import FaceResolvedEnergeticFlux, SurfaceFluxes


@dataclass(frozen=True)
class OutgoingNeutralParticleEvents3D:
    name: str
    face_count: int
    source_face: np.ndarray
    event_rate_s: np.ndarray
    event_position: np.ndarray
    event_velocity_sqrt_eV: np.ndarray

    def __post_init__(self):
        face = np.asarray(self.source_face, dtype=int).copy()
        rate = np.asarray(self.event_rate_s, dtype=float).copy()
        position = np.asarray(self.event_position, dtype=float).copy()
        velocity = np.asarray(self.event_velocity_sqrt_eV, dtype=float).copy()
        if (not self.name or int(self.face_count) != self.face_count
                or self.face_count <= 0 or face.ndim != 1
                or rate.shape != face.shape or position.shape != (face.size, 3)
                or velocity.shape != (face.size, 3)
                or np.any(face < 0) or np.any(face >= int(self.face_count))
                or np.any(~np.isfinite(rate)) or np.any(rate < 0.0)
                or np.any(~np.isfinite(position)) or np.any(~np.isfinite(velocity))
                or (velocity.size and np.any(np.linalg.norm(velocity, axis=1) <= 0.0))):
            raise ValueError("invalid outgoing neutral event measure")
        for value in (face, rate, position, velocity):
            value.setflags(write=False)
        object.__setattr__(self, "face_count", int(self.face_count))
        object.__setattr__(self, "source_face", face)
        object.__setattr__(self, "event_rate_s", rate)
        object.__setattr__(self, "event_position", position)
        object.__setattr__(self, "event_velocity_sqrt_eV", velocity)


@dataclass(frozen=True)
class NeutralSurfaceFlight3D:
    emitted: OutgoingNeutralParticleEvents3D
    incident: FaceResolvedEnergeticFlux
    termination: np.ndarray
    hit_face: np.ndarray
    emitted_rate_s: float
    landed_rate_s: float
    escaped_rate_s: float
    relative_particle_balance_error: float
    edge_launch_inset_count: int

    def __post_init__(self):
        termination = np.asarray(self.termination, dtype=np.int8).copy()
        hit_face = np.asarray(self.hit_face, dtype=int).copy()
        if (not isinstance(self.emitted, OutgoingNeutralParticleEvents3D)
                or not isinstance(self.incident, FaceResolvedEnergeticFlux)
                or self.incident.name != self.emitted.name
                or termination.shape != self.emitted.event_rate_s.shape
                or hit_face.shape != termination.shape
                or np.any(~np.isin(termination, (1, 2)))
                or np.any((termination == 1) & (hit_face < 0))
                or np.any((termination == 2) & (hit_face != -1))):
            raise ValueError("invalid neutral surface-flight lineage")
        classified = (
            float(np.sum(self.emitted.event_rate_s[termination == 1])),
            float(np.sum(self.emitted.event_rate_s[termination == 2])))
        rates = (float(self.emitted_rate_s), float(self.landed_rate_s),
                 float(self.escaped_rate_s))
        residual = rates[0] - rates[1] - rates[2]
        if (any(not np.isfinite(value) or value < 0.0 for value in rates)
                or not np.allclose(classified, rates[1:], rtol=5e-15, atol=0.0)
                or abs(residual) > 5e-15 * max(rates[0], np.finfo(float).tiny)
                or not np.isclose(
                    self.relative_particle_balance_error,
                    abs(residual) / max(rates[0], np.finfo(float).tiny),
                    rtol=0.0, atol=5e-16)
                or int(self.edge_launch_inset_count) != self.edge_launch_inset_count
                or not 0 <= self.edge_launch_inset_count <= len(termination)):
            raise ValueError("neutral surface flight does not close its particle ledger")
        termination.setflags(write=False); hit_face.setflags(write=False)
        object.__setattr__(self, "termination", termination)
        object.__setattr__(self, "hit_face", hit_face)
        object.__setattr__(self, "edge_launch_inset_count", int(self.edge_launch_inset_count))


@dataclass(frozen=True)
class HwangGiapisForwardScatter3DResult:
    emitted: OutgoingNeutralParticleEvents3D
    flight: NeutralSurfaceFlight3D
    eligible_incident_rate_s: float
    scattered_rate_s: float
    nonscattered_rate_s: float
    scattered_incident_energy_rate_eV_s: float
    retained_energy_rate_eV_s: float
    deposited_energy_rate_eV_s: float
    relative_surface_particle_balance_error: float
    relative_surface_energy_balance_error: float
    provenance: MappingProxyType


@dataclass(frozen=True)
class HwangGiapisSiO2ForwardScatter3D:
    """Eq. (4.2) probability and Eq. (4.3) hard-sphere energy loss."""

    material_id: int
    ion_species_name: str = "Cl+"
    neutral_species_name: str = "Cl_fast_neutral"
    critical_angle_deg: float = 45.0
    gas_to_effective_surface_mass_ratio: float = 1.0

    def __post_init__(self):
        if (int(self.material_id) != self.material_id or self.material_id <= 0
                or not self.ion_species_name or not self.neutral_species_name
                or self.ion_species_name == self.neutral_species_name
                or not np.isfinite(self.critical_angle_deg)
                or not 0.0 < self.critical_angle_deg < 90.0
                or not np.isfinite(self.gas_to_effective_surface_mass_ratio)
                or not 0.0 < self.gas_to_effective_surface_mass_ratio <= 1.0):
            raise ValueError("invalid Hwang--Giapis SiO2 forward-scatter inputs")
        object.__setattr__(self, "material_id", int(self.material_id))

    @property
    def provenance(self):
        return MappingProxyType({
            "model": "Hwang--Giapis neutralized SiO2 forward scattering",
            "material_id": self.material_id,
            "ion_species_name": self.ion_species_name,
            "neutral_species_name": self.neutral_species_name,
            "parameters": {
                "critical_angle_deg": self.critical_angle_deg,
                "gas_to_effective_surface_mass_ratio": (
                    self.gas_to_effective_surface_mass_ratio),
            },
            "bounds": {
                "critical_angle_deg": (30.0, 60.0),
                "gas_to_effective_surface_mass_ratio": (0.5, 1.0),
            },
            "evidence": {
                "critical_angle_deg": (
                    "Hwang & Giapis, JVST B 15, 70 (1997), DOI 10.1116/1.589258, "
                    "Eq. (4.2): 45 deg central value, 30--60 deg study"),
                "gas_to_effective_surface_mass_ratio": (
                    "Hwang & Giapis, JVST B 15, 70 (1997), DOI 10.1116/1.589258, "
                    "Eq. (4.3): m=1 central value; 0.5--1 described as realistic"),
            },
            "mass_ratio_note": "paper central value m=1.0 for Cl-covered SiO2",
            "charge_closure": "incident ion charge remains at SiO2; scattered projectile is neutral",
            "claim": "published benchmark model; parameters remain declared assumptions",
        })

    def scattering_probability(self, cosine_incidence):
        cosine = np.asarray(cosine_incidence, dtype=float)
        if np.any(~np.isfinite(cosine)) or np.any((cosine < 0.0) | (cosine > 1.0)):
            raise ValueError("incidence cosine must lie in [0,1]")
        angle = np.arccos(cosine)
        critical = np.deg2rad(self.critical_angle_deg)
        return np.clip((angle - critical) / (0.5 * np.pi - critical), 0.0, 1.0)

    def energy_retention_fraction(self, cosine_incidence):
        cosine = np.asarray(cosine_incidence, dtype=float)
        if np.any(~np.isfinite(cosine)) or np.any((cosine < 0.0) | (cosine > 1.0)):
            raise ValueError("incidence cosine must lie in [0,1]")
        angle = np.arccos(cosine)
        mass_ratio = float(self.gas_to_effective_surface_mass_ratio)
        sine_double = np.sin(2.0 * angle)
        cosine_double = np.cos(2.0 * angle)
        radical = np.sqrt(np.maximum(1.0 - mass_ratio ** 2 * sine_double ** 2, 0.0))
        loss = (2.0 * mass_ratio / (1.0 + mass_ratio) ** 2) * (
            1.0 + cosine_double * radical + mass_ratio * sine_double ** 2)
        return np.clip(1.0 - loss, 0.0, 1.0)

    def emit(self, incident_population, context: ChargedSurfaceContext3D):
        if (not isinstance(incident_population, FaceResolvedEnergeticFlux)
                or incident_population.name != self.ion_species_name
                or not isinstance(context, ChargedSurfaceContext3D)
                or incident_population.face_count != len(context.face_area_m2)
                or incident_population.event_position is None
                or incident_population.event_incident_direction is None):
            raise ValueError("forward scattering requires matching face-resolved Cl+ impact lineage")
        face = incident_population.event_face
        selected = context.face_material_id[face] == self.material_id
        probability = self.scattering_probability(
            incident_population.event_cosine_incidence)
        selected &= probability > 0.0
        source_face = face[selected]
        if not np.any(selected):
            return OutgoingNeutralParticleEvents3D(
                self.neutral_species_name, incident_population.face_count,
                np.empty(0, dtype=int), np.empty(0), np.empty((0, 3)), np.empty((0, 3)))
        normal = context.face_gas_normal[source_face]
        direction = incident_population.event_incident_direction[selected]
        incidence_cosine = incident_population.event_cosine_incidence[selected]
        specular = direction + 2.0 * incidence_cosine[:, None] * normal
        specular /= np.linalg.norm(specular, axis=1)[:, None]
        if np.any(np.einsum("rc,rc->r", specular, normal) < -2e-12):
            raise RuntimeError("neutralized forward-scatter launch points into the solid")
        retention = self.energy_retention_fraction(incidence_cosine)
        energy = retention * incident_population.event_energy_eV[selected]
        positive_energy = energy > 0.0
        selected_index = np.flatnonzero(selected)[positive_energy]
        source_face = face[selected_index]
        event_rate = (
            incident_population.event_flux_m2_s[selected_index]
            * context.face_area_m2[source_face] * probability[selected_index])
        return OutgoingNeutralParticleEvents3D(
            self.neutral_species_name, incident_population.face_count, source_face,
            event_rate, incident_population.event_position[selected_index],
            np.sqrt(energy[positive_energy])[:, None] * specular[positive_energy])


def trace_neutral_surface_events_3d(
        emitted: OutgoingNeutralParticleEvents3D, verts, faces, areas,
        face_gas_normals, *, domain_minimum, domain_maximum,
        mesh_length_unit_m=1e-6, launch_offset=1e-5,
        periodic_lateral=False, maximum_periodic_wraps=10000):
    """Trace exact straight neutral flights with hard float64 visibility.

    Lateral periodicity is handled by flight segments between cell boundaries;
    the ray wraps without changing direction or energy.  No finite time horizon
    is guessed.  ``maximum_periodic_wraps`` is only a non-escaping-cavity guard.
    """
    verts = np.asarray(verts, dtype=float)
    faces = np.asarray(faces, dtype=np.int64)
    areas = np.asarray(areas, dtype=float)
    normals = np.asarray(face_gas_normals, dtype=float)
    lower = np.asarray(domain_minimum, dtype=float)
    upper = np.asarray(domain_maximum, dtype=float)
    if (not isinstance(emitted, OutgoingNeutralParticleEvents3D)
            or verts.ndim != 2 or verts.shape[1] != 3
            or faces.ndim != 2 or faces.shape[1] != 3
            or areas.shape != (len(faces),) or normals.shape != (len(faces), 3)
            or emitted.face_count != len(faces)
            or lower.shape != (3,) or upper.shape != (3,) or np.any(upper <= lower)
            or np.any(~np.isfinite(verts)) or np.any(~np.isfinite(areas))
            or np.any(areas <= 0.0) or np.any(~np.isfinite(normals))
            or not np.allclose(np.linalg.norm(normals, axis=1), 1.0, rtol=0.0, atol=2e-6)
            or not np.isfinite(mesh_length_unit_m) or mesh_length_unit_m <= 0.0
            or not np.isfinite(launch_offset) or launch_offset <= 0.0
            or int(maximum_periodic_wraps) != maximum_periodic_wraps
            or maximum_periodic_wraps <= 0):
        raise ValueError("invalid neutral surface-flight mesh or controls")
    inset_position, inset_count = _inset_surface_launch_positions_3d(
        emitted.event_position, emitted.source_face, verts, faces, float(launch_offset))
    origin = inset_position + float(launch_offset) * normals[emitted.source_face]
    velocity = emitted.event_velocity_sqrt_eV
    direction = velocity / np.linalg.norm(velocity, axis=1)[:, None]
    if np.any(np.einsum("rc,rc->r", direction, normals[emitted.source_face]) < -2e-12):
        raise ValueError("outgoing neutral contains a solid-facing launch")
    hit_face = np.full(len(origin), -1, dtype=int)
    termination = np.full(len(origin), 2, dtype=np.int8)
    hit_position = np.zeros((len(origin), 3))
    epsilon = 2e-10 * float(np.max(upper - lower))
    for ray in range(len(origin)):
        position = origin[ray].copy()
        for _ in range(int(maximum_periodic_wraps) + 1):
            travel = np.full(3, np.inf)
            for axis in range(3):
                if direction[ray, axis] > 1e-15:
                    travel[axis] = (upper[axis] - position[axis]) / direction[ray, axis]
                elif direction[ray, axis] < -1e-15:
                    travel[axis] = (lower[axis] - position[axis]) / direction[ray, axis]
            positive = travel[travel > 1e-14]
            if not positive.size:
                raise RuntimeError("neutral ray has no forward domain exit")
            distance = float(np.min(positive))
            segment = direction[ray] * distance
            face, _, impact = _first_segment_triangle_hit_float64_3d(
                position, segment, verts, faces, upper - lower, False)
            if face >= 0:
                hit_face[ray] = int(face)
                hit_position[ray] = impact
                termination[ray] = 1
                break
            exit_axes = np.where(np.isclose(travel, distance, rtol=1e-11, atol=1e-13))[0]
            if (2 in exit_axes or not periodic_lateral
                    or any(axis not in (0, 1) for axis in exit_axes)):
                break
            position = position + direction[ray] * (distance + epsilon)
            for axis in (0, 1):
                length = upper[axis] - lower[axis]
                position[axis] = lower[axis] + (position[axis] - lower[axis]) % length
        else:
            raise RuntimeError(
                "neutral ray exceeded the periodic-wrap guard; possible non-escaping cavity")
    hit = termination == 1
    physical_area = areas * float(mesh_length_unit_m) ** 2
    hit_direction = direction[hit]
    cosine = -np.einsum("rc,rc->r", hit_direction, normals[hit_face[hit]])
    if np.any(cosine < -2e-8):
        raise RuntimeError("neutral flight hit a surface from its solid side")
    energy = np.sum(velocity[hit] ** 2, axis=1)
    incident = FaceResolvedEnergeticFlux(
        emitted.name, emitted.face_count, hit_face[hit],
        emitted.event_rate_s[hit] / physical_area[hit_face[hit]],
        energy, np.clip(cosine, 0.0, 1.0),
        event_position=hit_position[hit], event_incident_direction=hit_direction)
    emitted_rate = float(np.sum(emitted.event_rate_s))
    landed_rate = float(np.sum(emitted.event_rate_s[hit]))
    escaped_rate = float(np.sum(emitted.event_rate_s[~hit]))
    residual = emitted_rate - landed_rate - escaped_rate
    return NeutralSurfaceFlight3D(
        emitted, incident, termination, hit_face, emitted_rate, landed_rate,
        escaped_rate, abs(residual) / max(emitted_rate, np.finfo(float).tiny),
        inset_count)


def apply_hwang_giapis_forward_scatter_to_transport_3d(
        transport: BoundaryTransport3DResult, model: HwangGiapisSiO2ForwardScatter3D,
        context: ChargedSurfaceContext3D, verts, faces, areas, *,
        domain_minimum, domain_maximum, mesh_length_unit_m=1e-6,
        launch_offset=1e-5, periodic_lateral=False, maximum_periodic_wraps=10000):
    """Append the neutral second-flight measure without altering deposited ion charge."""
    if not isinstance(transport, BoundaryTransport3DResult):
        raise TypeError("transport must be BoundaryTransport3DResult")
    populations = {item.name: item for item in transport.surface_fluxes.energetic_fluxes}
    if model.ion_species_name not in populations:
        raise ValueError(f"transport has no {model.ion_species_name!r} impact population")
    if model.neutral_species_name in populations:
        raise ValueError("transport already contains the forward-scattered neutral species")
    incident = populations[model.ion_species_name]
    emitted = model.emit(incident, context)
    flight = trace_neutral_surface_events_3d(
        emitted, verts, faces, areas, context.face_gas_normal,
        domain_minimum=domain_minimum, domain_maximum=domain_maximum,
        mesh_length_unit_m=mesh_length_unit_m, launch_offset=launch_offset,
        periodic_lateral=periodic_lateral,
        maximum_periodic_wraps=maximum_periodic_wraps)
    face = incident.event_face
    material = context.face_material_id[face] == model.material_id
    probability = model.scattering_probability(incident.event_cosine_incidence)
    eligible_rate = float(np.dot(
        incident.event_flux_m2_s[material], context.face_area_m2[face[material]]))
    scattered_rate = float(np.sum(emitted.event_rate_s))
    nonscattered_rate = eligible_rate - scattered_rate
    weighted_incident_energy = float(np.dot(
        incident.event_flux_m2_s[material] * probability[material]
        * incident.event_energy_eV[material], context.face_area_m2[face[material]]))
    retained_energy = float(np.dot(
        emitted.event_rate_s, np.sum(emitted.event_velocity_sqrt_eV ** 2, axis=1)))
    deposited_energy = weighted_incident_energy - retained_energy
    particle_residual = eligible_rate - scattered_rate - nonscattered_rate
    energy_residual = weighted_incident_energy - retained_energy - deposited_energy
    particle_error = abs(particle_residual) / max(eligible_rate, np.finfo(float).tiny)
    energy_error = abs(energy_residual) / max(weighted_incident_energy, np.finfo(float).tiny)
    if particle_error > 5e-14 or energy_error > 5e-14 or nonscattered_rate < -1e-12 * eligible_rate:
        raise RuntimeError("forward-scatter surface ledger did not close")
    augmented = BoundaryTransport3DResult(
        SurfaceFluxes(
            transport.surface_fluxes.neutral_flux_m2_s,
            transport.surface_fluxes.energetic_fluxes + (flight.incident,)),
        transport.hit_probability, transport.escape_probability,
        transport.truncation_probability,
        transport.transport_model + "+hwang_giapis_neutralized_sio2_forward_scatter",
        transport.known_limitations + (
            "SiO2 forward scattering uses the Hwang--Giapis single-bounce specular hard-sphere model",
            "forward-scattered Cl is neutral and its incident ion charge remains deposited on SiO2",
        ), transport.lineage_replay_count, transport.lineage_replay_eligible_count,
        transport.edge_launch_inset_count + flight.edge_launch_inset_count,
        transport.trajectory_horizon_extension_count,
        transport.trajectory_initial_max_steps, transport.trajectory_final_max_steps,
        transport.trajectory_emergency_max_steps)
    result = HwangGiapisForwardScatter3DResult(
        emitted, flight, eligible_rate, scattered_rate, max(nonscattered_rate, 0.0),
        weighted_incident_energy, retained_energy, max(deposited_energy, 0.0),
        particle_error, energy_error, model.provenance)
    return augmented, result
