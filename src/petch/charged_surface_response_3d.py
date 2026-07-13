"""Conservative accounting contract for charged-particle surface response in 3-D.

This module does not choose reflection or emission physics.  It defines the rate measure emitted by
such a model and converts incident/outgoing charged events into signed deposited surface current.
Keeping this identity centralized prevents reflection, true-secondary emission, and future material
models from each inventing their own charge convention.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .sheath import ECHARGE
from .surface_kinetics import FaceResolvedEnergeticFlux


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

    outgoing_charge_rate = 0.0
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

    signed = positive - negative
    deposited_charge_rate = float(np.dot(signed, area))
    residual = incident_charge_rate - outgoing_charge_rate - deposited_charge_rate
    scale = max(
        abs(incident_charge_rate), abs(outgoing_charge_rate), abs(deposited_charge_rate),
        np.finfo(float).tiny)
    return ChargedSurfaceTransfer3D(
        positive, negative, signed, outgoing,
        incident_charge_rate_c_s=incident_charge_rate,
        outgoing_charge_rate_c_s=outgoing_charge_rate,
        deposited_charge_rate_c_s=deposited_charge_rate,
        charge_balance_residual_c_s=float(residual),
        relative_charge_balance_error=float(abs(residual) / scale))


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
    return ChargedSurfaceTransfer3D(
        positive, negative, signed, (),
        incident_charge_rate_c_s=deposited_charge_rate,
        outgoing_charge_rate_c_s=0.0,
        deposited_charge_rate_c_s=deposited_charge_rate,
        charge_balance_residual_c_s=0.0,
        relative_charge_balance_error=0.0)
