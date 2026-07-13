"""Conservative charged surface-response and full-field re-impact cascades."""
from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .boundary_transport_3d import (
    BoundaryTransport3DResult,
    ChargedSurfaceReimpactPopulation3D,
    trace_charged_surface_events_field_3d,
)
from .charged_surface_response_3d import (
    ChargedSurfaceContext3D,
    ChargedSurfaceResponse3D,
    ChargedSurfaceTransfer3D,
    perfect_absorber_surface_transfer_3d,
)
from .sheath import ECHARGE
from .surface_kinetics import FaceResolvedEnergeticFlux, SurfaceFluxes


@dataclass(frozen=True)
class ChargedSurfaceCascade3DResult:
    """Complete charge ledger for a finite response/re-impact cascade."""

    positive_deposition_current_density_a_m2: np.ndarray
    negative_deposition_current_density_a_m2: np.ndarray
    face_current_density_a_m2: np.ndarray
    transfers: tuple[ChargedSurfaceTransfer3D, ...]
    flights_by_bounce: tuple[tuple[ChargedSurfaceReimpactPopulation3D, ...], ...]
    unresolved_incident: tuple[FaceResolvedEnergeticFlux, ...]
    unresolved_charge_number_by_species: Mapping[str, int]
    initial_incident_charge_rate_c_s: float
    deposited_charge_rate_c_s: float
    escaped_charge_rate_c_s: float
    unresolved_charge_rate_c_s: float
    tail_closure_absolute_charge_rate_c_s: float
    tail_closure_relative_absolute_charge_rate: float
    tail_closure_l1_current_error_bound_relative: float
    charge_balance_residual_c_s: float
    relative_charge_balance_error: float
    completed: bool

    def __post_init__(self):
        positive = np.asarray(
            self.positive_deposition_current_density_a_m2, dtype=float).copy()
        negative = np.asarray(
            self.negative_deposition_current_density_a_m2, dtype=float).copy()
        signed = np.asarray(self.face_current_density_a_m2, dtype=float).copy()
        if (positive.ndim != 1 or negative.shape != positive.shape or signed.shape != positive.shape
                or np.any(~np.isfinite(positive)) or np.any(positive < 0.0)
                or np.any(~np.isfinite(negative)) or np.any(negative < 0.0)
                or not np.allclose(signed, positive - negative, rtol=1e-14, atol=0.0)):
            raise ValueError("invalid charged-cascade deposition currents")
        transfers = tuple(self.transfers)
        flights = tuple(tuple(items) for items in self.flights_by_bounce)
        unresolved = tuple(self.unresolved_incident)
        charge = {name: int(value) for name, value in self.unresolved_charge_number_by_species.items()}
        if (not transfers
                or any(not isinstance(item, ChargedSurfaceTransfer3D) for item in transfers)
                or any(any(not isinstance(item, ChargedSurfaceReimpactPopulation3D)
                            for item in items) for items in flights)
                or any(not isinstance(item, FaceResolvedEnergeticFlux) for item in unresolved)
                or set(charge) != {item.name for item in unresolved}
                or any(value == 0 for value in charge.values())
                or bool(unresolved) == bool(self.completed)):
            raise ValueError("invalid charged-cascade history or completion state")
        rates = np.asarray([
            self.initial_incident_charge_rate_c_s, self.deposited_charge_rate_c_s,
            self.escaped_charge_rate_c_s, self.unresolved_charge_rate_c_s,
            self.tail_closure_absolute_charge_rate_c_s,
            self.tail_closure_relative_absolute_charge_rate,
            self.tail_closure_l1_current_error_bound_relative,
            self.charge_balance_residual_c_s, self.relative_charge_balance_error], dtype=float)
        if (np.any(~np.isfinite(rates)) or self.relative_charge_balance_error < 0.0
                or self.tail_closure_absolute_charge_rate_c_s < 0.0
                or self.tail_closure_relative_absolute_charge_rate < 0.0
                or self.tail_closure_l1_current_error_bound_relative < 0.0):
            raise ValueError("invalid charged-cascade charge ledger")
        if not np.isclose(
                self.tail_closure_l1_current_error_bound_relative,
                2.0 * self.tail_closure_relative_absolute_charge_rate,
                rtol=0.0, atol=8.0 * np.finfo(float).eps):
            raise ValueError("charged-cascade tail error bound is inconsistent")
        if self.relative_charge_balance_error > 5e-13:
            raise ValueError("charged surface cascade does not conserve signed charge")
        for value in (positive, negative, signed):
            value.setflags(write=False)
        object.__setattr__(self, "positive_deposition_current_density_a_m2", positive)
        object.__setattr__(self, "negative_deposition_current_density_a_m2", negative)
        object.__setattr__(self, "face_current_density_a_m2", signed)
        object.__setattr__(self, "transfers", transfers)
        object.__setattr__(self, "flights_by_bounce", flights)
        object.__setattr__(self, "unresolved_incident", unresolved)
        object.__setattr__(
            self, "unresolved_charge_number_by_species", MappingProxyType(charge))


def _incident_charge_rate(populations, charge_number_by_species, face_area_m2):
    charge = dict(charge_number_by_species)
    area = np.asarray(face_area_m2, dtype=float)
    result = 0.0
    for population in populations:
        result += float(
            ECHARGE * int(charge[population.name])
            * np.dot(population.event_flux_m2_s, area[population.event_face]))
    return result


def _incident_absolute_charge_rate(populations, charge_number_by_species, face_area_m2):
    charge = dict(charge_number_by_species)
    area = np.asarray(face_area_m2, dtype=float)
    result = 0.0
    for population in populations:
        result += float(
            ECHARGE * abs(int(charge[population.name]))
            * np.dot(population.event_flux_m2_s, area[population.event_face]))
    return result


def solve_charged_surface_cascade_3d(
        incident_populations, charge_number_by_species, response: ChargedSurfaceResponse3D,
        context: ChargedSurfaceContext3D, verts, faces, areas, *,
        nodal_potential_v, potential_origin, potential_spacing,
        mesh_length_unit_m=1e-6, launch_offset=1e-5, fixed_dt=0.01,
        max_steps=10000, max_bounces=16, relative_tail_tolerance=0.0,
        periodic_lateral=False, device=None):
    """Alternate material response and full-field flight without losing capped charge.

    ``max_bounces`` caps response evaluations, not charge accounting.  If the cap is reached, the
    landed-but-unprocessed population is returned in ``unresolved_incident`` and included explicitly
    in the global charge ledger.  Production callers may refine the cap or apply an independently
    justified unbiased roulette; they may not treat an incomplete cascade as closed.

    A positive ``relative_tail_tolerance`` enables a deterministic conservative tail closure once
    the absolute charge rate remaining after at least one flight is below that fraction of the
    primary absolute charge rate. The remaining landed population is absorbed on its current faces.
    Global charge remains exact; the normalized L1 error of the spatial current distribution is
    bounded by twice the reported tail fraction. Zero retains the strict no-truncation behavior.
    """
    incident = tuple(incident_populations)
    supplied_charge = dict(charge_number_by_species)
    charge = {name: int(value) for name, value in supplied_charge.items()}
    if (not incident or any(not isinstance(item, FaceResolvedEnergeticFlux) for item in incident)
            or set(charge) != {item.name for item in incident}
            or any(int(value) != value or int(value) == 0 for value in supplied_charge.values())
            or not isinstance(context, ChargedSurfaceContext3D)
            or not hasattr(response, "evaluate")
            or int(max_bounces) != max_bounces or max_bounces <= 0
            or not np.isfinite(relative_tail_tolerance)
            or not 0.0 <= relative_tail_tolerance < 1.0):
        raise ValueError("invalid charged surface-cascade inputs")

    positive = np.zeros_like(context.face_area_m2)
    negative = np.zeros_like(context.face_area_m2)
    transfers = []
    flights_by_bounce = []
    escaped_charge_rate = 0.0
    current = incident
    current_charge = charge
    initial_charge_rate = _incident_charge_rate(
        incident, charge, context.face_area_m2)
    initial_absolute_charge_rate = _incident_absolute_charge_rate(
        incident, charge, context.face_area_m2)
    tail_closure_absolute_charge_rate = 0.0
    for _bounce in range(int(max_bounces)):
        expected_incident_charge_rate = _incident_charge_rate(
            current, current_charge, context.face_area_m2)
        expected_absolute_charge_rate = _incident_absolute_charge_rate(
            current, current_charge, context.face_area_m2)
        tail_close = bool(
            _bounce > 0 and relative_tail_tolerance > 0.0
            and expected_absolute_charge_rate
            <= float(relative_tail_tolerance) * initial_absolute_charge_rate)
        transfer = (
            perfect_absorber_surface_transfer_3d(
                current, current_charge, context.face_area_m2)
            if tail_close else response.evaluate(current, current_charge, context))
        if (not isinstance(transfer, ChargedSurfaceTransfer3D)
                or transfer.face_current_density_a_m2.shape != context.face_area_m2.shape):
            raise TypeError("charged surface response returned an incompatible transfer")
        local_scale = max(
            expected_absolute_charge_rate, abs(transfer.incident_charge_rate_c_s),
            np.finfo(float).tiny)
        if (abs(transfer.incident_charge_rate_c_s - expected_incident_charge_rate)
                    > 5e-13 * local_scale
                or transfer.relative_charge_balance_error > 5e-13):
            raise ValueError("charged surface response violated its local charge ledger")
        transfers.append(transfer)
        positive += transfer.positive_deposition_current_density_a_m2
        negative += transfer.negative_deposition_current_density_a_m2
        if tail_close:
            tail_closure_absolute_charge_rate = expected_absolute_charge_rate
            current = ()
            current_charge = {}
            break
        if not transfer.outgoing:
            current = ()
            current_charge = {}
            break
        outgoing_names = [item.name for item in transfer.outgoing]
        if len(set(outgoing_names)) != len(outgoing_names):
            raise ValueError("one response evaluation must use unique outgoing population names")
        flights = trace_charged_surface_events_field_3d(
            transfer.outgoing, verts, faces, areas, context.face_gas_normal,
            nodal_potential_v=nodal_potential_v, potential_origin=potential_origin,
            potential_spacing=potential_spacing, mesh_length_unit_m=mesh_length_unit_m,
            launch_offset=launch_offset, fixed_dt=fixed_dt, max_steps=max_steps,
            periodic_lateral=periodic_lateral, allow_truncation=False, device=device)
        flights_by_bounce.append(flights)
        for flight in flights:
            escaped_charge_rate += float(
                ECHARGE * flight.emitted.charge_number * flight.escaped_rate_s)
        landed = tuple(flight for flight in flights if flight.landed_rate_s > 0.0)
        current = tuple(flight.incident for flight in landed)
        current_charge = {
            flight.incident.name: flight.emitted.charge_number for flight in landed}
        if not current:
            break

    signed = positive - negative
    deposited_charge_rate = float(sum(item.deposited_charge_rate_c_s for item in transfers))
    unresolved_charge_rate = (
        _incident_charge_rate(current, current_charge, context.face_area_m2)
        if current else 0.0)
    residual = (
        initial_charge_rate - deposited_charge_rate
        - escaped_charge_rate - unresolved_charge_rate)
    scale = max(
        initial_absolute_charge_rate,
        abs(deposited_charge_rate), abs(escaped_charge_rate), abs(unresolved_charge_rate),
        np.finfo(float).tiny)
    tail_relative = float(
        tail_closure_absolute_charge_rate
        / max(initial_absolute_charge_rate, np.finfo(float).tiny))
    return ChargedSurfaceCascade3DResult(
        positive, negative, signed, tuple(transfers), tuple(flights_by_bounce),
        current, current_charge,
        initial_incident_charge_rate_c_s=initial_charge_rate,
        deposited_charge_rate_c_s=deposited_charge_rate,
        escaped_charge_rate_c_s=escaped_charge_rate,
        unresolved_charge_rate_c_s=unresolved_charge_rate,
        tail_closure_absolute_charge_rate_c_s=tail_closure_absolute_charge_rate,
        tail_closure_relative_absolute_charge_rate=tail_relative,
        tail_closure_l1_current_error_bound_relative=2.0 * tail_relative,
        charge_balance_residual_c_s=float(residual),
        relative_charge_balance_error=float(abs(residual) / scale),
        completed=not bool(current))


def _merge_face_resolved_populations_3d(populations):
    """Concatenate one species' exact event measures without histogramming lineage."""
    populations = tuple(populations)
    if (not populations
            or any(not isinstance(item, FaceResolvedEnergeticFlux) for item in populations)
            or len({item.name for item in populations}) != 1
            or len({item.face_count for item in populations}) != 1):
        raise ValueError("face-resolved populations must share one species and mesh")
    positions = [item.event_position for item in populations]
    directions = [item.event_incident_direction for item in populations]
    if any(value is None for value in positions) and not all(value is None for value in positions):
        raise ValueError("cannot merge partially preserved impact positions")
    if any(value is None for value in directions) and not all(value is None for value in directions):
        raise ValueError("cannot merge partially preserved incident directions")
    first = populations[0]
    return FaceResolvedEnergeticFlux(
        first.name, first.face_count,
        np.concatenate([item.event_face for item in populations]),
        np.concatenate([item.event_flux_m2_s for item in populations]),
        np.concatenate([item.event_energy_eV for item in populations]),
        np.concatenate([item.event_cosine_incidence for item in populations]),
        event_position=(None if positions[0] is None else np.concatenate(positions)),
        event_incident_direction=(
            None if directions[0] is None else np.concatenate(directions)))


def augment_transport_with_charged_reimpacts_3d(
        transport: BoundaryTransport3DResult, cascade: ChargedSurfaceCascade3DResult):
    """Return the chemistry-facing impact measure including every landed charged re-impact.

    Boundary hit/escape probabilities continue to describe particles launched at the plasma
    boundary. The energetic surface measure instead contains the primary hits plus all landed
    cascade hits. Events are concatenated by species with their face, position, energy, angle,
    and direction lineage intact. An incomplete cascade is refused because omitting its remaining
    impacts would silently violate the energetic-flux contract.
    """
    if not isinstance(transport, BoundaryTransport3DResult):
        raise TypeError("transport must be BoundaryTransport3DResult")
    if not isinstance(cascade, ChargedSurfaceCascade3DResult):
        raise TypeError("cascade must be ChargedSurfaceCascade3DResult")
    if not cascade.completed:
        raise ValueError("cannot construct a chemistry flux from an incomplete charged cascade")

    energetic = list(transport.surface_fluxes.energetic_fluxes)
    reimpacts = [
        flight.incident
        for bounce in cascade.flights_by_bounce
        for flight in bounce
        if flight.landed_rate_s > 0.0]
    if not reimpacts:
        return transport

    ordered_names = []
    grouped = {}
    passthrough = []
    for population in energetic + reimpacts:
        if isinstance(population, FaceResolvedEnergeticFlux):
            if population.name not in grouped:
                ordered_names.append(population.name)
                grouped[population.name] = []
            grouped[population.name].append(population)
        else:
            passthrough.append(population)
    merged = tuple(
        _merge_face_resolved_populations_3d(grouped[name]) for name in ordered_names)
    tail_limitation = (() if cascade.tail_closure_absolute_charge_rate_c_s == 0.0 else (
        "charged cascade tail was absorbed on its current faces with normalized spatial-current "
        f"L1 error bounded by {cascade.tail_closure_l1_current_error_bound_relative:.6g}",
    ))
    limitations = tuple(transport.known_limitations) + (
        "energetic surface flux includes conservative charged-particle re-impacts",
        "boundary hit/escape probabilities exclude subsequent surface-response flights",
    ) + tail_limitation
    return BoundaryTransport3DResult(
        SurfaceFluxes(
            transport.surface_fluxes.neutral_flux_m2_s,
            tuple(passthrough) + merged),
        transport.hit_probability, transport.escape_probability,
        transport.truncation_probability,
        transport.transport_model + " + charged_surface_reimpact_cascade",
        tuple(dict.fromkeys(limitations)))
