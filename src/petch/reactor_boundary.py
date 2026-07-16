"""Diagnostic-conditioned reactor/sheath boundary providers.

This module is the deliberately small bridge between reactor diagnostics (or a later equipment
model) and petch's authoritative :class:`~petch.boundary_state.PlasmaBoundaryState`.  It does not
invent missing reactor information: a self-bias scalar is not a sheath-voltage waveform, and an
assumed waveform can only produce a development/sensitivity boundary.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .boundary_state import (
    MaxwellianFluxVelocityDensity,
    PlasmaBoundaryState,
    SpeciesBoundaryState,
    collisionless_sheath_boundary_state,
    maxwellian_electron_boundary_state,
)
from .experimental_data import load_krueger_2024_evidence
from .sheath import (
    CollisionlessWaveformSheath,
    PeriodicSheathVoltage,
    bohm_speed,
)


_EVIDENCE_KINDS = {
    "measured",
    "published_distribution",
    "validated_reactor_model",
    "assumed",
}
_PREDICTIVE_EVIDENCE_KINDS = {
    "measured",
    "published_distribution",
    "validated_reactor_model",
}
_REACTOR_FLUX_EVIDENCE_KINDS = _EVIDENCE_KINDS | {
    "HPEM_simulation",
    "published_reactor_model_output",
}
_REACTOR_SPECIES_ROLES = {
    "neutral", "positive_ion", "negative_ion", "electron",
    "positive_ion_mixture", "negative_ion_mixture",
}
_BOLTZMANN_EV_PER_K = 8.617333262145e-5


def _is_sha256(value):
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


@dataclass(frozen=True)
class ReactorSpeciesFlux:
    """One reactor-to-wafer species flux fact, without an invented phase-space law.

    ``charge_number`` and ``mass_amu`` may be absent only for an explicitly unresolved mixture
    reported by a source as an aggregate (for example, Krüger Table I's single ``ions`` row).
    Such a row remains useful evidence, but it cannot be converted into a kinetic feature boundary
    until a species-resolved reactor output or an independently declared mixture closure is supplied.
    """

    name: str
    flux_m2_s: float
    role: str
    evidence_kind: str
    source_location: str
    charge_number: int | None = None
    mass_amu: float | None = None

    def __post_init__(self):
        if (
            not str(self.name).strip()
            or not np.isfinite(self.flux_m2_s)
            or self.flux_m2_s < 0.0
            or self.role not in _REACTOR_SPECIES_ROLES
            or self.evidence_kind not in _REACTOR_FLUX_EVIDENCE_KINDS
            or not str(self.source_location).strip()
        ):
            raise ValueError("invalid reactor species-flux record")
        unresolved_mixture = self.role in {
            "positive_ion_mixture", "negative_ion_mixture"}
        if unresolved_mixture:
            if self.charge_number is not None or self.mass_amu is not None:
                raise ValueError(
                    "an unresolved reactor mixture cannot declare a representative charge or mass")
            return
        if (
            self.charge_number is None
            or int(self.charge_number) != self.charge_number
            or self.mass_amu is None
            or not np.isfinite(self.mass_amu)
            or self.mass_amu <= 0.0
        ):
            raise ValueError("resolved reactor species require integer charge and positive mass")
        expected_sign = {
            "neutral": 0,
            "positive_ion": 1,
            "negative_ion": -1,
            "electron": -1,
        }[self.role]
        if (
            (expected_sign == 0 and self.charge_number != 0)
            or (expected_sign > 0 and self.charge_number <= 0)
            or (expected_sign < 0 and self.charge_number >= 0)
        ):
            raise ValueError("reactor species role and charge number disagree")

    @property
    def resolved(self):
        return self.charge_number is not None and self.mass_amu is not None

    @property
    def supports_predictive_boundary(self):
        return self.evidence_kind in _PREDICTIVE_EVIDENCE_KINDS


@dataclass(frozen=True)
class TabulatedReactorFluxDeck:
    """Provenance-bound wall-flux vector from measurements or a reactor calculation."""

    species_fluxes: tuple[ReactorSpeciesFlux, ...]
    source: str
    source_sha256: str
    reactor_model_validation_reference: str | None = None
    provenance: Mapping[str, object] = None

    def __post_init__(self):
        records = tuple(self.species_fluxes)
        if (
            not records
            or len({item.name for item in records}) != len(records)
            or any(not isinstance(item, ReactorSpeciesFlux) for item in records)
            or not str(self.source).strip()
            or not _is_sha256(self.source_sha256)
            or (
                self.reactor_model_validation_reference is not None
                and not str(self.reactor_model_validation_reference).strip()
            )
        ):
            raise ValueError("invalid tabulated reactor flux deck")
        object.__setattr__(self, "species_fluxes", records)
        object.__setattr__(
            self, "provenance",
            MappingProxyType({} if self.provenance is None else dict(self.provenance)),
        )

    @property
    def unresolved_species(self):
        return tuple(item.name for item in self.species_fluxes if not item.resolved)

    @property
    def supports_predictive_boundary(self):
        reactor_output = any(
            item.evidence_kind in {"HPEM_simulation", "published_reactor_model_output"}
            for item in self.species_fluxes
        )
        return (
            not self.unresolved_species
            and all(item.supports_predictive_boundary for item in self.species_fluxes)
            and (
                not reactor_output
                or bool(str(self.reactor_model_validation_reference or "").strip())
            )
        )

    def get(self, name):
        for item in self.species_fluxes:
            if item.name == name:
                return item
        raise KeyError(name)


def _thermal_flux_species(record, temperature_K, n_transverse, n_normal):
    temperature_eV = float(temperature_K) * _BOLTZMANN_EV_PER_K
    node, node_weight = np.polynomial.hermite.hermgauss(int(n_transverse))
    normal_node, normal_weight = np.polynomial.laguerre.laggauss(int(n_normal))
    ix, iy, iz = np.meshgrid(
        np.arange(node.size), np.arange(node.size), np.arange(normal_node.size),
        indexing="ij",
    )
    velocity = np.column_stack((
        np.sqrt(temperature_eV) * node[ix.ravel()],
        np.sqrt(temperature_eV) * node[iy.ravel()],
        np.sqrt(temperature_eV * normal_node[iz.ravel()]),
    ))
    weight = (
        node_weight[ix.ravel()]
        * node_weight[iy.ravel()]
        * normal_weight[iz.ravel()]
        / np.pi
    )
    return SpeciesBoundaryState(
        name=record.name,
        charge_number=int(record.charge_number),
        mass_amu=float(record.mass_amu),
        flux_m2_s=float(record.flux_m2_s),
        velocity_sqrt_eV=velocity,
        weight=weight,
        density_model=MaxwellianFluxVelocityDensity(temperature_eV),
        provenance={
            "provider": "tabulated_reactor_flux_deck",
            "role": record.role,
            "flux_evidence_kind": record.evidence_kind,
            "flux_source_location": record.source_location,
            "thermal_temperature_K": float(temperature_K),
            "thermal_distribution_closure": "half_maxwellian_at_declared_gas_temperature",
        },
    )


def _replace_species_flux(template, record, deck):
    if (
        template.name != record.name
        or template.charge_number != record.charge_number
        or not np.isclose(template.mass_amu, record.mass_amu, rtol=0.0, atol=1e-12)
    ):
        raise ValueError(
            f"kinetic template for {record.name!r} disagrees with reactor charge or mass")
    return SpeciesBoundaryState(
        name=template.name,
        charge_number=template.charge_number,
        mass_amu=template.mass_amu,
        flux_m2_s=float(record.flux_m2_s),
        velocity_sqrt_eV=template.velocity_sqrt_eV,
        weight=template.weight,
        phase_rad=template.phase_rad,
        position_m=template.position_m,
        density_model=template.density_model,
        density_model_2d=template.density_model_2d,
        provenance=dict(
            template.provenance,
            provider="tabulated_reactor_flux_deck",
            role=record.role,
            flux_evidence_kind=record.evidence_kind,
            flux_source_location=record.source_location,
            reactor_flux_deck_source=deck.source,
            reactor_flux_deck_sha256=deck.source_sha256,
        ),
    )


def build_tabulated_reactor_boundary(
        deck: TabulatedReactorFluxDeck, *,
        reference_plane_m: float,
        kinetic_templates: Mapping[str, SpeciesBoundaryState] | None = None,
        included_species: tuple[str, ...] | None = None,
        neutral_temperature_K: float = 300.0,
        n_transverse_neutral: int = 5,
        n_normal_neutral: int = 8,
        claim_mode: str = "development"):
    """Convert a reactor flux deck into the common kinetic boundary without filling gaps.

    Neutral rows may use an explicitly declared thermal half-Maxwellian closure. Charged rows
    always require a kinetic template carrying their IEAD/EEDF. ``included_species`` is explicit:
    selecting a development subset is legal and recorded, while the default requests the complete
    deck and therefore refuses unresolved aggregate rows. Predictive mode additionally requires the
    complete deck and predictive evidence for every supplied flux and distribution.
    """
    if not isinstance(deck, TabulatedReactorFluxDeck):
        raise TypeError("deck must be a TabulatedReactorFluxDeck")
    values = np.asarray([reference_plane_m, neutral_temperature_K], dtype=float)
    if (
        np.any(~np.isfinite(values))
        or reference_plane_m < 0.0
        or neutral_temperature_K <= 0.0
        or int(n_transverse_neutral) != n_transverse_neutral
        or n_transverse_neutral <= 0
        or int(n_normal_neutral) != n_normal_neutral
        or n_normal_neutral <= 0
        or claim_mode not in {"development", "predictive"}
    ):
        raise ValueError("invalid tabulated reactor-boundary controls")
    templates = {} if kinetic_templates is None else dict(kinetic_templates)
    if any(not isinstance(value, SpeciesBoundaryState) for value in templates.values()):
        raise TypeError("kinetic_templates must contain SpeciesBoundaryState values")
    names = tuple(item.name for item in deck.species_fluxes)
    selected_names = names if included_species is None else tuple(included_species)
    if (
        not selected_names
        or len(set(selected_names)) != len(selected_names)
        or not set(selected_names).issubset(names)
    ):
        raise ValueError("included_species must be a nonempty unique subset of the flux deck")
    if claim_mode == "predictive":
        if selected_names != names or not deck.supports_predictive_boundary:
            raise ValueError(
                "predictive reactor boundary requires the complete resolved deck and "
                "measurement/validated-model evidence")

    species = []
    for name in selected_names:
        record = deck.get(name)
        if not record.resolved:
            raise ValueError(
                f"reactor flux {name!r} is an unresolved mixture; supply species-resolved "
                "fluxes and kinetic distributions before building a feature boundary")
        if record.role == "neutral" and name not in templates:
            species.append(_thermal_flux_species(
                record, neutral_temperature_K,
                int(n_transverse_neutral), int(n_normal_neutral)))
        elif name in templates:
            species.append(_replace_species_flux(templates[name], record, deck))
        else:
            raise ValueError(
                f"charged reactor species {name!r} requires an explicit kinetic template")
    omitted = tuple(name for name in names if name not in selected_names)
    return PlasmaBoundaryState(
        species=tuple(species),
        reference_plane_m=float(reference_plane_m),
        provenance={
            "provider": "tabulated_reactor_flux_deck",
            "claim_mode": claim_mode,
            "supports_prediction": claim_mode == "predictive",
            "source": deck.source,
            "source_sha256": deck.source_sha256,
            "reactor_model_validation_reference": deck.reactor_model_validation_reference,
            "selected_species": selected_names,
            "omitted_species": omitted,
            "unresolved_species_in_complete_deck": deck.unresolved_species,
            "complete_flux_deck_used": not omitted,
            "neutral_temperature_K": float(neutral_temperature_K),
            "deck_provenance": dict(deck.provenance),
        },
    )


_KRUEGER_NEUTRAL_MASS_AMU = {
    "C3F4": 112.0263,
    "C2F3": 81.0178,
    "CF": 31.0094,
    "CF2": 50.0078,
    "CF3": 69.0062,
    "O": 15.999,
}


def load_krueger_2024_reactor_flux_deck(directory):
    """Load Krüger Table-I HPEM wall fluxes while preserving its unresolved ion mixture."""
    directory = Path(directory)
    evidence = load_krueger_2024_evidence(directory)
    source_path = directory / "base_case_boundary_fluxes.csv"
    digest = sha256(source_path.read_bytes()).hexdigest()
    records = []
    for item in evidence.boundary_fluxes:
        flux_m2_s = float(item.value_cm2_s) * 1.0e4
        if item.species == "ions":
            records.append(ReactorSpeciesFlux(
                name=item.species,
                flux_m2_s=flux_m2_s,
                role="positive_ion_mixture",
                evidence_kind=item.evidence_type,
                source_location=item.source_location,
            ))
        else:
            records.append(ReactorSpeciesFlux(
                name=item.species,
                flux_m2_s=flux_m2_s,
                role="neutral",
                evidence_kind=item.evidence_type,
                source_location=item.source_location,
                charge_number=0,
                mass_amu=_KRUEGER_NEUTRAL_MASS_AMU[item.species],
            ))
    return TabulatedReactorFluxDeck(
        species_fluxes=tuple(records),
        source="Krüger et al., JVST A 42, 043008 (2024), HPEM Table I",
        source_sha256=digest,
        reactor_model_validation_reference=None,
        provenance={
            "doi": "10.1116/6.0003554",
            "evidence_status": "published HPEM outputs; not measurements",
            "process": "10 mTorr C4F6/Ar/O2, 140/100/105 sccm, 1/40 MHz CCP",
            "ion_distribution_status": (
                "aggregate ion flux only; species-resolved ion flux and IEAD are absent "
                "from the bundled evidence"),
        },
    )


@dataclass(frozen=True)
class PlasmaDiagnosticState:
    """Minimal plasma-side state needed by the collisionless virtual sheath.

    If ``ion_flux_m2_s`` is omitted, ``electropositive_bohm_flux_closure`` must be explicitly true.
    That closure is intentionally not implicit: Bohm injection can be inaccurate in electronegative
    plasmas, where a measured ion flux or a richer presheath model is required.
    """
    electron_density_m3: float
    electron_temperature_eV: float
    ion_name: str
    ion_mass_amu: float
    source: str
    density_evidence_kind: str = "assumed"
    temperature_evidence_kind: str = "assumed"
    ion_flux_m2_s: float | None = None
    ion_flux_evidence_kind: str | None = None
    electropositive_bohm_flux_closure: bool = False

    def __post_init__(self):
        values = np.asarray([
            self.electron_density_m3,
            self.electron_temperature_eV,
            self.ion_mass_amu,
        ], dtype=float)
        if (np.any(~np.isfinite(values)) or np.any(values <= 0.0)
                or not str(self.ion_name).strip() or not str(self.source).strip()
                or self.density_evidence_kind not in _EVIDENCE_KINDS
                or self.temperature_evidence_kind not in _EVIDENCE_KINDS):
            raise ValueError("invalid diagnostic plasma state")
        if self.ion_flux_m2_s is None:
            if self.ion_flux_evidence_kind is not None:
                raise ValueError("ion-flux evidence requires an explicit ion flux")
        elif (not np.isfinite(self.ion_flux_m2_s) or self.ion_flux_m2_s <= 0.0
              or self.ion_flux_evidence_kind not in _EVIDENCE_KINDS):
            raise ValueError("invalid explicit ion-flux diagnostic")

    @property
    def ion_flux(self):
        if self.ion_flux_m2_s is not None:
            return float(self.ion_flux_m2_s), "explicit_ion_flux"
        if not self.electropositive_bohm_flux_closure:
            raise ValueError(
                "ion flux is missing; explicitly authorize the electropositive Bohm closure")
        return (float(self.electron_density_m3 * bohm_speed(
            self.electron_temperature_eV, self.ion_mass_amu)),
            "electropositive_bohm_flux")

    @property
    def supports_predictive_boundary(self):
        if (self.density_evidence_kind not in _PREDICTIVE_EVIDENCE_KINDS
                or self.temperature_evidence_kind not in _PREDICTIVE_EVIDENCE_KINDS):
            return False
        if self.ion_flux_m2_s is None:
            return self.electropositive_bohm_flux_closure
        return self.ion_flux_evidence_kind in _PREDICTIVE_EVIDENCE_KINDS


def _with_provenance(species, provenance):
    return SpeciesBoundaryState(
        name=species.name,
        charge_number=species.charge_number,
        mass_amu=species.mass_amu,
        flux_m2_s=species.flux_m2_s,
        velocity_sqrt_eV=species.velocity_sqrt_eV,
        weight=species.weight,
        phase_rad=species.phase_rad,
        position_m=species.position_m,
        density_model=species.density_model,
        provenance=provenance,
    )


def build_diagnostic_virtual_sheath_boundary(
        diagnostic: PlasmaDiagnosticState,
        waveform: PeriodicSheathVoltage,
        *,
        reference_plane_m: float,
        collisionless_justification: str,
        claim_mode: str = "development",
        model_validation_reference: str | None = None,
        ion_tangential_temperature_eV: float = 0.026,
        electron_flux_m2_s: float | None = None,
        n_phase: int = 256,
        n_transverse_ion: int = 3,
        n_transverse_electron: int = 5,
        n_normal_electron: int = 8,
        normal_energy_bins: int = 64,
        density_phase_count: int | None = None,
):
    """Build a common ion/electron boundary through the finite-transit sheath.

    ``claim_mode='predictive'`` is an evidence gate.  It requires diagnostic quantities and the
    complete sheath waveform to be measured/published or supplied by an independently validated
    reactor model, a collisionless-regime justification, and a validation reference for this reduced
    sheath closure.  Development mode retains every assumption in provenance and never upgrades the
    result into a prediction.

    The default electron closure enforces equal time-averaged ion/electron particle flux at the
    feature reference plane.  Kinetic feature charging still filters the Maxwellian electron phase
    space self-consistently; this boundary closure does not insert Boltzmann volume charge.
    """
    if not isinstance(diagnostic, PlasmaDiagnosticState):
        raise TypeError("diagnostic must be a PlasmaDiagnosticState")
    if not isinstance(waveform, PeriodicSheathVoltage):
        raise TypeError("waveform must be a PeriodicSheathVoltage")
    values = np.asarray([reference_plane_m, ion_tangential_temperature_eV], dtype=float)
    if (np.any(~np.isfinite(values)) or reference_plane_m < 0.0
            or ion_tangential_temperature_eV <= 0.0
            or claim_mode not in {"development", "predictive"}
            or not str(collisionless_justification).strip()):
        raise ValueError("invalid virtual-sheath boundary controls")

    probe_time = waveform.period_s * np.arange(4096, dtype=float) / 4096.0
    probe_voltage = waveform.voltage(probe_time)
    negative_fraction = float(np.mean(probe_voltage < 0.0))
    if claim_mode == "predictive":
        if (not diagnostic.supports_predictive_boundary
                or not waveform.supports_predictive_boundary
                or not str(model_validation_reference or "").strip()):
            raise ValueError(
                "predictive mode requires evidenced diagnostics, a measured/validated full "
                "sheath waveform, and a reduced-model validation reference")
        if np.min(probe_voltage) < -1e-9:
            raise ValueError("predictive sheath-voltage waveform cannot reverse sign")

    ion_flux, ion_flux_closure = diagnostic.ion_flux
    if electron_flux_m2_s is None:
        electron_flux = ion_flux
        electron_flux_closure = "ambipolar_time_average_equal_to_ion_flux"
    else:
        electron_flux = float(electron_flux_m2_s)
        if not np.isfinite(electron_flux) or electron_flux <= 0.0:
            raise ValueError("electron_flux_m2_s must be positive")
        electron_flux_closure = "explicit"

    sheath = CollisionlessWaveformSheath(
        waveform=waveform,
        Te_eV=diagnostic.electron_temperature_eV,
        ion_mass_amu=diagnostic.ion_mass_amu,
        density_m3=diagnostic.electron_density_m3,
    )
    ion_state = collisionless_sheath_boundary_state(
        sheath,
        ion_flux,
        n_phase=int(n_phase),
        ion_name=diagnostic.ion_name,
        reference_plane_m=float(reference_plane_m),
        tangential_temperature_eV=float(ion_tangential_temperature_eV),
        n_transverse=int(n_transverse_ion),
        normal_energy_bins=int(normal_energy_bins),
        density_phase_count=density_phase_count,
    )
    ion_base = ion_state.get(diagnostic.ion_name)
    shared = {
        "provider": "diagnostic_virtual_sheath",
        "claim_mode": claim_mode,
        "supports_prediction": claim_mode == "predictive",
        "diagnostic_source": diagnostic.source,
        "density_evidence_kind": diagnostic.density_evidence_kind,
        "temperature_evidence_kind": diagnostic.temperature_evidence_kind,
        "waveform_source": waveform.source,
        "waveform_evidence_kind": waveform.evidence_kind,
        "waveform_fundamental_frequency_hz": waveform.fundamental_frequency_hz,
        "waveform_harmonics": waveform.harmonic_number.tolist(),
        "waveform_negative_fraction_clipped_to_sheath_collapse": negative_fraction,
        "collisionless_justification": collisionless_justification,
        "model_validation_reference": model_validation_reference,
        "sheath_thickness_m": sheath.thickness,
    }
    ion = _with_provenance(ion_base, dict(
        ion_base.provenance,
        **shared,
        role="finite_transit_ion_iedf_iadf",
        ion_flux_closure=ion_flux_closure,
        ion_flux_evidence_kind=diagnostic.ion_flux_evidence_kind,
        electropositive_bohm_flux_closure=diagnostic.electropositive_bohm_flux_closure,
        ion_tangential_temperature_eV=float(ion_tangential_temperature_eV),
    ))

    electron_base = maxwellian_electron_boundary_state(
        diagnostic.electron_temperature_eV,
        electron_flux,
        n_transverse=int(n_transverse_electron),
        n_normal=int(n_normal_electron),
        reference_plane_m=float(reference_plane_m),
    ).get("electron")
    electron = _with_provenance(electron_base, dict(
        electron_base.provenance,
        **shared,
        role="analytic_half_maxwellian_electron_source",
        electron_flux_closure=electron_flux_closure,
    ))
    return PlasmaBoundaryState(
        species=(ion, electron),
        reference_plane_m=float(reference_plane_m),
        provenance=dict(
            shared,
            source="diagnostic_conditioned_reactor_to_feature_boundary",
            ion_flux_m2_s=ion_flux,
            electron_flux_m2_s=electron_flux,
            current_density_A_m2=float(
                1.602176634e-19 * (ion_flux - electron_flux)),
            volume_boltzmann_electron_term=False,
        ),
    )
