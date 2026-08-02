"""Diagnostic-conditioned reactor/sheath boundary providers.

This module is the deliberately small bridge between reactor diagnostics (or a later equipment
model) and petch's authoritative :class:`~petch.boundary_state.PlasmaBoundaryState`.  It does not
invent missing reactor information: a self-bias scalar is not a sheath-voltage waveform, and an
assumed waveform can only produce a development/sensitivity boundary.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
from hashlib import sha256
import json
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .boundary_state import (
    DiscreteEnergyAngleDensity2D,
    DiscreteEnergyPolarAzimuthDensity3D,
    MaxwellianFluxVelocityDensity,
    PlasmaBoundaryState,
    SpeciesBoundaryState,
    collisionless_sheath_boundary_state,
    maxwellian_electron_boundary_state,
)
from .angular_lift import axisymmetric_polar_weights
from .experimental_data import load_krueger_2024_base_boundary_fluxes
from .surface_kinetics import EnergeticFlux
from .sheath import (
    CollisionlessWaveformSheath,
    PeriodicSheathVoltage,
    bohm_speed,
)

# Angular resolution the Krüger Figure-4 digitization actually carries; the
# axisymmetric polar inversion is peeled on this grid when the caller has not
# compressed the quadrature to a coarser one.
_AXISYMMETRIC_LIFT_BIN_DEG = 0.25


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
KRUEGER_2024_IEAD_CSV_SHA256 = (
    "913d31be623ec5d52d226c8cea499e7f014cf4f5a27e017b519633c96e5e3ee3")
KRUEGER_2024_IEAD_METADATA_SHA256 = (
    "7904a700afdcd116c6f57ef35aeb5555661ffed0d004304d815d20f79840ca55")
KRUEGER_2024_TRANSFER_FLUX_SHA256 = (
    "716a416335a10a2b8a80ea971be2f1a99af41b09eca973d9f9549c60c6fdf9f3")
KRUEGER_2024_TRANSFER_IEAD_SHA256 = (
    "0b5e665031bb3c41fb35ff73a0c1213ba6c355f69f7730c54a9d02a1dc915e35")
KRUEGER_2024_TRANSFER_METADATA_SHA256 = (
    "44e85ee74d11c9cfe95ceb3781bfe122b3d6b5992bf4739dcfaf0b24561476b2")


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


@dataclass(frozen=True)
class Krueger2024DigitizedIEAD:
    """Digitized combined positive-ion IEAD from Krüger Figure 4(a)."""

    energy_eV: np.ndarray
    signed_angle_deg: np.ndarray
    probability_weight: np.ndarray
    source_pdf_sha256: str
    table_sha256: str
    metadata_sha256: str
    metadata: Mapping[str, object]

    def __post_init__(self):
        energy = np.asarray(self.energy_eV, dtype=float).copy()
        angle = np.asarray(self.signed_angle_deg, dtype=float).copy()
        weight = np.asarray(self.probability_weight, dtype=float).copy()
        if (energy.ndim != 1 or angle.shape != energy.shape or weight.shape != energy.shape
                or energy.size < 64 or np.any(~np.isfinite(energy))
                or np.any(energy < 0.0) or np.any(~np.isfinite(angle))
                or np.any(np.abs(angle) > 10.0) or np.any(~np.isfinite(weight))
                or np.any(weight < 0.0) or weight.sum() <= 0.0
                or not _is_sha256(self.source_pdf_sha256)
                or not _is_sha256(self.table_sha256)
                or not _is_sha256(self.metadata_sha256)):
            raise ValueError("invalid Krüger digitized IEAD")
        weight /= weight.sum()
        for value in (energy, angle, weight):
            value.setflags(write=False)
        object.__setattr__(self, "energy_eV", energy)
        object.__setattr__(self, "signed_angle_deg", angle)
        object.__setattr__(self, "probability_weight", weight)
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    @property
    def mean_energy_eV(self):
        return float(np.dot(self.probability_weight, self.energy_eV))

    def development_quadrature(self, *, energy_bin_eV=None, angle_bin_deg=None):
        """Return either the exact digitization or a provenance-priced centroid quadrature.

        Compression is a numerical acceleration for deterministic feature transport, not a new
        physical IEAD.  It bins the *joint* energy/signed-angle measure and replaces every occupied
        bin by its probability-weighted energy centroid and circular angle centroid.  Total
        probability and mean energy are therefore retained (up to roundoff), while the loss of
        higher moments is measured and returned in the manifest.  Omitting both widths returns the
        exact published digitization; supplying only one width is refused so the approximation can
        never be implicit.
        """
        if energy_bin_eV is None and angle_bin_deg is None:
            return (
                self.energy_eV,
                self.signed_angle_deg,
                self.probability_weight,
                MappingProxyType({
                    "mode": "exact_digitized_quadrature",
                    "source_node_count": int(self.energy_eV.size),
                    "node_count": int(self.energy_eV.size),
                    "probability_error": 0.0,
                    "mean_energy_relative_error": 0.0,
                    "second_energy_moment_relative_error": 0.0,
                    "mean_direction_maximum_absolute_error": 0.0,
                }),
            )
        if energy_bin_eV is None or angle_bin_deg is None:
            raise ValueError(
                "Krüger IEAD compression requires both energy and angle bin widths")
        widths = np.asarray([energy_bin_eV, angle_bin_deg], dtype=float)
        if np.any(~np.isfinite(widths)) or np.any(widths <= 0.0):
            raise ValueError("Krüger IEAD compression widths must be positive and finite")

        # Anchor bins to physical zero energy and -90 degrees rather than to the sampled extrema.
        # This makes the reduction invariant to adding/removing a zero-weight edge sample.
        key = np.column_stack((
            np.floor(self.energy_eV / widths[0]).astype(np.int64),
            np.floor((self.signed_angle_deg + 90.0) / widths[1]).astype(np.int64),
        ))
        _, inverse = np.unique(key, axis=0, return_inverse=True)
        node_count = int(inverse.max()) + 1
        weight = np.bincount(
            inverse, weights=self.probability_weight, minlength=node_count)
        energy = np.bincount(
            inverse,
            weights=self.probability_weight * self.energy_eV,
            minlength=node_count,
        ) / weight
        angle_rad = np.deg2rad(self.signed_angle_deg)
        mean_sine = np.bincount(
            inverse,
            weights=self.probability_weight * np.sin(angle_rad),
            minlength=node_count,
        ) / weight
        mean_cosine = np.bincount(
            inverse,
            weights=self.probability_weight * np.cos(angle_rad),
            minlength=node_count,
        ) / weight
        angle = np.rad2deg(np.arctan2(mean_sine, mean_cosine))
        weight /= weight.sum()

        source_second = float(np.dot(
            self.probability_weight, self.energy_eV * self.energy_eV))
        reduced_second = float(np.dot(weight, energy * energy))
        source_direction = np.array([
            np.dot(self.probability_weight, np.sin(angle_rad)),
            np.dot(self.probability_weight, np.cos(angle_rad)),
        ])
        reduced_angle_rad = np.deg2rad(angle)
        reduced_direction = np.array([
            np.dot(weight, np.sin(reduced_angle_rad)),
            np.dot(weight, np.cos(reduced_angle_rad)),
        ])
        packed = np.column_stack((energy, angle, weight)).astype("<f8", copy=False)
        diagnostics = MappingProxyType({
            "mode": "joint_probability_weighted_centroid_quadrature",
            "source_node_count": int(self.energy_eV.size),
            "node_count": int(node_count),
            "energy_bin_eV": float(widths[0]),
            "angle_bin_deg": float(widths[1]),
            "probability_error": float(abs(weight.sum() - 1.0)),
            "mean_energy_relative_error": float(
                abs(np.dot(weight, energy) - self.mean_energy_eV)
                / self.mean_energy_eV),
            "second_energy_moment_relative_error": float(
                abs(reduced_second - source_second) / source_second),
            "mean_direction_maximum_absolute_error": float(
                np.max(np.abs(reduced_direction - source_direction))),
            "compressed_quadrature_sha256": sha256(packed.tobytes()).hexdigest(),
        })
        return energy, angle, weight, diagnostics

    def energetic_flux(self, flux_m2_s, *, name="ions"):
        return EnergeticFlux(
            name, float(flux_m2_s), self.energy_eV,
            np.cos(np.deg2rad(self.signed_angle_deg)),
            self.probability_weight)

    def development_species(
            self, flux_m2_s, *, effective_mass_amu, mixture_closure, name="ions",
            energy_bin_eV=None, angle_bin_deg=None, azimuthal_closure=None,
            azimuthal_order=16):
        if (not np.isfinite(effective_mass_amu) or effective_mass_amu <= 0.0
                or not str(mixture_closure).strip()):
            raise ValueError("aggregate-ion development species requires mass and closure")
        if azimuthal_closure not in {None, "axisymmetric_uniform"}:
            raise ValueError(
                "azimuthal_closure must be None or 'axisymmetric_uniform'")
        if (int(azimuthal_order) != azimuthal_order
                or int(azimuthal_order) <= 0):
            raise ValueError("azimuthal_order must be a positive integer")
        energy, signed_angle, weight, quadrature = self.development_quadrature(
            energy_bin_eV=energy_bin_eV, angle_bin_deg=angle_bin_deg)
        speed = np.sqrt(energy)
        if azimuthal_closure is None:
            angle = np.deg2rad(signed_angle)
            velocity = np.column_stack((
                speed * np.sin(angle),
                np.zeros_like(speed),
                speed * np.cos(angle),
            ))
            quadrature_weight = weight
        else:
            # The published signed angle is the PROJECTION of the ion direction
            # into one plane, so the axisymmetric lift must invert that
            # marginalization rather than relabel it as the polar angle.  See
            # petch.angular_lift: identifying the two discards exactly sqrt(2)
            # of the angular width (RESULTS_ANGULAR_CONVERGENCE_P0_2026-08-02).
            lift_bin_deg = (
                _AXISYMMETRIC_LIFT_BIN_DEG if angle_bin_deg is None
                else float(angle_bin_deg))
            polar_deg, polar_weight, lift = axisymmetric_polar_weights(
                signed_angle, weight, bin_deg=lift_bin_deg)
            polar = np.deg2rad(polar_deg)
            azimuth = 2.0 * np.pi * (
                np.arange(int(azimuthal_order), dtype=float) + 0.5
            ) / int(azimuthal_order)
            transverse = (speed * np.sin(polar))[:, None]
            velocity = np.column_stack((
                (transverse * np.cos(azimuth)[None, :]).ravel(),
                (transverse * np.sin(azimuth)[None, :]).ravel(),
                np.repeat(speed * np.cos(polar), int(azimuthal_order)),
            ))
            quadrature_weight = np.repeat(
                polar_weight / int(azimuthal_order), int(azimuthal_order))
        return SpeciesBoundaryState(
            name=name, charge_number=1, mass_amu=float(effective_mass_amu),
            flux_m2_s=float(flux_m2_s), velocity_sqrt_eV=velocity,
            weight=quadrature_weight,
            density_model=(
                None if azimuthal_closure is None else
                DiscreteEnergyPolarAzimuthDensity3D(
                    energy, polar_deg, polar_weight)),
            density_model_2d=DiscreteEnergyAngleDensity2D(
                energy, signed_angle, weight),
            provenance={
                "provider": "krueger_2024_digitized_combined_iead",
                "claim_mode": "development",
                "supports_prediction": False,
                "aggregate_positive_ion_mixture": True,
                "effective_mass_amu": float(effective_mass_amu),
                "mixture_closure": str(mixture_closure),
                "source_pdf_sha256": self.source_pdf_sha256,
                "digitized_table_sha256": self.table_sha256,
                "digitization_metadata_sha256": self.metadata_sha256,
                "numerical_quadrature": dict(quadrature),
                "three_dimensional_azimuthal_closure": (
                    "none; quadrature lies in the published signed-angle plane"
                    if azimuthal_closure is None else
                    "axisymmetric uniform azimuth; explicit development closure because the "
                    "published IEAD supplies only one signed-angle plane"),
                "three_dimensional_azimuthal_order": (
                    None if azimuthal_closure is None else int(azimuthal_order)),
                "three_dimensional_polar_inversion": (
                    None if azimuthal_closure is None else
                    "onion-peel Abel inversion of the published planar marginal "
                    "(petch.angular_lift); the polar rms exceeds the planar rms "
                    "by sqrt(2) for any axisymmetric measure"),
                "three_dimensional_polar_inversion_diagnostics": (
                    None if azimuthal_closure is None else dict(lift)),
            })


@dataclass(frozen=True)
class Krueger2024TransferBoundaryData:
    """Checksum-bound Figure-16 HPEM boundary inputs for held-out transfer cases."""

    flux_m2_s_by_oxygen_ratio: Mapping[float, Mapping[str, float]]
    flux_interval_m2_s_by_oxygen_ratio: Mapping[float, Mapping[str, tuple[float, float]]]
    iead_by_low_frequency_power_kw: Mapping[float, Krueger2024DigitizedIEAD]
    flux_table_sha256: str
    iead_table_sha256: str
    metadata_sha256: str
    metadata: Mapping[str, object]

    def __post_init__(self):
        flux = {
            float(ratio): MappingProxyType({
                str(name): float(value) for name, value in values.items()
            })
            for ratio, values in self.flux_m2_s_by_oxygen_ratio.items()}
        interval = {
            float(ratio): MappingProxyType({
                str(name): (float(value[0]), float(value[1]))
                for name, value in values.items()
            })
            for ratio, values in self.flux_interval_m2_s_by_oxygen_ratio.items()}
        iead = {
            float(power): value
            for power, value in self.iead_by_low_frequency_power_kw.items()}
        expected_ratios = {0.5, 1.0, 1.5, 2.5}
        expected_powers = {0.0, 4.0, 6.0, 8.0}
        if (set(flux) != expected_ratios or set(interval) != expected_ratios
                or set(iead) != expected_powers
                or any(set(flux[key]) != set(interval[key]) for key in flux)
                or any(not isinstance(value, Krueger2024DigitizedIEAD)
                       for value in iead.values())
                or any(not _is_sha256(value) for value in (
                    self.flux_table_sha256, self.iead_table_sha256,
                    self.metadata_sha256))):
            raise ValueError("invalid Krüger Figure-16 transfer boundary data")
        for ratio in flux:
            for name, value in flux[ratio].items():
                lower, upper = interval[ratio][name]
                if (not np.isfinite(value) or value <= 0.0
                        or not np.isfinite(lower) or not np.isfinite(upper)
                        or lower <= 0.0 or not lower <= value <= upper):
                    raise ValueError("invalid Krüger Figure-16 flux interval")
        object.__setattr__(
            self, "flux_m2_s_by_oxygen_ratio", MappingProxyType(flux))
        object.__setattr__(
            self, "flux_interval_m2_s_by_oxygen_ratio", MappingProxyType(interval))
        object.__setattr__(
            self, "iead_by_low_frequency_power_kw", MappingProxyType(iead))
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))


def load_krueger_2024_digitized_iead(directory):
    """Load the checksum-bound Figure-4 IEAD digitization and its uncertainty metadata."""
    directory = Path(directory)
    table_path = directory / "digitized_figure4_iead.csv"
    metadata_path = directory / "digitized_figure4_iead_metadata.json"
    table_sha = sha256(table_path.read_bytes()).hexdigest()
    metadata_sha = sha256(metadata_path.read_bytes()).hexdigest()
    if (table_sha != KRUEGER_2024_IEAD_CSV_SHA256
            or metadata_sha != KRUEGER_2024_IEAD_METADATA_SHA256):
        raise ValueError("Krüger IEAD digitization checksum mismatch")
    table = np.genfromtxt(table_path, delimiter=",", names=True)
    if table.dtype.names != (
            "energy_eV", "signed_angle_deg", "probability_weight"):
        raise ValueError("Krüger IEAD table schema mismatch")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    source_pdf_sha = metadata.get("source_pdf_sha256", "")
    if source_pdf_sha != (
            "65b7750b2b773c3725d8f09f778b5b728ce9974a4548a5d522d19256f6bf9a51"):
        raise ValueError("Krüger IEAD metadata source checksum mismatch")
    output = Krueger2024DigitizedIEAD(
        table["energy_eV"], table["signed_angle_deg"],
        table["probability_weight"], source_pdf_sha,
        table_sha, metadata_sha, metadata)
    expected_mean = metadata["resampled_summary"]["mean_energy_eV"]
    if not np.isclose(output.mean_energy_eV, expected_mean, rtol=0.0, atol=1e-9):
        raise ValueError("Krüger IEAD table and metadata summaries disagree")
    return output


def load_krueger_2024_transfer_boundary_data(directory):
    """Load Figure-16 transfer flux/IEAD tables without exposing held-out profiles."""
    directory = Path(directory)
    flux_path = directory / "digitized_figure16a_transfer_fluxes.csv"
    iead_path = directory / "digitized_figure16b_power_ieads.csv"
    metadata_path = directory / "digitized_figure16_metadata.json"
    checksums = tuple(sha256(path.read_bytes()).hexdigest() for path in (
        flux_path, iead_path, metadata_path))
    if checksums != (
            KRUEGER_2024_TRANSFER_FLUX_SHA256,
            KRUEGER_2024_TRANSFER_IEAD_SHA256,
            KRUEGER_2024_TRANSFER_METADATA_SHA256):
        raise ValueError("Krüger Figure-16 transfer-boundary checksum mismatch")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("source_pdf_sha256") != (
            "65b7750b2b773c3725d8f09f778b5b728ce9974a4548a5d522d19256f6bf9a51"):
        raise ValueError("Krüger Figure-16 source checksum mismatch")

    flux = {}
    interval = {}
    with flux_path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        expected = {
            "oxygen_to_fluorocarbon_ratio", "species", "flux_cm2_s",
            "figure_digitized_flux_cm2_s", "selected_value_source",
            "digitization_lower_cm2_s", "digitization_upper_cm2_s",
        }
        if set(reader.fieldnames or ()) != expected:
            raise ValueError("Krüger Figure-16 flux-table schema mismatch")
        for row in reader:
            ratio = float(row["oxygen_to_fluorocarbon_ratio"])
            name = str(row["species"])
            flux.setdefault(ratio, {})[name] = float(row["flux_cm2_s"]) * 1e4
            interval.setdefault(ratio, {})[name] = (
                float(row["digitization_lower_cm2_s"]) * 1e4,
                float(row["digitization_upper_cm2_s"]) * 1e4,
            )

    table = np.genfromtxt(iead_path, delimiter=",", names=True)
    if table.dtype.names != (
            "low_frequency_power_kw", "energy_eV", "signed_angle_deg",
            "probability_weight"):
        raise ValueError("Krüger Figure-16 IEAD-table schema mismatch")
    iead = {}
    for power in np.unique(table["low_frequency_power_kw"]):
        selected = table["low_frequency_power_kw"] == power
        iead[float(power)] = Krueger2024DigitizedIEAD(
            table["energy_eV"][selected],
            table["signed_angle_deg"][selected],
            table["probability_weight"][selected],
            metadata["source_pdf_sha256"], checksums[1], checksums[2],
            dict(metadata, selected_low_frequency_power_kw=float(power)))
    return Krueger2024TransferBoundaryData(
        flux, interval, iead, checksums[0], checksums[1], checksums[2], metadata)


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


def _direction_marginalized_thermal_flux_species(
        record, temperature_K, polar_order, azimuthal_order):
    """Integrate neutral speed out of a flux-weighted half Maxwellian exactly.

    A Maxwellian incident *flux* factorizes in spherical velocity coordinates into a speed
    density and the angular density ``p(mu, phi) = 2*mu/(2*pi)``, where ``mu`` is the incident
    normal direction cosine.  Krüger's current reduced surface mechanism is energy independent
    for neutral reactants and its field-free ballistic transport consumes only ray direction.
    Under those two declared restrictions, replacing the speed distribution by its mean energy
    and integrating only the angular marginal is an exact model reduction, not a colder gas or a
    narrowed IADF.  The continuous Maxwellian density is retained for any sampling-based path.

    The polar integral uses Gauss--Legendre nodes on ``mu in [0, 1]`` with the physical ``2*mu``
    weight; the periodic azimuth uses an equally spaced trapezoidal rule.  Endpoint operator and
    angular-order refinement audits remain required because hard visibility is discontinuous.
    """
    values = np.asarray(
        [temperature_K, polar_order, azimuthal_order], dtype=float)
    if (np.any(~np.isfinite(values)) or temperature_K <= 0.0
            or int(polar_order) != polar_order or int(azimuthal_order) != azimuthal_order
            or int(polar_order) < 2 or int(azimuthal_order) < 4):
        raise ValueError(
            "direction-marginalized neutral quadrature requires polar order >= 2 "
            "and azimuthal order >= 4")
    polar_order = int(polar_order)
    azimuthal_order = int(azimuthal_order)
    temperature_eV = float(temperature_K) * _BOLTZMANN_EV_PER_K

    legendre_node, legendre_weight = np.polynomial.legendre.leggauss(polar_order)
    mu = 0.5 * (legendre_node + 1.0)
    mu_weight = legendre_weight * mu  # (w/2) * 2*mu
    phi = 2.0 * np.pi * np.arange(azimuthal_order, dtype=float) / azimuthal_order
    mu_grid, phi_grid = np.meshgrid(mu, phi, indexing="ij")
    tangent = np.sqrt(np.maximum(1.0 - mu_grid ** 2, 0.0))
    direction = np.column_stack((
        (tangent * np.cos(phi_grid)).ravel(),
        (tangent * np.sin(phi_grid)).ravel(),
        mu_grid.ravel(),
    ))
    weight = np.repeat(mu_weight / azimuthal_order, azimuthal_order)

    # In the engine's velocity convention |v|^2 is energy in eV.  The incident-flux
    # Maxwellian has <E> = 2*T, so this representative speed preserves its first energy moment.
    velocity = np.sqrt(2.0 * temperature_eV) * direction
    first_direction = np.einsum("s,sc->c", weight, direction)
    second_direction = np.einsum("s,si,sj->ij", weight, direction, direction)
    expected_first = np.array([0.0, 0.0, 2.0 / 3.0])
    expected_second = np.diag([0.25, 0.25, 0.5])
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
            "numerical_quadrature": {
                "method": "analytic_speed_marginal_plus_angular_quadrature",
                "polar_rule": "gauss_legendre_in_mu_with_2mu_flux_weight",
                "polar_order": polar_order,
                "azimuthal_rule": "periodic_trapezoidal",
                "azimuthal_order": azimuthal_order,
                "node_count": int(weight.size),
                "analytically_marginalized_speed": True,
                "representative_energy_eV": float(2.0 * temperature_eV),
                "maximum_first_direction_moment_error": float(
                    np.max(np.abs(first_direction - expected_first))),
                "maximum_second_direction_moment_error": float(
                    np.max(np.abs(second_direction - expected_second))),
                "validity_domain": (
                    "field-free neutral transport with energy-independent neutral surface laws"),
                "hard_visibility_refinement_required": True,
            },
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
    boundary_fluxes = load_krueger_2024_base_boundary_fluxes(directory)
    source_path = directory / "base_case_boundary_fluxes.csv"
    digest = sha256(source_path.read_bytes()).hexdigest()
    records = []
    for item in boundary_fluxes:
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


def build_krueger_2024_development_boundary(
        directory, *, reference_plane_m, effective_ion_mass_amu=39.948,
        ion_mixture_closure=(
            "singly charged aggregate uses Ar mass for trajectory time scaling; "
            "electrostatic path at fixed energy/charge is mass independent"),
        neutral_temperature_K=300.0,
        n_transverse_neutral=5, n_normal_neutral=8,
        neutral_direction_polar_order=None, neutral_direction_azimuthal_order=None,
        ion_energy_bin_eV=None, ion_angle_bin_deg=None,
        ion_azimuthal_closure=None, ion_azimuthal_order=16):
    """Build the complete published Krüger flux/IEAD boundary as an explicit development closure.

    The HPEM table publishes only an aggregate positive-ion flux and a combined IEAD. This helper
    therefore cannot produce a predictive boundary: it retains the aggregate ``ions`` population,
    requires an explicit effective-mass closure, and provides only the source-plane 2-D sampler.
    """
    values = np.asarray(
        [reference_plane_m, effective_ion_mass_amu, neutral_temperature_K],
        dtype=float)
    if (np.any(~np.isfinite(values)) or reference_plane_m < 0.0
            or effective_ion_mass_amu <= 0.0 or neutral_temperature_K <= 0.0
            or not str(ion_mixture_closure).strip()
            or int(n_transverse_neutral) != n_transverse_neutral
            or int(n_normal_neutral) != n_normal_neutral
            or n_transverse_neutral <= 0 or n_normal_neutral <= 0):
        raise ValueError("invalid Krüger development-boundary controls")
    deck = load_krueger_2024_reactor_flux_deck(directory)
    iead = load_krueger_2024_digitized_iead(directory)
    directional = (
        neutral_direction_polar_order is not None
        or neutral_direction_azimuthal_order is not None)
    if directional and (
            neutral_direction_polar_order is None
            or neutral_direction_azimuthal_order is None):
        raise ValueError(
            "direction-marginalized neutral quadrature requires both angular orders")
    if directional:
        species = [
            _direction_marginalized_thermal_flux_species(
                record, neutral_temperature_K,
                neutral_direction_polar_order, neutral_direction_azimuthal_order)
            for record in deck.species_fluxes if record.role == "neutral"
        ]
        neutral_quadrature = {
            "method": "analytic_speed_marginal_plus_angular_quadrature",
            "polar_order": int(neutral_direction_polar_order),
            "azimuthal_order": int(neutral_direction_azimuthal_order),
            "nodes_per_species": int(
                neutral_direction_polar_order * neutral_direction_azimuthal_order),
            "validity_domain": (
                "field-free neutral transport with energy-independent neutral surface laws"),
        }
    else:
        species = [
            _thermal_flux_species(
                record, neutral_temperature_K,
                int(n_transverse_neutral), int(n_normal_neutral))
            for record in deck.species_fluxes if record.role == "neutral"
        ]
        neutral_quadrature = {
            "method": "tensor_hermite_laguerre_velocity_quadrature",
            "transverse_order": int(n_transverse_neutral),
            "normal_order": int(n_normal_neutral),
            "nodes_per_species": int(
                n_transverse_neutral * n_transverse_neutral * n_normal_neutral),
        }
    ion = deck.get("ions")
    species.append(iead.development_species(
        ion.flux_m2_s, effective_mass_amu=effective_ion_mass_amu,
        mixture_closure=ion_mixture_closure, name=ion.name,
        energy_bin_eV=ion_energy_bin_eV, angle_bin_deg=ion_angle_bin_deg,
        azimuthal_closure=ion_azimuthal_closure,
        azimuthal_order=ion_azimuthal_order))
    return PlasmaBoundaryState(
        tuple(species), reference_plane_m=float(reference_plane_m),
        provenance={
            "provider": "krueger_2024_published_flux_and_digitized_iead",
            "claim_mode": "development",
            "supports_prediction": False,
            "source": deck.source,
            "source_sha256": deck.source_sha256,
            "selected_species": tuple(item.name for item in species),
            "complete_published_flux_table_used": True,
            "aggregate_ion_mixture_unresolved": True,
            "effective_ion_mass_amu": float(effective_ion_mass_amu),
            "ion_mixture_closure": str(ion_mixture_closure),
            "ion_azimuthal_closure": ion_azimuthal_closure,
            "ion_azimuthal_order": (
                None if ion_azimuthal_closure is None else int(ion_azimuthal_order)),
            "neutral_temperature_K": float(neutral_temperature_K),
            "neutral_quadrature": neutral_quadrature,
            "ion_quadrature": dict(species[-1].provenance["numerical_quadrature"]),
            "total_boundary_quadrature_nodes": int(sum(
                item.velocity_sqrt_eV.shape[0] for item in species)),
            "digitized_iead_table_sha256": iead.table_sha256,
            "digitization_metadata_sha256": iead.metadata_sha256,
            "predictive_blockers": (
                "species-resolved positive-ion/hot-neutral composition is unpublished",
                "the IEAD is HPEM output rather than a measured wafer distribution",
                "3-D ion azimuth is a symmetry closure rather than a measured distribution",
            ),
        })


def build_krueger_2024_transfer_boundary(
        directory, *, reference_plane_m, low_frequency_power_kw,
        oxygen_to_fluorocarbon_ratio=None, effective_ion_mass_amu=39.948,
        ion_mixture_closure=(
            "singly charged aggregate uses Ar mass for trajectory time scaling; "
            "electrostatic path at fixed energy/charge is mass independent"),
        neutral_temperature_K=300.0,
        neutral_direction_polar_order=8, neutral_direction_azimuthal_order=16,
        ion_energy_bin_eV=None, ion_angle_bin_deg=None,
        ion_azimuthal_closure=None, ion_azimuthal_order=16):
    """Build one Figure-16 transfer boundary while retaining its HPEM evidence grade.

    Oxygen-ratio cases use the Figure-16(a) species flux vector and the 6 kW IEAD.  Power-sweep
    cases (``oxygen_to_fluorocarbon_ratio=None``) use Table-I base fluxes because Figure 16 publishes
    only the power-dependent IEADs; that constant-flux closure is explicit in provenance and prevents
    a quantitative reactor-boundary claim.  In either mode the held-out *profiles* remain unread.
    """
    power = float(low_frequency_power_kw)
    if power not in {0.0, 4.0, 6.0, 8.0}:
        raise ValueError("Krüger transfer power must be one of 0, 4, 6, or 8 kW")
    transfer = load_krueger_2024_transfer_boundary_data(directory)
    base = load_krueger_2024_reactor_flux_deck(directory)
    base_by_name = {item.name: item for item in base.species_fluxes}
    expected_names = set(base_by_name)
    if oxygen_to_fluorocarbon_ratio is None:
        selected_flux = {
            item.name: item.flux_m2_s for item in base.species_fluxes}
        flux_source = "Table I base fluxes held constant across the power sweep"
        flux_interval = None
    else:
        ratio = float(oxygen_to_fluorocarbon_ratio)
        if ratio not in transfer.flux_m2_s_by_oxygen_ratio:
            raise ValueError("Krüger oxygen ratio must be one of 0.5, 1, 1.5, or 2.5")
        if power != 6.0:
            raise ValueError("Figure-16 oxygen-ratio fluxes belong to the 6 kW process")
        selected_flux = dict(transfer.flux_m2_s_by_oxygen_ratio[ratio])
        flux_interval = dict(
            transfer.flux_interval_m2_s_by_oxygen_ratio[ratio])
        flux_source = f"Figure 16(a), O2/C4F6={ratio:g}"
    if set(selected_flux) != expected_names:
        raise RuntimeError("Figure-16 boundary species disagree with the base reactor deck")

    neutral_records = []
    for name, template in base_by_name.items():
        if template.role != "neutral":
            continue
        neutral_records.append(ReactorSpeciesFlux(
            name, selected_flux[name], template.role,
            "published_reactor_model_output", flux_source,
            charge_number=template.charge_number, mass_amu=template.mass_amu))
    species = [
        _direction_marginalized_thermal_flux_species(
            record, neutral_temperature_K,
            neutral_direction_polar_order, neutral_direction_azimuthal_order)
        for record in neutral_records
    ]
    ion_template = base_by_name["ions"]
    species.append(transfer.iead_by_low_frequency_power_kw[power].development_species(
        selected_flux["ions"], effective_mass_amu=effective_ion_mass_amu,
        mixture_closure=ion_mixture_closure, name=ion_template.name,
        energy_bin_eV=ion_energy_bin_eV, angle_bin_deg=ion_angle_bin_deg,
        azimuthal_closure=ion_azimuthal_closure,
        azimuthal_order=ion_azimuthal_order))
    return PlasmaBoundaryState(
        tuple(species), reference_plane_m=float(reference_plane_m),
        provenance={
            "provider": "krueger_2024_figure16_transfer_boundary",
            "claim_mode": "held_out_development_transfer",
            "supports_prediction": False,
            "held_out_profiles_used_to_construct_boundary": False,
            "low_frequency_power_kw": power,
            "oxygen_to_fluorocarbon_ratio": (
                None if oxygen_to_fluorocarbon_ratio is None
                else float(oxygen_to_fluorocarbon_ratio)),
            "flux_source": flux_source,
            "flux_digitization_interval_m2_s": flux_interval,
            "power_sweep_constant_flux_closure": (
                oxygen_to_fluorocarbon_ratio is None),
            "effective_ion_mass_amu": float(effective_ion_mass_amu),
            "ion_mixture_closure": str(ion_mixture_closure),
            "ion_azimuthal_closure": ion_azimuthal_closure,
            "ion_azimuthal_order": (
                None if ion_azimuthal_closure is None else int(ion_azimuthal_order)),
            "neutral_temperature_K": float(neutral_temperature_K),
            "neutral_quadrature": {
                "method": "analytic_speed_marginal_plus_angular_quadrature",
                "polar_order": int(neutral_direction_polar_order),
                "azimuthal_order": int(neutral_direction_azimuthal_order),
                "nodes_per_species": int(
                    neutral_direction_polar_order
                    * neutral_direction_azimuthal_order),
            },
            "ion_quadrature": dict(
                species[-1].provenance["numerical_quadrature"]),
            "total_boundary_quadrature_nodes": int(sum(
                item.weight.size for item in species)),
            "figure16_flux_table_sha256": transfer.flux_table_sha256,
            "figure16_iead_table_sha256": transfer.iead_table_sha256,
            "figure16_metadata_sha256": transfer.metadata_sha256,
            "predictive_blockers": (
                "Figure-16 boundary values are HPEM outputs, not measurements",
                "species-resolved positive-ion/hot-neutral composition is unpublished",
                "3-D ion azimuth is a symmetry closure rather than a measured distribution",
                ("power-dependent wafer fluxes are unpublished and held at Table-I values"
                 if oxygen_to_fluorocarbon_ratio is None else
                 "Figure-16 flux digitization and HPEM model uncertainty remain"),
            ),
        })


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


def append_global_current_balance_maxwellian_electrons(
        boundary: PlasmaBoundaryState, *, electron_temperature_eV: float,
        temperature_source: str, temperature_evidence_kind: str,
        electron_name: str = "electron", n_transverse: int = 5,
        n_normal: int = 8):
    """Append the thermal-electron closure used by MCFPM-style charging models.

    The incident electron *number* flux is chosen so the charge-weighted, time-averaged flux at
    the feature reference plane is globally neutral.  The electron energy and angle law is the
    flux-weighted half Maxwellian produced by :func:`maxwellian_electron_boundary_state`; its
    cosine angular marginal is the Lambertian closure described by Wang and Kushner (2010).
    Local currents are deliberately not balanced here.  They remain kinetic outputs of trajectory
    transport through the self-consistent feature field.

    This helper adds a boundary current only.  It never inserts a Boltzmann electron density into
    the feature-volume Poisson equation.
    """
    if not isinstance(boundary, PlasmaBoundaryState):
        raise TypeError("boundary must be a PlasmaBoundaryState")
    if (
        not np.isfinite(electron_temperature_eV)
        or electron_temperature_eV <= 0.0
        or not str(temperature_source).strip()
        or not str(temperature_evidence_kind).strip()
        or not str(electron_name).strip()
        or int(n_transverse) != n_transverse
        or int(n_normal) != n_normal
        or n_transverse <= 0
        or n_normal <= 0
    ):
        raise ValueError("invalid global-current-balance electron closure")
    if electron_name in {item.name for item in boundary.species}:
        raise ValueError(f"boundary already contains species {electron_name!r}")

    positive_charge_flux = float(sum(
        max(int(item.charge_number), 0) * float(item.flux_m2_s)
        for item in boundary.species))
    negative_charge_flux = float(sum(
        max(-int(item.charge_number), 0) * float(item.flux_m2_s)
        for item in boundary.species))
    electron_flux = positive_charge_flux - negative_charge_flux
    if not np.isfinite(electron_flux) or electron_flux <= 0.0:
        raise ValueError(
            "existing boundary has no positive net charge flux for a balancing electron source")

    base = maxwellian_electron_boundary_state(
        float(electron_temperature_eV), electron_flux,
        n_transverse=int(n_transverse), n_normal=int(n_normal),
        reference_plane_m=boundary.reference_plane_m,
        electron_name=str(electron_name),
    ).get(str(electron_name))
    electron = _with_provenance(base, dict(
        base.provenance,
        provider="global_current_balance_maxwellian_electrons",
        role="charge_carrier",
        flux_closure="time_averaged_global_charge_neutrality",
        positive_charge_flux_per_m2_s=positive_charge_flux,
        preexisting_negative_charge_flux_per_m2_s=negative_charge_flux,
        electron_flux_m2_s=electron_flux,
        electron_temperature_eV=float(electron_temperature_eV),
        temperature_source=str(temperature_source),
        temperature_evidence_kind=str(temperature_evidence_kind),
        angular_distribution="Lambertian_cosine_flux_marginal",
        local_current_balance="self_consistent_kinetic_feature_transport",
        volume_boltzmann_electron_term=False,
    ))
    return PlasmaBoundaryState(
        species=tuple(boundary.species) + (electron,),
        reference_plane_m=boundary.reference_plane_m,
        provenance=dict(
            boundary.provenance,
            global_current_balance_electron_closure={
                "electron_name": str(electron_name),
                "electron_flux_m2_s": electron_flux,
                "electron_temperature_eV": float(electron_temperature_eV),
                "temperature_source": str(temperature_source),
                "temperature_evidence_kind": str(temperature_evidence_kind),
                "angular_distribution": "Lambertian_cosine_flux_marginal",
                "local_balance_is_not_imposed": True,
                "volume_boltzmann_electron_term": False,
            },
        ),
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
