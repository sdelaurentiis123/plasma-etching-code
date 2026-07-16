"""Adapters from experimental diagnostics to the common plasma-boundary contract.

Experimental measurements and missing closures remain separate. In particular, a self-bias voltage
is not silently converted into an IEDF, and an integrated radical/ion ratio is not silently promoted
to a species-resolved flux vector.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .boundary_state import (
    EnergyCosineAngleDensity2D, IonEnergyTransverseDensity2D,
    IonEnergyTransverseMaxwellianDensity, MaxwellianFluxVelocityDensity,
    PlasmaBoundaryState, SpeciesBoundaryState,
    maxwellian_electron_boundary_state,
)
from .experimental_data import (
    Jeon2022ElectronBiasControl, Jeon2022PlasmaControl,
    jeon_2022_bohm_ion_flux_m2_s,
)
from .sheath import bohm_speed


HWANG_GIAPIS_1997_IEDF_SHA256 = (
    "540601fc95bc85e5c906d9d3e5d566f966e2761d4c65fc8a8167aaaf4c28adea")
HWANG_GIAPIS_1997_EEDF_SHA256 = (
    "17ae2728a0e3d5561fdd7b898d1d69a1b9a7267fa1c65a645daac25916119af6")
HWANG_GIAPIS_1997_PDF_SHA256 = (
    "30a6871d6416f27e8dbbb45e9eabbca79cddf7632872f8ed185a9e193832f63d")


def _load_hwang_giapis_1997_iedf(path, *, verify_checksum=True):
    path = Path(path)
    payload = path.read_bytes()
    digest = sha256(payload).hexdigest()
    if verify_checksum and digest != HWANG_GIAPIS_1997_IEDF_SHA256:
        raise ValueError("Hwang--Giapis Fig. 4(a) IEDF checksum mismatch")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    required = {
        "normal_energy_lower_eV", "normal_energy_upper_eV", "probability_mass",
        "digitized_curve_height_px", "source_pdf_sha256", "source_pdf_page",
        "source_figure"}
    if not rows or set(rows[0]) != required:
        raise ValueError("unexpected Hwang--Giapis IEDF schema")
    lower = np.asarray([float(row["normal_energy_lower_eV"]) for row in rows])
    upper = np.asarray([float(row["normal_energy_upper_eV"]) for row in rows])
    mass = np.asarray([float(row["probability_mass"]) for row in rows])
    if (np.any(~np.isfinite(lower)) or np.any(~np.isfinite(upper))
            or np.any(~np.isfinite(mass)) or np.any(mass < 0.0)
            or not np.allclose(lower[1:], upper[:-1], rtol=0.0, atol=1e-12)
            or np.any(upper <= lower) or not np.isclose(mass.sum(), 1.0, atol=2e-10)
            or {row["source_pdf_sha256"] for row in rows}
            != {HWANG_GIAPIS_1997_PDF_SHA256}
            or {row["source_pdf_page"] for row in rows} != {"4"}
            or {row["source_figure"] for row in rows} != {"Fig. 4(a)"}):
        raise ValueError("invalid Hwang--Giapis digitized IEDF rows")
    return np.concatenate((lower[:1], upper)), mass / mass.sum(), digest


def _load_hwang_giapis_1997_eedf(path, *, verify_checksum=True):
    path = Path(path)
    payload = path.read_bytes()
    digest = sha256(payload).hexdigest()
    if verify_checksum and digest != HWANG_GIAPIS_1997_EEDF_SHA256:
        raise ValueError("Hwang--Giapis Fig. 4(b) EEDF checksum mismatch")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    required = {
        "energy_lower_eV", "energy_upper_eV", "probability_mass",
        "digitized_curve_height_px", "source_pdf_sha256", "source_pdf_page",
        "source_figure"}
    if not rows or set(rows[0]) != required:
        raise ValueError("unexpected Hwang--Giapis EEDF schema")
    lower = np.asarray([float(row["energy_lower_eV"]) for row in rows])
    upper = np.asarray([float(row["energy_upper_eV"]) for row in rows])
    mass = np.asarray([float(row["probability_mass"]) for row in rows])
    if (np.any(~np.isfinite(lower)) or np.any(~np.isfinite(upper))
            or np.any(~np.isfinite(mass)) or np.any(mass < 0.0)
            or not np.allclose(lower[1:], upper[:-1], rtol=0.0, atol=1e-12)
            or np.any(upper <= lower) or not np.isclose(mass.sum(), 1.0, atol=2e-10)
            or {row["source_pdf_sha256"] for row in rows}
            != {HWANG_GIAPIS_1997_PDF_SHA256}
            or {row["source_pdf_page"] for row in rows} != {"4"}
            or {row["source_figure"] for row in rows} != {"Fig. 4(b)"}):
        raise ValueError("invalid Hwang--Giapis digitized EEDF rows")
    return np.concatenate((lower[:1], upper)), mass / mass.sum(), digest


def build_hwang_giapis_1997_boundary_state(
        iedf_csv_path, eedf_csv_path=None, *, reference_plane_m,
        plasma_density_m3=1.0e18,
        electron_temperature_eV=4.0, ion_tangential_temperature_eV=0.5,
        ion_mass_amu=35.45, n_transverse_ion=3,
        n_transverse_electron=5, n_normal_electron=8,
        verify_checksum=True):
    """Build the primary-source boundary used for the Nozawa notch replay.

    The ion normal-energy mass is digitized from Hwang--Giapis Fig. 4(a),
    rather than regenerated with the adjacent symmetric Child-sheath model.
    A 0.5 eV transverse Maxwellian makes the 2-D IADF energy-dependent (the
    low-energy wings are broader), as stated in Sec. III B.  When
    ``eedf_csv_path`` is supplied, the paper's Fig. 4(b) EEDF and explicit
    Fig. 5(b) ``cos(theta)**0.6`` EADF fit define the 2-D electron source.

    The ordinary three-dimensional density models remain explicit closures for
    three-dimensional consumers; a 2-D solver consumes ``density_model_2d``
    directly and never folds an unmodeled out-of-plane energy into its plane.
    Ion and electron particle fluxes are equal, with their absolute value
    derived from the declared density and Bohm speed.
    """
    values = np.asarray([
        reference_plane_m, plasma_density_m3, electron_temperature_eV,
        ion_tangential_temperature_eV, ion_mass_amu], dtype=float)
    if (np.any(~np.isfinite(values)) or reference_plane_m < 0.0
            or plasma_density_m3 <= 0.0 or electron_temperature_eV <= 0.0
            or ion_tangential_temperature_eV <= 0.0 or ion_mass_amu <= 0.0
            or int(n_transverse_ion) != n_transverse_ion or n_transverse_ion <= 0
            or int(n_transverse_electron) != n_transverse_electron
            or n_transverse_electron <= 0
            or int(n_normal_electron) != n_normal_electron
            or n_normal_electron <= 0):
        raise ValueError("invalid Hwang--Giapis boundary controls")
    edges, energy_mass, iedf_digest = _load_hwang_giapis_1997_iedf(
        iedf_csv_path, verify_checksum=verify_checksum)
    eedf_edges = None
    eedf_mass = None
    eedf_digest = None
    if eedf_csv_path is not None:
        eedf_edges, eedf_mass, eedf_digest = _load_hwang_giapis_1997_eedf(
            eedf_csv_path, verify_checksum=verify_checksum)
    node, node_weight = np.polynomial.hermite.hermgauss(int(n_transverse_ion))
    transverse = np.sqrt(float(ion_tangential_temperature_eV)) * node
    transverse_weight = node_weight / np.sqrt(np.pi)
    energy = 0.5 * (edges[:-1] + edges[1:])
    ix, iy, iz = np.meshgrid(
        np.arange(node.size), np.arange(node.size), np.arange(energy.size),
        indexing="ij")
    velocity = np.column_stack((
        transverse[ix.ravel()], transverse[iy.ravel()],
        np.sqrt(energy[iz.ravel()])))
    weight = (transverse_weight[ix.ravel()] * transverse_weight[iy.ravel()]
              * energy_mass[iz.ravel()])
    flux = float(plasma_density_m3) * bohm_speed(
        electron_temperature_eV, ion_mass_amu)
    shared_source = {
        "source": "Hwang & Giapis, JVST B 15, 70 (1997)",
        "doi": "10.1116/1.589258",
        "pressure_mTorr": 3.0,
        "rf_frequency_hz": 4.0e5,
        "rf_bias_peak_to_peak_v": 60.0,
        "mean_sheath_voltage_v": 37.0,
        "plasma_density_m3": float(plasma_density_m3),
        "electron_temperature_eV": float(electron_temperature_eV),
    }
    ion = SpeciesBoundaryState(
        "Cl+", 1, float(ion_mass_amu), flux, velocity, weight,
        density_model=IonEnergyTransverseMaxwellianDensity(
            edges, energy_mass, float(ion_tangential_temperature_eV)),
        provenance=dict(
            shared_source, role="digitized_nonlinear_sheath_iedf",
            iedf_csv_sha256=iedf_digest,
            source_pdf_sha256=HWANG_GIAPIS_1997_PDF_SHA256,
            source_figure="Fig. 4(a)", source_pdf_page=4,
            ion_tangential_temperature_eV=float(ion_tangential_temperature_eV),
            reported_iadf_hwhm_deg=4.3,
            energy_dependent_iadf=True,
            two_dimensional_projection=(
                "one in-plane Maxwellian tangent plus digitized normal energy"),
            supports_prediction_within_declared_benchmark=True),
        density_model_2d=IonEnergyTransverseDensity2D(
            edges, energy_mass, float(ion_tangential_temperature_eV)))
    electron = maxwellian_electron_boundary_state(
        electron_temperature_eV, flux,
        n_transverse=int(n_transverse_electron), n_normal=int(n_normal_electron),
        electron_name="electron", reference_plane_m=reference_plane_m).get("electron")
    electron_density_2d = (
        None if eedf_edges is None
        else EnergyCosineAngleDensity2D(eedf_edges, eedf_mass, 0.6))
    electron = SpeciesBoundaryState(
        electron.name, electron.charge_number, electron.mass_amu, electron.flux_m2_s,
        electron.velocity_sqrt_eV, electron.weight,
        density_model=electron.density_model,
        provenance=dict(
            shared_source,
            role=(
                "digitized_eedf_and_analytic_eadf_in_2d"
                if electron_density_2d is not None
                else "analytic_half_maxwellian_closure"),
            source_figure=(
                ["Fig. 4(b)", "Fig. 5(b)"]
                if electron_density_2d is not None else None),
            eedf_csv_sha256=eedf_digest,
            eadf_cosine_power=(0.6 if electron_density_2d is not None else None),
            three_dimensional_density_model="analytic_half_maxwellian_closure",
            two_dimensional_density_model=(
                "digitized_energy_times_cosine_power_angle"
                if electron_density_2d is not None else "legacy_3d_projection"),
            equal_particle_flux_to_ions=True,
            supports_prediction_within_declared_benchmark=(
                electron_density_2d is not None)),
        density_model_2d=electron_density_2d)
    return PlasmaBoundaryState(
        (ion, electron), float(reference_plane_m),
        provenance=dict(
            shared_source, model="Hwang--Giapis 1997 plasma-to-feature boundary",
            equal_ion_electron_particle_flux=True,
            ion_flux_m2_s=flux, iedf_csv_sha256=iedf_digest,
            eedf_csv_sha256=eedf_digest,
            two_dimensional_source_is_source_faithful=(
                electron_density_2d is not None),
            reference_plane_height_above_sio2_m=float(reference_plane_m),
            declared_reference_height_in_paper_m=3.7e-6))


@dataclass(frozen=True)
class Jeon2022BoundaryClosure:
    ion_name: str
    ion_mass_amu: float
    ion_normal_energy_eV: np.ndarray
    ion_normal_energy_weight: np.ndarray
    ion_tangential_temperature_eV: float
    neutral_flux_fraction: Mapping[str, float]
    neutral_mass_amu: Mapping[str, float]
    neutral_temperature_K: float
    provenance: Mapping[str, object]
    supports_prediction_within_declared_domain: bool = False

    def __post_init__(self):
        energy = np.asarray(self.ion_normal_energy_eV, dtype=float).copy()
        weight = np.asarray(self.ion_normal_energy_weight, dtype=float).copy()
        fraction = dict(self.neutral_flux_fraction)
        mass = dict(self.neutral_mass_amu)
        if (not self.ion_name or self.ion_mass_amu <= 0.0 or energy.ndim != 1
                or weight.shape != energy.shape or energy.size == 0
                or np.any(~np.isfinite(energy)) or np.any(energy < 0.0)
                or np.any(~np.isfinite(weight)) or np.any(weight < 0.0) or weight.sum() <= 0.0
                or self.ion_tangential_temperature_eV <= 0.0
                or self.neutral_temperature_K <= 0.0 or not fraction
                or set(fraction) != set(mass) or any(not name for name in fraction)
                or any(not np.isfinite(value) or value < 0.0 for value in fraction.values())
                or any(not np.isfinite(value) or value <= 0.0 for value in mass.values())
                or not np.isclose(sum(fraction.values()), 1.0, rtol=0.0, atol=2e-13)
                or not isinstance(self.supports_prediction_within_declared_domain, bool)):
            raise ValueError("invalid Jeon boundary closure")
        weight /= weight.sum()
        energy.setflags(write=False); weight.setflags(write=False)
        object.__setattr__(self, "ion_normal_energy_eV", energy)
        object.__setattr__(self, "ion_normal_energy_weight", weight)
        object.__setattr__(self, "neutral_flux_fraction", MappingProxyType(fraction))
        object.__setattr__(self, "neutral_mass_amu", MappingProxyType(mass))
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))


def _same_jeon_condition(left, right):
    return (left.condition_family == right.condition_family
            and left.c4f8_fraction == right.c4f8_fraction
            and left.pulse_off_ms == right.pulse_off_ms)


def build_jeon_2022_boundary_state(
        plasma_control: Jeon2022PlasmaControl,
        electron_bias_control: Jeon2022ElectronBiasControl,
        closure: Jeon2022BoundaryClosure, *, reference_plane_m, n_transverse_ion=3,
        n_transverse_neutral=5, n_normal_neutral=8):
    """Build one dimensional boundary state from measurements plus an explicit closure."""
    if (not isinstance(plasma_control, Jeon2022PlasmaControl)
            or not isinstance(electron_bias_control, Jeon2022ElectronBiasControl)
            or not isinstance(closure, Jeon2022BoundaryClosure)):
        raise TypeError("Jeon boundary construction requires typed evidence and closure")
    if not _same_jeon_condition(plasma_control, electron_bias_control):
        raise ValueError("Jeon plasma controls refer to different experimental conditions")
    if not np.isfinite(reference_plane_m):
        raise ValueError("reference_plane_m must be finite")
    if int(n_transverse_ion) <= 0 or int(n_transverse_neutral) <= 0 or int(n_normal_neutral) <= 0:
        raise ValueError("boundary quadrature orders must be positive")

    ion_flux = jeon_2022_bohm_ion_flux_m2_s(electron_bias_control)
    nodes, gh_weight = np.polynomial.hermite.hermgauss(int(n_transverse_ion))
    transverse = np.sqrt(closure.ion_tangential_temperature_eV) * nodes
    transverse_weight = gh_weight / np.sqrt(np.pi)
    ix, iy, iz = np.meshgrid(
        np.arange(nodes.size), np.arange(nodes.size),
        np.arange(closure.ion_normal_energy_eV.size), indexing="ij")
    ion_velocity = np.column_stack((
        transverse[ix.ravel()], transverse[iy.ravel()],
        np.sqrt(closure.ion_normal_energy_eV[iz.ravel()])))
    ion_weight = (transverse_weight[ix.ravel()] * transverse_weight[iy.ravel()]
                  * closure.ion_normal_energy_weight[iz.ravel()])
    ion = SpeciesBoundaryState(
        closure.ion_name, 1, closure.ion_mass_amu, ion_flux,
        ion_velocity, ion_weight,
        provenance={
            "role": "explicit_iedf_closure",
            "closure": dict(closure.provenance),
            "supports_prediction": closure.supports_prediction_within_declared_domain,
        })

    temperature_eV = closure.neutral_temperature_K * 8.617333262145e-5
    hermite_node, hermite_weight = np.polynomial.hermite.hermgauss(int(n_transverse_neutral))
    laguerre_node, laguerre_weight = np.polynomial.laguerre.laggauss(int(n_normal_neutral))
    nx, ny, nz = np.meshgrid(
        np.arange(hermite_node.size), np.arange(hermite_node.size),
        np.arange(laguerre_node.size), indexing="ij")
    neutral_velocity = np.column_stack((
        np.sqrt(temperature_eV) * hermite_node[nx.ravel()],
        np.sqrt(temperature_eV) * hermite_node[ny.ravel()],
        np.sqrt(temperature_eV * laguerre_node[nz.ravel()])))
    neutral_weight = (hermite_weight[nx.ravel()] * hermite_weight[ny.ravel()]
                      * laguerre_weight[nz.ravel()] / np.pi)
    total_neutral_flux = plasma_control.neutral_to_ion_flux_ratio * ion_flux
    neutral_species = tuple(SpeciesBoundaryState(
        name, 0, closure.neutral_mass_amu[name], total_neutral_flux * fraction,
        neutral_velocity, neutral_weight,
        density_model=MaxwellianFluxVelocityDensity(temperature_eV),
        provenance={
            "role": "species_composition_closure",
            "integrated_flux_source": plasma_control.source_location,
            "closure": dict(closure.provenance),
            "supports_prediction": closure.supports_prediction_within_declared_domain,
        }) for name, fraction in closure.neutral_flux_fraction.items())
    return PlasmaBoundaryState(
        (ion,) + neutral_species, float(reference_plane_m),
        provenance={
            "source": "Jeong_2022_diagnostics_plus_explicit_closure",
            "electron_density": electron_bias_control.source_location,
            "self_bias_energy_scale_v": electron_bias_control.self_bias_magnitude_v,
            "self_bias_is_not_iedf": True,
            "integrated_neutral_to_ion_ratio": plasma_control.neutral_to_ion_flux_ratio,
            "closure_supports_prediction": closure.supports_prediction_within_declared_domain,
        })
