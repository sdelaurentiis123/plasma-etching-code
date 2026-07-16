"""User-facing Nozawa/Hwang notching replay through the common 3-D engine.

This module converts one checksum-verified experimental row into an explicit
geometry, plasma boundary, electrostatic topology, surface mechanism, and run
manifest.  It deliberately supports the Hwang--Giapis grating-edge/open-area
cell first.  The shared-pad and individual-pad curves require the physical pad
collector geometry (or a separately evidenced terminal-current model), which
the published numeric curves do not specify; those families are refused rather
than replaced by a width-only fit.

The ``smoke`` mode executes the exact hard-visibility operator for one short
physical-time update.  It proves installation and operator wiring, not charge
saturation or experimental validation.  The ``experiment`` mode uses the
signed R2 quasi-static saturation gates and therefore stops with a replayable
checkpoint when the field is not stationary.  Both modes use the same public
``PhysicalChargingProcess`` and the same physics objects.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import sysconfig
import tempfile
from types import MappingProxyType
from typing import Mapping

import numpy as np
from scipy.ndimage import label
from scipy.spatial import cKDTree
from skimage.measure import marching_cubes

from .charging_coevolution_3d import SurfaceChargingSaturationError
from .charging_poisson_3d import NodalPoissonSystem3D
from .chlorine_poly_si import HwangGiapisClSiMechanism
from .experimental_boundary import build_hwang_giapis_1997_boundary_state
from .feature_step_3d import FeatureGeometry3D
from .hwang_giapis_scatter_3d import HwangGiapisSiO2ForwardScatter3D
from .notching_validation_3d import (
    NOZAWA_1995_NOTCH_CURVES_SHA256,
    NotchingBenchmarkProtocol3D,
    load_nozawa_1995_notch_observations,
)
from .physical_api import PhysicalChargingProcess
from .profile_observables_3d import measure_trench_profile_observables_3d
from .threed import extract_mesh_3d, reinit_narrow


NOZAWA_POLY_SI_MATERIAL_ID = 1
NOZAWA_PHOTORESIST_MATERIAL_ID = 2
NOZAWA_SIO2_MATERIAL_ID = 3
NOZAWA_REPLAY_SCHEMA = "petch-nozawa-1995-replay-v1"
NOZAWA_STATIONARITY_CONTRACT = "CCA-2026-07-13-R2"
NOZAWA_REPORTED_POLY_ETCH_RATE_M_S = 2325.0e-10 / 60.0

_SERIES_NUMBER = re.compile(r"(?:^|;)([LS])=([0-9.]+)um(?=;|$)")


def _sha256_path(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _charging_surface_mesh_fingerprint_3d(geometry: FeatureGeometry3D) -> str:
    """Bind a face-charge vector to the exact geometry and marching-cubes ordering."""
    verts, faces, _centroids, _areas = extract_mesh_3d(geometry.phi, geometry.dx)
    digest = sha256()
    for array, dtype in (
            (geometry.phi, "<f8"), (geometry.material_id, "<i8"),
            (verts, "<f8"), (faces, "<i8")):
        digest.update(np.ascontiguousarray(array, dtype=dtype).tobytes())
    digest.update(np.asarray(
        [geometry.dx, geometry.mesh_length_unit_m, *geometry.mesh_origin_m],
        dtype="<f8").tobytes())
    return digest.hexdigest()


def _json_value(value):
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, (float, np.floating)):
        value = float(value)
        if not np.isfinite(value):
            raise ValueError("run artifact contains a non-finite float")
        return value
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    raise TypeError(f"cannot serialize {type(value).__name__} into a replay artifact")


def _write_json(path: Path, payload) -> str:
    encoded = (json.dumps(_json_value(payload), indent=2, sort_keys=True) + "\n").encode()
    with tempfile.NamedTemporaryFile(
            dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", delete=False) as stream:
        temporary = Path(stream.name)
        stream.write(encoded)
    os.replace(temporary, path)
    return sha256(encoded).hexdigest()


def _atomic_npz(path: Path, **arrays):
    """Replace one restart artifact only after its compressed payload is complete."""
    path = Path(path)
    with tempfile.NamedTemporaryFile(
            dir=path.parent, prefix=f".{path.name}.", suffix=".npz", delete=False) as stream:
        temporary = Path(stream.name)
        np.savez_compressed(stream, **arrays)
    os.replace(temporary, path)


def _failure_diagnostic_value(value):
    """Serialize refusal diagnostics without laundering undefined ratios into finite values."""
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return "nan" if np.isnan(value) else ("+inf" if value > 0.0 else "-inf")
    if isinstance(value, np.ndarray):
        return _failure_diagnostic_value(value.tolist())
    if isinstance(value, Mapping):
        return {
            str(key): _failure_diagnostic_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_failure_diagnostic_value(item) for item in value]
    return value


def default_nozawa_data_path(filename: str) -> Path:
    """Locate bundled evidence in a checkout or a wheel ``data-files`` install."""
    if filename not in {
            "digitized_notch_curves.csv", "fig4a_ion_energy_distribution.csv"}:
        raise ValueError("unknown Nozawa replay evidence file")
    root = Path(__file__).resolve().parents[2]
    source = (
        root / "data/experimental/nozawa_1995" / filename
        if filename.startswith("digitized")
        else root / "data/experimental/hwang_giapis_1997" / filename)
    installed = (Path(sysconfig.get_path("data")) / "share/petch/data/experimental"
                 / ("nozawa_1995" if filename.startswith("digitized")
                    else "hwang_giapis_1997") / filename)
    for candidate in (source, installed):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        f"{filename} is unavailable; reinstall petch with bundled experimental data "
        "or pass an explicit evidence path")


@dataclass(frozen=True)
class Nozawa1995ReplayCondition3D:
    condition_id: str
    target_family: str
    line_width_um: float
    space_width_um: float
    open_area_width_um: float
    measured_notch_depth_um: float
    digitization_uncertainty_um: float
    split: str
    electrical_topology: str
    source_csv_sha256: str = NOZAWA_1995_NOTCH_CURVES_SHA256

    def __post_init__(self):
        values = np.asarray([
            self.line_width_um, self.space_width_um, self.open_area_width_um,
            self.measured_notch_depth_um, self.digitization_uncertainty_um], dtype=float)
        if (not self.condition_id or self.target_family != "open_area_width"
                or np.any(~np.isfinite(values)) or np.any(values[:3] <= 0.0)
                or self.measured_notch_depth_um < 0.0
                or self.digitization_uncertainty_um <= 0.0
                or self.split != "calibration"
                or self.electrical_topology != "four_line_full_mirror_cell"
                or len(self.source_csv_sha256) != 64):
            raise ValueError("invalid supported Nozawa replay condition")


def load_nozawa_1995_replay_condition(
        condition_id: str, *, observations_path=None) -> Nozawa1995ReplayCondition3D:
    """Load one supported open-area row without weakening the frozen evidence split."""
    path = (default_nozawa_data_path("digitized_notch_curves.csv")
            if observations_path is None else Path(observations_path))
    observations = load_nozawa_1995_notch_observations(
        path, expected_sha256=NOZAWA_1995_NOTCH_CURVES_SHA256)
    matches = [item for item in observations if item.condition_id == condition_id]
    if len(matches) != 1:
        raise KeyError(f"unknown Nozawa condition {condition_id!r}")
    observation = matches[0]
    if observation.target_family != "open_area_width":
        raise NotImplementedError(
            "the shared-pad and individual-pad rows need the physical pad collector "
            "geometry or an independently evidenced terminal-current model; a local "
            "line/space cell would erase the experiment's controlling topology")
    parsed = {name: float(value) for name, value in _SERIES_NUMBER.findall(
        observation.series_label + ";")}
    if set(parsed) != {"L", "S"}:
        raise ValueError("Nozawa series label does not declare line and space widths")
    return Nozawa1995ReplayCondition3D(
        observation.condition_id, observation.target_family,
        parsed["L"], parsed["S"], observation.control_value,
        observation.notch_depth_um, observation.digitization_uncertainty_um,
        observation.split, "four_line_full_mirror_cell",
        NOZAWA_1995_NOTCH_CURVES_SHA256)


@dataclass(frozen=True)
class Nozawa1995GeometryContract3D:
    geometry: FeatureGeometry3D
    requested_line_width_um: float
    requested_space_width_um: float
    requested_open_area_width_um: float
    realized_cell_width_um: float
    realized_open_area_width_um: float
    source_plane_um: float
    line_centers_um: tuple[float, float, float, float]
    conductor_component_count: int
    topology: str = "full_periodic_reflection_of_X-Y_half_cell"

    def __post_init__(self):
        if (not isinstance(self.geometry, FeatureGeometry3D)
                or self.conductor_component_count != 4
                or len(self.line_centers_um) != 4
                or self.realized_open_area_width_um <= 0.0):
            raise ValueError("invalid Nozawa geometry contract")

    @property
    def manifest(self):
        return MappingProxyType({
            "topology": self.topology,
            "requested_line_width_um": self.requested_line_width_um,
            "requested_space_width_um": self.requested_space_width_um,
            "requested_open_area_width_um": self.requested_open_area_width_um,
            "realized_cell_width_um": self.realized_cell_width_um,
            "realized_open_area_width_um": self.realized_open_area_width_um,
            "open_area_discretization_error_um": (
                self.realized_open_area_width_um - self.requested_open_area_width_um),
            "source_plane_um": self.source_plane_um,
            "line_centers_um": self.line_centers_um,
            "line_count": 4,
            "floating_conductor_components": self.conductor_component_count,
            "material_ids": {
                "poly_si": NOZAWA_POLY_SI_MATERIAL_ID,
                "photoresist": NOZAWA_PHOTORESIST_MATERIAL_ID,
                "sio2": NOZAWA_SIO2_MATERIAL_ID,
            },
            "vertical_stack_um": {
                "sio2": 0.1, "poly_si": 0.3, "photoresist": 1.0,
                "source_height_above_sio2": 3.7,
            },
            "boundary_replication": (
                "Hwang--Giapis Fig. 3 contains X and Y between mirror planes at the "
                "centers of an ordinary space and the varied open area. Reflecting that "
                "half-cell creates four physical lines in one full periodic cell, so "
                "particle and Poisson boundaries use the published topology exactly"),
        })


def _horizontal_line_levelset(x, centers, width):
    return np.maximum.reduce([
        0.5 * float(width) - np.abs(x - float(center)) for center in centers])


def _project_periodic_endpoint_planes(field):
    output = np.asarray(field, dtype=float).copy()
    for axis in (0, 1):
        first = [slice(None)] * output.ndim
        last = [slice(None)] * output.ndim
        first[axis] = 0
        last[axis] = -1
        first = tuple(first)
        last = tuple(last)
        seam = 0.5 * (output[first] + output[last])
        output[first] = seam
        output[last] = seam
    return output


def make_nozawa_1995_open_area_geometry_3d(
        condition: Nozawa1995ReplayCondition3D, *, dx_um=0.1,
        mesh_length_unit_m=1e-6) -> Nozawa1995GeometryContract3D:
    """Build the full periodic counterpart of Hwang--Giapis Fig. 3.

    The source paper places lines X and Y between mirror planes through the
    center of the ordinary space before X and the varied open area after Y.
    Reflecting that half-cell therefore produces four lines in the full periodic
    cell.  Periodicizing the unreflected half-cell would instead join an ordinary
    space directly to the open area and change both conductor collection and
    trajectory topology.
    """
    if not isinstance(condition, Nozawa1995ReplayCondition3D):
        raise TypeError("condition must be a supported Nozawa replay condition")
    dx = float(dx_um)
    if (not np.isfinite(dx) or dx <= 0.0 or dx > 0.1
            or not np.isfinite(mesh_length_unit_m) or mesh_length_unit_m <= 0.0):
        raise ValueError("Nozawa geometry requires 0 < dx <= 0.1 um")
    line_width = condition.line_width_um
    space_width = condition.space_width_um
    # Published half-domain:
    #   S/2 | X(L) | S | Y(L) | W/2
    # Its reflected periodic counterpart is 4L + 3S + W and contains
    # X, Y, mirror(Y), mirror(X).
    requested_cell = (
        4.0 * line_width + 3.0 * space_width + condition.open_area_width_um)
    x_intervals = max(4, int(round(requested_cell / dx)))
    cell_width = x_intervals * dx
    open_width = cell_width - 4.0 * line_width - 3.0 * space_width
    if open_width <= 0.5 * dx:
        raise ValueError("grid snapping eliminated the declared open area")
    source_plane = 0.1 + 3.7
    z_intervals = int(round(source_plane / dx))
    if not np.isclose(z_intervals * dx, source_plane, atol=1e-12, rtol=0.0):
        raise ValueError("dx must resolve the 0.1/0.3/1.0 um vertical benchmark stack")
    shape = (x_intervals + 1, 3, z_intervals + 1)
    x = np.arange(shape[0]) * dx
    y = np.arange(shape[1]) * dx
    z = np.arange(shape[2]) * dx
    X, _Y, Z = np.meshgrid(x, y, z, indexing="ij")
    half_cell = 0.5 * cell_width
    first_center = 0.5 * (space_width + line_width)
    second_center = first_center + line_width + space_width
    centers = (
        first_center,
        second_center,
        2.0 * half_cell - second_center,
        2.0 * half_cell - first_center)
    if (not np.all(np.diff(centers) > line_width)
            or not np.isclose(
                centers[2] - centers[1] - line_width,
                open_width, atol=5e-13, rtol=0.0)):
        raise RuntimeError("reflected Nozawa line centers do not realize the declared gaps")
    horizontal = _horizontal_line_levelset(X, centers, line_width)
    oxide_analytic = 0.1 - Z
    poly_analytic = np.minimum.reduce((horizontal, Z - 0.1, 0.4 - Z))
    resist_analytic = np.minimum.reduce((horizontal, Z - 0.4, 1.4 - Z))
    bandwidth = cell_width + source_plane
    oxide_phi = _project_periodic_endpoint_planes(
        reinit_narrow(oxide_analytic, dx, bandwidth))
    poly_phi = _project_periodic_endpoint_planes(
        reinit_narrow(poly_analytic, dx, bandwidth))
    resist_phi = _project_periodic_endpoint_planes(
        reinit_narrow(resist_analytic, dx, bandwidth))
    material_levelsets = {
        NOZAWA_POLY_SI_MATERIAL_ID: poly_phi,
        NOZAWA_PHOTORESIST_MATERIAL_ID: resist_phi,
        NOZAWA_SIO2_MATERIAL_ID: oxide_phi,
    }
    analytic = np.maximum.reduce((poly_phi, resist_phi, oxide_phi))
    phi = _project_periodic_endpoint_planes(
        reinit_narrow(analytic, dx, bandwidth))
    ids = np.asarray(tuple(material_levelsets), dtype=int)
    stack = np.stack([material_levelsets[int(item)] for item in ids])
    owner = ids[np.argmax(stack, axis=0)]
    material = np.where(phi >= 0.0, owner, 0)
    geometry = FeatureGeometry3D(
        phi, material, dx, mesh_length_unit_m,
        material_levelsets=material_levelsets)
    conductor, count = label(
        (geometry.material_id == NOZAWA_POLY_SI_MATERIAL_ID) & (geometry.phi >= 0.0))
    if count != 4:
        raise RuntimeError(
            f"Nozawa geometry resolved {count} poly components instead of four")
    # Both duplicated endpoint planes must describe one periodic physical plane.
    if (not np.array_equal(geometry.material_id[0], geometry.material_id[-1])
            or not np.allclose(geometry.phi[0], geometry.phi[-1], atol=5e-12, rtol=0.0)
            or not np.array_equal(geometry.material_id[:, 0], geometry.material_id[:, -1])
            or not np.allclose(
                geometry.phi[:, 0], geometry.phi[:, -1], atol=5e-12, rtol=0.0)):
        raise RuntimeError("Nozawa full-cell geometry is not periodic at duplicated endpoints")
    return Nozawa1995GeometryContract3D(
        geometry, line_width, space_width, condition.open_area_width_um,
        cell_width, open_width, source_plane, centers, count)


def _cell_center_average(field):
    field = np.asarray(field, dtype=float)
    return sum(
        field[i:i + field.shape[0] - 1,
              j:j + field.shape[1] - 1,
              k:k + field.shape[2] - 1]
        for i in (0, 1) for j in (0, 1) for k in (0, 1)) / 8.0


def make_nozawa_1995_poisson_system_3d(geometry: FeatureGeometry3D):
    """Build Q1 electrostatics with four independently floating poly lines."""
    if (not isinstance(geometry, FeatureGeometry3D)
            or geometry.material_levelsets is None
            or set(geometry.material_levelsets) != {
                NOZAWA_POLY_SI_MATERIAL_ID,
                NOZAWA_PHOTORESIST_MATERIAL_ID,
                NOZAWA_SIO2_MATERIAL_ID}):
        raise ValueError("geometry is not a Nozawa three-material replay cell")
    ids = np.asarray(sorted(geometry.material_levelsets), dtype=int)
    center_fields = np.stack([
        _cell_center_average(geometry.material_levelsets[int(item)]) for item in ids])
    owner = ids[np.argmax(center_fields, axis=0)]
    solid = np.max(center_fields, axis=0) >= 0.0
    epsilon_r = np.ones(owner.shape)
    epsilon_r[solid & (owner == NOZAWA_PHOTORESIST_MATERIAL_ID)] = 1.6
    epsilon_r[solid & (owner == NOZAWA_SIO2_MATERIAL_ID)] = 3.9
    # The ideal-conductor constraint makes the interior poly permittivity immaterial because
    # grad(V)=0 there; 11.7 is retained only as a transparent numerical material label.
    epsilon_r[solid & (owner == NOZAWA_POLY_SI_MATERIAL_ID)] = 11.7
    conductor_mask = (
        (geometry.material_id == NOZAWA_POLY_SI_MATERIAL_ID) & (geometry.phi >= 0.0))
    components, count = label(conductor_mask)
    if count != 4:
        raise RuntimeError("Nozawa Poisson builder requires four resolved poly lines")
    # Stable ids follow x position, not scipy's traversal details.
    centers = []
    for old_id in range(1, count + 1):
        centers.append((float(np.mean(np.where(components == old_id)[0])), old_id))
    conductor_ids = np.zeros_like(components)
    for new_id, (_center, old_id) in enumerate(sorted(centers), start=1):
        conductor_ids[components == old_id] = new_id
    fixed = np.zeros(geometry.phi.shape, dtype=bool)
    fixed[:, :, -1] = True
    return NodalPoissonSystem3D(
        epsilon_r, geometry.dx * geometry.mesh_length_unit_m, fixed,
        periodic_axes=(0, 1), floating_conductor_node_ids=conductor_ids)


@dataclass(frozen=True)
class Nozawa1995ReplaySetup3D:
    condition: Nozawa1995ReplayCondition3D
    geometry_contract: Nozawa1995GeometryContract3D
    process: PhysicalChargingProcess
    protocol: NotchingBenchmarkProtocol3D
    mode: str
    preflight_manifest: Mapping[str, object]

    def __post_init__(self):
        if (self.mode not in {"smoke", "charge_audit", "experiment"}
                or not isinstance(self.process, PhysicalChargingProcess)
                or not isinstance(self.protocol, NotchingBenchmarkProtocol3D)):
            raise ValueError("invalid Nozawa replay setup")
        object.__setattr__(
            self, "preflight_manifest", MappingProxyType(dict(self.preflight_manifest)))


def make_nozawa_1995_replay_setup(
        condition_id="fig10_l06s06_04", *, mode="smoke", dx_um=0.1,
        n_position=16, seed=1701, charging_timestep_s=1.0e-9,
        maximum_charging_steps=4000, terminal_window_s=2.0e-6,
        charging_timestep_policy="fixed", stochastic_gain_exponent=0.75,
        stochastic_gain_offset_steps=16,
        profile_steps=64, trajectory_emergency_max_steps=65536,
        observations_path=None, iedf_path=None, transport_device="cpu"):
    """Assemble, but do not run, one source-backed common-engine experiment."""
    if mode not in {"smoke", "charge_audit", "experiment"}:
        raise ValueError("mode must be 'smoke', 'charge_audit', or 'experiment'")
    if (int(n_position) != n_position or n_position < 4
            or int(seed) != seed or seed < 0
            or int(maximum_charging_steps) != maximum_charging_steps
            or maximum_charging_steps <= 0
            or int(profile_steps) != profile_steps or profile_steps <= 0
            or charging_timestep_policy not in {"fixed", "decreasing_gain"}
            or (mode != "smoke" and charging_timestep_policy == "fixed"
                and (not np.isfinite(terminal_window_s) or terminal_window_s <= 0.0))
            or not np.isfinite(stochastic_gain_exponent)
            or not 0.5 < stochastic_gain_exponent <= 1.0
            or int(stochastic_gain_offset_steps) != stochastic_gain_offset_steps
            or stochastic_gain_offset_steps <= 0
            or int(trajectory_emergency_max_steps) != trajectory_emergency_max_steps
            or trajectory_emergency_max_steps < 1024
            or not np.isfinite(charging_timestep_s) or charging_timestep_s <= 0.0):
        raise ValueError("invalid Nozawa replay numerical controls (n_position must be >= 4)")
    observation_path = (default_nozawa_data_path("digitized_notch_curves.csv")
                        if observations_path is None else Path(observations_path))
    ion_path = (default_nozawa_data_path("fig4a_ion_energy_distribution.csv")
                if iedf_path is None else Path(iedf_path))
    condition = load_nozawa_1995_replay_condition(
        condition_id, observations_path=observation_path)
    geometry_contract = make_nozawa_1995_open_area_geometry_3d(
        condition, dx_um=dx_um)
    geometry = geometry_contract.geometry
    boundary = build_hwang_giapis_1997_boundary_state(
        ion_path, reference_plane_m=geometry_contract.source_plane_um * 1e-6)
    observations = load_nozawa_1995_notch_observations(
        observation_path, expected_sha256=NOZAWA_1995_NOTCH_CURVES_SHA256)
    protocol = NotchingBenchmarkProtocol3D(
        observations, NOZAWA_1995_NOTCH_CURVES_SHA256,
        calibration_parameter_bounds={
            "overetch_fluence_scale": (0.5, 2.0),
            "sio2_scatter_critical_angle_deg": (30.0, 60.0),
        },
        stationarity_contract_revision=NOZAWA_STATIONARITY_CONTRACT,
        stationarity_contract_approved=True)
    charging_options = dict(
        patch_scales_m=(max(2.0 * float(dx_um), 0.2) * 1e-6, 0.8e-6),
        potential_rate_tolerance_v_s=1.0e3,
        timestep_s=float(charging_timestep_s),
        maximum_steps=(1 if mode == "smoke" else int(maximum_charging_steps)),
        current_balance_tolerance=0.08,
        timestep_policy=str(charging_timestep_policy),
        stochastic_gain_exponent=float(stochastic_gain_exponent),
        stochastic_gain_offset_steps=int(stochastic_gain_offset_steps),
        scramble_mode="fresh", compatible_q1_charge_state=True)
    if charging_timestep_policy != "decreasing_gain" and mode != "smoke":
        charging_options["terminal_window_s"] = float(terminal_window_s)
    profile_duration = (
        float(charging_timestep_s) if mode == "smoke"
        else (0.0 if mode == "charge_audit"
              else 2.0 * 0.3e-6 / NOZAWA_REPORTED_POLY_ETCH_RATE_M_S))
    n_steps = 1 if mode in {"smoke", "charge_audit"} else int(profile_steps)
    solver_options = dict(
        n_position=int(n_position), seed=int(seed), trajectory_fixed_dt=0.005,
        trajectory_max_steps=1024, trajectory_adaptive_horizon=True,
        trajectory_emergency_max_steps=int(trajectory_emergency_max_steps),
        periodic_lateral=True,
        neutral_forward_scatter=HwangGiapisSiO2ForwardScatter3D(
            NOZAWA_SIO2_MATERIAL_ID),
        neutral_forward_scatter_options={
            "periodic_lateral": True, "launch_offset": 1e-5,
            "maximum_periodic_wraps": 10000},
        reinitialize=(mode == "experiment"), transport_device=transport_device,
        profile_motion_enabled=(mode != "smoke"),
        bias_mode=("physical_time_resolved" if mode == "smoke" else "quasi_static"),
        experimental_claim=False)
    process = PhysicalChargingProcess(
        geometry=geometry, boundary=boundary,
        species_role={"Cl+": "energetic_bombardment", "electron": "charge_carrier"},
        mechanism=HwangGiapisClSiMechanism(),
        charging_system_builder=make_nozawa_1995_poisson_system_3d,
        etchable_material_ids=(NOZAWA_POLY_SI_MATERIAL_ID,),
        duration_s=profile_duration, n_steps=n_steps,
        source_bounds=(
            0.0, geometry_contract.realized_cell_width_um,
            0.0, (geometry.phi.shape[1] - 1) * geometry.dx),
        source_z=geometry_contract.source_plane_um,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        charging_options=charging_options, solver_options=solver_options)
    preflight = {
        "schema": NOZAWA_REPLAY_SCHEMA,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_status": (
            "operational hard-visibility smoke; not a saturation or validation claim"
            if mode == "smoke" else
            (("bounded decreasing-gain charge warm start; zero physical-time credit, "
              "no saturation or profile claim")
             if mode == "charge_audit"
             and charging_timestep_policy == "decreasing_gain" else
             "fixed-geometry signed-R2 charging audit; no profile claim"
             if mode == "charge_audit" else
             "development experimental replay; validation requires all C4 uncertainty gates")),
        "mode": mode,
        "condition": _json_value(condition.__dict__),
        "geometry": _json_value(geometry_contract.manifest),
        "evidence": {
            "nozawa_observations_path_name": observation_path.name,
            "nozawa_observations_sha256": _sha256_path(observation_path),
            "hwang_iedf_path_name": ion_path.name,
            "hwang_iedf_sha256": _sha256_path(ion_path),
        },
        "protocol_commit_sha256": protocol.commit_sha256,
        "calibration_split_frozen_before_run": True,
        "numerics": {
            "dx_um": float(dx_um), "n_position": int(n_position), "seed": int(seed),
            "charging_timestep_s": float(charging_timestep_s),
            "charging_timestep_policy": str(charging_timestep_policy),
            "maximum_charging_steps": int(maximum_charging_steps),
            "terminal_window_s": (
                None if mode == "smoke" or charging_timestep_policy == "decreasing_gain"
                else float(terminal_window_s)),
            "stochastic_gain_exponent": (
                float(stochastic_gain_exponent)
                if charging_timestep_policy == "decreasing_gain" else None),
            "stochastic_gain_offset_steps": (
                int(stochastic_gain_offset_steps)
                if charging_timestep_policy == "decreasing_gain" else None),
            "profile_steps": n_steps, "transport_device": transport_device,
            "profile_motion_enabled": bool(mode != "smoke"),
            "trajectory_emergency_max_steps": int(trajectory_emergency_max_steps),
        },
        "operators": {
            "visibility": "exact hard visibility with certified float64 replay",
            "charge": "fresh-scramble physical surface-current ODE",
            "electrostatics": (
                "compatible mixed-surface Q1 charge; four floating equipotentials with "
                "separately conserved conductor totals"),
            "poly_si_removal": "Hwang--Giapis Eq. (4.1)",
            "sio2_scatter": "neutralized Hwang--Giapis Eqs. (4.2)--(4.3)",
            "see": "off, as in the source model",
            "surface_conduction": "poly equipotential only; dielectric leakage off",
        },
        "claim_blockers": [
            "published Nozawa measurement uncertainty is not quantified",
            "the yield prefactor fixes computational etch time and is nonpredictive",
            "a result must pass charging, grid, timestep, sample, and held-out gates",
        ],
    }
    return Nozawa1995ReplaySetup3D(
        condition, geometry_contract, process, protocol, mode, preflight)


def measure_nozawa_1995_edge_notch_3d(setup, geometry):
    """Measure the Hwang--Giapis target sidewall, not generic lateral erosion.

    In Fig. 3 of Hwang & Giapis (1997), X and Y are the final ordinary line/space pair and the
    varied open area lies between Y and Z.  Charging of that open area creates a notch on the
    *inner* (X-facing) foot of line Y.  In the full periodic cell used here, that is the right
    boundary of the X--Y opening.  Depth is referenced to the 0.4 um poly-Si top; the bottom
    0.1 um band is the notch band immediately above the 0.1 um SiO2 interface.
    """
    if not isinstance(setup, Nozawa1995ReplaySetup3D):
        raise TypeError("setup must be a Nozawa1995ReplaySetup3D")
    centers = setup.geometry_contract.line_centers_um
    left_center, middle_center = centers[:2]
    observables = measure_trench_profile_observables_3d(
        geometry,
        lateral_bounds_m=(left_center * 1e-6, middle_center * 1e-6, 0.0, 0.2e-6),
        opening_center_x_m=0.5 * (left_center + middle_center) * 1e-6,
        feature_top_z_m=0.4e-6,
        reference_depth_interval_m=(0.0, 0.1e-6),
        bow_depth_interval_m=(0.1e-6, 0.2e-6),
        notch_depth_interval_m=(0.2e-6, 0.3e-6),
        minimum_longitudinal_rows=3)
    return observables


def _surface_displacement_summary(initial, final):
    """Measure geometric motion at the zero level set, not off-interface phi drift."""
    spacing = (float(initial.dx),) * 3
    initial_vertices = marching_cubes(initial.phi, 0.0, spacing=spacing)[0]
    final_vertices = marching_cubes(final.phi, 0.0, spacing=spacing)[0]
    initial_to_final = cKDTree(final_vertices).query(initial_vertices)[0]
    final_to_initial = cKDTree(initial_vertices).query(final_vertices)[0]
    return {
        "profile_surface_mean_displacement_um": float(
            0.5 * (np.mean(initial_to_final) + np.mean(final_to_initial))),
        "profile_surface_rms_displacement_um": float(np.sqrt(0.5 * (
            np.mean(initial_to_final ** 2) + np.mean(final_to_initial ** 2)))),
        "profile_surface_max_displacement_um": float(max(
            np.max(initial_to_final), np.max(final_to_initial))),
        "profile_grid_spacing_um": float(initial.dx),
    }


def _plot_replay(setup, result, path, displacement):
    # Keep one writable cache across replay directories.  A per-run cache made
    # every successful CLI invocation rebuild Matplotlib's font index, which is
    # both slow and noisy in unattended jobs and read-only container installs.
    cache = Path(tempfile.gettempdir()) / "petch-matplotlib-cache"
    cache.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch

    initial = setup.geometry_contract.geometry
    final = result.geometry
    step = result.steps[-1]
    potential = step.charging.potential_v
    notch = measure_nozawa_1995_edge_notch_3d(setup, final)
    predicted_notch_um = 1e6 * notch.maximum_right_notch_depth_m
    j = initial.phi.shape[1] // 2
    x = np.arange(initial.phi.shape[0]) * initial.dx
    z = np.arange(initial.phi.shape[2]) * initial.dx
    figure, axes = plt.subplots(1, 3, figsize=(15.5, 4.5), constrained_layout=True)
    material = np.ma.masked_where(
        initial.material_id[:, j, :].T == 0, initial.material_id[:, j, :].T)
    material_cmap = ListedColormap(("#4c78a8", "#f2cf5b", "#9d755d"))
    material_cmap.set_bad(alpha=0.0)
    axes[0].pcolormesh(
        x, z, material, shading="nearest", cmap=material_cmap,
        norm=BoundaryNorm((0.5, 1.5, 2.5, 3.5), material_cmap.N))
    axes[0].contour(
        x, z, initial.phi[:, j, :].T, levels=[0.0], colors="black", linewidths=.8)
    axes[0].set(
        title="Starting material stack (before etch)", xlabel="x (um)", ylabel="z (um)",
        ylim=(-0.02, 1.65))
    axes[0].legend(
        handles=(Patch(color="#4c78a8", label="poly-Si"),
                 Patch(color="#f2cf5b", label="photoresist"),
                 Patch(color="#9d755d", label="SiO2")),
        loc="upper right", fontsize=8)

    axes[1].contour(
        x, z, initial.phi[:, j, :].T, levels=[0.0], colors="#777777",
        linewidths=3.0, linestyles="dashed")
    axes[1].contour(
        x, z, final.phi[:, j, :].T, levels=[0.0], colors="#e45756",
        linewidths=1.5)
    shift_nm = 1.0e3 * displacement["profile_surface_max_displacement_um"]
    shift_cells = (displacement["profile_surface_max_displacement_um"]
                   / displacement["profile_grid_spacing_um"])
    profile_title = (
        (f"Smoke displacement: {shift_nm:.2f} nm ({shift_cells:.3f} cell)"
         if setup.mode == "smoke" else
         "Fixed geometry: charging only, no profile motion")
        if setup.mode != "experiment" else
        f"Target Y-foot notch: {predicted_notch_um:.3f} um "
        f"(experiment {setup.condition.measured_notch_depth_um:.3f} um)")
    axes[1].set(title=profile_title, xlabel="x (um)", ylabel="z (um)")
    if setup.mode != "experiment":
        axes[1].set_ylim(-0.02, 1.65)
    else:
        target_x = setup.geometry_contract.line_centers_um[1] - (
            0.5 * setup.condition.line_width_um)
        axes[1].set(
            xlim=(target_x - 0.45, target_x + 0.35),
            ylim=(0.02, 0.48))
    axes[1].plot([], [], color="#777777", linewidth=3, linestyle="--", label="initial")
    axes[1].plot(
        [], [], color="#e45756", linewidth=1.5,
        label=("after smoke" if setup.mode == "smoke" else "final profile"))
    axes[1].legend(loc="upper right", fontsize=8)

    image = axes[2].pcolormesh(
        x, z, potential[:, j, :].T, shading="auto", cmap="coolwarm")
    axes[2].contour(
        x, z, final.phi[:, j, :].T, levels=[0.0], colors="black", linewidths=.7)
    axes[2].set(
        title=(
            f"Potential after {step.charging.physical_time_s * 1e6:g} us "
            f"({'R2 converged' if step.charging.converged else 'not equilibrium'})"),
        xlabel="x (um)", ylabel="z (um)")
    figure.colorbar(image, ax=axes[2], label="V")
    figure.suptitle(
        f"{setup.condition.condition_id} — "
        + (("hard-visibility operator smoke; experimental notch not scored"
            if setup.mode == "smoke"
            else "fixed-geometry charging audit; experimental notch not scored")
           if setup.mode != "experiment"
           else "Nozawa/Hwang open-area replay; exact target sidewall scored"),
        fontsize=11)
    figure.savefig(path, dpi=170)
    plt.close(figure)


def run_nozawa_1995_replay(setup: Nozawa1995ReplaySetup3D, output_directory):
    """Run one setup and always leave a machine-readable success or refusal artifact."""
    if not isinstance(setup, Nozawa1995ReplaySetup3D):
        raise TypeError("setup must be Nozawa1995ReplaySetup3D")
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=False)
    preflight_hash = _write_json(output / "preflight_manifest.json", setup.preflight_manifest)
    progress_checkpoint = output / "charging_progress_checkpoint.npz"
    heartbeat_path = output / "heartbeat.json"
    last_checkpointed_step = None
    surface_mesh_fingerprint = _charging_surface_mesh_fingerprint_3d(
        setup.geometry_contract.geometry)

    def persist_progress(*, sigma_c_per_m2, charge_node_c, potential_v,
                         history_item, accepted_steps, rejected_steps,
                         physical_time_s, pseudo_time_s, resume_sampling_epoch):
        nonlocal last_checkpointed_step
        item = dict(history_item)
        checkpoint_due = bool(
            accepted_steps != last_checkpointed_step
            and (accepted_steps % 25 == 0 or item["saturation_gates_satisfied"]))
        if checkpoint_due:
            _atomic_npz(
                progress_checkpoint,
                sigma_c_per_m2=sigma_c_per_m2,
                charge_node_c=charge_node_c,
                potential_v=potential_v,
                accepted_steps=np.asarray(accepted_steps),
                rejected_steps=np.asarray(rejected_steps),
                physical_time_s=np.asarray(physical_time_s),
                pseudo_time_s=np.asarray(pseudo_time_s),
                resume_sampling_epoch=np.asarray(resume_sampling_epoch),
                resume_stochastic_gain_age_steps=np.asarray(
                    item["stochastic_gain_age_steps"]),
                condition_id=np.asarray(setup.condition.condition_id),
                surface_mesh_fingerprint=np.asarray(surface_mesh_fingerprint))
            last_checkpointed_step = accepted_steps
        if checkpoint_due or item["evaluation"] % 10 == 0:
            _write_json(heartbeat_path, _failure_diagnostic_value({
                "schema": "petch-nozawa-1995-heartbeat-v1",
                "status": "running",
                "updated_utc": datetime.now(timezone.utc).isoformat(),
                "condition_id": setup.condition.condition_id,
                "accepted_steps": accepted_steps,
                "rejected_steps": rejected_steps,
                "physical_time_s": physical_time_s,
                "pseudo_time_s": pseudo_time_s,
                "resume_sampling_epoch": resume_sampling_epoch,
                "resume_stochastic_gain_age_steps": item[
                    "stochastic_gain_age_steps"],
                "potential_rate_max_v_s": item["gate_potential_rate_max_v_s"],
                "patch_b2_rms": item["gate_patch_rms_relative_imbalance"],
                "patch_b2_max": item["gate_patch_max_relative_imbalance"],
                "node_rms": item["rms_relative_current_imbalance_node"],
                "node_worst": item["max_relative_current_imbalance_node"],
                "terminal_window_ready": item["terminal_window_ready"],
                "b1_satisfied": item["b1_potential_saturation_satisfied"],
                "b2_satisfied": item["b2_patch_balance_satisfied"],
                "checkpoint": (
                    progress_checkpoint.name if progress_checkpoint.is_file() else None),
            }))

    charging_options = dict(setup.process.charging_options)
    if "progress_callback" in charging_options:
        raise ValueError("Nozawa replay reserves progress_callback for its restart artifacts")
    charging_options["progress_callback"] = persist_progress
    process = replace(setup.process, charging_options=charging_options)
    try:
        result = process.run()
    except SurfaceChargingSaturationError as error:
        _atomic_npz(
            output / "charging_failure_checkpoint.npz",
            sigma_c_per_m2=error.sigma_c_per_m2,
            accepted_steps=np.asarray(error.accepted_steps),
            rejected_steps=np.asarray(error.rejected_steps),
            physical_time_s=np.asarray(error.physical_time_s),
            pseudo_time_s=np.asarray(error.pseudo_time_s),
            resume_sampling_epoch=np.asarray(error.resume_sampling_epoch),
            resume_stochastic_gain_age_steps=np.asarray(
                error.resume_stochastic_gain_age_steps),
            condition_id=np.asarray(setup.condition.condition_id),
            surface_mesh_fingerprint=np.asarray(surface_mesh_fingerprint))
        history_hash = _write_json(
            output / "charging_history.json",
            [_failure_diagnostic_value(dict(item)) for item in error.history])
        summary = {
            "schema": NOZAWA_REPLAY_SCHEMA,
            "status": "scientific_refusal",
            "reason": str(error),
            "accepted_steps": error.accepted_steps,
            "rejected_steps": error.rejected_steps,
            "physical_time_s": error.physical_time_s,
            "pseudo_time_s": error.pseudo_time_s,
            "resume_sampling_epoch": error.resume_sampling_epoch,
            "resume_stochastic_gain_age_steps": (
                error.resume_stochastic_gain_age_steps),
            "last_diagnostics": (
                _failure_diagnostic_value(dict(error.history[-1]))
                if error.history else None),
            "preflight_manifest_sha256": preflight_hash,
            "charging_history_sha256": history_hash,
            "restart_note": (
                "the face-charge checkpoint is preserved; the engine refused profile motion "
                "because the declared stationary-field gate did not pass"),
        }
        _write_json(output / "summary.json", summary)
        _write_json(heartbeat_path, {
            "schema": "petch-nozawa-1995-heartbeat-v1",
            "status": "scientific_refusal",
            "updated_utc": datetime.now(timezone.utc).isoformat(),
            "condition_id": setup.condition.condition_id,
            "accepted_steps": error.accepted_steps,
            "physical_time_s": error.physical_time_s,
            "pseudo_time_s": error.pseudo_time_s,
            "resume_sampling_epoch": error.resume_sampling_epoch,
            "resume_stochastic_gain_age_steps": (
                error.resume_stochastic_gain_age_steps),
            "reason": str(error),
        })
        return None, summary
    step = result.steps[-1]
    # ``final_step`` evaluates current at the accepted endpoint and also contains one proposed
    # Euler update.  The reported state/potential therefore pair with ``poisson_before``.
    poisson = step.charging.final_step.poisson_before
    np.savez_compressed(
        output / "final_state.npz", phi=result.geometry.phi,
        material_id=result.geometry.material_id,
        sigma_c_per_m2=result.surface_charge_c_per_m2,
        potential_v=step.charging.potential_v,
        resume_sampling_epoch=np.asarray(
            step.charging.diagnostics["resume_sampling_epoch"]),
        resume_stochastic_gain_age_steps=np.asarray(
            step.charging.diagnostics["resume_stochastic_gain_age_steps"]),
        condition_id=np.asarray(setup.condition.condition_id),
        surface_mesh_fingerprint=np.asarray(surface_mesh_fingerprint))
    engine_manifest_hash = _write_json(
        output / "engine_run_manifest.json", result.run_manifest)
    history_hash = _write_json(
        output / "charging_history.json",
        [_failure_diagnostic_value(dict(item)) for item in step.charging.history])
    displacement = _surface_displacement_summary(
        setup.geometry_contract.geometry, result.geometry)
    notch = measure_nozawa_1995_edge_notch_3d(setup, result.geometry)
    predicted_notch_um = 1e6 * notch.maximum_right_notch_depth_m
    experimental_error_um = predicted_notch_um - setup.condition.measured_notch_depth_um
    summary = {
        "schema": NOZAWA_REPLAY_SCHEMA,
        "status": "pass",
        "scientific_status": setup.preflight_manifest["scientific_status"],
        "condition_id": setup.condition.condition_id,
        "measured_notch_depth_um": setup.condition.measured_notch_depth_um,
        "charging_converged_signed_r2": step.charging.converged,
        "charging_accepted_steps": step.charging.accepted_steps,
        "charging_physical_time_s": step.charging.physical_time_s,
        "charging_pseudo_time_s": step.charging.pseudo_time_s,
        "charging_timestep_policy": step.charging.timestep_policy,
        "resume_sampling_epoch": step.charging.diagnostics[
            "resume_sampling_epoch"],
        "resume_stochastic_gain_age_steps": step.charging.diagnostics[
            "resume_stochastic_gain_age_steps"],
        "potential_min_v": float(np.min(step.charging.potential_v)),
        "potential_max_v": float(np.max(step.charging.potential_v)),
        "floating_conductor_ids": poisson.floating_conductor_ids,
        "floating_conductor_voltage_v": poisson.floating_conductor_voltage_v,
        "maximum_floating_conductor_voltage_spread_v": (
            poisson.maximum_floating_conductor_voltage_spread_v),
        "charge_conservation_relative_error": step.charging.history[-1][
            "charge_conservation_relative_error"],
        "neutral_forward_scatter_rate_s": step.diagnostics[
            "neutral_forward_scatter_rate_s"],
        "neutral_forward_scatter_particle_balance_error": step.diagnostics[
            "neutral_forward_scatter_particle_balance_error"],
        "neutral_forward_scatter_energy_balance_error": step.diagnostics[
            "neutral_forward_scatter_energy_balance_error"],
        "predicted_target_notch_depth_um": predicted_notch_um,
        "predicted_control_side_notch_depth_um": (
            1e6 * notch.maximum_left_notch_depth_m),
        "notch_depth_error_um": experimental_error_um,
        "absolute_notch_depth_error_um": abs(experimental_error_um),
        "digitization_uncertainty_um": setup.condition.digitization_uncertainty_um,
        "preflight_manifest_sha256": preflight_hash,
        "engine_run_manifest_sha256": engine_manifest_hash,
        "charging_history_sha256": history_hash,
        "validated_notch_prediction": False,
        **displacement,
    }
    _write_json(output / "summary.json", summary)
    _plot_replay(setup, result, output / "overview.png", displacement)
    _write_json(heartbeat_path, {
        "schema": "petch-nozawa-1995-heartbeat-v1",
        "status": "complete",
        "updated_utc": datetime.now(timezone.utc).isoformat(),
        "condition_id": setup.condition.condition_id,
        "charging_converged_signed_r2": step.charging.converged,
        "predicted_target_notch_depth_um": predicted_notch_um,
    })
    return result, summary


def _default_output(condition_id, mode):
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return Path("results/nozawa_1995_user_replay") / f"{condition_id}-{mode}-{stamp}"


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Run the checksum-verified Nozawa/Hwang notching replay")
    parser.add_argument("--condition", default="fig10_l06s06_04")
    parser.add_argument(
        "--mode", choices=("preflight", "smoke", "charge_audit", "experiment"),
                        default="smoke")
    parser.add_argument("--dx-um", type=float, default=0.1)
    parser.add_argument("--n-position", type=int, default=16)
    parser.add_argument("--seed", type=int, default=1701)
    parser.add_argument("--charging-timestep-s", type=float, default=1e-9)
    parser.add_argument("--maximum-charging-steps", type=int, default=4000)
    parser.add_argument("--terminal-window-s", type=float, default=2e-6)
    parser.add_argument(
        "--charging-timestep-policy",
        choices=("fixed", "decreasing_gain"), default="fixed")
    parser.add_argument("--stochastic-gain-exponent", type=float, default=0.75)
    parser.add_argument("--stochastic-gain-offset-steps", type=int, default=16)
    parser.add_argument("--profile-steps", type=int, default=64)
    parser.add_argument("--trajectory-emergency-max-steps", type=int, default=65536)
    parser.add_argument("--transport-device", default="cpu")
    parser.add_argument("--restart-checkpoint")
    parser.add_argument("--restart-sampling-epoch", type=int)
    parser.add_argument("--restart-stochastic-gain-age-steps", type=int)
    parser.add_argument("--observations")
    parser.add_argument("--iedf")
    parser.add_argument("--output")
    args = parser.parse_args(argv)
    setup_mode = "smoke" if args.mode == "preflight" else args.mode
    setup = make_nozawa_1995_replay_setup(
        args.condition, mode=setup_mode, dx_um=args.dx_um,
        n_position=args.n_position, seed=args.seed,
        charging_timestep_s=args.charging_timestep_s,
        maximum_charging_steps=args.maximum_charging_steps,
        terminal_window_s=args.terminal_window_s,
        charging_timestep_policy=args.charging_timestep_policy,
        stochastic_gain_exponent=args.stochastic_gain_exponent,
        stochastic_gain_offset_steps=args.stochastic_gain_offset_steps,
        profile_steps=args.profile_steps,
        trajectory_emergency_max_steps=args.trajectory_emergency_max_steps,
        observations_path=args.observations,
        iedf_path=args.iedf, transport_device=args.transport_device)
    if args.restart_checkpoint is not None:
        if args.restart_sampling_epoch is None or args.restart_sampling_epoch < 0:
            parser.error(
                "--restart-checkpoint requires an explicit nonnegative "
                "--restart-sampling-epoch")
        checkpoint_path = Path(args.restart_checkpoint).resolve()
        expected_surface_mesh_fingerprint = _charging_surface_mesh_fingerprint_3d(
            setup.geometry_contract.geometry)
        expected_surface_face_count = len(extract_mesh_3d(
            setup.geometry_contract.geometry.phi,
            setup.geometry_contract.geometry.dx)[1])
        expected_geometry_manifest = _json_value(
            setup.geometry_contract.manifest)
        with np.load(checkpoint_path, allow_pickle=False) as checkpoint:
            if "sigma_c_per_m2" not in checkpoint.files:
                parser.error("restart checkpoint has no sigma_c_per_m2")
            sigma = np.asarray(checkpoint["sigma_c_per_m2"], dtype=float).copy()
            if sigma.ndim != 1 or np.any(~np.isfinite(sigma)):
                parser.error("restart surface charge must be one finite face vector")
            if sigma.shape != (expected_surface_face_count,):
                parser.error(
                    "restart surface charge face count disagrees with the "
                    "reconstructed condition geometry")
            restart_metadata = {
                name: _json_value(checkpoint[name])
                for name in (
                    "accepted_steps", "physical_time_s", "pseudo_time_s",
                    "resume_sampling_epoch", "resume_stochastic_gain_age_steps",
                    "condition_id", "surface_mesh_fingerprint")
                if name in checkpoint.files
            }
        recorded_condition_id = restart_metadata.get("condition_id")
        if (recorded_condition_id is not None
                and str(recorded_condition_id) != setup.condition.condition_id):
            parser.error(
                "restart checkpoint condition_id disagrees with --condition")
        recorded_mesh_fingerprint = restart_metadata.get(
            "surface_mesh_fingerprint")
        if (recorded_mesh_fingerprint is not None
                and str(recorded_mesh_fingerprint)
                != expected_surface_mesh_fingerprint):
            parser.error(
                "restart checkpoint surface mesh fingerprint disagrees with "
                "the reconstructed condition")
        source_preflight_path = checkpoint_path.parent / "preflight_manifest.json"
        source_preflight_hash = None
        source_preflight_condition_id = None
        source_preflight_geometry = None
        if source_preflight_path.is_file():
            try:
                source_preflight = json.loads(source_preflight_path.read_text())
                source_preflight_condition_id = source_preflight.get(
                    "condition", {}).get("condition_id")
                source_preflight_geometry = source_preflight.get("geometry")
            except (OSError, ValueError, TypeError) as error:
                parser.error(f"could not verify restart source preflight: {error}")
            if (source_preflight_condition_id is not None
                    and source_preflight_condition_id != setup.condition.condition_id):
                parser.error(
                    "restart source preflight condition disagrees with --condition")
            if (source_preflight_geometry is not None
                    and source_preflight_geometry != expected_geometry_manifest):
                parser.error(
                    "restart source geometry disagrees with the reconstructed "
                    "condition geometry")
            source_preflight_hash = _sha256_path(source_preflight_path)
        charging_options = dict(setup.process.charging_options)
        charging_options["initial_sampling_epoch"] = int(
            args.restart_sampling_epoch)
        if args.charging_timestep_policy == "decreasing_gain":
            recorded_gain_age = restart_metadata.get(
                "resume_stochastic_gain_age_steps")
            declared_gain_age = args.restart_stochastic_gain_age_steps
            if declared_gain_age is None:
                if recorded_gain_age is None:
                    parser.error(
                        "a decreasing-gain restart requires "
                        "--restart-stochastic-gain-age-steps when the checkpoint "
                        "does not record one")
                declared_gain_age = int(recorded_gain_age)
            if declared_gain_age < 0:
                parser.error(
                    "--restart-stochastic-gain-age-steps must be nonnegative")
            if (recorded_gain_age is not None
                    and int(recorded_gain_age) != int(declared_gain_age)):
                parser.error(
                    "declared stochastic gain age disagrees with the checkpoint")
            charging_options["initial_stochastic_gain_age_steps"] = int(
                declared_gain_age)
        elif args.restart_stochastic_gain_age_steps is not None:
            parser.error(
                "--restart-stochastic-gain-age-steps requires "
                "--charging-timestep-policy decreasing_gain")
        solver_options = dict(setup.process.solver_options)
        solver_options["initial_sigma_c_per_m2"] = sigma
        process = replace(
            setup.process, charging_options=charging_options,
            solver_options=solver_options)
        preflight = dict(setup.preflight_manifest)
        preflight["restart"] = {
            "checkpoint_name": checkpoint_path.name,
            "checkpoint_sha256": _sha256_path(checkpoint_path),
            "declared_first_unused_sampling_epoch": int(
                args.restart_sampling_epoch),
            "declared_stochastic_gain_age_steps": (
                charging_options.get("initial_stochastic_gain_age_steps")),
            "recorded_checkpoint_metadata": restart_metadata,
            "reconstructed_surface_mesh_fingerprint": (
                expected_surface_mesh_fingerprint),
            "source_preflight_manifest_sha256": source_preflight_hash,
            "source_preflight_condition_id": source_preflight_condition_id,
            "source_preflight_geometry_verified": bool(
                source_preflight_geometry is not None),
            "initial_state_policy": (
                "project into the compatible mixed conductor/dielectric state while "
                "preserving the exact Poisson field and each floating-conductor total"),
        }
        setup = replace(setup, process=process, preflight_manifest=preflight)
    elif args.restart_sampling_epoch is not None:
        parser.error("--restart-sampling-epoch requires --restart-checkpoint")
    elif args.restart_stochastic_gain_age_steps is not None:
        parser.error(
            "--restart-stochastic-gain-age-steps requires --restart-checkpoint")
    output = Path(args.output) if args.output else _default_output(args.condition, args.mode)
    if args.mode == "preflight":
        output.mkdir(parents=True, exist_ok=False)
        digest = _write_json(output / "preflight_manifest.json", setup.preflight_manifest)
        _write_json(output / "summary.json", {
            "schema": NOZAWA_REPLAY_SCHEMA, "status": "preflight_pass",
            "preflight_manifest_sha256": digest,
            "next_command": "repeat with --mode smoke or --mode experiment",
        })
        print(output)
        return 0
    _result, summary = run_nozawa_1995_replay(setup, output)
    print(json.dumps(_json_value(summary), indent=2, sort_keys=True))
    print(f"artifacts: {output}")
    return 0 if summary["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
