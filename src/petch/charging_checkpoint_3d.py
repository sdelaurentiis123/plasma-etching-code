"""Safe, versioned step-boundary checkpoints for charged profile co-evolution."""
from __future__ import annotations

from dataclasses import dataclass, field
from hashlib import sha256
import json
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

import numpy as np

from .charging_coevolution_3d import ChargingCoevolution3DResult
from .feature_step_3d import FeatureGeometry3D
from .fluorocarbon_lamagna import LaMagnaFluorocarbonState
from .material_mechanism_3d import MaterialSurfaceState3D
from .physical_sputtering import PhysicalSputterState
from .surface_kinetics import SiO2SurfaceState
from .tabulated_chemistry import TabulatedSiSurfaceState


CHARGING_CHECKPOINT_SCHEMA = "petch-charging-checkpoint-3d-v1"

_STATE_TYPES = {
    "none": None,
    "petch.PhysicalSputterState": PhysicalSputterState,
    "petch.SiO2SurfaceState": SiO2SurfaceState,
    "petch.TabulatedSiSurfaceState": TabulatedSiSurfaceState,
    "petch.MaterialSurfaceState3D": MaterialSurfaceState3D,
    "petch.LaMagnaFluorocarbonState": LaMagnaFluorocarbonState,
}
_STATE_IDS = {value: key for key, value in _STATE_TYPES.items() if value is not None}


def _manifest_sha256(manifest):
    encoded = json.dumps(
        dict(manifest), sort_keys=True, separators=(",", ":"),
        ensure_ascii=True).encode("utf-8")
    return sha256(encoded).hexdigest()


@dataclass(frozen=True)
class PhysicalChargingCheckpoint3D:
    """A profile-step boundary with geometry, charge, and matching surface state.

    The disk representation is NPZ plus an embedded JSON metadata scalar, loaded with
    ``allow_pickle=False``. Boundary conditions, mechanisms, and solver settings deliberately stay
    in the caller's versioned process config; their prior manifest is bound by SHA-256 so a resume
    can report whether the operator/config changed.
    """

    geometry: FeatureGeometry3D
    sigma_c_per_m2: np.ndarray
    surface_state_type: str
    surface_state_fields: Mapping[str, np.ndarray]
    surface_state_mesh_fingerprint: str
    completed_duration_s: float
    completed_steps: int
    source_manifest_sha256: str
    surface_state_upper_bounds: Mapping[str, float | None] = field(default_factory=dict)
    surface_state_remap_modes: Mapping[str, str] = field(default_factory=dict)
    schema: str = CHARGING_CHECKPOINT_SCHEMA

    def __post_init__(self):
        sigma = np.asarray(self.sigma_c_per_m2, dtype=float).copy()
        fields = {}
        for name, supplied in dict(self.surface_state_fields).items():
            if (not isinstance(name, str) or not name
                    or any(character not in "abcdefghijklmnopqrstuvwxyz0123456789_"
                           for character in name)):
                raise ValueError("checkpoint surface-state field names must be safe identifiers")
            value = np.asarray(supplied, dtype=float).copy()
            if np.any(~np.isfinite(value)):
                raise ValueError("checkpoint surface-state fields must be finite")
            value.setflags(write=False)
            fields[name] = value
        upper = (dict(self.surface_state_upper_bounds)
                 if self.surface_state_upper_bounds
                 else {name: None for name in fields})
        modes = (dict(self.surface_state_remap_modes)
                 if self.surface_state_remap_modes
                 else {name: "conservative" for name in fields})
        if (not isinstance(self.geometry, FeatureGeometry3D)
                or sigma.ndim != 1 or np.any(~np.isfinite(sigma))
                or self.surface_state_type not in _STATE_TYPES
                or (self.surface_state_type == "none") != (not fields)
                or set(upper) != set(fields)
                or any(value is not None and (
                    not np.isfinite(value) or value <= 0.0) for value in upper.values())
                or set(modes) != set(fields)
                or any(value not in {"conservative", "intensive"}
                       for value in modes.values())
                or not isinstance(self.surface_state_mesh_fingerprint, str)
                or (self.surface_state_type != "none"
                    and not self.surface_state_mesh_fingerprint)
                or not np.isfinite(self.completed_duration_s)
                or self.completed_duration_s < 0.0
                or int(self.completed_steps) != self.completed_steps
                or self.completed_steps < 0
                or len(self.source_manifest_sha256) != 64
                or any(character not in "0123456789abcdef"
                       for character in self.source_manifest_sha256)
                or self.schema != CHARGING_CHECKPOINT_SCHEMA):
            raise ValueError("invalid physical charging checkpoint")
        sigma.setflags(write=False)
        object.__setattr__(self, "sigma_c_per_m2", sigma)
        object.__setattr__(self, "surface_state_fields", MappingProxyType(fields))
        object.__setattr__(self, "surface_state_upper_bounds", MappingProxyType(upper))
        object.__setattr__(self, "surface_state_remap_modes", MappingProxyType(modes))
        object.__setattr__(self, "completed_duration_s", float(self.completed_duration_s))
        object.__setattr__(self, "completed_steps", int(self.completed_steps))

    @classmethod
    def from_result(cls, result: ChargingCoevolution3DResult):
        if not isinstance(result, ChargingCoevolution3DResult):
            raise TypeError("result must be ChargingCoevolution3DResult")
        state = result.surface_state
        if state is None:
            state_type = "none"
            fields = {}
        else:
            state_type = _STATE_IDS.get(type(state))
            if state_type is None or not hasattr(state, "conservative_surface_fields"):
                raise TypeError(
                    f"surface state {type(state).__name__} has no registered safe checkpoint codec")
            fields = state.conservative_surface_fields()
        upper = ({} if state is None else state.conservative_surface_upper_bounds())
        modes = ({} if state is None else (
            state.surface_field_remap_modes()
            if hasattr(state, "surface_field_remap_modes")
            else {name: "conservative" for name in fields}))
        return cls(
            geometry=result.geometry,
            sigma_c_per_m2=result.sigma_c_per_m2,
            surface_state_type=state_type,
            surface_state_fields=fields,
            surface_state_mesh_fingerprint=result.surface_state_mesh_fingerprint,
            completed_duration_s=result.duration_s,
            completed_steps=len(result.steps),
            source_manifest_sha256=_manifest_sha256(result.run_manifest),
            surface_state_upper_bounds=upper, surface_state_remap_modes=modes)

    def restore_surface_state(self):
        state_class = _STATE_TYPES[self.surface_state_type]
        if state_class is None:
            return None
        if state_class is MaterialSurfaceState3D:
            return state_class(
                self.surface_state_fields, self.surface_state_upper_bounds,
                self.surface_state_remap_modes)
        bare = state_class.bare()
        return bare.with_conservative_surface_fields(self.surface_state_fields)

    def save(self, path):
        target = Path(path)
        if target.suffix != ".npz":
            raise ValueError("charging checkpoint path must end in .npz")
        material_layers = self.geometry.material_levelsets
        metadata = dict(
            schema=self.schema,
            dx=self.geometry.dx,
            mesh_length_unit_m=self.geometry.mesh_length_unit_m,
            mesh_origin_m=self.geometry.mesh_origin_m,
            material_levelset_ids=([] if material_layers is None else sorted(material_layers)),
            surface_state_type=self.surface_state_type,
            surface_state_field_names=sorted(self.surface_state_fields),
            surface_state_upper_bounds=dict(self.surface_state_upper_bounds),
            surface_state_remap_modes=dict(self.surface_state_remap_modes),
            surface_state_mesh_fingerprint=self.surface_state_mesh_fingerprint,
            completed_duration_s=self.completed_duration_s,
            completed_steps=self.completed_steps,
            source_manifest_sha256=self.source_manifest_sha256)
        arrays = dict(
            metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)),
            geometry_phi=self.geometry.phi,
            geometry_material_id=self.geometry.material_id,
            sigma_c_per_m2=self.sigma_c_per_m2)
        if material_layers is not None:
            arrays.update({
                f"material_levelset_{material_id}": material_layers[material_id]
                for material_id in sorted(material_layers)})
        arrays.update({
            f"surface_state_{name}": value
            for name, value in self.surface_state_fields.items()})
        np.savez_compressed(target, **arrays)

    @classmethod
    def load(cls, path):
        source = Path(path)
        with np.load(source, allow_pickle=False) as archive:
            metadata = json.loads(str(archive["metadata_json"]))
            if metadata.get("schema") != CHARGING_CHECKPOINT_SCHEMA:
                raise ValueError("unsupported charging checkpoint schema")
            layer_ids = tuple(int(value) for value in metadata["material_levelset_ids"])
            fields = {
                name: archive[f"surface_state_{name}"].copy()
                for name in metadata["surface_state_field_names"]}
            expected_keys = {
                "metadata_json", "geometry_phi", "geometry_material_id", "sigma_c_per_m2",
                *(f"material_levelset_{value}" for value in layer_ids),
                *(f"surface_state_{name}" for name in fields),
            }
            if set(archive.files) != expected_keys:
                raise ValueError("charging checkpoint contains missing or undeclared arrays")
            geometry = FeatureGeometry3D(
                archive["geometry_phi"].copy(),
                archive["geometry_material_id"].copy(),
                metadata["dx"], metadata["mesh_length_unit_m"],
                tuple(metadata["mesh_origin_m"]),
                material_levelsets=(None if not layer_ids else {
                    value: archive[f"material_levelset_{value}"].copy()
                    for value in layer_ids}))
            sigma = archive["sigma_c_per_m2"].copy()
        return cls(
            geometry=geometry, sigma_c_per_m2=sigma,
            surface_state_type=metadata["surface_state_type"],
            surface_state_fields=fields,
            surface_state_mesh_fingerprint=metadata[
                "surface_state_mesh_fingerprint"],
            completed_duration_s=metadata["completed_duration_s"],
            completed_steps=metadata["completed_steps"],
            source_manifest_sha256=metadata["source_manifest_sha256"],
            surface_state_upper_bounds=metadata.get("surface_state_upper_bounds", {}),
            surface_state_remap_modes=metadata.get("surface_state_remap_modes", {}),
            schema=metadata["schema"])
