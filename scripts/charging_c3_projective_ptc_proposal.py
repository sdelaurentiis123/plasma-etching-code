#!/usr/bin/env python3
"""Create a safeguarded compatible-Q1 projective/PTC checkpoint proposal."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from charging_task1_physical_time_3d import _geometry_and_poisson  # noqa: E402
from petch.charging_coevolution_3d import (  # noqa: E402
    propose_compatible_q1_pseudo_time_step_3d,
)


def _hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _atomic_json(path: Path, value) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    os.replace(temporary, path)


def _atomic_npz(path: Path, **arrays) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
    os.replace(temporary, path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--current-audit", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--pseudo-timestep-s", type=float, required=True)
    parser.add_argument("--maximum-potential-jump-v", type=float, default=5.0)
    parser.add_argument("--grid-dx-um", type=float, default=0.25)
    parser.add_argument("--sampling-epoch-offset", type=int, default=1)
    args = parser.parse_args()
    if args.sampling_epoch_offset < 1:
        parser.error("--sampling-epoch-offset must reserve at least one unused audit epoch")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    with np.load(args.checkpoint) as checkpoint:
        payload = {name: np.asarray(checkpoint[name]).copy() for name in checkpoint.files}
    required = {
        "sigma_c_per_m2", "face_charge_c", "charge_node_c", "potential_v",
        "vertices", "faces", "areas", "resume_sampling_epoch", "scramble_mode",
        "scramble_base_seed", "sampling_seed_stride"}
    if not required.issubset(payload):
        parser.error("checkpoint lacks required C3 replay state")
    if str(np.asarray(payload["scramble_mode"]).item()) != "fresh":
        parser.error("projective stochastic PTC requires a fresh-scramble checkpoint")
    if ("compatible_q1_charge_state" in payload
            and not bool(np.asarray(payload["compatible_q1_charge_state"]).item())):
        parser.error("projective PTC requires a compatible-Q1 checkpoint")

    with np.load(args.current_audit) as audit:
        if "ensemble_mean_positive_face_current_density_a_m2" in audit:
            positive = np.asarray(
                audit["ensemble_mean_positive_face_current_density_a_m2"], dtype=float)
            negative = np.asarray(
                audit["ensemble_mean_negative_face_current_density_a_m2"], dtype=float)
            current_estimator = "independent fixed-state ensemble mean"
            current_replicates = int(np.asarray(audit["replicate_count"]).item())
        else:
            if ("terminal_window_ready" not in audit
                    or not bool(np.asarray(audit["terminal_window_ready"]).item())):
                parser.error("current audit lacks a complete terminal-window estimate")
            positive = np.asarray(
                audit["terminal_window_positive_face_current_density_a_m2"], dtype=float)
            negative = np.asarray(
                audit["terminal_window_negative_face_current_density_a_m2"], dtype=float)
            current_estimator = "physical-time terminal-window mean"
            current_replicates = None
    if (positive.shape != np.asarray(payload["sigma_c_per_m2"]).shape
            or negative.shape != positive.shape
            or np.any(~np.isfinite(positive)) or np.any(~np.isfinite(negative))
            or np.any(positive < 0.0) or np.any(negative < 0.0)):
        parser.error("invalid terminal-window face currents")

    geometry, poisson = _geometry_and_poisson(args.grid_dx_um)
    proposal = propose_compatible_q1_pseudo_time_step_3d(
        poisson, payload["vertices"], payload["faces"], payload["areas"],
        payload["sigma_c_per_m2"], positive - negative,
        args.pseudo_timestep_s,
        potential_origin=(0.0, 0.0, 0.0), potential_spacing=geometry.dx,
        mesh_length_unit_m=geometry.mesh_length_unit_m,
        maximum_potential_jump_v=args.maximum_potential_jump_v)
    source_epoch = int(np.asarray(payload["resume_sampling_epoch"]).item())
    audit_epoch = source_epoch + int(args.sampling_epoch_offset)
    payload.update(
        sigma_c_per_m2=proposal.sigma_c_per_m2,
        face_charge_c=proposal.face_charge_c,
        charge_node_c=proposal.charge_node_c,
        potential_v=proposal.potential_v,
        resume_sampling_epoch=np.asarray(audit_epoch),
        compatible_q1_charge_state=np.asarray(True),
        projective_source_sampling_epoch=np.asarray(source_epoch),
        projective_audit_sampling_epoch=np.asarray(audit_epoch),
        projective_pseudo_timestep_s=np.asarray(args.pseudo_timestep_s))
    checkpoint_path = output / "face_checkpoint.npz"
    _atomic_npz(checkpoint_path, **payload)
    manifest = dict(
        schema="petch.charging.c3.compatible-projective-ptc-proposal.v1",
        status="proposal only; not accepted and not a convergence claim",
        pseudo_timestep_s=float(args.pseudo_timestep_s),
        source_physical_time_advanced_s=0.0,
        source_sampling_epoch=source_epoch,
        reserved_fresh_audit_sampling_epoch=audit_epoch,
        current_estimator=current_estimator,
        current_replicates=current_replicates,
        maximum_potential_jump_v=float(args.maximum_potential_jump_v),
        proposal_diagnostics=dict(proposal.diagnostics),
        audit_contract=(
            "score candidates with paired fresh exact hard-visibility kinetic evaluations; "
            "select on Q1-resolved residual/potential rate; confirm the selected candidate at "
            "an unused epoch and with a fixed-physical-time continuation"),
        provenance=dict(
            source_checkpoint_name=args.checkpoint.name,
            source_checkpoint_sha256=_hash(args.checkpoint),
            source_current_audit_name=args.current_audit.name,
            source_current_audit_sha256=_hash(args.current_audit),
            proposal_script_sha256=_hash(Path(__file__).resolve()),
            charging_coevolution_sha256=_hash(
                ROOT / "src/petch/charging_coevolution_3d.py"),
            charging_poisson_sha256=_hash(ROOT / "src/petch/charging_poisson_3d.py")),
        artifact=dict(name=checkpoint_path.name, sha256=_hash(checkpoint_path)))
    manifest_path = output / "proposal.json"
    _atomic_json(manifest_path, manifest)
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
