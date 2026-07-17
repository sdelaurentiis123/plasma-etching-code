#!/usr/bin/env python3
"""Run the frozen Krüger oxygen/power transfer inputs without reading their outcomes.

This is an execution supervisor, not a scorer.  It consumes a checksum-bound physics reveal and
condition-specific published HPEM boundary inputs already handled by the common pilot.  Completed
runs and physical clogging events are preserved; wall-budget checkpoints resume automatically.
The experimental transfer-observation table is not imported or opened here.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
PILOT = ROOT / "scripts" / "krueger_2024_trench_pilot.py"
CASES = (
    ("oxygen_o05", ("--boundary-case", "oxygen_ratio", "--oxygen-ratio", "0.5",
                    "--low-frequency-power-kw", "6")),
    ("oxygen_o10", ("--boundary-case", "oxygen_ratio", "--oxygen-ratio", "1.0",
                    "--low-frequency-power-kw", "6")),
    ("oxygen_o15", ("--boundary-case", "oxygen_ratio", "--oxygen-ratio", "1.5",
                    "--low-frequency-power-kw", "6")),
    ("oxygen_o25", ("--boundary-case", "oxygen_ratio", "--oxygen-ratio", "2.5",
                    "--low-frequency-power-kw", "6")),
    ("power_0kw", ("--boundary-case", "power_sweep", "--low-frequency-power-kw", "0")),
    ("power_4kw", ("--boundary-case", "power_sweep", "--low-frequency-power-kw", "4")),
    ("power_6kw", ("--boundary-case", "power_sweep", "--low-frequency-power-kw", "6")),
    ("power_8kw", ("--boundary-case", "power_sweep", "--low-frequency-power-kw", "8")),
)
SUCCESS_STATUSES = frozenset({"complete"})


def _sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _write(path, payload):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _selected(case_set):
    if case_set == "all":
        return CASES
    prefix = "oxygen_" if case_set == "oxygen" else "power_"
    return tuple(item for item in CASES if item[0].startswith(prefix))


def _verify_frozen_files(freeze):
    """Refuse execution if code or reactor-boundary inputs changed after the reveal."""
    for section in ("source_sha256", "boundary_data_sha256"):
        mapping = freeze.get(section)
        if not isinstance(mapping, dict) or not mapping:
            raise ValueError(f"frozen reveal has no {section} manifest")
        for relative, expected in mapping.items():
            path = ROOT / relative
            if not path.is_file() or _sha(path) != expected:
                raise ValueError(f"frozen file checksum mismatch: {relative}")


def run(args):
    freeze_path = Path(args.freeze)
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    if (freeze.get("schema") != "petch.krueger-2024.frozen-physics-reveal.v2"
            or freeze.get("held_out_profile_data_read") is not False
            or not isinstance(freeze.get("reveal_sha256"), str)):
        raise ValueError("transfer campaign requires a sealed Krueger physics reveal")
    canonical = dict(freeze)
    claimed = canonical.pop("reveal_sha256")
    actual = sha256(json.dumps(
        canonical, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()
    if actual != claimed:
        raise ValueError("frozen physics reveal checksum does not match its content")
    _verify_frozen_files(freeze)
    physics = freeze["frozen_physics"]
    numerics = freeze["authority_numerics"]
    if (not abs(float(numerics.get("dx_um", 0.0)) - 0.005) <= 1e-15
            or numerics.get("topology_change_policy") != "continue_gas_cavity"
            or numerics.get("surface_state_remap_backend") not in (
                "indexed_knn", "common_refinement")):
        raise ValueError("held-out execution requires the sealed R1.9 authority operator")
    output_root = Path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    campaign_path = output_root / "campaign.json"
    campaign = {
        "schema": "petch.krueger-2024.held-out-transfer-execution.v1",
        "scientific_status": "blind transfer execution; experimental outcomes not yet revealed",
        "frozen_reveal_sha256": claimed,
        "freeze_file_sha256": _sha(freeze_path),
        "case_set": args.case_set,
        "transport_device": args.transport_device,
        "cases": {},
        "held_out_profile_data_read": False,
    }
    if campaign_path.exists():
        existing = json.loads(campaign_path.read_text(encoding="utf-8"))
        if existing.get("frozen_reveal_sha256") != claimed:
            raise ValueError("output directory contains a campaign for another reveal")
        campaign["cases"].update(existing.get("cases", {}))

    common = [
        sys.executable, str(PILOT),
        "--duration-s", str(numerics["duration_s"]),
        "--n-steps", str(numerics["n_steps"]),
        "--dx-um", str(numerics["dx_um"]),
        "--substrate-top-um", str(numerics["geometry"]["substrate_top_um"]),
        "--domain-height-um", str(numerics["geometry"]["domain_height_um"]),
        "--n-position", str(numerics["n_position"]),
        "--compress-boundary-quadrature",
        "--neutral-transverse-order", str(numerics["neutral_transverse_order"]),
        "--neutral-normal-order", str(numerics["neutral_normal_order"]),
        "--neutral-direction-polar-order",
        str(numerics["neutral_direction_polar_order"]),
        "--neutral-direction-azimuthal-order",
        str(numerics["neutral_direction_azimuthal_order"]),
        "--ion-energy-bin-eV", str(numerics["ion_energy_bin_eV"]),
        "--ion-angle-bin-deg", str(numerics["ion_angle_bin_deg"]),
        "--ion-azimuthal-order", str(numerics["ion_azimuthal_order"]),
        "--effective-mask-crosslinked-growth-fraction",
        str(physics["effective_mask_crosslinked_growth_fraction"]),
        "--oxide-etch-yield-scale", str(physics["oxide_etch_yield_scale"]),
        "--ballistic-transport", "face_gather",
        "--transport-device", str(args.transport_device),
        "--face-quadrature-points", str(numerics["ballistic_face_quadrature_points"]),
        "--radiosity-rays", str(numerics["radiosity_rays_per_face"]),
        "--radiosity-tolerance", str(numerics["radiosity_relative_tolerance"]),
        "--radiosity-max-iterations", str(numerics["radiosity_maximum_iterations"]),
        "--seed", str(numerics["seed"]),
        "--adaptive-profile-timestep",
        "--minimum-step-s", str(numerics["minimum_step_duration_s"]),
        "--target-displacement-cells", str(numerics["target_displacement_cells"]),
        "--maximum-displacement-cells", str(numerics["maximum_displacement_cells"]),
        "--adaptive-shrink-factor", str(numerics["adaptive_shrink_factor"]),
        "--adaptive-growth-factor", str(numerics["adaptive_growth_factor"]),
        "--adaptive-safety-factor", str(numerics["adaptive_safety_factor"]),
        "--maximum-accepted-steps", str(numerics["maximum_accepted_steps"]),
        "--topology-change-policy", str(numerics["topology_change_policy"]),
        "--surface-state-remap-backend",
        str(numerics["surface_state_remap_backend"]),
        "--max-wall-s", str(args.max_wall_s),
    ]
    for case_name, case_args in _selected(args.case_set):
        output = output_root / case_name
        audit_path = output / "audit.json"
        attempts = 0
        while True:
            status = None
            if audit_path.exists():
                status = json.loads(audit_path.read_text(encoding="utf-8")).get("status")
            if status in SUCCESS_STATUSES:
                break
            if attempts >= int(args.maximum_resume_count) + 1:
                raise RuntimeError(f"{case_name} exhausted its bounded resume budget")
            command = common + ["--output", str(output), *case_args]
            if audit_path.exists() and (output / "checkpoint.npz").exists():
                command.append("--resume")
            subprocess.run(command, cwd=ROOT, check=True)
            attempts += 1

        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        campaign["cases"][case_name] = {
            "status": audit["status"],
            "audit_path_name": f"{case_name}/audit.json",
            "audit_sha256": _sha(audit_path),
            "config_hash": audit["config_hash"],
            "final_metrics": audit.get("final_metrics"),
            "terminal_event": audit.get("terminal_event"),
            "resume_invocations_this_call": attempts,
            "boundary_provenance": audit.get("boundary_provenance"),
        }
        _write(campaign_path, campaign)
    canonical_campaign = json.dumps(campaign, sort_keys=True, separators=(",", ":"))
    campaign["execution_sha256"] = sha256(canonical_campaign.encode("utf-8")).hexdigest()
    _write(campaign_path, campaign)
    print(json.dumps({
        name: value["status"] for name, value in campaign["cases"].items()
    }, indent=2, sort_keys=True))
    return campaign


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--case-set", choices=("oxygen", "power", "all"), default="all")
    parser.add_argument("--transport-device", default="cpu")
    parser.add_argument("--max-wall-s", type=float, default=1800.0)
    parser.add_argument("--maximum-resume-count", type=int, default=3)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
