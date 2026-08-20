#!/usr/bin/env python3
"""Run one preregistered Guo/Kwon -> Krueger forecast without target observables.

The numerical operator is selected by the committed deterministic-prefix time
and space gates.  This supervisor never imports the Krueger target profile or
depth.  It resumes bounded wall-time checkpoints until the requested horizon
or an explicitly declared physical clogging event is reached.
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
SUCCESS = frozenset({"complete", "terminal_feature_clogged"})
ION_CLOSURES = {
    "nominal_unresolved": None,
    "all_cf2": "CF2",
    "all_cf3": "CF3",
}


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _write(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def pilot_command(args: argparse.Namespace, *, resume: bool) -> list[str]:
    command = [
        sys.executable,
        str(PILOT),
        "--output", str(Path(args.output)),
        "--duration-s", "60",
        "--n-steps", "3840",
        "--dx-um", "0.01",
        "--substrate-top-um", "1.8",
        "--domain-height-um", "2.8",
        "--n-position", "16",
        "--compress-boundary-quadrature",
        "--neutral-transverse-order", "5",
        "--neutral-normal-order", "2",
        "--neutral-direction-polar-order", "8",
        "--neutral-direction-azimuthal-order", "16",
        "--ion-energy-bin-eV", "250",
        "--ion-angle-bin-deg", "0.25",
        "--ion-azimuthal-order", "16",
        "--ion-flux-normalization", "1",
        "--surface-model", "guo_tml",
        "--guo-translating-layer-thickness-nm", "2.5",
        "--ballistic-transport", "face_gather",
        "--transport-device", str(args.transport_device),
        "--face-quadrature-points", "3",
        "--radiosity-backend", "deterministic_extruded_2d",
        "--radiosity-tolerance", "1e-12",
        "--radiosity-max-iterations", "2000",
        "--exchange-method", "analytic_occlusion",
        "--exchange-geometry-tolerance", "1e-9",
        "--exchange-relative-tolerance", "1e-5",
        "--adaptive-profile-timestep",
        "--minimum-step-s", "1e-5",
        "--target-displacement-cells", "0.35",
        "--maximum-displacement-cells", "0.75",
        "--adaptive-shrink-factor", "0.5",
        "--adaptive-growth-factor", "1.5",
        "--adaptive-safety-factor", "0.9",
        "--maximum-accepted-steps", "10000",
        "--topology-change-policy", "continue_gas_cavity",
        "--surface-state-remap-backend", "common_refinement",
        "--max-wall-s", str(float(args.max_wall_s)),
        "--blind-execution",
    ]
    formula = ION_CLOSURES[str(args.case)]
    if formula is not None:
        command.extend(("--guo-aggregate-ion-formula", formula))
    if resume:
        command.append("--resume")
    return command


def run(args: argparse.Namespace) -> dict:
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    audit_path = output / "audit.json"
    checkpoint_path = output / "checkpoint.npz"
    attempts = 0
    while True:
        status = None
        if audit_path.exists():
            audit = json.loads(audit_path.read_text(encoding="utf-8"))
            status = audit.get("status")
            if audit.get("experimental_outcomes_read") is not False:
                raise ValueError(
                    "output directory contains a non-blind Krueger execution")
            if status in SUCCESS:
                break
        if attempts > int(args.maximum_resume_count):
            raise RuntimeError("forecast exhausted its bounded resume budget")
        resume = audit_path.exists() and checkpoint_path.exists()
        subprocess.run(pilot_command(args, resume=resume), cwd=ROOT, check=True)
        attempts += 1

    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    receipt = {
        "schema": "petch.guo-krueger.no-fit-forecast-execution.v1",
        "scientific_status": (
            "blind published-boundary transfer sensitivity; no Krueger target "
            "observable was read during execution"
        ),
        "case": str(args.case),
        "aggregate_ion_identity_closure": ION_CLOSURES[str(args.case)],
        "ion_flux_normalization": 1.0,
        "guo_translating_layer_thickness_nm": 2.5,
        "nominal_time_step_s": 0.015625,
        "dx_nm": 10.0,
        "status": audit["status"],
        "config_hash": audit["config_hash"],
        "audit_sha256": _sha256(audit_path),
        "checkpoint_sha256": _sha256(checkpoint_path),
        "final_metrics": audit.get("final_metrics"),
        "terminal_event": audit.get("terminal_event"),
        "experimental_outcomes_read": audit.get("experimental_outcomes_read"),
        "invocations_this_call": attempts,
        "source_sha256": {
            "pilot": _sha256(PILOT),
            "guo_mechanism": _sha256(ROOT / "src" / "petch" / "guo_c4f8_sio2_feature.py"),
            "material_router": _sha256(ROOT / "src" / "petch" / "amorphous_carbon_mask.py"),
        },
    }
    _write(output / "forecast_execution.json", receipt)
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", choices=tuple(ION_CLOSURES), required=True)
    parser.add_argument("--transport-device", default="cpu")
    parser.add_argument("--max-wall-s", type=float, default=1800.0)
    parser.add_argument("--maximum-resume-count", type=int, default=200)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
