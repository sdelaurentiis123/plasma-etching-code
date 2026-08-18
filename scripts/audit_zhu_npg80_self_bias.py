#!/usr/bin/env python3
"""Build/check the target-free Oxford-80 self-bias sensitivity receipt."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from petch.reactor_global.oxford80_self_bias import (
    Oxford80RIECondition,
    build_oxford80_self_bias_transfer,
    load_oxford80_self_bias_evidence,
)


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "experimental" / "zhu_2026_tio2_npg80"
DEFAULT_EVIDENCE = DATA_DIR / "oxford80_self_bias_evidence.csv"
DEFAULT_OUTPUT = DATA_DIR / "oxford80_self_bias_transfer.json"
DEFAULT_MANIFEST = DATA_DIR / "recipe_manifest.json"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest(path: Path = DEFAULT_MANIFEST) -> dict:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("schema") != "petch.experimental-recipe.v1":
        raise ValueError("unexpected recipe-manifest schema")
    if payload.get("measurement_state") != "pre_sem_specific_condition":
        raise ValueError("condition is not marked as a pre-SEM commitment")
    if payload["outcomes"]["specific_condition_sem_received"]:
        raise ValueError("self-bias transfer cannot consume a revealed SEM")
    image = ROOT / payload["source_image"]["path"]
    if _sha256(image) != payload["source_image"]["sha256"]:
        raise ValueError("recipe screenshot checksum mismatch")
    return payload


def build_target_condition(manifest: dict) -> Oxford80RIECondition:
    process = manifest["process"]
    return Oxford80RIECondition(
        tool_model=manifest["tool"]["model"],
        rf_power_W=process["table_rf_forward_power_setpoint_W"],
        pressure_mTorr=1000.0 * process["pressure_Torr"],
        gas_flows_sccm=process["gases_sccm"],
        electrode_temperature_C=process["table_temperature_C"],
        duration_s=process["etch_time_s"],
    )


def build_self_bias_receipt(
    *,
    evidence_path: Path = DEFAULT_EVIDENCE,
) -> dict:
    manifest = load_manifest()
    observations = load_oxford80_self_bias_evidence(evidence_path)
    transfer = build_oxford80_self_bias_transfer(
        build_target_condition(manifest), observations)

    return {
        "schema": "petch.oxford80-self-bias-transfer.v1",
        "condition_id": manifest["condition_id"],
        "frozen_before_specific_condition_sem": True,
        "sem_target_used": transfer.sem_target_used,
        "measured_target_bias_used": transfer.measured_target_bias_used,
        "voltage_convention": (
            "magnitude of negative powered-electrode DC self-bias"
        ),
        "evidence_table": {
            "path": str(evidence_path.relative_to(ROOT)),
            "sha256": _sha256(evidence_path),
            "observation_count": len(observations),
        },
        "target": {
            "tool_model": transfer.target.tool_model,
            "rf_power_W": transfer.target.rf_power_W,
            "pressure_mTorr": transfer.target.pressure_mTorr,
            "reduced_drive_W_per_mTorr": (
                transfer.target.reduced_drive_W_per_mTorr
            ),
            "gas_flows_sccm": dict(transfer.target.gas_flows_sccm),
            "duration_s": transfer.target.duration_s,
        },
        "mechanical_anchor_selection": {
            "rule": (
                "same active-gas set, uncensored voltage, then minimum "
                "absolute log separation in reported forward-power/pressure"
            ),
            "source_id": transfer.matched_chemistry_reduced_drive_source_id,
            "anchor_V": transfer.matched_chemistry_reduced_drive_anchor_V,
            "is_target_measurement": False,
            "supports_unique_target_bias": False,
        },
        "printed_reference_window": {
            "lower_V": transfer.printed_reference_window_V[0],
            "upper_V": transfer.printed_reference_window_V[1],
            "is_probability_interval": (
                transfer.printed_window_is_probability_interval
            ),
            "contains_all_physically_allowed_target_values": False,
            "warning": (
                "Exact-NGP80 endpoint statements are censored (>300 V at "
                "start and <~200 V at end); this printed-source window is "
                "not a credible interval."
            ),
        },
        "observations": [
            {
                "source_id": item.source_id,
                "tool_model": item.tool_model,
                "tool_relation": item.tool_relation,
                "rf_power_W": item.rf_power_W,
                "pressure_mTorr": item.pressure_mTorr,
                "reduced_drive_W_per_mTorr": (
                    item.reduced_drive_W_per_mTorr
                ),
                "gas_flows_sccm": dict(item.gas_flows_sccm),
                "bias_relation": item.bias_relation,
                "bias_lower_V": item.bias_lower_V,
                "bias_upper_V": item.bias_upper_V,
                "run_phase": item.run_phase,
                "source_page": item.source_page,
                "source_pdf_sha256": item.source_pdf_sha256,
                "source_locator": item.source_locator,
                "loading_context": item.loading_context,
            }
            for item in observations
        ],
        "sensitivity_histories": [
            {
                "name": history.name,
                "time_s": history.time_s.tolist(),
                "bias_magnitude_V": history.bias_magnitude_V.tolist(),
                "source_ids": list(history.source_ids),
                "interpretation": history.interpretation,
                "endpoints_are_censor_thresholds": (
                    history.endpoints_are_censor_thresholds
                ),
                "measured_on_target_condition": (
                    history.measured_on_target_condition
                ),
                "supports_absolute_depth_prediction": (
                    history.supports_absolute_depth_prediction
                ),
            }
            for history in transfer.histories
        ],
        "certification": {
            "claim_class": (
                "target-free machine-family self-bias sensitivity transfer"
            ),
            "supports_unique_target_bias": transfer.supports_unique_target_bias,
            "supports_absolute_depth_prediction": (
                transfer.supports_absolute_depth_prediction
            ),
            "permitted_use": (
                "propagate each deterministic voltage history through the "
                "collisional sheath and report depth/profile sensitivity"
            ),
            "forbidden_use": (
                "select a history after viewing the held-out SEM or call the "
                "276 V family anchor a measured target-machine boundary"
            ),
        },
    }


def _canonical_json(payload: dict) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if the committed transfer differs from a clean rebuild",
    )
    args = parser.parse_args()
    rendered = _canonical_json(build_self_bias_receipt(
        evidence_path=args.evidence))
    if args.check:
        if not args.output.exists():
            raise SystemExit(f"missing committed transfer: {args.output}")
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("committed Oxford-80 self-bias transfer is stale")
        print(rendered, end="")
        return
    args.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
