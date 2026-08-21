#!/usr/bin/env python3
"""Reveal only the Bosch heldout rows after verifying the committed v8 seal."""
from __future__ import annotations

import argparse
from hashlib import md5, sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from petch.bosch_process_data import WAFER_MEASUREMENT_89_POINT_MD5
from scripts.extract_bosch_calibration_measurements import EXPECTED_HEADER
from scripts.seal_bosch_wafer_boundary_map_heldout_prediction import (
    OUTPUT as PREDICTION,
    SEAL as PREDICTION_SEAL,
    V1_PREREGISTRATION,
)

DATA = ROOT / "data" / "experimental" / "zenodo_17122442"
SOURCE = DATA / "Si_Oxide_etch_89_points.csv"
OUTPUT = DATA / "revealed_heldout_Si_Oxide_etch_89_points.csv"
MANIFEST = DATA / "revealed_heldout_Si_Oxide_etch_89_points_manifest.json"


def _hash_bytes(payload):
    return sha256(payload).hexdigest()


def _hash(path):
    return _hash_bytes(Path(path).read_bytes())


def build_assets(
        source_path=SOURCE, preregistration_path=V1_PREREGISTRATION,
        prediction_path=PREDICTION, seal_path=PREDICTION_SEAL):
    preregistration_bytes = Path(preregistration_path).read_bytes()
    preregistration = json.loads(preregistration_bytes)
    prediction_bytes = Path(prediction_path).read_bytes()
    prediction = json.loads(prediction_bytes)
    seal_bytes = Path(seal_path).read_bytes()
    seal = json.loads(seal_bytes)
    if (
        seal["prediction_sha256"] != _hash_bytes(prediction_bytes)
        or seal["heldout_outcomes_read"] is not False
        or seal["heldout_prediction_written"] is not True
        or seal["eligible_for_separate_outcome_score_after_commit"] is not True
        or prediction["target_firewall"]["heldout_outcomes_read"] is not False
        or prediction["target_firewall"]["heldout_prediction_written"] is not True
    ):
        raise RuntimeError("Bosch heldout prediction is not hash sealed")
    allowed = frozenset(
        preregistration["split_rule"]["heldout_experiment_keys"])
    forbidden = frozenset(
        preregistration["split_rule"]["calibration_experiment_keys"])
    if len(allowed) != 20 or allowed & forbidden:
        raise RuntimeError("Bosch chronological split is stale")

    source = Path(source_path).read_bytes()
    if md5(source).hexdigest() != WAFER_MEASUREMENT_89_POINT_MD5:
        raise ValueError("Bosch mixed 89-point source checksum mismatch")
    lines = source.splitlines(keepends=True)
    if not lines or lines[0].rstrip(b"\r\n") != EXPECTED_HEADER:
        raise ValueError("unexpected Bosch mixed 89-point source header")
    output = [lines[0]]
    copied_keys = set()
    calibration_keys = set()
    unknown_keys = set()
    copied_rows = 0
    for line in lines[1:]:
        if not line.strip():
            continue
        key_bytes, separator, _numeric_payload = line.partition(b",")
        if not separator:
            raise ValueError("invalid Bosch mixed outcome row")
        key = key_bytes.decode("ascii")
        if key in allowed:
            output.append(line)
            copied_keys.add(key)
            copied_rows += 1
        elif key in forbidden:
            calibration_keys.add(key)
        else:
            unknown_keys.add(key)
    if unknown_keys:
        raise ValueError(
            f"mixed source contains keys outside the frozen split: {unknown_keys}")
    if copied_rows != 89 * len(copied_keys):
        raise ValueError("revealed Bosch heldout asset has an incomplete wafer")
    if len(copied_keys) != 13 or len(calibration_keys) != 75:
        raise ValueError("unexpected measured Bosch heldout/calibration split")

    output_bytes = b"".join(output)
    manifest = {
        "schema": "petch-bosch-heldout-measurement-reveal-v1",
        "source_file": Path(source_path).name,
        "source_md5": md5(source).hexdigest(),
        "source_row_count": len(lines) - 1,
        "output_file": OUTPUT.name,
        "output_sha256": _hash_bytes(output_bytes),
        "output_row_count": copied_rows,
        "output_experiment_key_count": len(copied_keys),
        "revealed_experiment_keys": sorted(copied_keys),
        "missing_heldout_process_keys": sorted(allowed - copied_keys),
        "calibration_rows_copied": False,
        "numeric_outcome_fields_parsed_by_broker": False,
        "prediction_committed_before_reveal": True,
        "prediction_sha256": _hash_bytes(prediction_bytes),
        "prediction_seal_sha256": _hash_bytes(seal_bytes),
        "preregistration_sha256": _hash_bytes(preregistration_bytes),
        "heldout_outcomes_revealed_after_seal": True,
    }
    return output_bytes, (
        json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    if args.write == args.check:
        parser.error("select exactly one of --write or --check")
    output, manifest = build_assets()
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_bytes() != output:
            raise SystemExit("revealed Bosch heldout outcome asset is stale")
        if not MANIFEST.exists() or MANIFEST.read_bytes() != manifest:
            raise SystemExit("revealed Bosch heldout outcome manifest is stale")
        print("Bosch heldout outcome reveal is current")
        return 0
    OUTPUT.write_bytes(output)
    MANIFEST.write_bytes(manifest)
    print(OUTPUT.relative_to(ROOT))
    print(MANIFEST.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
