#!/usr/bin/env python3
"""Extract the preregistered Bosch calibration rows without parsing targets.

The official 89-point CSV mixes calibration and held-out wafers.  This broker
uses only the first comma-delimited field as a key, copies allowed calibration
rows byte-for-byte, and never converts any outcome column.  The downstream
fitter reads only the extracted asset and its manifest.
"""
from __future__ import annotations

import argparse
from hashlib import md5, sha256
import json
from pathlib import Path

from petch.bosch_process_data import WAFER_MEASUREMENT_89_POINT_MD5


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "zenodo_17122442"
SOURCE = DATA / "Si_Oxide_etch_89_points.csv"
OUTPUT = DATA / "calibration_Si_Oxide_etch_89_points.csv"
MANIFEST = DATA / "calibration_Si_Oxide_etch_89_points_manifest.json"
PREREGISTRATION = (
    ROOT / "results" / "curated" / "zenodo_bosch_reactor_depth_holdout_v1"
    / "preregistration.json"
)
EXPECTED_HEADER = (
    b"experiment_key,lot_number,wafer_number,X,Y,preox_thickness,"
    b"postox_thickness,postox_thickness_nan,stepheight,oxide_etch,si_etch"
)


def _sha256(payload: bytes) -> str:
    return sha256(payload).hexdigest()


def build_assets(source_path=SOURCE, preregistration_path=PREREGISTRATION):
    source_path = Path(source_path)
    preregistration_path = Path(preregistration_path)
    source = source_path.read_bytes()
    if md5(source).hexdigest() != WAFER_MEASUREMENT_89_POINT_MD5:
        raise ValueError("Bosch 89-point source checksum mismatch")
    preregistration_bytes = preregistration_path.read_bytes()
    preregistration = json.loads(preregistration_bytes)
    split = preregistration["split_rule"]
    allowed = frozenset(split["calibration_experiment_keys"])
    forbidden = frozenset(split["heldout_experiment_keys"])
    if not allowed or allowed & forbidden:
        raise ValueError("invalid Bosch preregistration key split")

    lines = source.splitlines(keepends=True)
    if not lines or lines[0].rstrip(b"\r\n") != EXPECTED_HEADER:
        raise ValueError("unexpected Bosch 89-point source header")
    output = [lines[0]]
    copied_keys = set()
    heldout_keys = set()
    unknown_keys = set()
    copied_rows = 0
    heldout_rows = 0
    for line in lines[1:]:
        if not line.strip():
            continue
        key_bytes, separator, _numeric_payload = line.partition(b",")
        if not separator:
            raise ValueError("invalid Bosch source row")
        key = key_bytes.decode("ascii")
        if key in allowed:
            output.append(line)
            copied_keys.add(key)
            copied_rows += 1
        elif key in forbidden:
            heldout_keys.add(key)
            heldout_rows += 1
        else:
            unknown_keys.add(key)
    if unknown_keys:
        raise ValueError(f"source contains keys outside the frozen split: {unknown_keys}")
    if copied_rows != 89 * len(copied_keys):
        raise ValueError("calibration extraction contains incomplete 89-point wafers")
    if heldout_rows != 89 * len(heldout_keys):
        raise ValueError("heldout source contains incomplete 89-point wafers")

    output_bytes = b"".join(output)
    missing_calibration = sorted(allowed - copied_keys)
    missing_heldout = sorted(forbidden - heldout_keys)
    manifest = {
        "schema": "petch-bosch-calibration-measurement-extraction-v1",
        "source_record": "https://zenodo.org/records/17122442",
        "source_file": source_path.name,
        "source_md5": md5(source).hexdigest(),
        "source_row_count": len(lines) - 1,
        "preregistration_file": str(preregistration_path.relative_to(ROOT)),
        "preregistration_sha256": _sha256(preregistration_bytes),
        "calibration_allowlist_sha256": _sha256(
            ("\n".join(sorted(allowed)) + "\n").encode("ascii")),
        "output_file": OUTPUT.name,
        "output_sha256": _sha256(output_bytes),
        "output_row_count": copied_rows,
        "output_experiment_key_count": len(copied_keys),
        "missing_calibration_process_keys": missing_calibration,
        "excluded_heldout_row_count": heldout_rows,
        "excluded_heldout_experiment_key_count": len(heldout_keys),
        "missing_heldout_process_keys": missing_heldout,
        "splitter_numeric_outcome_fields_parsed": False,
        "mixed_source_bytes_opened_by_splitter": True,
        "heldout_rows_copied_to_fit_asset": False,
        "fit_process_must_open_only_output_file": True,
    }
    manifest_bytes = (
        json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8")
    return output_bytes, manifest_bytes


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    output, manifest = build_assets()
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_bytes() != output:
            raise SystemExit("calibration measurement asset is stale")
        if not MANIFEST.exists() or MANIFEST.read_bytes() != manifest:
            raise SystemExit("calibration measurement manifest is stale")
        return 0
    OUTPUT.write_bytes(output)
    MANIFEST.write_bytes(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
