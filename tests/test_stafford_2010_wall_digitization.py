import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "digitize_stafford_2010_figure8.py"
DATA = ROOT / "data" / "experimental" / "stafford_2010"
CSV_PATH = DATA / "figure8_chlorine_wall_recombination.csv"
MANIFEST_PATH = DATA / "digitization_manifest.json"


def _module():
    spec = importlib.util.spec_from_file_location(
        "digitize_stafford_2010_figure8", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_stafford_committed_digitization_is_exactly_reproducible():
    module = _module()
    expected_csv = module.csv_text()
    assert CSV_PATH.read_text(encoding="utf-8") == expected_csv
    assert (
        MANIFEST_PATH.read_text(encoding="utf-8")
        == module.manifest_text(expected_csv)
    )


def test_stafford_figure8_marker_identity_and_domain_are_preserved():
    records = list(csv.DictReader(io.StringIO(
        CSV_PATH.read_text(encoding="utf-8"))))
    assert len(records) == 39
    assert sum(
        row["material"] == "anodized_aluminum" for row in records) == 23
    assert sum(row["material"] == "stainless_steel" for row in records) == 16
    assert {float(row["pressure_mTorr"]) for row in records} == {
        1.25, 5.0, 10.0, 20.0}
    assert all(
        row["reported_icp_power_range_W"] == "100-600"
        for row in records
    )
    assert all(row["observable_basis"] == "spinning_wall_LH_recombination"
               for row in records)


def test_stafford_digitization_matches_source_text_range_without_false_error():
    records = list(csv.DictReader(CSV_PATH.open(encoding="utf-8")))
    stainless = [
        float(row["cl_recombination_probability"])
        for row in records
        if row["material"] == "stainless_steel"
    ]
    assert min(stainless) == 0.0039811
    assert max(stainless) == 0.0315555

    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    assert manifest["digitization"]["marker_count"] == 39
    assert manifest["digitization"]["measurement_uncertainty"].startswith(
        "not assigned")
    assert any(
        "one constant gamma_Cl" in item
        for item in manifest["claim_boundary"]["not_valid"]
    )
    csv_hash = hashlib.sha256(
        CSV_PATH.read_bytes()).hexdigest()
    assert manifest["output"]["sha256"] == csv_hash
