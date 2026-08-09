import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from petch import Hirsch2020PulsedDCAntiSynergySensitivity


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "experimental" / "hirsch_2020_pae_iae"


def test_digitization_manifest_hash_and_curve_match_runtime_card():
    manifest = json.loads(
        (DATA / "digitization_manifest.json").read_text(encoding="utf-8"))
    csv_path = DATA / "figure8_relative_pae_yield.csv"
    assert sha256(csv_path.read_bytes()).hexdigest() == manifest["output"][
        "csv_sha256"]
    with csv_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    duty = np.asarray([
        float(row["dc_bias_duty_cycle_percent"]) for row in rows])
    relative = np.asarray([float(row["relative_pae_yield"]) for row in rows])
    card = Hirsch2020PulsedDCAntiSynergySensitivity()
    np.testing.assert_allclose(duty, card.duty_cycle_percent, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(relative, card.relative_pae_yield, rtol=0.0, atol=0.0)
    assert manifest["extraction"]["feature_depth_used"] is False
    assert manifest["evidence_class"]["not"].startswith(
        "a direct spectral photon-flux")
