import base64
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "partner-private" / "arun_resona_clockgate_2026"


def _payload():
    return json.loads(
        (CASE / "results" / "etch_explorer_data.json").read_text(
            encoding="utf-8"))


def _decode(value, shape):
    return np.frombuffer(base64.b64decode(value), dtype=np.uint8).reshape(shape)


def test_arun_explorer_payload_binds_exact_geometry_without_target_sem():
    payload = _payload()
    source = json.loads(
        (CASE / "results" / "source_scale_audit.json").read_text(
            encoding="utf-8"))
    assert payload["schema"] == "petch.partner.clockgate-etch-explorer.v1"
    assert payload["target_sem_used"] is False
    assert payload["refusal"]["unique_absolute_profile_ready"] is False
    assert (
        payload["inputs"]["polygon_receipt"]["source_gds_sha256"]
        == source["reference_gds"]["local_sha256"]
    )
    assert payload["inputs"]["footprint_um"] == [98.8, 62.8]
    assert payload["inputs"]["mask_height_um"] == 30.0


def test_arun_explorer_maps_have_one_common_periodic_grid():
    payload = _payload()
    grid = payload["geometry"]["grid"]
    shape = tuple(grid["shape_xy"])
    opening = _decode(grid["gas_opening"], shape)
    neutral = _decode(
        payload["transport"]["neutral_direct"]["transmission_map"], shape)
    assert shape == (247, 157)
    assert np.any(opening == 0) and np.any(opening == 255)
    assert np.all(neutral[opening == 0] == 0)
    for encoded in payload["transport"]["ion"]["transmission_maps"].values():
        ion = _decode(encoded, shape)
        assert np.all(ion[opening == 0] == 0)
        assert np.all(ion[opening > 0] >= neutral[opening > 0])


def test_arun_explorer_records_radical_entrance_as_dominant_gap():
    payload = _payload()
    ion = payload["transport"]["ion"]["statistics_over_opening"]["1.5"]
    neutral = payload["transport"]["neutral_direct"][
        "statistics_over_opening"]
    assert ion["mean"] > 0.9
    assert neutral["mean"] < 0.08
    assert ion["mean"] / neutral["mean"] > 10.0
    assert "printed-polymer F/O wall return" in payload[
        "surface_transfer"]["nonpredictive_inputs"]


def test_arun_explorer_is_a_self_contained_html_fragment():
    page = (CASE / "explorer" / "arun_etch_explorer.html").read_text(
        encoding="utf-8")
    assert "__ARUN_EXPLORER_DATA__" not in page
    assert "petch.partner.clockgate-etch-explorer.v1" in page
    assert "<style>" in page and "<script>" in page
    assert "<!DOCTYPE" not in page and "<html" not in page and "<body" not in page
