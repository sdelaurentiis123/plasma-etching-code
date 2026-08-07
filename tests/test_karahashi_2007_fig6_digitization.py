import json
from pathlib import Path

import pytest

from scripts.digitize_karahashi_2007_fig6 import (
    CSV_PATH,
    DEFAULT_RENDER,
    MANIFEST_PATH,
    csv_text,
    manifest,
    verify_render,
)


ROOT = Path(__file__).parents[1]


@pytest.fixture(scope="module")
def record():
    payload = csv_text()
    return payload, manifest(__import__("hashlib").sha256(
        payload.encode("utf-8")).hexdigest())


def test_angular_board_has_all_markers_and_unknown_energy(record):
    payload, audit = record
    assert audit["derived_checks"]["point_count"] == 20
    assert audit["experiment_scope"]["ion_energy_eV"] is None
    assert audit["experiment_scope"][
        "ion_energy_status"] == "not_reported_by_source"
    assert payload.count("not_reported_by_source") == 20


def test_pixel_setpoints_and_source_receipt_are_bounded(record):
    _, audit = record
    assert audit["source"]["pdf_sha256"] == (
        "093b18b91b0a6d910fc414779ee8320b7a046ac4cad38ef5de0b7f2dd25a2d79")
    assert audit["source"]["rendered_page_sha256"] == (
        "075cd97e5d7b5238f028f546d5c63cfc0c96cd6ff33a769b640f130d73b5fe1d")
    assert audit["derived_checks"][
        "maximum_marker_x_angle_offset_deg"] < 0.80


def test_pil_render_and_axes_are_verified():
    crop = verify_render(DEFAULT_RENDER)
    assert crop.size == (2050, 1100)


def test_source_text_ratio_cross_checks(record):
    _, audit = record
    ratios = audit["derived_checks"]["yield_ratio_60_to_0"]
    assert ratios["CF+"] == pytest.approx(2.2518, abs=2e-4)
    assert ratios["CF3+"] == pytest.approx(1.2849, abs=2e-4)
    assert ratios["CF+"] > ratios["CF2+"] > ratios["CF3+"]


def test_cross_figure_energy_is_inference_not_reported_fact(record):
    _, audit = record
    rmse = audit["derived_checks"][
        "normal_incidence_cross_figure_candidate_rmse"]
    assert rmse["1000"] < 0.02
    assert min(rmse[key] for key in ("750", "1500", "2000")) > 0.15
    assert audit["derived_checks"][
        "best_candidate_status"] == (
            "strong_inference_only_not_source_reported")
    assert not audit["claim_boundary"]["production_surface_model_use"]


def test_committed_digitization_is_current(record):
    payload, audit = record
    assert CSV_PATH.read_text(encoding="utf-8") == payload
    assert json.loads(MANIFEST_PATH.read_text(encoding="utf-8")) == audit
