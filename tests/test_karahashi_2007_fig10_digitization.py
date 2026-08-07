import json
from pathlib import Path

import pytest

from scripts.digitize_karahashi_2007_fig10 import (
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


def test_product_board_is_complete_and_condition_specific(record):
    _, audit = record
    assert audit["derived_checks"]["point_count"] == 9
    assert audit["experiment_scope"]["incident_ion"] == "CF3+"
    assert audit["experiment_scope"]["energies_eV"] == [500, 1000, 2000]
    assert audit["experiment_scope"]["detected_products"] == [
        "SiF", "SiF2", "SiF4"]


def test_pil_render_axes_and_source_are_pinned(record):
    _, audit = record
    assert audit["source"]["pdf_sha256"] == (
        "093b18b91b0a6d910fc414779ee8320b7a046ac4cad38ef5de0b7f2dd25a2d79")
    assert audit["source"]["rendered_page_sha256"] == (
        "e9e49375c53aa95e819302f434b68068f784219d118e2cf2fa3bb78c9c19a352")
    crop = verify_render(DEFAULT_RENDER)
    assert crop.size == (1550, 1400)
    assert audit["derived_checks"][
        "maximum_marker_x_energy_offset_eV"] < 1.5
    visual = audit["digitization"]["full_resolution_visual_inspection"]
    assert visual["passed"]
    assert visual["overlay_sha256"] == (
        "71b556bf7f71afe18c0680ae678ae5f4ce8de082e62bedca778354d964c239fe")


def test_product_fractions_close_the_plotted_normalization(record):
    _, audit = record
    sums = audit["derived_checks"]["fraction_sum_percent_by_energy"]
    assert sums["500"] == pytest.approx(99.5, abs=0.1)
    assert sums["1000"] == pytest.approx(99.9, abs=0.1)
    assert sums["2000"] == pytest.approx(98.9, abs=0.1)
    assert all(abs(total - 100.0) < 1.2 for total in sums.values())


def test_source_trends_and_claim_boundary_are_preserved(record):
    _, audit = record
    assert audit["derived_checks"][
        "dominant_product_at_or_below_1000_eV"] == "SiF2"
    assert audit["derived_checks"]["sif_fraction_increases_with_energy"]
    assert audit["derived_checks"]["sif4_fraction_decreases_with_energy"]
    assert not audit["claim_boundary"]["production_escape_parameter_use"]
    assert "a product escape probability or diffusion length" in (
        audit["claim_boundary"]["not_valid"])


def test_committed_digitization_is_current(record):
    payload, audit = record
    assert CSV_PATH.read_text(encoding="utf-8") == payload
    assert json.loads(MANIFEST_PATH.read_text(encoding="utf-8")) == audit
