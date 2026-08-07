import hashlib
import json

import pytest

from scripts.digitize_karahashi_2007_fig9 import (
    CSV_PATH,
    MANIFEST_PATH,
    PEAK_MARKERS,
    csv_text,
    manifest,
    rows,
)


@pytest.fixture(scope="module")
def record():
    payload = csv_text()
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    return payload, manifest(digest)


def test_digitization_is_exactly_replayable(record):
    payload, audit = record
    assert CSV_PATH.read_text(encoding="utf-8") == payload
    assert json.loads(MANIFEST_PATH.read_text(encoding="utf-8")) == audit


def test_all_product_peak_markers_and_timing_order_are_retained(record):
    _, audit = record
    assert [point.product for point in PEAK_MARKERS] == [
        "SiF",
        "SiF2",
        "SiF4",
    ]
    peak = {
        row["desorbed_product"]: float(
            row["peak_measured_sample_time_ms"]
        )
        for row in rows()
    }
    assert peak["SiF"] == pytest.approx(0.2391, abs=1.0e-4)
    assert peak["SiF2"] == pytest.approx(0.2993, abs=1.0e-4)
    assert peak["SiF4"] == pytest.approx(0.5403, abs=1.0e-4)
    assert peak["SiF"] < peak["SiF2"] < peak["SiF4"]
    assert audit["derived_checks"]["sif4_peak_is_near_half_millisecond"]


def test_beam_condition_and_instrument_convolution_are_explicit(record):
    _, audit = record
    scope = audit["experiment_scope"]
    assert scope["ion_energy_eV"] == 1000
    assert scope["incidence_angle_deg"] == 30
    assert scope["beam_pulse_fwhm_us"] == 100
    boundaries = " ".join(audit["claim_boundary"]["not_valid"])
    assert "instrument-response-deconvolved" in boundaries
    assert "diffusion coefficient" in boundaries
    assert not audit["claim_boundary"]["production_escape_parameter_use"]


def test_visual_audit_and_arbitrary_ordinate_boundary_are_pinned(record):
    _, audit = record
    digitization = audit["digitization"]
    assert digitization["full_resolution_visual_inspection"]["passed"]
    assert len(
        digitization["full_resolution_visual_inspection"]["overlay_sha256"]
    ) == 64
    assert digitization["time_digitization_bound_ms"] < 0.006
    assert "do not digitize intensity" in digitization["ordinate_policy"]
