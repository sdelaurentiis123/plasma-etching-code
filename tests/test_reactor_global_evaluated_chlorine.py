import csv
from pathlib import Path

import numpy as np
import pytest

from petch.reactor_global import (
    ATOMIC_CHLORINE_IONIZATION_THRESHOLD_EV,
    RateContext,
    nist_hayes_atomic_chlorine_ionization_rate,
)

ROOT = Path(__file__).resolve().parents[1]
TABLE25 = (
    ROOT / "research_sources" / "digitized"
    / "christophorou_olthoff_1999_table25_atomic_cl_ionization.csv"
)
TABLE25_MANIFEST = TABLE25.with_name(
    "christophorou_olthoff_1999_table25_manifest.md")


def test_nist_hayes_table25_transcription_and_evidence():
    coefficient = nist_hayes_atomic_chlorine_ionization_rate()
    assert len(coefficient.electron_energy_eV) == 48
    assert coefficient.electron_energy_eV[:4] == (11.0, 12.0, 13.0, 14.0)
    assert coefficient.electron_energy_eV[-4:] == (
        170.0, 180.0, 190.0, 200.0)
    np.testing.assert_allclose(
        np.asarray(coefficient.cross_section_m2[:4]) / 1.0e-20,
        [0.00, 0.01, 0.02, 0.24],
        rtol=0.0,
        atol=2.0e-16,
    )
    np.testing.assert_allclose(
        np.asarray(coefficient.cross_section_m2[-4:]) / 1.0e-20,
        [2.81, 2.72, 2.68, 2.63],
        rtol=0.0,
        atol=5.0e-16,
    )
    assert coefficient.threshold_eV == (
        ATOMIC_CHLORINE_IONIZATION_THRESHOLD_EV)
    assert coefficient.relative_uncertainty == 0.14
    assert coefficient.evidence_kind == "measured"


def test_nist_hayes_executable_table_matches_pixel_audited_csv():
    coefficient = nist_hayes_atomic_chlorine_ionization_rate()
    with TABLE25.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    np.testing.assert_allclose(
        coefficient.electron_energy_eV,
        [float(row["electron_energy_eV"]) for row in rows],
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_allclose(
        np.asarray(coefficient.cross_section_m2) / 1.0e-20,
        [float(row["cross_section_1e_minus_20_m2"]) for row in rows],
        rtol=0.0,
        atol=5.0e-16,
    )
    manifest = TABLE25_MANIFEST.read_text(encoding="utf-8")
    assert "original-resolution visual review of all 48" in manifest
    assert (
        "6a01d03172e2d49619998e0593d14e9b547ad01803f45a6677068768cf599c25"
        in manifest
    )


@pytest.mark.parametrize("temperature", [2.0, 3.0, 5.0, 8.0, 10.0])
def test_nist_hayes_rate_is_positive_on_reactor_temperature_domain(
        temperature):
    coefficient = nist_hayes_atomic_chlorine_ionization_rate()
    assert coefficient.coefficient_si(RateContext(temperature)) > 0.0
    assert coefficient.maxwellian_kernel_tail_fraction(temperature) <= 1.0e-6


def test_nist_hayes_rate_rejects_temperature_with_material_unknown_tail():
    coefficient = nist_hayes_atomic_chlorine_ionization_rate()
    with pytest.raises(ValueError, match="unmeasured cross-section tail"):
        coefficient.coefficient_si(RateContext(15.0))
