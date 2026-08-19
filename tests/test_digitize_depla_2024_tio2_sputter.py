import csv
import json

from scripts.digitize_depla_2024_tio2_sputter import (
    CSV_PATH,
    MANIFEST_PATH,
    csv_text,
    manifest_text,
)


def _rows():
    with CSV_PATH.open(newline="", encoding="utf-8") as handle:
        return {int(row["argon_ion_energy_eV"]): row for row in csv.DictReader(handle)}


def test_committed_curve_and_receipt_replay_exactly():
    payload = csv_text()
    assert CSV_PATH.read_text(encoding="utf-8") == payload
    assert MANIFEST_PATH.read_text(encoding="utf-8") == manifest_text(payload)


def test_total_atom_axis_is_not_confused_with_formula_unit_yield():
    rows = _rows()
    for row in rows.values():
        atoms = float(row["fitted_total_atom_yield_per_ion"])
        units = float(row["stoichiometric_tio2_formula_units_per_ion"])
        assert abs(atoms / 3.0 - units) < 1.0e-6
        assert row["evidence_class"] == "digitized_semiempirical_fit_curve"
    assert abs(
        float(rows[276]["stoichiometric_tio2_formula_units_per_ion"])
        - 0.192143
    ) < 1.0e-6


def test_manifest_refuses_reactive_or_direct_measurement_promotion():
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    invalid = " ".join(manifest["claim_boundary"]["not_valid"])
    assert "direct low-energy" in invalid
    assert "reactive CHF3/SF6/O2" in invalid
    assert "Oxford NPG80 target coefficient" in invalid
    assert manifest["source_evidence"]["tio2_oxygen_surface_binding_energy_eV"] == 3.52
    assert manifest["source_evidence"]["tio2_oxygen_surface_binding_energy_95pct_eV"] == 1.14
