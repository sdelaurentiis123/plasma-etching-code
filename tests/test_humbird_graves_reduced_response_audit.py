import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT
    / "results"
    / "curated"
    / "humbird_graves_reduced_response"
    / "audit.json"
)


def test_calibration_excluded_surface_gates_pass_without_feature_depth():
    audit = json.loads(AUDIT.read_text())
    assert audit["all_held_out_gates_pass"]
    assert audit["feature_depth_values_loaded"] == []
    assert audit["reactor_flux_normalizations_loaded"] == []
    assert "not human-blind" in audit["calibration_exclusion_semantics"]
    gates = audit["held_out_gate_results"]
    assert gates[
        "carbon_and_fluorine_inventory_normalized_rmse"]["value"] < 0.2
    assert gates["instantaneous_yield_normalized_rmse"]["value"] < 0.2
    assert gates[
        "cumulative_silicon_etch_normalized_rmse"]["value"] < 0.15


def test_every_surface_coefficient_has_declared_provenance():
    audit = json.loads(AUDIT.read_text())
    provenance = audit["parameter_values_and_provenance"]["evidence"]
    assert "regressed" in provenance["carbon_kinetics"]
    assert "regressed" in provenance["fluorine_kinetics"]
    assert "regressed" in provenance["silicon_yield"]
    assert "one ion-renewed site" in provenance[
        "atomic_f_renewal_response"]
