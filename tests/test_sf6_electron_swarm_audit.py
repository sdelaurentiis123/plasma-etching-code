import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "results" / "curated" / "sf6_electron_swarm_v1" / "audit.json"


def test_sf6_audit_keeps_source_replay_transport_and_depth_claims_separate():
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    assert audit["schema"] == "petch.sf6_electron_swarm_audit.v1"
    certification = audit["certification"]
    assert certification["supports_use_as_sf6_component_collision_input"] is True
    assert certification["primary_positive_and_negative_branching_resolved"] is False
    assert certification["bulk_flux_nonconservative_transport_resolved"] is False
    assert certification["independently_validated_in_target_mixture_band"] is False
    assert certification["supports_unique_target_reactor_state"] is False
    assert certification["supports_wafer_flux"] is False
    assert certification["supports_feature_depth"] is False
    assert certification["feature_depth_used"] is False


def test_sf6_low_energy_attachment_boundary_layer_is_numerically_converged():
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    convergence = audit["numerical_convergence"]
    assert convergence["coarse_actual_cell_count"] == 892
    assert convergence["fine_actual_cell_count"] == 966
    assert convergence[
        "maximum_absolute_flux_drift_relative_change"
    ] < 1.0e-6
    assert convergence[
        "maximum_absolute_attachment_rate_relative_change"
    ] < 1.0e-6


def test_sf6_transport_definition_mismatch_is_exposed_not_fitted_away():
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    diagnostics = audit["transport_definition_diagnostics"]
    assert diagnostics["drift"]["measurement_equivalent_grade"] is False
    assert diagnostics["effective_ionization"][
        "measurement_equivalent_grade"
    ] is False
    critical = diagnostics["critical_field"]
    assert 430.0 < critical["predicted_temporal_zero_Td"] < 450.0
    assert critical["source_spatial_townsend_zero_Td"] == 359.3
    assert critical["measurement_equivalent_grade"] is False


def test_sf6_aggregate_vibration_closure_is_subpercent_on_converged_grid():
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    sensitivity = audit["closure_sensitivity"]
    assert max(abs(item) for item in sensitivity["flux_drift_relative_change"]) < .004
