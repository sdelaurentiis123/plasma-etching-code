import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT / "results" / "curated" / "chf3_electron_swarm_v1"
    / "audit.json"
)


def test_chf3_swarm_receipt_keeps_source_replay_and_depth_claims_separate():
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    assert audit["schema"] == "petch.chf3_electron_swarm_audit.v1"
    certification = audit["certification"]
    assert certification["supports_use_as_electron_transport_input"] is True
    assert certification[
        "supports_independent_chemical_branch_validation"
    ] is False
    assert certification["supports_target_reactor_state_prediction"] is False
    assert certification["supports_wafer_flux"] is False
    assert certification["supports_feature_depth"] is False
    assert certification["feature_depth_used"] is False
    assert audit["comparison_quantity"][
        "bulk_flux_distinction_resolved"
    ] is False


def test_nist_backbone_improves_the_declared_transport_band_and_is_converged():
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    models = audit["models"]
    working = models["kushner_zhang_working_set"]["comparison"]
    evaluated = models[
        "nist_evaluated_constant_join_ratio"
    ]["comparison"]
    key = "declared_reactor_development_band"
    assert evaluated[key]["maximum_absolute_relative_residual"] < 0.10
    assert (
        evaluated[key]["median_absolute_relative_residual"]
        < working[key]["median_absolute_relative_residual"]
    )
    assert audit["numerical_convergence"][
        "maximum_absolute_flux_drift_relative_change"
    ] < 0.002


def test_high_energy_transport_closure_is_irrelevant_in_40_to_100_td_band():
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    models = audit["models"]
    constant = models[
        "nist_evaluated_constant_join_ratio"
    ]["comparison"]["declared_reactor_development_band"]
    tapered = models[
        "nist_evaluated_linear_return_to_working_set_at_120eV"
    ]["comparison"]["declared_reactor_development_band"]
    assert abs(
        constant["maximum_absolute_relative_residual"]
        - tapered["maximum_absolute_relative_residual"]
    ) < 1.0e-9
