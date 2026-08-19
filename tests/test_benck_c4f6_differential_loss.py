import json

from scripts.audit_benck_c4f6_differential_loss import DEFAULT_OUTPUT, audit


def test_differential_loss_repair_is_nonnegative_and_target_free():
    payload = audit()
    diagnostic = payload["diagnostics"]

    assert not any(payload["calibration_firewall"].values())
    assert diagnostic["condition_count"] == 4
    assert diagnostic["temperature_count_per_condition"] == 5
    assert diagnostic["all_Ar_mixture_rows_have_nonnegative_sources"] is True
    assert (
        diagnostic["all_Ar_mixture_rows_admit_positive_CF3plus_selective_loss"]
        is True
    )
    assert 0.10 < diagnostic[
        "minimum_required_Ar_mixture_CF3plus_selective_loss_over_common_loss"
    ] < 0.12
    assert 0.60 < diagnostic[
        "maximum_required_Ar_mixture_CF3plus_selective_loss_over_common_loss"
    ] < 0.62
    assert diagnostic["maximum_CF2_ratio_replay_absolute_error"] < 2.0e-16


def test_pure_feed_rejects_one_fixed_differential_loss():
    payload = audit()
    diagnostic = payload["diagnostics"]
    decision = payload["physics_decision"]

    assert diagnostic["pure_C4F6_rows_requiring_positive_loss"] == 1
    assert (
        diagnostic[
            "pure_C4F6_rows_requiring_missing_CF3plus_source_or_changed_branching"
        ]
        == 4
    )
    assert decision["negative_CF3_density_is_required"] is False
    assert decision[
        "Ar_mixture_model_form_can_be_repaired_by_positive_differential_loss"
    ] is True
    assert decision["single_condition_independent_CF3plus_loss_is_certified"] is False
    assert payload["certification"]["supports_krueger_boundary"] is False
    assert payload["certification"]["supports_feature_depth"] is False


def test_committed_differential_loss_board_is_exact_replay():
    committed = json.loads(
        (DEFAULT_OUTPUT / "audit.json").read_text(encoding="utf-8")
    )
    assert committed == audit()
