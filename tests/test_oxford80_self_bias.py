from pathlib import Path

import numpy as np
import pytest

from petch.reactor_global.oxford80_self_bias import (
    Oxford80RIECondition,
    build_oxford80_self_bias_transfer,
    load_oxford80_self_bias_evidence,
)
from scripts.audit_zhu_npg80_self_bias import (
    DEFAULT_EVIDENCE,
    build_self_bias_receipt,
)


def _target() -> Oxford80RIECondition:
    return Oxford80RIECondition(
        tool_model="Oxford PlasmaPro NPG80 RIE",
        rf_power_W=150.0,
        pressure_mTorr=30.0,
        gas_flows_sccm={"CHF3": 55.0, "SF6": 5.0, "O2": 1.0, "Ar": 0.0},
        electrode_temperature_C=20.0,
        duration_s=1200.0,
    )


def test_primary_rows_preserve_exact_values_checksums_and_censoring():
    observations = load_oxford80_self_bias_evidence(DEFAULT_EVIDENCE)
    by_id = {item.source_id: item for item in observations}

    assert len(observations) == 5
    start = by_id["harmon-2019-ngp80-chf3-start"]
    end = by_id["harmon-2019-ngp80-chf3-end"]
    ternary = by_id["penaud-2006-plasmalab80-ternary"]
    same_pressure = by_id["zou-2004-plasmalab80-sf6-o2-30mtorr"]

    assert start.bias_relation == "greater_than"
    assert start.bias_lower_V == 300.0
    assert start.bias_upper_V is None
    assert end.bias_relation == "less_than_approx"
    assert end.bias_lower_V is None
    assert end.bias_upper_V == 200.0
    assert ternary.midpoint_V == 276.0
    assert ternary.source_pdf_sha256 == (
        "fceed7b5bd6d2ed43a7fe5fe3b1f1b8eeb8ae338ad3d1144326cbc41d2ba3e46"
    )
    assert same_pressure.bias_lower_V == 360.0
    assert same_pressure.bias_upper_V == 387.0


def test_anchor_is_same_chemistry_and_exact_reduced_drive_match():
    transfer = build_oxford80_self_bias_transfer(
        _target(), load_oxford80_self_bias_evidence(DEFAULT_EVIDENCE))
    by_id = {item.source_id: item for item in transfer.observations}
    anchor = by_id[transfer.matched_chemistry_reduced_drive_source_id]

    assert transfer.target.active_gases == frozenset({"CHF3", "SF6", "O2"})
    assert anchor.active_gases == transfer.target.active_gases
    assert anchor.reduced_drive_W_per_mTorr == 5.0
    assert transfer.target.reduced_drive_W_per_mTorr == 5.0
    assert transfer.matched_chemistry_reduced_drive_anchor_V == 276.0


def test_histories_are_deterministic_sensitivity_not_prediction():
    transfer = build_oxford80_self_bias_transfer(
        _target(), load_oxford80_self_bias_evidence(DEFAULT_EVIDENCE))
    histories = {item.name: item for item in transfer.histories}
    drift = histories["exact-NGP80 conditioning thresholds"]

    assert drift.endpoints_are_censor_thresholds is True
    np.testing.assert_allclose(drift.bias_magnitude_V, [300.0, 200.0])
    assert drift.at(600.0) == 250.0
    with pytest.raises(ValueError, match="outside"):
        drift.at(1200.1)
    assert transfer.printed_reference_window_V == (200.0, 400.0)
    assert transfer.printed_window_is_probability_interval is False
    assert transfer.censored_data_extend_outside_printed_window is True
    assert transfer.supports_unique_target_bias is False
    assert transfer.supports_absolute_depth_prediction is False
    assert all(
        history.measured_on_target_condition is False
        and history.supports_absolute_depth_prediction is False
        for history in transfer.histories
    )


def test_receipt_is_target_free_and_binds_evidence_table_checksum():
    receipt = build_self_bias_receipt(evidence_path=Path(DEFAULT_EVIDENCE))

    assert receipt["sem_target_used"] is False
    assert receipt["measured_target_bias_used"] is False
    assert receipt["evidence_table"]["observation_count"] == 5
    assert receipt["mechanical_anchor_selection"]["anchor_V"] == 276.0
    assert receipt["certification"]["supports_unique_target_bias"] is False
    assert receipt["certification"]["supports_absolute_depth_prediction"] is False
    assert "held-out SEM" in receipt["certification"]["forbidden_use"]
