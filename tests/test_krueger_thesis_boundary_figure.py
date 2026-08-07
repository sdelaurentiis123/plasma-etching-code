from scripts.audit_krueger_thesis_boundary_figure import build_audit


def test_figure_6_17_does_not_promote_missing_curves_to_boundary_data():
    audit = build_audit()
    transcription = audit["visual_transcription"]
    assert transcription["curve_labels"] == [
        "CF2", "C3F4", "O", "C2F3", "CF", "Ions", "CF3", "CO", "C",
    ]
    assert transcription["stable_c4f6_parent_curve_present"] is False
    assert transcription["positive_ion_species_resolved"] is False
    assert transcription["positive_ion_label"] == "Ions"

    boundary = audit["claim_boundary"]
    assert "stable C4F6 wafer flux is zero" in boundary["does_not_prove"]
    assert "independently predictive" in boundary["depth_consequence"]


def test_figure_6_17_scope_is_not_the_base_boundary():
    audit = build_audit()
    assert audit["figure_condition"]["low_frequency_power_kW"] == 6.0
    assert audit["figure_condition"]["high_frequency_power_kW"] == 2.5
    assert "not the 8 kW base-case" in audit["figure_condition"]["scope"]
    assert audit["method"]["numerical_curve_digitization_performed"] is False
