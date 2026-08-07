import pytest

from scripts.audit_krueger_iead_energy_support import build_audit


def test_energy_support_audit_pins_the_uncovered_krueger_tail():
    audit = build_audit()
    domains = {row["domain_id"]: row for row in audit["domains"]}

    assert audit["source"]["node_count"] == 878
    assert audit["source"]["probability_sum"] == pytest.approx(1.0)
    assert audit["source"]["probability_weighted_mean_energy_eV"] == pytest.approx(
        3465.1107871409254
    )
    assert domains["guo_yin_regression"][
        "iead_probability_above_maximum"
    ] == pytest.approx(1.0)
    assert domains["an_2026_released_nnp_outputs"][
        "iead_probability_above_maximum"
    ] == pytest.approx(0.9489846506185552)
    assert domains["karahashi_mass_selected_beam"][
        "iead_probability_above_maximum"
    ] == pytest.approx(0.8681695284947692)
    assert domains["tachi_1982_si_target_lead"][
        "iead_probability_above_maximum"
    ] == pytest.approx(0.7285794839278725)


def test_energy_overlap_never_promotes_missing_species_or_chemistry():
    audit = build_audit()

    assert not audit["atomic_accuracy_verdict"]["granted"]
    assert all(
        not row["grants_krueger_prediction_authority"]
        for row in audit["domains"]
    )
    assert all(not row["species_match"] for row in audit["domains"])
    assert all(not row["chemistry_match"] for row in audit["domains"])
    tachi = next(
        row
        for row in audit["domains"]
        if row["domain_id"] == "tachi_1982_si_target_lead"
    )
    assert not tachi["target_match"]
    assert tachi["authority"] == "target_mismatched_lead_not_surface_support"
