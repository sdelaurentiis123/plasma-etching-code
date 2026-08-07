import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = (
    ROOT
    / "results"
    / "curated"
    / "kounis_melas_ale_products"
    / "audit.json"
)


def test_ale_product_audit_is_atom_balanced_and_does_not_overclaim_depth():
    audit = json.loads(AUDIT.read_text())

    assert "not experimental validation" in audit["claim"]
    assert "not" in audit["claim"] and "held-out depth" in audit["claim"]
    assert audit["source"]["ion_energy_eV"] == 215.0
    assert not audit["dose_windows"]["interpolation_used"]
    sequence = audit["integrated_sequence"]
    assert sequence["product_routing_complete"]
    assert sequence["exchange_residual_atoms_m2"] == {
        "Cl_atom": 0.0,
        "Si_atom": 0.0,
    }
    assert sequence["remaining_cl_atoms_m2"] > 0.0
    assert sequence["equivalent_removed_si_depth_nm"] > 0.0
    assert sequence["all_products_lack_launch_distribution"]


def test_ale_product_source_condition_and_csv_values_pass_vision_gate():
    audit = json.loads(AUDIT.read_text())
    vision = audit["vision_condition_and_value_audit"]

    assert vision["status"] == "passed"
    assert "215 eV" in vision["condition_disambiguation"]
    assert "80 eV" in vision["condition_disambiguation"]
    assert vision["paper_sha256"] == (
        "d8d374d412d99e625b2b989399d564976332f0c342b37b0249d1039bca4b5bb1")
    for panel in vision["panels"].values():
        assert panel["caption_ion_energy_eV"] == 215.0
        assert panel["csv_nodes_match_within_1p25_pixels"]
        assert len(panel["markers"]) == 5


def test_independent_chang_comparison_remains_a_diagnostic_not_a_fit():
    audit = json.loads(AUDIT.read_text())
    diagnostic = audit["independent_continuous_rie_diagnostic"]

    assert not diagnostic["fit_to_comparison_used"]
    assert diagnostic["maximum_relative_error"] > 0.20
    assert "not an exact" in diagnostic["verdict"]
