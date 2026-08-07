import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
AUDIT_PATH = (
    ROOT / "results" / "curated" / "guo_krueger_surface_transfer"
    / "audit.json"
)


def _load_auditor():
    path = ROOT / "scripts" / "audit_guo_krueger_surface_transfer.py"
    spec = importlib.util.spec_from_file_location(
        "guo_krueger_surface_transfer_audit", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_frozen_guo_krueger_surface_audit_rebuilds_without_target_fit():
    frozen = json.loads(AUDIT_PATH.read_text())
    rebuilt = _load_auditor().build_audit(ROOT)

    assert frozen == rebuilt
    assert rebuilt["status"] == "passed"
    firewall = rebuilt["calibration_firewall"]
    assert firewall["feature_depth_used_by_surface_solver"] is False
    assert firewall["surface_or_boundary_parameters_adjusted"] == []
    assert firewall["optimization_performed"] is False
    assert firewall[
        "target_used_only_after_all_surface_solves_for_scoring"] is True


def test_nominal_no_fit_surface_match_is_real_but_not_a_depth_prediction():
    audit = json.loads(AUDIT_PATH.read_text())

    nominal = audit["nominal_no_fit_surface_result"]
    assert nominal["sio2_yield_per_wafer_ion"] == pytest.approx(
        2.613053833986497)
    assert abs(nominal["atom_ledger_residual"]) < 1.0e-12
    assert nominal["steady_state_residual"] < 2.0e-8
    assert audit["score"]["within_five_percent"] is True
    assert audit["score"]["absolute_percentage_error"] == pytest.approx(
        0.036583339102081336)
    assert audit["linear_depth_diagnostic_only"][
        "surface_yield_ratio_scaled_depth_nm"] == pytest.approx(
            855.1812547592172)
    assert audit["linear_depth_diagnostic_only"][
        "authoritative_feature_prediction"] is False
    assert audit["verdict"]["feature_depth_match_earned"] is False
    assert audit["verdict"]["exact_825nm_prediction_authorized"] is False


def test_energy_and_species_sensitivities_block_an_accidental_promotion():
    audit = json.loads(AUDIT_PATH.read_text())

    support = audit["energy_support"]
    assert support[
        "krueger_iead_probability_within_guo_fit_support"] == 0.0
    independent = support["independent_sqrt_form_support"]
    assert independent[
        "chemical_probability_weight_within_support"] == pytest.approx(
            0.13183047150523108)
    assert independent[
        "physical_probability_weight_within_support"] == pytest.approx(
            0.06569097356735597)
    capped = audit["energy_law_sensitivity"][
        "energy_censored_at_source_maximum_370_eV"]
    assert capped["sio2_yield_per_wafer_ion"] == pytest.approx(
        0.8749408224693355)

    neutral = audit["neutral_mapping_sensitivity"]
    assert neutral["omit_C3F4_outside_source_species_list"][
        "sio2_yield_per_wafer_ion"] == pytest.approx(2.0462375961407697)
    assert neutral["source_species_intersection_only"][
        "sio2_yield_per_wafer_ion"] == pytest.approx(1.7428114061753623)

    ions = audit["unpublished_ion_composition_sensitivity"]
    assert ions[
        "converged_yield_range_including_deposition"] == pytest.approx([
            -0.7199876957278631, 3.66123223577177])
    assert ions["deposition_cases"] == ["all_C3F3_positive_ions"]
    assert ions["nonconverged_cases"] == []


def test_source_typo_and_atomicity_boundaries_remain_explicit():
    audit = json.loads(AUDIT_PATH.read_text())

    angular = audit["angular_source_defect"]
    assert angular["literal_negative_sample_count"] == 5
    assert angular[
        "repair_independently_traced_to_original_fit"] is False
    atomicity = audit["atomicity_statement"]
    assert atomicity["atom_balanced"] is True
    assert atomicity["atomistic_trajectory_or_interatomic_potential"] is False
    assert atomicity["atomic_level_accuracy_claimed"] is False
