import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "compare_krueger_2024_radiosity_ray_refinement",
    ROOT / "scripts" / "compare_krueger_2024_radiosity_ray_refinement.py",
)
COMPARE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMPARE)


def _integration(*, ray_shift, oxide):
    state = {
        "polymer_units_m2": [1.0 + ray_shift, 3.0],
        "coverage": [0.25, 0.75],
    }
    exchange = {
        "removed": {
            "SiO2_formula_unit": [2.0 + ray_shift, 4.0],
            "polymer_unit": [0.0, 0.5],
        },
        "outgoing": {"SiO2_formula_unit": [2.0 + ray_shift, 4.0]},
        "unresolved": {"SiO2_formula_unit": [0.0, 0.0]},
        "deposited": {"polymer_unit": [0.0, 0.5]},
    }
    return {
        "final_state_fields": state,
        "final_state_fields_sha256": f"state-{ray_shift}-{oxide}",
        "per_face_integrated_exchange_units_m2": exchange,
        "per_face_integrated_exchange_sha256": f"exchange-{ray_shift}-{oxide}",
        "integrated_exchange": {
            "removed": {
                "SiO2_formula_unit": oxide,
                "polymer_unit": 0.5,
            },
            "outgoing": {"SiO2_formula_unit": oxide},
            "unresolved": {"SiO2_formula_unit": 0.0},
            "deposited": {"polymer_unit": 0.5},
        },
        "oxide_removal": {
            "integrated_formula_units": oxide,
            "integrated_volume_m3": oxide * 1.0e-28,
            "effective_mean_normal_velocity_m_s": oxide * 1.0e-12,
        },
        "displacement": {
            "maximum_gross_displacement_m": 2.0e-9,
            "maximum_gross_displacement_dx": 0.02,
            "per_face_integrated_recession_m": [1.0e-9 + ray_shift * 1.0e-10, 2.0e-9],
            "per_face_integrated_growth_m": [0.0, 0.5e-9],
        },
    }


def _audit(rays, *, ray_shift=0.0, face_area=True, checkpoint="checkpoint"):
    r17_nominal = _integration(ray_shift=ray_shift, oxide=10.0 + ray_shift)
    r17_tight = _integration(ray_shift=ray_shift, oxide=10.5 + ray_shift)
    r19_nominal = _integration(ray_shift=ray_shift, oxide=9.0 + 0.5 * ray_shift)
    r19_tight = _integration(ray_shift=ray_shift, oxide=9.5 + 0.5 * ray_shift)
    difference = (
        r19_tight["oxide_removal"]["integrated_formula_units"]
        - r17_tight["oxide_removal"]["integrated_formula_units"])
    payload = {
        "schema": COMPARE.EXPECTED_AUDIT_SCHEMA,
        "status": "pass",
        "scientific_scope": "manufactured frozen checkpoint",
        "data_firewall": {
            "boundary_case": "base",
            "held_out_observations_loaded": False,
        },
        "checkpoint": {
            "audit_sha256": "audit-source",
            "checkpoint_sha256": checkpoint,
            "metadata": {"physical_time_s": 60.0},
        },
        "parameter_provenance": {
            "r17_audit_sha256": "r17",
            "r19_audit_sha256": "r19",
            "parameter_pairs": {"r17": {"a": 1.0}, "r19": {"a": 2.0}},
        },
        "transport_operator": {
            "transport_seed": 241,
            "neutral_radiosity_seed": 10241,
        },
        "radiosity_operator": {
            "rays_per_face": rays,
            "maximum_rays_per_face": 8 * rays,
            "seed_offset": 10000,
            "relative_tolerance": 1.0e-12,
        },
        "direct_transport": {
            "direct_surface_flux_sha256": "same-direct-flux",
            "boundary_provenance": {"case": "base"},
        },
        "form_factors": {
            "face_count": 2,
            "requested_rays_per_face": rays,
            "rays_per_face": rays,
        },
        "execution_budget": {"next_profile_step_s": 0.125},
        "provenance": {
            "source": {"files": {"engine.py": "engine-hash"}},
            "base_inputs": {"files": {"input.json": "input-hash"}},
            "runtime_selection": {"requested_rays_per_face": rays},
        },
        "horizons": [{
            "fraction_of_next_profile_step": 0.0625,
            "horizon_s": 0.0078125,
            "common_pass": True,
            "parameter_results": {
                "r17": {
                    "all_gates_pass": True,
                    "nominal": r17_nominal,
                    "tight": r17_tight,
                },
                "r19": {
                    "all_gates_pass": True,
                    "nominal": r19_nominal,
                    "tight": r19_tight,
                },
            },
            "paired_oxide_removal_direction": {
                "r19_minus_r17_integrated_formula_units": difference,
                "r19_to_r17_ratio": (
                    r19_tight["oxide_removal"]["integrated_formula_units"]
                    / r17_tight["oxide_removal"]["integrated_formula_units"]),
                "direction": "r19_lower",
            },
        }],
    }
    if face_area:
        payload["face_area_m2"] = [1.0, 9.0]
    return payload


def test_manufactured_receipt_compares_every_collection_and_uses_face_area():
    lower = _audit(8, ray_shift=0.0)
    higher = _audit(16, ray_shift=1.0)

    receipt = COMPARE.compare_audits(lower, higher)

    assert receipt["status"] == "complete"
    assert receipt["compatibility"]["nested_extension_factor"] == 2
    assert receipt["decision_contract"]["post_hoc_pass_tolerance_applied"] is False
    assert receipt["array_norm_contract"]["relative_l1"]["kind"] == (
        "physical_area_weighted_relative_l1")

    comparison = receipt["comparisons"]["r17"]["nominal"]
    assert comparison["final_state_fields"]["item_count"] == 2
    assert comparison["per_face_integrated_exchange_units_m2"]["item_count"] == 5
    assert comparison["per_face_displacement"]["item_count"] == 2
    assert comparison["integrated_exchange"]["item_count"] == 5
    assert comparison["oxide_removal"]["item_count"] == 3

    polymer = comparison["final_state_fields"]["by_item"]["polymer_units_m2"]
    assert polymer["symmetric_relative_linf_error"] == pytest.approx(1.0 / 3.0)
    assert polymer["symmetric_relative_l1_error"] == pytest.approx(1.0 / 29.0)
    assert polymer["l1_norm"] == "physical_area_weighted_relative_l1"
    assert receipt["paired_r19_minus_r17_oxide_removal"]["nominal"][
        "direction_preserved"]
    assert receipt["paired_r19_minus_r17_oxide_removal"]["tight"][
        "lower_ray"]["direction"] == "r19_lower"


def test_receipt_honestly_labels_normalized_l1_without_face_area():
    lower = _audit(8, ray_shift=0.0, face_area=False)
    higher = _audit(16, ray_shift=1.0, face_area=False)

    receipt = COMPARE.compare_audits(lower, higher)
    polymer = receipt["comparisons"]["r17"]["nominal"][
        "final_state_fields"]["by_item"]["polymer_units_m2"]

    assert receipt["array_norm_contract"]["relative_l1"]["kind"] == (
        "normalized_unweighted_relative_l1")
    assert polymer["l1_norm"] == "normalized_unweighted_relative_l1"
    assert polymer["symmetric_relative_l1_error"] == pytest.approx(1.0 / 5.0)


@pytest.mark.parametrize(
    "mutate, match",
    [
        (lambda audit: audit["checkpoint"].update(checkpoint_sha256="other"),
         "checkpoint_sha256"),
        (lambda audit: audit["direct_transport"].update(
            direct_surface_flux_sha256="other"), "direct_surface_flux_sha256"),
        (lambda audit: audit["form_factors"].update(
            rays_per_face=8, requested_rays_per_face=8), "same actual ray count"),
        (lambda audit: audit["horizons"][0].update(horizon_s=0.01),
         "different physical durations"),
    ],
)
def test_refuses_incompatible_provenance_and_horizon(mutate, match):
    lower = _audit(8)
    higher = _audit(16, ray_shift=0.1)
    mutate(higher)
    if "same actual ray count" in match:
        higher["provenance"]["runtime_selection"]["requested_rays_per_face"] = 8
        higher["radiosity_operator"]["rays_per_face"] = 8
        higher["radiosity_operator"]["maximum_rays_per_face"] = 64

    with pytest.raises(ValueError, match=match):
        COMPARE.compare_audits(lower, higher)


def test_run_writes_source_hashes_and_selects_explicit_horizon(tmp_path):
    lower_path = tmp_path / "rays8" / "audit.json"
    higher_path = tmp_path / "rays16" / "audit.json"
    lower_path.parent.mkdir()
    higher_path.parent.mkdir()
    lower_path.write_text(json.dumps(_audit(8)), encoding="utf-8")
    higher_path.write_text(json.dumps(_audit(16, ray_shift=0.25)), encoding="utf-8")
    output = tmp_path / "comparison" / "audit.json"

    args = COMPARE.parse_args([
        "--audit-a", str(higher_path),
        "--audit-b", str(lower_path),
        "--horizon-fraction", "0.0625",
        "--output", str(output),
    ])
    receipt = COMPARE.run(args)
    loaded = json.loads(output.read_text(encoding="utf-8"))

    assert loaded == receipt
    assert loaded["selected_horizon"]["selection_rule"] == (
        "explicit_horizon_fraction")
    assert loaded["sources"]["lower_ray"]["input_sha256"] == (
        COMPARE._sha256(lower_path))
    assert loaded["sources"]["higher_ray"]["input_sha256"] == (
        COMPARE._sha256(higher_path))
    assert loaded["sources"]["lower_ray"]["actual_rays_per_face"] == 8
    assert loaded["sources"]["higher_ray"]["actual_rays_per_face"] == 16
