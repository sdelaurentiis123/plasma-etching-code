#!/usr/bin/env python3
"""Audit what the Levinson 1997 article can identify without target fitting."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT / "results" / "curated"
    / "levinson_1997_feature_identifiability" / "audit.json"
)


def build_audit() -> dict:
    """Return the claim-scoped controlled-beam feature ledger."""
    cases = [
        {
            "figure_panel": "11a",
            "trench_width_um": 2.67,
            "measured_depth_um": 0.38,
            "unshadowed_ion_to_neutral_flux_ratio": 0.004,
            "ion_energy_eV": 100.0,
            "initial_mask_height_um": 0.97,
            "mask_removed_in_image": False,
        },
        {
            "figure_panel": "11b",
            "trench_width_um": 5.0,
            "measured_depth_um": 1.2,
            "unshadowed_ion_to_neutral_flux_ratio": 0.008,
            "ion_energy_eV": 100.0,
            "initial_mask_height_um": 0.8,
            "mask_removed_in_image": True,
        },
        {
            "figure_panel": "11c",
            "trench_width_um": 1.9,
            "measured_depth_um": 1.18,
            "unshadowed_ion_to_neutral_flux_ratio": 0.008,
            "ion_energy_eV": 100.0,
            "initial_mask_height_um": 0.75,
            "mask_removed_in_image": True,
        },
    ]
    for case in cases:
        case["measured_depth_to_width_ratio"] = (
            case["measured_depth_um"] / case["trench_width_um"]
        )
        case["initial_mask_aspect_ratio"] = (
            case["initial_mask_height_um"] / case["trench_width_um"]
        )

    return {
        "schema": "petch.levinson-1997-feature-identifiability.v1",
        "audit_id": "LEVINSON-1997-CL2-SI-FEATURE-2026-08-06-R1",
        "source": {
            "citation": (
                "J. A. Levinson, E. S. G. Shaqfeh, M. Balooch, "
                "and A. V. Hamza, JVST A 15, 1902-1912 (1997)"
            ),
            "doi": "10.1116/1.580658",
            "full_text_route": (
                "author-uploaded ResearchGate full-text transcription"
            ),
            "source_pixels_archived": False,
            "source_pdf_archived": False,
        },
        "apparatus": {
            "ion_species": "Ar+",
            "neutral_species": "Cl2",
            "ion_direction": "normal to wafer",
            "neutral_distribution": "isotropic effusive background",
            "source_energy_range_eV": [10.0, 1200.0],
            "source_design_current_density_mA_cm2_at_100eV": 0.45,
            "fraction_within_plus_minus_10eV": "at least 0.70",
            "energy_distribution_shape": "near Gaussian",
            "sample_position_current_measured": True,
            "sample_position_current_values_published_for_figure11": False,
        },
        "surface_closure": {
            "coverage_equation": (
                "Q = 1 / (1 + A*Ychem*GI/(2*S0*GN))"
            ),
            "yield_equation": "Ytot = Ysput*(1-Q) + Ychem*Q",
            "initial_sticking_probability_S0": 0.75,
            "average_product_chlorine_stoichiometry_A": 2.0,
            "parameter_class": (
                "S0 and A regressed to planar beam-yield data; not "
                "first-principles constants"
            ),
            "independent_of_feature_profiles": True,
        },
        "figure11_cases": cases,
        "absolute_depth_identifiability": {
            "case_specific_ion_flux_published": False,
            "case_specific_exposure_time_published": False,
            "case_specific_ion_fluence_published": False,
            "typical_saturated_rate_A_per_min": [15.0, 20.0],
            "apparatus_design_current_is_case_boundary": False,
            "target_depth_may_select_simulation_time": False,
            "article_identifies_absolute_feature_depth_prediction": False,
            "missing_minimum": (
                "case-specific ion fluence, or measured sample-position "
                "ion current plus exposure time"
            ),
        },
        "allowed_now": {
            "surface_yield_curve_replay": True,
            "dimensionless_neutral_transport_replay": True,
            "normalized_arde_trend_test": True,
            "profile_shape_test_after_original_pixel_archive": True,
            "absolute_feature_depth_test": False,
        },
        "observed_model_discrepancies": {
            "microtrenching": True,
            "sidewall_slope": True,
            "source_simulation_reproduced_them": False,
            "source_candidates": [
                "ion reflection from sloped sidewalls",
                "sputtered oxide-mask redeposition",
                "oxide-mask charging and ion focusing",
            ],
        },
        "verdict": (
            "Excellent boundary-identified Cl2/Si surface and neutral-transport "
            "board; the 1997 article alone does not identify absolute feature "
            "depth because Figure 11 omits exposure time and case-specific "
            "ion fluence. Do not fit time to the reported depth."
        ),
    }


def audit_text() -> str:
    return json.dumps(build_audit(), indent=2) + "\n"


def main() -> None:
    payload = audit_text()
    if OUTPUT.read_text(encoding="utf-8") != payload:
        raise RuntimeError("committed Levinson feature audit is stale")
    print(payload, end="")


if __name__ == "__main__":
    main()
