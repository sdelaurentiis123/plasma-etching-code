"""Evidence arithmetic for the Krueger 2024 absolute-depth question.

This module does not predict an etch profile.  It separates three quantities
that were previously conflated:

* the run-average removal normalized by the published wafer-plane ion flux;
* the instantaneous yield that would be required at one assumed final-geometry
  floor-delivery fraction; and
* measured surface-physics examples whose reactor boundary is not Krueger's.

The result is an identifiability audit: it can invalidate an impossibility
proof without manufacturing the missing species-resolved boundary.
"""

from __future__ import annotations

import math


ELEMENTARY_MASS_UNIT_KG = 1.66053906660e-27
BOLTZMANN_J_K = 1.380649e-23
TORR_PA = 133.32236842105263


def effective_yield_from_depth(
        depth_nm, duration_s, formula_density_m3, wafer_ion_flux_m2_s):
    """Return removed formula units per wafer-plane incident ion."""
    values = (
        float(depth_nm), float(duration_s), float(formula_density_m3),
        float(wafer_ion_flux_m2_s))
    if (any(not math.isfinite(value) for value in values)
            or values[0] < 0.0
            or values[1] <= 0.0
            or values[2] <= 0.0
            or values[3] <= 0.0):
        raise ValueError("depth-yield inputs must be finite and physical")
    return values[0] * 1e-9 * values[2] / (values[1] * values[3])


def ideal_gas_impingement_flux(pressure_pa, molecule_mass_amu, temperature_k):
    """Maxwellian one-sided molecular impingement flux, particles m^-2 s^-1."""
    pressure = float(pressure_pa)
    mass = float(molecule_mass_amu) * ELEMENTARY_MASS_UNIT_KG
    temperature = float(temperature_k)
    if (not math.isfinite(pressure) or pressure < 0.0
            or not math.isfinite(mass) or mass <= 0.0
            or not math.isfinite(temperature) or temperature <= 0.0):
        raise ValueError("impingement inputs must be finite and physical")
    return pressure / math.sqrt(
        2.0 * math.pi * mass * BOLTZMANN_J_K * temperature)


def parent_fraction_for_flux_ratio(
        desired_floor_parent_to_ion_ratio,
        undissociated_parent_impingement_flux_m2_s,
        wafer_ion_flux_m2_s,
        parent_floor_delivery=1.0,
        ion_floor_delivery=1.0):
    """Fraction of a declared parent-flux reference needed for a ratio.

    The returned quantity is a boundary sensitivity, not a dissociation
    measurement.  It lumps plasma destruction, chamber composition, and
    wafer coupling into one survival/availability fraction.  In the Krueger
    audit the reference uses the inlet flow fraction at the total pressure;
    that is a no-dissociation scale estimate, not a strict pressure bound.
    """
    ratio = float(desired_floor_parent_to_ion_ratio)
    parent_flux = float(undissociated_parent_impingement_flux_m2_s)
    ion_flux = float(wafer_ion_flux_m2_s)
    parent_delivery = float(parent_floor_delivery)
    ion_delivery = float(ion_floor_delivery)
    values = (ratio, parent_flux, ion_flux, parent_delivery, ion_delivery)
    if (any(not math.isfinite(value) for value in values)
            or ratio < 0.0 or parent_flux <= 0.0 or ion_flux <= 0.0
            or parent_delivery <= 0.0 or ion_delivery <= 0.0):
        raise ValueError("flux-ratio inputs must be finite and physical")
    return (
        ratio * ion_flux * ion_delivery / (parent_flux * parent_delivery))


def krueger_2024_depth_identifiability(
        *, simulated_depth_nm, takada_400eV_peak_yield,
        karahashi_pure_ion_peak_yield):
    """Build the deterministic depth-evidence audit.

    `takada_400eV_peak_yield` and `karahashi_pure_ion_peak_yield` come from
    checksum-bound digitizations.  The separate 900 eV Takada value is a
    verbatim text datum and is intentionally labeled as such below.
    """
    formula_density = 2.2e28
    target_depth_nm = 825.0
    duration_s = 60.0
    wafer_ion_flux = 1.2e20
    ion_floor_delivery = 0.70
    parent_floor_delivery = 0.10

    target_effective_yield = effective_yield_from_depth(
        target_depth_nm, duration_s, formula_density, wafer_ion_flux)
    simulated_effective_yield = effective_yield_from_depth(
        simulated_depth_nm, duration_s, formula_density, wafer_ion_flux)

    total_pressure_pa = 10e-3 * TORR_PA
    c4f6_feed_fraction = 140.0 / (140.0 + 100.0 + 105.0)
    c4f6_partial_pressure_reference_pa = (
        total_pressure_pa * c4f6_feed_fraction)
    c4f6_undissociated_flux_reference = ideal_gas_impingement_flux(
        c4f6_partial_pressure_reference_pa, 162.03, 300.0)

    wafer_ratio_reference = c4f6_undissociated_flux_reference / wafer_ion_flux
    floor_ratio_reference = (
        wafer_ratio_reference * parent_floor_delivery / ion_floor_delivery)

    # The open TMRSJ paper reports 2.5; the related JAP article's publisher
    # abstract reports 2.4 for the same nominal 900 eV, ratio-1 condition.
    # Preserve the source discrepancy as a range rather than choosing the
    # value that happens to sit closest to the feature target.
    takada_tmrsj_900eV_ratio1_text_yield = 2.5
    takada_jap_900eV_ratio1_abstract_yield = 2.4
    final_floor_required = target_effective_yield / ion_floor_delivery

    return {
        "audit_id": "KRUEGER-2024-DEPTH-IDENTIFIABILITY-R1",
        "inputs": {
            "target_depth_nm": target_depth_nm,
            "simulated_depth_nm": float(simulated_depth_nm),
            "duration_s": duration_s,
            "sio2_formula_density_m3": formula_density,
            "published_aggregate_wafer_ion_flux_m2_s": wafer_ion_flux,
            "published_positive_ion_species_composition": None,
            "published_stable_c4f6_wafer_flux_m2_s": None,
            "diagnostic_final_geometry_ion_floor_delivery": ion_floor_delivery,
            "diagnostic_final_geometry_neutral_floor_delivery": (
                parent_floor_delivery),
        },
        "run_average_wafer_ion_normalization": {
            "target_sio2_per_wafer_ion": target_effective_yield,
            "simulation_sio2_per_wafer_ion": simulated_effective_yield,
            "gap_sio2_per_wafer_ion": (
                target_effective_yield - simulated_effective_yield),
            "target_to_simulation_ratio": (
                target_effective_yield / simulated_effective_yield),
            "meaning": (
                "depth-integrated lower-bound normalization; it does not "
                "assert that every wafer-plane ion reaches the evolving floor"
            ),
        },
        "final_geometry_diagnostic": {
            "target_sio2_per_delivered_floor_ion_if_delivery_were_constant": (
                final_floor_required),
            "meaning": (
                "counterfactual instantaneous normalization using the final "
                "0.70 delivery for all 60 s; not the evolving-feature history"
            ),
        },
        "direct_surface_evidence": {
            "karahashi_pure_cf3_ion_peak_digitized_sio2_per_ion": float(
                karahashi_pure_ion_peak_yield),
            "karahashi_scope": (
                "CF3+, normal incidence, radical-free, digitized within "
                "250-2000 eV support"
            ),
            "takada_c5f8_ar_400eV_ratio1_digitized_sio2_per_ar_ion": float(
                takada_400eV_peak_yield),
            "takada_c5f8_ar_900eV_ratio1_sio2_per_ar_ion": {
                "tmrsj_open_full_text": (
                    takada_tmrsj_900eV_ratio1_text_yield),
                "jap_publisher_abstract": (
                    takada_jap_900eV_ratio1_abstract_yield),
                "policy": "retain 2.4--2.5 source discrepancy as a range",
            },
            "takada_900eV_to_target_run_average_ratio_range": [
                takada_jap_900eV_ratio1_abstract_yield
                / target_effective_yield,
                takada_tmrsj_900eV_ratio1_text_yield
                / target_effective_yield,
            ],
            "interpretation": (
                "a measured stable-parent/ion regime reaches 95--99% of the "
                "wafer-ion-normalized target, so 1.5 is not a universal "
                "surface ceiling; the 2.4--2.5 source discrepancy is retained "
                "and C5F8 cannot be transplanted to C4F6"
            ),
        },
        "external_c4f6_boundary_evidence": {
            "source_bibkey": "kim-2021-coatings",
            "measured_condition": {
                "feed_sccm": {"C4F6": 40.0, "Ar": 80.0},
                "pressure_mTorr": 20.0,
                "rf_power_W": 300.0,
                "oxygen_present": False,
            },
            "neutral_mass_spectrum_contains_c4f6_parent_signal": True,
            "reported_positive_ion_density_order": [
                "CF+", "C3F3+", "CF3+", "CF2+", "C3F5+",
            ],
            "species_energy_distributions_differ": True,
            "absolute_wafer_flux_calibrated": False,
            "transferable_to_krueger_boundary": False,
            "meaning": (
                "direct qualitative evidence that both omitted boundary "
                "classes exist in a C4F6 CCP; detector count rates and the "
                "different reactor condition cannot normalize Krueger"
            ),
        },
        "c4f6_parent_boundary_sensitivity": {
            "assumptions": {
                "total_pressure_mTorr": 10.0,
                "feed_sccm": {"C4F6": 140.0, "Ar": 100.0, "O2": 105.0},
                "gas_temperature_k": 300.0,
                "c4f6_mass_amu": 162.03,
                "partial_pressure_policy": (
                    "total pressure times inlet flow fraction; "
                    "no-dissociation reference, not a strict bound"
                ),
            },
            "c4f6_partial_pressure_feed_fraction_reference_pa": (
                c4f6_partial_pressure_reference_pa),
            "undissociated_c4f6_impingement_feed_fraction_reference_m2_s": (
                c4f6_undissociated_flux_reference),
            "parent_to_published_ion_ratio_wafer_reference": (
                wafer_ratio_reference),
            "parent_to_ion_ratio_floor_reference_with_diagnostic_delivery": (
                floor_ratio_reference),
            "availability_fraction_of_reference_for_wafer_ratio_0p25": (
                parent_fraction_for_flux_ratio(
                    0.25, c4f6_undissociated_flux_reference, wafer_ion_flux)),
            "availability_fraction_of_reference_for_floor_ratio_0p25": (
                parent_fraction_for_flux_ratio(
                    0.25, c4f6_undissociated_flux_reference, wafer_ion_flux,
                    parent_floor_delivery, ion_floor_delivery)),
            "availability_fraction_of_reference_for_floor_ratio_1": (
                parent_fraction_for_flux_ratio(
                    1.0, c4f6_undissociated_flux_reference, wafer_ion_flux,
                    parent_floor_delivery, ion_floor_delivery)),
            "meaning": (
                "order-of-magnitude feed-fraction reference only; "
                "species-dependent residence time can shift the chamber "
                "fraction, and the source publishes neither parent survival "
                "nor parent wafer flux"
            ),
        },
        "verdict": {
            "former_universal_1p5_ceiling_proof_valid": False,
            "published_inputs_identify_absolute_depth": False,
            "exact_825nm_prediction_authorized": False,
            "measured_missing_mechanism_is_plausible": True,
            "reason": (
                "species-resolved positive-ion composition and stable C4F6 "
                "wafer flux are unpublished, while direct beam evidence shows "
                "that both ion identity and parent-molecule co-incidence can "
                "change the yield by order unity"
            ),
            "required_new_evidence": [
                "species-resolved positive-ion flux and IEAD at the wafer",
                "stable C4F6 molecular flux at the wafer",
                "direct C4F6/ion co-incidence yield and deposition curve",
            ],
        },
    }
