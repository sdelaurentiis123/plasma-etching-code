import subprocess
import sys
from pathlib import Path

import numpy as np

from petch.reactor_global.sf6_electron_collisions import (
    SF6_MASS_AMU,
    derive_nist_evaluated_sf6_replay,
    derive_nist_product_resolved_sf6_replay,
    load_nist_2000_sf6_attachment_rate_curve,
    load_nist_2000_sf6_drift_curve,
    load_nist_2000_sf6_effective_ionization_curve,
    load_nist_2000_sf6_momentum_transfer,
    load_nist_2000_sf6_partial_attachment_curves,
    load_nist_2000_sf6_partial_ionization_anchors,
    load_nist_2000_sf6_total_attachment,
    load_nist_2000_sf6_total_scattering,
)


ROOT = Path(__file__).resolve().parents[1]


def test_committed_sf6_source_tables_replay_without_source_pdf():
    subprocess.run([
        sys.executable,
        str(ROOT / "scripts" / "extract_sf6_electron_evidence.py"),
        "--check",
    ], check=True, cwd=ROOT)


def test_exact_nist_sf6_rows_keep_units_and_evidence_classes_distinct():
    total = load_nist_2000_sf6_total_scattering()
    momentum = load_nist_2000_sf6_momentum_transfer()
    attachment = load_nist_2000_sf6_total_attachment()
    drift = load_nist_2000_sf6_drift_curve()
    effective = load_nist_2000_sf6_effective_ionization_curve()
    rate = load_nist_2000_sf6_attachment_rate_curve()
    assert total.electron_energy_eV[[0, -1]].tolist() == [.035, 4000.0]
    assert total.cross_section_m2[[0, -1]].tolist() == [
        379.8e-20, 2.64e-20,
    ]
    assert momentum.electron_energy_eV[[0, -1]].tolist() == [2.75, 700.0]
    assert momentum.cross_section_m2[[0, -1]].tolist() == [16e-20, .66e-20]
    assert attachment.cross_section_m2[[0, -1]].tolist() == [
        7617e-20, .003e-20,
    ]
    assert drift.recommended_mask.sum() == 16
    assert drift.drift_velocity_m_s[
        np.where(drift.reduced_electric_field_Td == 275.0)[0][0]
    ] == 170000.0
    np.testing.assert_allclose(
        effective.effective_ionization_coefficient_m2[
            np.where(effective.reduced_electric_field_Td == 350.0)[0][0]
        ],
        -2.43e-22,
        rtol=0.0,
        atol=1.0e-36,
    )
    np.testing.assert_allclose(
        rate.attachment_rate_coefficient_m3_s[
            np.where(rate.reduced_electric_field_Td == 100.0)[0][0]
        ],
        887e-18,
        rtol=0.0,
        atol=1.0e-30,
    )


def test_evaluated_deck_has_one_momentum_and_aggregate_balanced_branches():
    replay = derive_nist_evaluated_sf6_replay()
    kinds = [process.kind for process in replay.derived_deck.processes]
    assert kinds == [
        "MOMENTUM", "EXCITATION", "EXCITATION", "IONIZATION", "ATTACHMENT",
    ]
    momentum = replay.derived_deck.processes[0]
    assert momentum.mass_ratio == 5.48579909065e-4 / SF6_MASS_AMU
    assert [
        process.electron_number_change
        for process in replay.derived_deck.processes
    ] == [0, 0, 0, 1, -1]
    assert all(
        process.electron_energy_eV[-1] >= replay.maximum_energy_eV
        for process in replay.derived_deck.processes
    )
    assert replay.supports_resolved_primary_chemistry is False
    assert replay.supports_feature_depth is False


def test_total_scattering_is_deconvolved_and_reconstructed_exactly_once():
    replay = derive_nist_evaluated_sf6_replay()
    total = load_nist_2000_sf6_total_scattering()
    probe = total.electron_energy_eV[total.electron_energy_eV <= 2.5]
    assembled = np.zeros(probe.shape)
    for process in replay.derived_deck.processes:
        assembled += np.interp(
            probe,
            process.electron_energy_eV,
            process.cross_section_m2,
        )
    np.testing.assert_allclose(
        assembled,
        total.cross_section_m2[:probe.size],
        rtol=2.0e-15,
        atol=1.0e-30,
    )


def test_vibrational_loss_is_an_explicit_small_transport_sensitivity():
    nominal = derive_nist_evaluated_sf6_replay()
    alternate = derive_nist_evaluated_sf6_replay(
        vibrational_energy_loss_eV=.117
    )
    assert nominal.derived_deck.payload_sha256 != alternate.derived_deck.payload_sha256
    assert nominal.derived_deck.processes[1].energy_loss_eV == .095
    assert alternate.derived_deck.processes[1].energy_loss_eV == .117


def test_product_tables_replay_all_reported_charged_fragments():
    positive = load_nist_2000_sf6_partial_ionization_anchors()
    negative = load_nist_2000_sf6_partial_attachment_curves()
    assert [curve.product for curve in positive] == [
        "SF5+", "SF4+", "SF3+", "SF2+", "SF+", "S+", "F+",
        "SF4++", "SF2++",
    ]
    assert [curve.product for curve in negative] == [
        "SF6-", "SF5-", "SF4-", "SF3-", "SF2-", "F2-", "F-",
    ]
    assert all(curve.electron_energy_eV.tolist() == [100.0] for curve in positive)
    np.testing.assert_allclose(
        sum(curve.cross_section_m2[0] for curve in positive),
        6.51e-20,
        rtol=0.0,
        atol=1.0e-34,
    )
    # Blank Table 27 cells are absent source rows, never manufactured zeros.
    assert negative[2].electron_energy_eV[[0, -1]].tolist() == [3.5, 8.5]
    assert negative[-1].electron_energy_eV[[0, -1]].tolist() == [2.0, 15.0]


def test_resolved_attachment_reconstructs_evaluated_total_exactly():
    replay = derive_nist_product_resolved_sf6_replay()
    aggregate = next(
        process for process in replay.aggregate_replay.derived_deck.processes
        if process.kind == "ATTACHMENT"
    )
    branches = tuple(
        process for process in replay.derived_deck.processes
        if process.kind == "ATTACHMENT"
    )
    assert [process.product for process in branches] == [
        "SF5-", "SF6-", "SF4-", "SF3-", "SF2-", "F2-", "F-",
    ]
    probe = np.unique(np.concatenate((
        aggregate.electron_energy_eV,
        *(process.electron_energy_eV for process in branches),
    )))
    expected = np.interp(
        probe, aggregate.electron_energy_eV, aggregate.cross_section_m2)
    actual = sum(
        np.interp(probe, process.electron_energy_eV, process.cross_section_m2)
        for process in branches
    )
    np.testing.assert_allclose(actual, expected, rtol=2e-15, atol=1e-35)
    assert replay.maximum_attachment_rounding_rescale_fraction < .012
    assert replay.maximum_sf5_source_peak_normalized_residual < .012
    by_product = {process.product: process for process in branches}
    dominant = lambda product, energy: np.interp(
        energy,
        by_product[product].electron_energy_eV,
        by_product[product].cross_section_m2,
    )
    assert dominant("SF6-", .001) > dominant("SF5-", .001)
    assert dominant("SF5-", .3) > dominant("SF6-", .3)
    assert dominant("F-", 5.0) > dominant("SF5-", 5.0)


def test_resolved_ionization_preserves_total_and_measured_100eV_fractions():
    replay = derive_nist_product_resolved_sf6_replay()
    aggregate = next(
        process for process in replay.aggregate_replay.derived_deck.processes
        if process.kind == "IONIZATION"
    )
    branches = tuple(
        process for process in replay.derived_deck.processes
        if process.kind == "IONIZATION"
    )
    probe = np.unique(np.concatenate((
        aggregate.electron_energy_eV,
        *(process.electron_energy_eV for process in branches),
    )))
    expected = np.interp(
        probe, aggregate.electron_energy_eV, aggregate.cross_section_m2)
    actual = sum(
        np.interp(probe, process.electron_energy_eV, process.cross_section_m2)
        for process in branches
    )
    np.testing.assert_allclose(actual, expected, rtol=2e-15, atol=1e-35)
    source = {
        curve.product: curve.cross_section_m2[0]
        for curve in replay.ionization_source_anchors
    }
    for process in branches:
        at_100 = np.interp(
            100.0, process.electron_energy_eV, process.cross_section_m2)
        np.testing.assert_allclose(
            at_100 / replay.evaluated_total_ionization_at_100eV_m2,
            source[process.product] / replay.ionization_anchor_sum_m2,
            rtol=2e-15,
            atol=1e-16,
        )
    assert {
        process.product: process.electron_number_change for process in branches
    }["SF4++"] == 2
    assert replay.supports_direct_positive_ion_curves is False
    assert replay.supports_feature_depth is False
