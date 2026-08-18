import subprocess
import sys
from pathlib import Path

import numpy as np

from petch.reactor_global.sf6_electron_collisions import (
    SF6_MASS_AMU,
    derive_nist_evaluated_sf6_replay,
    load_nist_2000_sf6_attachment_rate_curve,
    load_nist_2000_sf6_drift_curve,
    load_nist_2000_sf6_effective_ionization_curve,
    load_nist_2000_sf6_momentum_transfer,
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
