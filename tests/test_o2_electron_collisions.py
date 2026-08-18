import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from petch.reactor_global.electron_collision_deck import (
    ElectronCollisionDeck,
    ElectronCollisionProcess,
)
from petch.reactor_global.electron_collision_mixture import (
    compose_electron_collision_decks,
)
from petch.reactor_global.o2_electron_collisions import (
    O2_MASS_AMU,
    SONG_2026_O2_WORKBOOK_SHA256,
    load_song_2026_o2_replay,
)


ROOT = Path(__file__).resolve().parents[1]
LOCAL_SOURCE_CANDIDATES = (
    Path(os.environ.get("PETCH_SONG_2026_O2_WORKBOOK", "")),
    Path("/private/tmp/o2_song_2026_supplement.xlsx"),
)


def _local_source() -> Path:
    for path in LOCAL_SOURCE_CANDIDATES:
        if str(path) and path.is_file():
            return path
    pytest.skip("licensed Song 2026 O2 source workbook not supplied")


def _manufactured_deck(target: str, payload: str) -> ElectronCollisionDeck:
    return ElectronCollisionDeck(
        processes=(ElectronCollisionProcess(
            kind="MOMENTUM",
            target=target,
            product=target,
            electron_energy_eV=(0.0, 1.0),
            cross_section_m2=(1.0e-20, 1.0e-20),
            mass_ratio=1.0e-5,
        ),),
        payload_sha256=payload,
        source_database="manufactured",
        retrieved_at="2026-08-18",
        source_reference="manufactured unit test",
    )


def test_exact_external_workbook_replays_recommended_o2_curves_and_units():
    replay = load_song_2026_o2_replay(_local_source())
    assert replay.source_workbook_sha256 == SONG_2026_O2_WORKBOOK_SHA256
    assert [curve.label for curve in replay.source_curves] == [
        "momentum_transfer",
        "vibrational_v0_to_v1",
        "electronic_a1Delta_g",
        "electronic_b1Sigma_g_plus",
        "electronic_c1Sigma_u_minus",
        "electronic_A3Delta_u",
        "electronic_A3Sigma_u_plus",
        "neutral_dissociation_total",
        "positive_ionization_total",
        "dissociative_attachment_total",
    ]
    momentum = replay.curve("momentum_transfer")
    assert momentum.electron_energy_eV[[0, -1]].tolist() == [.001, 1000.0]
    np.testing.assert_allclose(
        momentum.cross_section_m2[[0, -1]], [.35e-20, .1e-20],
        rtol=0.0, atol=1.0e-36,
    )
    vibration = replay.curve("vibrational_v0_to_v1")
    assert vibration.electron_energy_eV.size == 1800
    np.testing.assert_allclose(vibration.cross_section_m2.max(), 6.5413e-20)
    attachment = replay.curve("dissociative_attachment_total")
    np.testing.assert_allclose(attachment.cross_section_m2.max(), .0141e-20)


def test_o2_kinetic_deck_keeps_charge_sign_and_use_boundary_explicit():
    replay = load_song_2026_o2_replay(_local_source())
    assert [process.kind for process in replay.derived_deck.processes] == [
        "MOMENTUM",
        "EXCITATION", "EXCITATION", "EXCITATION", "EXCITATION",
        "EXCITATION", "EXCITATION", "EXCITATION",
        "IONIZATION", "ATTACHMENT",
    ]
    momentum = replay.derived_deck.processes[0]
    assert momentum.mass_ratio == 5.48579909065e-4 / O2_MASS_AMU
    ionization = replay.derived_deck.processes[-2]
    assert ionization.electron_number_change == 1
    assert ionization.product == "O2+/O+ aggregate"
    assert replay.source_artifact_committed is False
    assert replay.supports_resolved_primary_chemistry is True
    assert replay.supports_direct_swarm_validation is False
    assert replay.supports_target_reactor_state_prediction is False
    assert replay.supports_feature_depth is False


def test_o2_unmeasured_closures_change_the_derived_hash():
    path = _local_source()
    nominal = load_song_2026_o2_replay(path)
    tail = load_song_2026_o2_replay(
        path, high_energy_tail_closure="linear_to_zero_at_30eV")
    onset = load_song_2026_o2_replay(
        path,
        dissociation_onset_closure="linear_from_physical_threshold",
    )
    assert nominal.derived_deck.payload_sha256 != tail.derived_deck.payload_sha256
    assert nominal.derived_deck.payload_sha256 != onset.derived_deck.payload_sha256


def test_external_o2_ingestion_fails_closed_on_wrong_bytes(tmp_path):
    source = tmp_path / "wrong.xlsx"
    source.write_bytes(b"not the licensed source workbook")
    with pytest.raises(RuntimeError, match="checksum changed"):
        load_song_2026_o2_replay(source)


def test_collision_deck_composition_preserves_processes_and_rejects_overlap():
    first = _manufactured_deck("A", "1" * 64)
    second = _manufactured_deck("B", "2" * 64)
    mixed = compose_electron_collision_decks(
        (first, second), retrieved_at="2026-08-18", mixture_name="A/B")
    assert mixed.targets == ("A", "B")
    assert mixed.processes == first.processes + second.processes
    assert "mole_fractions_embedded\": false" in mixed.source_reference
    with pytest.raises(ValueError, match="appears in multiple decks"):
        compose_electron_collision_decks(
            (first, first), retrieved_at="2026-08-18", mixture_name="bad")


def test_committed_npg80_feed_electron_audit_replays_without_licensed_source():
    subprocess.run([
        sys.executable,
        str(ROOT / "scripts" / "audit_zhu_npg80_feed_electron_kinetics.py"),
        "--check",
    ], check=True, cwd=ROOT)
