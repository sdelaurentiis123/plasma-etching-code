from dataclasses import replace

import pytest

from petch.reactor_global.chlorine_siglo_replay import (
    DeclaredLogLinearWallExtrapolation,
    LEGACY_SIGLO_CL2_2013_SHA256,
    _EXPECTED_ROWS,
    _derive_legacy_siglo_comsol_chlorine_replay,
    derive_legacy_siglo_cl2_replay,
)
from petch.reactor_global.chlorine_wall import (
    stafford_2010_conditioned_wall_recombination_provider,
)
from petch.reactor_global.electron_collision_deck import (
    ElectronCollisionDeck,
    ElectronCollisionProcess,
)


def _manufactured_topology_fixture():
    processes = []
    for kind, product, loss in _EXPECTED_ROWS:
        if kind == "ELASTIC":
            process = ElectronCollisionProcess(
                kind=kind,
                target="Cl2",
                product=product,
                electron_energy_eV=(0.02, 100.0),
                cross_section_m2=(2.0e-20, 2.0e-20),
                mass_ratio=7.68e-6,
            )
        elif kind == "ATTACHMENT":
            process = ElectronCollisionProcess(
                kind=kind,
                target="Cl2",
                product=product,
                electron_energy_eV=(0.0, 100.0),
                cross_section_m2=(1.0e-20, 0.0),
            )
        else:
            process = ElectronCollisionProcess(
                kind=kind,
                target="Cl2",
                product=product,
                electron_energy_eV=(0.0, 100.0),
                cross_section_m2=(0.0, 1.0e-20),
                energy_loss_eV=loss,
            )
        processes.append(process)
    return ElectronCollisionDeck(
        processes=tuple(processes),
        payload_sha256=LEGACY_SIGLO_CL2_2013_SHA256,
        source_database="manufactured topology fixture",
        retrieved_at="2026-08-08",
        source_reference="test only; contains no LXCat data",
    )


def test_legacy_replay_declares_padding_multiplicity_and_every_mapping():
    replay = derive_legacy_siglo_cl2_replay(
        _manufactured_topology_fixture(), maximum_energy_eV=200.0)
    assert replay.derived_deck.processes[1].electron_energy_eV[0] == 0.0
    assert replay.derived_deck.processes[0].electron_energy_eV[-1] == 200.0
    assert replay.derived_deck.processes[14].electron_number_change == 2
    assert replay.derived_deck.processes[15].electron_number_change == 2
    assert len(replay.collision_chemistry.mappings) == 15
    assert replay.missing_reactor_channels == (
        "atomic_chlorine_ionization",
        "electron_detachment_from_Clminus",
        "tracked_vibrational_and_electronic_state_kinetics",
    )
    assert not replay.supports_direct_swarm_grade
    assert not replay.supports_reactor_state_prediction
    assert not replay.supports_feature_depth


def test_legacy_replay_refuses_signature_or_hash_drift():
    deck = _manufactured_topology_fixture()
    damaged = replace(
        deck,
        processes=(
            replace(deck.processes[0], product="silently changed"),
        ) + deck.processes[1:],
    )
    with pytest.raises(RuntimeError, match="topology/signatures changed"):
        derive_legacy_siglo_cl2_replay(damaged)
    with pytest.raises(RuntimeError, match="raw hash mismatch"):
        derive_legacy_siglo_cl2_replay(replace(
            deck, payload_sha256="0" * 64))


def test_declared_wall_extension_is_sensitivity_only_and_keeps_other_domains():
    direct = stafford_2010_conditioned_wall_recombination_provider(
        "anodized_aluminum")
    extrapolated = DeclaredLogLinearWallExtrapolation(
        provider=direct,
        allowed_cl_to_cl2_ratio=(1.0e-5, 1.5),
        source="declared source-replay sensitivity; no reactor fit",
    )
    boundary = extrapolated.predict(
        cl_to_cl2_ratio=0.9,
        pressure_Pa=2.0 * 0.1333223684,
        icp_power_W=300.0,
        gas_temperature_K=300.0,
    )
    assert boundary.evidence_kind == "sensitivity"
    assert not boundary.supports_local_prediction
    assert boundary.provenance["direct_marker_ratio_domain"] == (
        0.105610, 0.779646)
    assert boundary.provenance["coefficient_selection_target"] is None
    with pytest.raises(ValueError, match="extrapolation interval"):
        extrapolated.predict(
            cl_to_cl2_ratio=1.6,
            pressure_Pa=2.0 * 0.1333223684,
            icp_power_W=300.0,
            gas_temperature_K=300.0,
        )
    with pytest.raises(ValueError, match="source wall-provider domain"):
        extrapolated.predict(
            cl_to_cl2_ratio=0.9,
            pressure_Pa=2.0 * 0.1333223684,
            icp_power_W=300.0,
            gas_temperature_K=333.0,
        )


def test_two_target_replay_adds_atomic_momentum_and_measured_ionization():
    molecular = derive_legacy_siglo_cl2_replay(
        _manufactured_topology_fixture(), maximum_energy_eV=200.0)
    energies = tuple(25.0 * index / 27.0 for index in range(28))
    replay = _derive_legacy_siglo_comsol_chlorine_replay(
        molecular,
        atomic_momentum_energy_eV=energies,
        atomic_momentum_cross_section_m2=(2.0e-20,) * 28,
        atomic_momentum_payload_sha256="a" * 64,
    )
    assert replay.derived_deck.targets == ("Cl", "Cl2")
    assert replay.derived_deck.processes[-2].kind == "MOMENTUM"
    assert replay.derived_deck.processes[-2].target == "Cl"
    assert replay.derived_deck.processes[-1].kind == "IONIZATION"
    assert replay.derived_deck.processes[-1].target == "Cl"
    assert replay.collision_chemistry.mappings[-1].reaction_name == (
        "atomic_chlorine_ionization")
    assert replay.missing_reactor_channels == (
        "electron_detachment_from_Clminus",
        "tracked_vibrational_and_electronic_state_kinetics",
    )
    assert not replay.supports_reactor_state_prediction
