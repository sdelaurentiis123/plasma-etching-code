import pytest

from scripts.audit_zhu_npg80_species_resolved_feature_boundary import (
    build_receipt,
)


@pytest.fixture(scope="module")
def receipt():
    return build_receipt()


def test_all_reactor_species_reach_the_feature_boundary(receipt):
    inventory = receipt["inventory"]
    assert inventory["positive_ion_species_count"] == 20
    assert inventory["thermal_neutral_species_count"] == 37
    assert receipt["certification"][
        "reactor_species_identity_preserved"
    ] is True
    assert receipt["certification"][
        "species_specific_iadf_interface_supported"
    ] is True


def test_flux_and_deterministic_replay_close(receipt):
    cases = receipt["boundary_cases"]
    assert len(cases) == 2
    assert cases[0]["boundary_sha256"] != cases[1]["boundary_sha256"]
    for row in cases:
        assert abs(row["ion_flux_relative_conservation_residual"]) < 2e-15
        assert abs(row["neutral_flux_relative_conservation_residual"]) < 2e-15
        assert row["deterministic_quadrature"] is True
        assert row["monte_carlo"] is False


def test_charge_energy_and_evidence_gates_are_explicit(receipt):
    multiply_charged = receipt["inventory"]["multiply_charged_ions"]
    assert set(multiply_charged) == {"SF2++", "SF4++"}
    assert all(row["charge_number"] == 2 for row in multiply_charged.values())
    singly_charged_energy = (
        receipt["inputs"]["powered_electrode_sheath_drop_V"]
        + 0.5 * receipt["inputs"]["electron_temperature_eV"]
    )
    assert multiply_charged["SF2++"]["impact_energy_eV"] > (
        1.99 * singly_charged_energy
    )
    certification = receipt["certification"]
    assert certification["conditional_feature_transport_boundary_executable"]
    assert certification["supports_unique_absolute_profile_prediction"] is False
    assert certification["target_machine_iead_measured"] is False
    assert certification["target_tio2_cr_surface_probabilities_validated"] is False
