import numpy as np

from scripts.audit_zhu_npg80_axisymmetric_ccp import (
    DEFAULT_STATE,
    _load_input,
    build_receipt,
)


def test_effective_bohm_mass_uses_inventory_not_graded_flux():
    state = _load_input(DEFAULT_STATE)
    original = state.effective_bohm_mass_amu
    # The withheld 0-D axial flux is a grade, not an input to the aggregate
    # multi-ion Bohm mass or spatial solve.
    modified = type(state)(
        condition_id=state.condition_id,
        geometry=state.geometry,
        positive_ion_density_m3=state.positive_ion_density_m3,
        global_axial_positive_ion_flux_m2_s={
            name: value * (2.0 if name == "CF3+" else 0.5)
            for name, value in state.global_axial_positive_ion_flux_m2_s.items()
        },
        electron_density_m3=state.electron_density_m3,
        electronegativity=state.electronegativity,
        mean_electron_energy_eV=state.mean_electron_energy_eV,
        total_neutral_density_m3=state.total_neutral_density_m3,
        ion_temperature_eV=state.ion_temperature_eV,
        ion_momentum_mean_free_path_m=state.ion_momentum_mean_free_path_m,
        source=state.source,
    )
    assert modified.effective_bohm_mass_amu == original


def test_axisymmetric_ccp_closes_global_flux_without_target_fit():
    receipt = build_receipt()
    central = receipt["central_48x16_result"]
    certification = receipt["certification"]

    assert receipt["sem_target_used"] is False
    assert receipt["measured_depth_target_used"] is False
    assert abs(central["global_to_spatial_relative_residual"]) < 0.003
    assert 1.005 < central["central_3mm_to_full_electrode_flux_ratio"] < 1.008
    assert 1.01 < central["first_to_last_annulus_flux_ratio"] < 1.02
    assert certification["global_flux_reproduced_within_1_percent"] is True
    assert certification["absolute_target_wafer_flux_supported"] is False


def test_axisymmetric_ccp_grid_convergence_and_particle_ledger():
    receipt = build_receipt()
    convergence = receipt["grid_convergence"]
    assert convergence["passed_0p1_percent"] is True
    assert convergence["central_to_fine_full_flux_relative_change"] < 3.0e-5
    assert convergence["central_to_fine_optic_flux_relative_change"] < 3.0e-5
    for row in receipt["resolution_board"]:
        assert row["maximum_species_ledger_relative_residual"] < 2.0e-13
        assert row["maximum_inventory_relative_residual"] < 1.0e-12


def test_committed_radial_profile_is_center_high_and_smooth():
    receipt = build_receipt()
    central = receipt["central_48x16_result"]
    radius = np.asarray(central["radial_center_m"])
    flux = np.asarray(central["total_lower_endcap_flux_m2_s"])
    assert np.all(np.diff(radius) > 0.0)
    assert np.all(np.diff(flux) < 0.0)
    assert np.max(np.abs(np.diff(flux, n=2))) / np.mean(flux) < 2.0e-4
    assert np.isclose(sum(central["species_flux_fraction"].values()), 1.0)
    local = receipt["central_3mm_smooth_radial_resolution"]
    assert local["center_to_next_annulus_flux_ratio"] > 1.0
    assert local["resolved_flux_change_percent"] < 0.001
    certification = receipt["certification"]
    assert (
        certification[
            "central_3mm_average_enhancement_vs_full_electrode_percent"
        ]
        > local["resolved_flux_change_percent"]
    )
