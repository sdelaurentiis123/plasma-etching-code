from pathlib import Path

import numpy as np
import pytest

from petch.charging2d import _build_edge_array_geometry, solve_edge_array_charging, solve_trench_charging
from petch.experimental_boundary import build_hwang_giapis_1997_boundary_state


IEDF = (
    Path(__file__).resolve().parents[1] / "data" / "experimental"
    / "hwang_giapis_1997" / "fig4a_ion_energy_distribution.csv")
EEDF = (
    Path(__file__).resolve().parents[1] / "data" / "experimental"
    / "hwang_giapis_1997" / "fig4b_electron_energy_distribution.csv")


def test_split_cell_stays_symmetric_without_open_area():
    r = solve_trench_charging(2.0, W=12, pad=6, mouth=50, n_per_iter=500, n_iter=8,
                              seed=11, poly_mode="split")
    assert abs(r["V_poly_left"] - r["V_poly_right"]) < 2.0


def test_edge_open_line_of_sight_reports_current_budget():
    r = solve_trench_charging(2.0, W=12, pad=6, mouth=50, n_per_iter=500, n_iter=8,
                              seed=12, poly_mode="edge_open",
                              edge_open_model="line_of_sight", edge_open_samples=4096)
    eo = r["diag"]["edge_open"]
    assert eo["model"] == "line_of_sight"
    assert eo["electron_gross"] >= 0.0
    assert eo["ion_gross"] >= 0.0
    assert np.isfinite(eo["net_electron"])
    assert "V_poly_edge" in r and "V_poly_neighbor" in r


@pytest.mark.parametrize("electron_model", ("mc", "viewfactor"))
def test_edge_array_line_of_sight_tops_up_gross_species_not_net_current(
        electron_model):
    r = solve_edge_array_charging(
        1.2, W=12, mouth=40, n_per_iter=320, n_iter=8, seed=13,
        open_width_um=1.0, rf_bursts=False,
        edge_open_model="line_of_sight", edge_open_samples=4096,
        electron_model=electron_model)
    eo = r["diag"]["edge_open"]
    assert eo["boundary_electron_gross"] >= 0.0
    assert eo["boundary_ion_gross"] >= 0.0
    assert eo["electron_gross"] == pytest.approx(max(
        eo["explicit_electron_gross"], eo["target_electron_gross"]))
    assert eo["ion_gross"] == pytest.approx(max(
        eo["explicit_ion_gross"], eo["target_ion_gross"]))
    assert eo["boundary_net_electron"] == pytest.approx(
        eo["boundary_electron_gross"] - eo["boundary_ion_gross"])


def test_edge_array_accepts_the_unified_hwang_boundary_state():
    boundary = build_hwang_giapis_1997_boundary_state(
        IEDF, EEDF, reference_plane_m=1.4e-6)
    coarse = solve_edge_array_charging(
        1.2, W=10, mouth=16, n_per_iter=160, n_iter=3, seed=14,
        feature_w_um=0.5, open_width_um=1.0, rf_bursts=False,
        plasma_boundary=boundary, source_model="boundary_state",
        return_final_ion_lineage=True, final_audit_samples=320,
        final_audit_seed=99)
    refined = solve_edge_array_charging(
        1.2, W=10, mouth=16, n_per_iter=160, n_iter=3, seed=14,
        feature_w_um=0.5, open_width_um=1.0, rf_bursts=False,
        plasma_boundary=boundary, source_model="boundary_state",
        return_final_ion_lineage=True, final_audit_samples=640,
        final_audit_seed=99)
    assert coarse["diag"]["trace"]["source_model"] == "PlasmaBoundaryState"
    assert coarse["diag"]["trace"]["source_launch_plane"] == (
        "sheath_lower_boundary")
    assert coarse["diag"]["trace"]["source_z_grid"] == pytest.approx(1.0)
    assert coarse["diag"]["trace"]["nominal_reference_height_m"] == (
        pytest.approx(1.4e-6))
    assert coarse["diag"]["trace"]["explicit_2d_source_projection"] is True
    assert coarse["diag"]["trace"]["final_sampling_protocol"] == (
        "dedicated_nested_boundary_units_v1")
    assert coarse["final_ion_lineage"]["source_particle_count"] == 320
    assert refined["final_ion_lineage"]["source_particle_count"] == 640
    for name in (
            "hit_type", "hit_ix", "hit_iz", "impact_energy_eV",
            "hit_x_grid", "hit_z_grid", "hit_vx_sqrt_eV",
            "hit_vz_sqrt_eV", "termination"):
        assert np.array_equal(
            coarse["final_ion_lineage"][name],
            refined["final_ion_lineage"][name][:320])
    assert np.all(
        refined["final_ion_lineage"]["impact_energy_eV"] >= 0.0)
    assert refined["final_ion_lineage"]["truncation_count"] == 0


def test_edge_array_legacy_mouth_launch_is_explicit_and_reported():
    result = solve_edge_array_charging(
        1.2, W=10, mouth=16, n_per_iter=80, n_iter=2, seed=15,
        feature_w_um=0.5, open_width_um=1.0, rf_bursts=False,
        source_launch_plane="feature_mouth_legacy")
    assert result["diag"]["trace"]["source_launch_plane"] == (
        "feature_mouth_legacy")
    assert result["diag"]["trace"]["source_z_grid"] == pytest.approx(15.5)
    assert result["geom"]["source_launch_plane"] == (
        "feature_mouth_legacy")


def test_edge_array_refuses_the_old_mouth_height_as_boundary_reference():
    boundary = build_hwang_giapis_1997_boundary_state(
        IEDF, EEDF, reference_plane_m=0.8e-6)
    with pytest.raises(
            ValueError,
            match="sheath-boundary-to-SiO2 height"):
        solve_edge_array_charging(
            1.2, W=10, mouth=16, n_per_iter=80, n_iter=2, seed=16,
            feature_w_um=0.5, open_width_um=1.0, rf_bursts=False,
            plasma_boundary=boundary, source_model="boundary_state")


def test_edge_array_geometry_labels_conductors_and_trench():
    g = _build_edge_array_geometry(2.0, W=12, mouth=40, open_width_um=1.0)
    assert g["edge0"] < g["edge1"] <= g["trench0"] < g["trench1"] <= g["neigh0"] < g["neigh1"]
    assert g["neigh1"] <= g["right_trench0"] < g["right_trench1"] <= g["next0"]
    assert int((g["cond"] == 1).sum()) > 0
    assert int((g["cond"] == 2).sum()) > 0
    assert int(g["floor_trench_mask"].sum()) == 12
    assert not g["solid"][g["right_trench0"]:g["right_trench1"], 40:-1].any()


def test_hwang_mirror_cell_has_source_domain_width_and_no_extra_line():
    g = _build_edge_array_geometry(
        2.6, W=20, mouth=50, poly_um=0.3, feature_w_um=0.5,
        open_width_um=2.0, domain_model="hwang_mirror_cell")
    # Fig. 3 width: half of the 2 um open area + edge line + trench
    # + neighboring line + half of the ordinary 0.5 um space.
    assert g["nx"] == 110
    assert g["mirror_x"]
    assert g["domain_model"] == "hwang_mirror_cell"
    assert g["next0"] == g["nx"]
    assert not g["solid"][g["right_trench0"]:g["right_trench1"], 50:-1].any()


def test_edge_array_tiny_solve_reports_edge_observables():
    r = solve_edge_array_charging(1.2, W=12, mouth=40, n_per_iter=160, n_iter=3,
                                  seed=21, open_width_um=1.0, rf_bursts=False,
                                  return_final_ion_lineage=True,
                                  final_audit_samples=777)
    assert r["floor_flux"] >= 0.0
    assert "V_poly_edge" in r and "V_poly_neighbor" in r
    assert r["diag"]["edge_open"]["model"] == "explicit_geometry"
    assert r["diag"]["trace"]["last_ion"]["survivor_frac"] < 0.05
    assert r["diag"]["trace"]["last_ion"]["truncation_frac"] == 0.0
    lineage = r["final_ion_lineage"]
    assert lineage["source_particle_count"] == 777
    assert lineage["hit_type"].shape == (777,)
    assert not lineage["hit_type"].flags.writeable
    assert np.allclose(
        lineage["impact_energy_eV"],
        lineage["hit_vx_sqrt_eV"] ** 2 + lineage["hit_vz_sqrt_eV"] ** 2)
    assert np.all(np.isin(lineage["termination"], (1, 2)))


def test_edge_array_adaptive_horizon_matches_one_shot_long_horizon():
    controls = dict(
        AR=1.2, W=8, mouth=20, n_per_iter=80, n_iter=2,
        seed=25, open_width_um=1.0, rf_bursts=False,
        return_final_ion_lineage=True, final_audit_samples=320,
        final_audit_seed=125)
    adaptive = solve_edge_array_charging(
        **controls, trace_step_cap_factor=1.0,
        trace_adaptive_horizon=True,
        trace_emergency_step_cap_factor=32.0)
    one_shot = solve_edge_array_charging(
        **controls, trace_step_cap_factor=32.0,
        trace_adaptive_horizon=False,
        trace_emergency_step_cap_factor=32.0)

    assert (
        adaptive["diag"]["trace"]["last_ion"]["horizon_extension_count"]
        > 0)
    assert adaptive["diag"]["trace"]["last_ion"]["truncated"] == 0
    assert adaptive["diag"]["trace"]["last_electron"]["truncated"] == 0
    for name in (
            "hit_type", "hit_ix", "hit_iz", "impact_energy_eV",
            "hit_x_grid", "hit_z_grid", "hit_vx_sqrt_eV",
            "hit_vz_sqrt_eV", "termination"):
        assert np.array_equal(
            adaptive["final_ion_lineage"][name],
            one_shot["final_ion_lineage"][name])
    assert adaptive["diag"]["residual_snapshot"] == pytest.approx(
        one_shot["diag"]["residual_snapshot"])


def test_edge_array_refuses_unresolved_trajectory_without_adaptive_horizon():
    with pytest.raises(RuntimeError, match="exhausted max_steps"):
        solve_edge_array_charging(
            1.2, W=8, mouth=20, n_per_iter=80, n_iter=2,
            seed=26, open_width_um=1.0, rf_bursts=False,
            trace_step_cap_factor=0.1,
            trace_adaptive_horizon=False,
            trace_emergency_step_cap_factor=0.1)


def test_edge_array_closes_only_a_declared_bounded_trajectory_tail():
    result = solve_edge_array_charging(
        1.2, W=8, mouth=20, n_per_iter=80, n_iter=2,
        seed=26, open_width_um=1.0, rf_bursts=False,
        trace_step_cap_factor=4.0,
        trace_adaptive_horizon=False,
        trace_emergency_step_cap_factor=4.0,
        trace_relative_tail_tolerance=0.05,
        return_final_ion_lineage=True,
        final_audit_samples=320, final_audit_seed=126)

    horizon = result["diag"]["trace"]["campaign_horizon"]
    assert horizon["ion"]["tail_closure_evaluation_count"] > 0
    assert horizon["electron"]["tail_closure_evaluation_count"] > 0
    for species in horizon.values():
        assert (
            species[
                "maximum_tail_closure_l1_surface_current_error_bound_relative"]
            <= 0.05)
    assert result["final_ion_lineage"]["tail_closure_count"] > 0


def test_edge_array_decreasing_gain_reports_state_stationarity():
    result = solve_edge_array_charging(
        1.2, W=10, mouth=30, n_per_iter=160, n_iter=12,
        seed=31, open_width_um=1.0, rf_bursts=False,
        stochastic_gain_exponent=0.75, stochastic_gain_offset=8.0)
    gain = result["diag"]["stochastic_gain"]
    stationarity = result["diag"]["potential_stationarity"]
    assert gain["mode"] == "robbins_monro"
    assert gain["final"] < gain["initial"]
    assert stationarity["tail_iterations"] >= 5
    assert stationarity["first_half_iterations"] > 0
    assert stationarity["second_half_iterations"] > 0
    assert stationarity["surface_rms_drift_v"] >= 0.0
    assert stationarity["floor_trench"]["sample_count"] == 10
    assert stationarity["target_notch_band_field"]["sample_count"] > 1
    assert stationarity["surface_maximum_drift_location"]["ix"] >= 0
    assert set(result["diag"]["region_current_balance"]) == {
        "floor", "pr", "poly_edge", "poly_neighbor"}
    continuation = result["continuation_state"]
    assert continuation["surface_potential_v"].shape == result["V"].shape
    assert continuation["stochastic_gain_age"] == 12


def test_edge_array_continuation_preserves_gain_age_and_accepts_saved_state():
    first = solve_edge_array_charging(
        1.2, W=10, mouth=30, n_per_iter=160, n_iter=8,
        seed=32, open_width_um=1.0, rf_bursts=False,
        stochastic_gain_exponent=0.75, stochastic_gain_offset=8.0)
    state = first["continuation_state"]
    second = solve_edge_array_charging(
        1.2, W=10, mouth=30, n_per_iter=160, n_iter=4,
        seed=33, open_width_um=1.0, rf_bursts=False,
        stochastic_gain_exponent=0.75, stochastic_gain_offset=8.0,
        initial_surface_potential_v=state["surface_potential_v"],
        initial_edge_potential_v=state["edge_potential_v"],
        initial_neighbor_potential_v=state["neighbor_potential_v"],
        stochastic_gain_age=state["stochastic_gain_age"])
    gain = second["diag"]["stochastic_gain"]
    assert gain["starting_age"] == 8
    assert gain["ending_age"] == 12
    assert second["diag"]["potential_history"][0]["gain_age"] == 9
    expected_first_gain = (8.0 / 16.0) ** 0.75
    assert gain["initial"] == pytest.approx(expected_first_gain)
    assert second["continuation_state"]["stochastic_gain_age"] == 12


def test_edge_array_progress_callback_reports_bounded_state_history():
    records = []
    result = solve_edge_array_charging(
        1.2, W=10, mouth=30, n_per_iter=120, n_iter=7,
        seed=34, open_width_um=1.0, rf_bursts=False,
        stochastic_gain_exponent=0.75,
        progress_callback=records.append)
    assert records == result["diag"]["potential_history"]
    assert records[-1]["iteration"] == 7
    assert records[-1]["gain_age"] == 7


@pytest.mark.parametrize("exponent", (0.5, 1.01, np.nan))
def test_edge_array_refuses_invalid_decreasing_gain(exponent):
    with pytest.raises(ValueError, match="stochastic_gain_exponent"):
        solve_edge_array_charging(
            1.2, W=8, mouth=20, n_per_iter=80, n_iter=3,
            stochastic_gain_exponent=exponent)


@pytest.mark.parametrize("gain_age", (-1, 1.5, True))
def test_edge_array_refuses_invalid_stochastic_gain_age(gain_age):
    with pytest.raises(ValueError, match="stochastic_gain_age"):
        solve_edge_array_charging(
            1.2, W=8, mouth=20, n_per_iter=80, n_iter=3,
            stochastic_gain_exponent=0.75,
            stochastic_gain_age=gain_age)


def test_edge_array_refuses_unknown_source_launch_plane():
    with pytest.raises(ValueError, match="source_launch_plane"):
        solve_edge_array_charging(
            1.2, W=8, mouth=20, n_per_iter=80, n_iter=3,
            source_launch_plane="near_the_feature")
