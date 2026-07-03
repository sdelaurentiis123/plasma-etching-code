import numpy as np

from petch.charging2d import _build_edge_array_geometry, solve_edge_array_charging, solve_trench_charging


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


def test_edge_array_geometry_labels_conductors_and_trench():
    g = _build_edge_array_geometry(2.0, W=12, mouth=40, open_width_um=1.0)
    assert g["edge0"] < g["edge1"] <= g["trench0"] < g["trench1"] <= g["neigh0"] < g["neigh1"]
    assert int((g["cond"] == 1).sum()) > 0
    assert int((g["cond"] == 2).sum()) > 0
    assert int(g["floor_trench_mask"].sum()) == 12


def test_edge_array_tiny_solve_reports_edge_observables():
    r = solve_edge_array_charging(1.2, W=12, mouth=40, n_per_iter=160, n_iter=3,
                                  seed=21, open_width_um=1.0, rf_bursts=False)
    assert r["floor_flux"] >= 0.0
    assert "V_poly_edge" in r and "V_poly_neighbor" in r
    assert r["diag"]["edge_open"]["model"] == "explicit_geometry"
    assert r["diag"]["trace"]["last_ion"]["survivor_frac"] < 0.05
