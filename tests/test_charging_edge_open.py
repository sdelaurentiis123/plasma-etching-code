import numpy as np

from petch.charging2d import solve_trench_charging


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
