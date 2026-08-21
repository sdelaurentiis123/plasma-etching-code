from hashlib import sha256
import json
from pathlib import Path

from scripts.audit_bosch_recipe_path_spatial_residual import (
    OUTPUT,
    V7_FIT,
    V7_PREREGISTRATION,
    V7_RESPONSE_TABLE,
    build,
)


def _hash(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def test_spatial_residual_discovery_is_reproducible_and_target_firewalled():
    written = json.loads(OUTPUT.read_text(encoding="utf-8"))
    rebuilt = build()

    assert rebuilt == written
    assert written["heldout_outcomes_read"] is False
    assert written["heldout_prediction_written"] is False
    assert written["eligible_for_prediction_seal"] is False
    assert written["surface_laws_changed"] is False
    assert written["positive_ion_boundary_changed"] is False
    assert written["input_hashes"]["v7_preregistration"] == _hash(
        V7_PREREGISTRATION)
    assert written["input_hashes"]["v7_fit"] == _hash(V7_FIT)
    assert written["input_hashes"]["v7_response_table"] == _hash(
        V7_RESPONSE_TABLE)


def test_spatial_residual_localizes_fixed_map_and_voltage_edge_mode():
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    decomposition = payload["decomposition"]
    proxies = payload["whole_lot_output_space_proxies"]
    diagnostic = payload["selected_physical_family_diagnostic"]

    assert decomposition[
        "fraction_of_squared_residual_in_fixed_mean_map"] > 0.85
    assert proxies["shared_89_point_intercept"]["metrics"][
        "normalized_shape_rmse_percent"] < 0.637
    assert proxies["shared_89_point_intercept"]["metrics"][
        "silicon_point_rmse_um"] < 0.4866
    assert diagnostic["standardized_slope_outer_mean_percent"] > 0.0
    assert diagnostic["standardized_slope_inner_mean_percent"] < 0.0
    assert diagnostic["standardized_slope_map_edge_basis_pearson"] > 0.7
