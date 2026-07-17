import importlib.util
from pathlib import Path

import numpy as np

from petch.amorphous_carbon_mask import build_krueger_2024_material_router_3d
from petch.material_mechanism_3d import MaterialSurfaceState3D
from petch.surface_kinetics import EnergeticFlux, SurfaceFluxes


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "krueger_2024_frozen_surface_chemistry",
    ROOT / "scripts" / "krueger_2024_frozen_surface_chemistry.py",
)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def test_runtime_provenance_binds_full_petch_tree():
    assert "src/petch/feature_geometry_state_3d.py" in AUDIT.SOURCE_PATHS
    expected = {
        str(path.relative_to(ROOT))
        for path in (ROOT / "src" / "petch").glob("*.py")
        if not path.name.startswith("._")
    }
    assert expected.issubset(AUDIT.SOURCE_PATHS)


def _manufactured_case(*, polymer_units=1.0e15, depositing=True, flux_scale=1.0):
    material = np.asarray([1, 2], dtype=int)
    mechanism = build_krueger_2024_material_router_3d(
        effective_mask_crosslinked_growth_fraction=0.9,
        oxide_etch_yield_scale=0.56,
    )
    initial = mechanism.initial_state_by_material(material)
    fields = {
        name: np.asarray(values).copy() for name, values in initial.fields.items()
    }
    fields["m1__polymer_units_m2"][0] = polymer_units
    fields["m1__complex_fraction"][0] = 0.5
    fields["m2__polymer_units_m2"][1] = polymer_units
    state = MaterialSurfaceState3D(
        fields, initial.upper_bounds, initial.remap_modes
    )
    polymer_flux = 1.0e19 * flux_scale if depositing else 0.0
    neutral = {
        "C3F4": np.zeros(2),
        "C2F3": np.full(2, polymer_flux),
        "CF": np.full(2, polymer_flux),
        "CF2": np.full(2, polymer_flux),
        "CF3": np.full(2, polymer_flux),
        "O": np.full(2, 1.0e20 * flux_scale),
    }
    energetic = EnergeticFlux(
        "ions", np.full(2, 1.0e20 * flux_scale),
        np.asarray([500.0]), np.asarray([1.0]), np.asarray([1.0]),
    )
    return (
        state,
        SurfaceFluxes(neutral, (energetic,)),
        np.full(2, 1.0e-12),
        material,
        mechanism,
    )


def test_microstep_closes_ledger_and_preserves_state_and_flux_inputs():
    state, flux, area, material, mechanism = _manufactured_case()
    state_before = {
        name: np.asarray(values).copy() for name, values in state.fields.items()
    }
    flux_hash = AUDIT.surface_flux_sha256(flux)

    result = AUDIT.advance_frozen_surface_chemistry_3d(
        state, flux, area, material, mechanism,
        horizon_s=1.0e-3, substeps=2, dx_m=1.0e-6,
    )

    assert result["integrated_exchange"][
        "maximum_step_ledger_residual_units_m2"] == 0.0
    assert result["integrated_exchange"][
        "maximum_cumulative_ledger_residual_units_m2"] == 0.0
    assert result["oxide_removal"]["integrated_formula_units"] > 0.0
    assert result["state"] is not state
    assert all(
        np.array_equal(values, state.fields[name])
        for name, values in state_before.items()
    )
    assert AUDIT.surface_flux_sha256(flux) == flux_hash


def test_n_2n_refinement_covers_removed_deposited_and_state_increments():
    state, flux, area, material, mechanism = _manufactured_case()

    result = AUDIT.evaluate_frozen_horizon(
        state, flux, area, material, mechanism,
        horizon_s=1.0e-3, coarse_substeps=1, dx_m=1.0e-6,
    )

    assert result["all_gates_pass"]
    assert result["refinement"]["maximum_relative_error"] < 0.01
    assert set(result["refinement"]["by_exchange_inventory"]) == {
        "removed_units", "deposited_units"
    }
    assert set(result["refinement"]["by_final_state_increment"]) == set(
        state.fields
    )


def test_displacement_and_reaction_probability_drift_each_refuse(monkeypatch):
    state, flux, area, material, mechanism = _manufactured_case()
    monkeypatch.setitem(AUDIT.GATES, "maximum_gross_displacement_dx", 0.0)
    displaced = AUDIT.evaluate_frozen_horizon(
        state, flux, area, material, mechanism,
        horizon_s=1.0e-3, coarse_substeps=1, dx_m=1.0e-6,
    )
    assert not displaced["gates"]["gross_displacement_within_frozen_geometry_limit"]
    assert not displaced["all_gates_pass"]

    monkeypatch.setitem(AUDIT.GATES, "maximum_gross_displacement_dx", 0.05)
    monkeypatch.setitem(
        AUDIT.GATES, "maximum_neutral_reaction_probability_absolute_drift", 0.0
    )
    drifted = AUDIT.evaluate_frozen_horizon(
        state, flux, area, material, mechanism,
        horizon_s=1.0e-2, coarse_substeps=1, dx_m=1.0e-6,
    )
    assert drifted["fine"]["neutral_reaction_probability_drift"][
        "maximum_absolute_probability_drift"] > 0.0
    assert not drifted["gates"]["neutral_reaction_probability_stable"]
    assert not drifted["all_gates_pass"]


def test_film_depletion_is_detected_without_mutating_checkpoint_state():
    state, flux, area, material, mechanism = _manufactured_case(
        polymer_units=1.0e15, depositing=False, flux_scale=10.0
    )
    before = np.asarray(state.fields["m1__polymer_units_m2"]).copy()

    result = AUDIT.advance_frozen_surface_chemistry_3d(
        state, flux, area, material, mechanism,
        horizon_s=1.0, substeps=2, dx_m=1.0e-6,
    )

    assert result["film_depletion"]["effectively_depleted_oxide_face_count"] == 1
    assert result["film_depletion"]["exactly_depleted_oxide_face_count"] == 1
    assert np.array_equal(before, state.fields["m1__polymer_units_m2"])
