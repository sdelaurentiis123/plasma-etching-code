import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from petch.amorphous_carbon_mask import build_krueger_2024_material_router_3d
from petch.material_mechanism_3d import MaterialSurfaceState3D
from petch.neutral_radiosity_3d import (
    DiffuseFormFactors3D,
    solve_diffuse_neutral_radiosity_3d,
)
from petch.surface_kinetics import EnergeticFlux, SurfaceFluxes


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "krueger_2024_frozen_radiosity_chemistry",
    ROOT / "scripts" / "krueger_2024_frozen_radiosity_chemistry.py",
)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def _manufactured_case(*, flux_scale=1.0):
    material = np.asarray([1, 2], dtype=int)
    mechanism = build_krueger_2024_material_router_3d(
        effective_mask_crosslinked_growth_fraction=0.9,
        oxide_etch_yield_scale=0.56,
    )
    initial = mechanism.initial_state_by_material(material)
    fields = {
        name: np.asarray(values).copy() for name, values in initial.fields.items()
    }
    fields["m1__polymer_units_m2"][0] = 1.0e15
    fields["m1__complex_fraction"][0] = 0.5
    fields["m2__polymer_units_m2"][1] = 1.0e15
    state = MaterialSurfaceState3D(
        fields, initial.upper_bounds, initial.remap_modes)
    neutral = {
        "C3F4": np.zeros(2),
        "C2F3": np.full(2, 1.0e19 * flux_scale),
        "CF": np.full(2, 1.0e19 * flux_scale),
        "CF2": np.full(2, 1.0e19 * flux_scale),
        "CF3": np.full(2, 1.0e19 * flux_scale),
        "O": np.full(2, 1.0e20 * flux_scale),
    }
    energetic = EnergeticFlux(
        "ions", np.full(2, 1.0e20 * flux_scale),
        np.asarray([500.0]), np.asarray([1.0]), np.asarray([1.0]))
    flux = SurfaceFluxes(neutral, (energetic,))
    factors = DiffuseFormFactors3D(
        face_count=2,
        source_face=np.asarray([0, 1]),
        target_face=np.asarray([1, 0]),
        transfer_fraction=np.asarray([0.25, 0.25]),
        escape_fraction=np.asarray([0.75, 0.75]),
        rays_per_face=4,
    )
    role = {name: "neutral_reactant" for name in neutral}
    role["ions"] = "energetic_bombardment"
    return (
        state, flux, np.full(2, 1.0e-12), factors,
        np.arange(2), material, role, mechanism,
    )


def test_runtime_provenance_binds_full_petch_tree():
    assert "src/petch/feature_geometry_state_3d.py" in AUDIT.SOURCE_PATHS
    expected = {
        str(path.relative_to(ROOT))
        for path in (ROOT / "src" / "petch").glob("*.py")
        if not path.name.startswith("._")
    }
    assert expected.issubset(AUDIT.SOURCE_PATHS)


def test_cli_defaults_preserve_archived_8ray_horizon_contract():
    args = AUDIT.parse_args([])

    assert args.rays_per_face == 8
    assert str(args.q3_reference).endswith(
        "frozen_checkpoint_q3_seed241_current/audit.json")
    assert args.requested_horizon_fractions == AUDIT.HORIZON_FRACTIONS
    assert args.horizon_fractions == AUDIT.HORIZON_FRACTIONS
    assert args.maximum_horizon_fraction == 1.0
    assert AUDIT._certification_mode_for_rays(args.rays_per_face) == (
        "archived_q3_8ray_plus_production_facewise_cache_exact")


def test_cli_selects_one_nondefault_short_horizon_and_refuses_invalid_values():
    args = AUDIT.parse_args([
        "--rays-per-face", "16",
        "--horizon-fractions", "0.0625", "0.125", "0.25",
        "--maximum-horizon-fraction", "0.0625",
    ])
    assert args.rays_per_face == 16
    assert args.requested_horizon_fractions == (0.0625, 0.125, 0.25)
    assert args.horizon_fractions == (0.0625,)
    assert AUDIT._certification_mode_for_rays(args.rays_per_face) == (
        "nondefault_rays_production_facewise_cache_exact")

    invalid = (
        ["--rays-per-face", "0"],
        ["--rays-per-face", "3"],
        ["--rays-per-face", "2048"],
        ["--horizon-fractions", "0"],
        ["--horizon-fractions", "1.1"],
        ["--horizon-fractions", "0.25", "0.125"],
        ["--horizon-fractions", "0.125", "0.125"],
        ["--horizon-fractions", "0.25", "--maximum-horizon-fraction", "0.125"],
    )
    for argv in invalid:
        with pytest.raises(SystemExit):
            AUDIT.parse_args(argv)


def test_production_capture_receives_requested_ray_level(monkeypatch):
    factors = DiffuseFormFactors3D(
        face_count=1,
        source_face=np.asarray([], dtype=int),
        target_face=np.asarray([], dtype=int),
        transfer_fraction=np.asarray([], dtype=float),
        escape_fraction=np.asarray([1.0]),
        rays_per_face=16,
    )
    base_transport = SimpleNamespace(surface_fluxes=SurfaceFluxes({"A": [1.0]}))
    observed = {}

    monkeypatch.setattr(
        AUDIT.feature_step_module, "estimate_diffuse_form_factors_3d",
        lambda *args, **kwargs: factors)
    monkeypatch.setattr(
        AUDIT.feature_step_module, "_apply_diffuse_neutral_transport",
        lambda transport, *args, **kwargs: transport)

    def fake_evaluate(*args, **kwargs):
        observed["rays_per_face"] = kwargs["radiosity_rays"]
        AUDIT.feature_step_module.estimate_diffuse_form_factors_3d()
        AUDIT.feature_step_module._apply_diffuse_neutral_transport(base_transport)
        result = SimpleNamespace(
            active_face_index=np.asarray([0]),
            active_face_area=np.asarray([1.0]))
        return result, SimpleNamespace(), 0.01

    monkeypatch.setattr(AUDIT, "_evaluate", fake_evaluate)
    monkeypatch.setattr(
        AUDIT, "extract_mesh_3d",
        lambda phi, dx: (
            np.zeros((3, 3)), np.asarray([[0, 1, 2]]),
            np.zeros((1, 3)), np.asarray([1.0])))
    source = {
        "geometry": SimpleNamespace(
            phi=np.zeros((2, 2, 2)), dx=1.0, mesh_length_unit_m=1.0),
        "state": object(),
        "fingerprint": "fingerprint",
    }

    captured = AUDIT._evaluate_with_captured_radiosity(
        source, {}, seed=241, rays_per_face=16)

    assert observed["rays_per_face"] == 16
    assert captured[4].rays_per_face == 16
    assert captured[6] == 1


def test_nondefault_production_facewise_certification_skips_archived_q3():
    state, direct, area, factors, active, material, _, mechanism = (
        _manufactured_case())
    cached, diagnostics, _ = AUDIT.solve_cached_radiosity(
        direct, area, factors, state, mechanism, active, material,
        relative_tolerance=1.0e-12, maximum_iterations=2000)

    certification = AUDIT.certify_production_facewise_cache_exact(
        cached, cached, area, diagnostics, diagnostics)

    assert certification["mode"] == (
        "nondefault_rays_production_facewise_cache_exact")
    assert certification["archived_q3_comparison_performed"] is False
    assert certification["facewise_production_replay"]["hash_equal"]
    assert certification["all_gates_pass"]

    changed_neutral = dict(cached.neutral_flux_m2_s)
    changed_neutral["CF"] = np.asarray(changed_neutral["CF"]).copy()
    changed_neutral["CF"][0] *= 1.001
    changed = SurfaceFluxes(changed_neutral, cached.energetic_fluxes)
    failed = AUDIT.certify_production_facewise_cache_exact(
        cached, changed, area, diagnostics, diagnostics)
    assert not failed["gates"]["facewise_hash_reproduction"]
    assert not failed["all_gates_pass"]


def test_cached_form_factors_reproduce_species_solves_and_balance():
    state, direct, area, factors, active, material, _, mechanism = (
        _manufactured_case())

    cached, diagnostics, probability = AUDIT.solve_cached_radiosity(
        direct, area, factors, state, mechanism, active, material,
        relative_tolerance=1.0e-12, maximum_iterations=2000)

    for name, direct_flux in direct.neutral_flux_m2_s.items():
        p = probability[name]
        if np.any(p > 0.0):
            reference = solve_diffuse_neutral_radiosity_3d(
                direct_flux, area, factors.source_face, factors.target_face,
                factors.transfer_fraction, factors.escape_fraction, p,
                relative_tolerance=1.0e-12, maximum_iterations=2000)
            assert np.allclose(
                cached.neutral_flux_m2_s[name], reference.incident_flux_m2_s,
                rtol=1.0e-14, atol=0.0)
        else:
            assert np.array_equal(
                cached.neutral_flux_m2_s[name], direct_flux)
        assert diagnostics[name]["relative_balance_error"] <= 5.0e-12
    facewise = AUDIT._facewise_flux_comparison(cached, cached, area)
    assert facewise["hash_equal"]
    assert facewise["maximum_relative_linf_error"] == 0.0
    assert facewise["maximum_area_weighted_relative_l1_error"] == 0.0
    changed_neutral = dict(cached.neutral_flux_m2_s)
    changed_neutral["CF"] = np.asarray(changed_neutral["CF"]).copy()
    changed_neutral["CF"][0] *= 1.001
    changed = SurfaceFluxes(changed_neutral, cached.energetic_fluxes)
    mismatch = AUDIT._facewise_flux_comparison(cached, changed, area)
    assert not mismatch["hash_equal"]
    assert mismatch["maximum_relative_linf_error"] > 0.0


def test_embedded_controller_rejects_halves_accepts_fine_and_closes_ledgers(
        monkeypatch):
    state, direct, area, factors, active, material, role, mechanism = (
        _manufactured_case())
    before = {
        name: np.asarray(values).copy() for name, values in state.fields.items()
    }
    flux_hash = AUDIT.surface_flux_sha256(direct)

    monkeypatch.setitem(AUDIT.GATES, "maximum_embedded_relative_error", 1.0e-5)
    result = AUDIT.co_integrate_frozen_radiosity_chemistry(
        state, direct, area, factors, active, material, role, mechanism,
        horizon_s=1.0e-3, maximum_local_dp=1.0e-2,
        minimum_substep_s=1.0e-6, relative_tolerance=1.0e-12,
        maximum_iterations=2000, dx_m=1.0e-6)

    assert result["rejected_substeps"] > 0
    assert result["accepted_substeps"] > 1
    assert result["rejection_history"][0]["kind"] == "embedded_step_error"
    assert result["minimum_accepted_substep_s"] >= 1.0e-6
    assert result["maximum_local_dp_seen"] <= 1.0e-2
    assert result["maximum_embedded_relative_error_seen"] <= 1.0e-5
    embedded = result["accepted_embedded_error_maxima"]
    assert set(embedded["by_state_field_increment"]) == set(state.fields)
    assert set(embedded["by_exchange_inventory"]) == {
        "removed", "outgoing", "unresolved", "deposited"}
    assert embedded["per_face_integrated_recession_relative_error"] >= 0.0
    assert embedded["per_face_integrated_growth_relative_error"] >= 0.0
    assert result["maximum_step_ledger_residual_units_m2"] == 0.0
    assert result["maximum_cumulative_ledger_residual_units_m2"] == 0.0
    assert result["oxide_removal"]["integrated_formula_units"] > 0.0
    assert all(
        np.array_equal(values, state.fields[name])
        for name, values in before.items())
    assert AUDIT.surface_flux_sha256(direct) == flux_hash


def test_public_receipt_persists_exact_final_state_and_per_face_exchange():
    state, direct, area, factors, active, material, role, mechanism = (
        _manufactured_case())
    result = AUDIT.co_integrate_frozen_radiosity_chemistry(
        state, direct, area, factors, active, material, role, mechanism,
        horizon_s=1.0e-3, maximum_local_dp=1.0e-2,
        minimum_substep_s=1.0e-6, relative_tolerance=1.0e-12,
        maximum_iterations=2000, dx_m=1.0e-6)

    public = AUDIT._public_integration(result)

    assert "state" not in public
    assert set(public["final_state_fields"]) == set(result["state"].fields)
    for name, values in result["state"].fields.items():
        assert np.array_equal(public["final_state_fields"][name], values)
        assert public["final_state_fields"][name] is not values
    assert public["final_state_fields_sha256"] == (
        AUDIT.surface_state_fields_sha256(result["state"]))
    per_face = public["per_face_integrated_exchange_units_m2"]
    assert set(per_face) == {"removed", "outgoing", "unresolved", "deposited"}
    for inventory in per_face.values():
        for values in inventory.values():
            assert np.asarray(values).shape == material.shape
    assert public["per_face_integrated_exchange_sha256"] == (
        AUDIT.per_face_exchange_sha256(per_face))


def test_array_manifest_hash_is_order_independent_and_value_sensitive():
    left = [
        ("b", np.asarray([2.0], dtype=np.float64)),
        ("a", np.asarray([1.0], dtype=np.float64)),
    ]
    right = list(reversed(left))
    changed = [
        ("a", np.asarray([1.0], dtype=np.float64)),
        ("b", np.asarray([3.0], dtype=np.float64)),
    ]

    assert AUDIT._array_manifest_sha256(left) == AUDIT._array_manifest_sha256(right)
    assert AUDIT._array_manifest_sha256(left) != AUDIT._array_manifest_sha256(changed)


def test_tolerance_halving_agrees_and_frozen_displacement_gate_refuses(monkeypatch):
    case = _manufactured_case()
    result = AUDIT.evaluate_co_integrated_horizon(
        *case,
        horizon_s=1.0e-3, minimum_substep_s=1.0e-6,
        relative_tolerance=1.0e-12, maximum_iterations=2000,
        dx_m=1.0e-6)
    assert result["all_gates_pass"]
    assert result["tolerance_halving_relative_oxide_error"] <= 0.01

    monkeypatch.setitem(AUDIT.GATES, "maximum_gross_displacement_dx", 0.0)
    stale = AUDIT.evaluate_co_integrated_horizon(
        *case,
        horizon_s=1.0e-3, minimum_substep_s=1.0e-6,
        relative_tolerance=1.0e-12, maximum_iterations=2000,
        dx_m=1.0e-6)
    assert not stale["gates"]["gross_displacement_within_frozen_geometry_limit"]
    assert not stale["all_gates_pass"]


def test_refuses_when_probability_gate_requires_subminimum_step():
    state, direct, area, factors, active, material, role, mechanism = (
        _manufactured_case())

    with pytest.raises(AUDIT.CoIntegrationRefusal, match="minimum dt"):
        AUDIT.co_integrate_frozen_radiosity_chemistry(
            state, direct, area, factors, active, material, role, mechanism,
            horizon_s=1.0e-3, maximum_local_dp=0.0,
            minimum_substep_s=5.0e-4, relative_tolerance=1.0e-12,
            maximum_iterations=2000, dx_m=1.0e-6)
