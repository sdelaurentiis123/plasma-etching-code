import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "krueger_2024_replicated_form_factor_closure",
    ROOT / "scripts" / "krueger_2024_replicated_form_factor_closure.py",
)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)
from petch.boundary_transport_3d import (  # noqa: E402
    _diffuse_form_factor_ray_samples_3d,
)


def _patch_geometry(*, x_nm=80.0, y_nm=20.0, z_nm=10.0, dx_nm=10.0):
    shape = tuple(int(round(value / dx_nm)) + 1 for value in (x_nm, y_nm, z_nm))
    return SimpleNamespace(
        phi=np.zeros(shape), dx=float(dx_nm), mesh_length_unit_m=1.0e-9,
        mesh_origin_m=np.zeros(3))


def _complete_patch_record(*, patch_scale_m=20.0e-9, passed=True):
    sensitivity = []
    for threshold in AUDIT.PATCH_SUPPORT_SENSITIVITY_THRESHOLDS:
        sensitivity.append({
            "minimum_mean_support_fraction": threshold,
            "eligible_mean_patch_count": 1,
            "excluded_mean_patch_count": 0,
            "gate_defined": True,
            "pass": bool(passed),
        })
    return {
        "patch_scale_m": patch_scale_m,
        "patch_count": 1,
        "minimum_mean_support_fraction": (
            AUDIT.DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION),
        "eligible_mean_patch_count": 1,
        "excluded_mean_patch_count": 0,
        "excluded_mean_surface_area_fraction": 0.0,
        "excluded_mean_projected_support_fraction": 0.0,
        "represented_nominal_projected_area_m2": {
            "minimum": patch_scale_m ** 2,
            "maximum": patch_scale_m ** 2,
            "total": patch_scale_m ** 2,
        },
        "support_threshold_sensitivity": {
            "thresholds": sensitivity,
            "primary_threshold": AUDIT.DEFAULT_MINIMUM_MEAN_SUPPORT_FRACTION,
            "primary_gate_conclusion": bool(passed),
            "gate_conclusion_stable_over_predeclared_thresholds": True,
        },
        "pass": bool(passed),
    }


def _source(path, *, boundary_case="base"):
    path.mkdir()
    (path / "audit.json").write_text(json.dumps({
        "status": "complete",
        "configuration": {
            "boundary_case": boundary_case,
            "effective_mask_crosslinked_growth_fraction": 0.9,
            "oxide_etch_yield_scale": 0.56,
        },
    }), encoding="utf-8")
    (path / "checkpoint.npz").write_bytes(b"checkpoint")


def test_cli_defaults_are_stage_a_only_and_bounded():
    args = AUDIT.parse_args([])

    assert args.stage == "stage_a"
    assert args.authorize_stage_b is False
    assert args.stage_a_audit is None
    assert args.replicate_seeds == AUDIT.REPLICATE_SEEDS
    assert len(args.replicate_seeds) == 8
    assert AUDIT.RAY_LEVELS == (8, 16, 32)
    assert AUDIT.DIAGNOSTIC_LEVEL_PAIR == (8, 16)
    assert AUDIT.AUTHORITATIVE_LEVEL_PAIR == (16, 32)
    assert AUDIT.FINAL_STAGE_A_LEVEL == 32
    assert AUDIT.FROZEN_RESPONSE_HORIZON_FRACTION == pytest.approx(1.0 / 1024.0)
    assert args.maximum_direct_transport_wall_s == 120.0
    assert args.maximum_form_factor_replicate_wall_s == 60.0
    assert args.maximum_endpoint_job_wall_s == 90.0
    assert args.maximum_total_wall_s == 300.0
    assert args.maximum_process_count == 4


def test_cli_refuses_implicit_stage_b_and_budget_or_seed_expansion(tmp_path):
    stage_a = tmp_path / "stage_a.json"
    invalid = (
        ["--stage", "stage_b", "--stage-a-audit", str(stage_a)],
        ["--stage", "stage_b", "--authorize-stage-b"],
        ["--authorize-stage-b"],
        ["--stage-a-audit", str(stage_a)],
        ["--replicate-seeds", "1", "2", "3", "4"],
        ["--replicate-seeds", "1", "2", "3", "4", "5", "6", "7", "7"],
        ["--maximum-direct-transport-wall-s", "121"],
        ["--maximum-form-factor-replicate-wall-s", "61"],
        ["--maximum-endpoint-job-wall-s", "91"],
        ["--maximum-total-wall-s", "901"],
        ["--maximum-total-wall-s", "301"],
        ["--maximum-process-count", "5"],
    )
    for argv in invalid:
        with pytest.raises(SystemExit):
            AUDIT.parse_args(argv)

def test_non_base_is_refused_before_checkpoint_read(tmp_path, monkeypatch):
    source = tmp_path / "heldout"
    _source(source, boundary_case="oxygen_ratio")
    monkeypatch.setattr(
        AUDIT, "_load_checkpoint",
        lambda _path: pytest.fail("non-base checkpoint must not be opened"),
    )

    with pytest.raises(ValueError, match="checkpoint was not opened"):
        AUDIT._load_sealed_base_source(source)


def test_stage_b_contract_is_structurally_held():
    with pytest.raises(RuntimeError, match="structurally held"):
        AUDIT._stage_b_job_specs()


def test_direct_run_refuses_stage_b_before_checkpoint_read(monkeypatch):
    monkeypatch.setattr(
        AUDIT, "_load_sealed_base_source",
        lambda _path: pytest.fail("Stage B must be refused before source I/O"))
    args = SimpleNamespace(
        stage="stage_b", authorize_stage_b=True, stage_a_audit="old.json",
        plan_only=False)

    with pytest.raises(ValueError, match="structurally held"):
        AUDIT.run(args)


def test_scalar_uncertainty_separates_confidence_from_authority_bias():
    unbiased = AUDIT._scalar_replicate_score(
        1.0, [0.99, 1.01, 1.0, 1.0, 1.01, 0.99, 1.0, 1.0],
        absolute_tolerance=0.01, relative_tolerance=0.05)
    biased = AUDIT._scalar_replicate_score(
        1.2, [0.99, 1.01, 1.0, 1.0, 1.01, 0.99, 1.0, 1.0],
        absolute_tolerance=0.01, relative_tolerance=0.05)

    assert unbiased["confidence_half_width"] > 0.0
    assert unbiased["authority_to_replicate_mean_bias"] == pytest.approx(0.0)
    assert unbiased["pass"]
    assert biased["authority_to_replicate_mean_bias"] == pytest.approx(0.2)
    assert not biased["pass"]


def test_stage_a_gate_holds_instead_of_automatically_escalating():
    level_summary = {
        8: {"replicate_count": 8, "parameter_scores": {
            "r17": {"all_gates_pass": True}, "r19": {"all_gates_pass": True}}},
        16: {"replicate_count": 8, "parameter_scores": {
            "r17": {"all_gates_pass": False}, "r19": {"all_gates_pass": False}}},
        32: {"replicate_count": 8, "parameter_scores": {
            "r17": {
                "all_gates_pass": True,
                "maximum_gross_displacement_dx": 0.01,
                "physical_patch_fields": {"field": [_complete_patch_record()]},
            },
            "r19": {
                "all_gates_pass": False,
                "maximum_gross_displacement_dx": 0.01,
                "physical_patch_fields": {"field": [_complete_patch_record()]},
            }}},
    }
    authoritative_nested = {
        "r17": {
            "all_gates_pass": True,
            "physical_patch_fields": {"field": [_complete_patch_record()]},
        },
        "r19": {
            "all_gates_pass": True,
            "physical_patch_fields": {"field": [_complete_patch_record()]},
        },
    }
    authoritative_paired_nested = {
        "r17": {
            "all_gates_pass": True,
            "physical_patch_fields": {"field": [_complete_patch_record()]},
        },
        "r19": {
            "all_gates_pass": True,
            "physical_patch_fields": {"field": [_complete_patch_record()]},
        },
    }
    paired = {"all_gates_pass": True}
    nested_sampling = {
        "8_to_16": {"all_gates_pass": False},
        "16_to_32": {"all_gates_pass": True},
    }

    gates = AUDIT._stage_a_gate(
        level_summary, authoritative_nested, authoritative_paired_nested,
        paired, {8: 0.0, 16: 0.0, 32: 0.0},
        nested_sampling, claimed_feature_extent_m=45.0e-9)

    assert gates["level32_replicated_uncertainty"] is False
    assert gates["exact_nested_sampling_extension_16_to_32"] is True
    assert AUDIT.RAY_LEVELS == (8, 16, 32)
    assert max(AUDIT.RAY_LEVELS) == 32


def test_stage_a_authorization_binds_checkpoint_levels_and_seeds(
        tmp_path, monkeypatch):
    checkpoint = tmp_path / "checkpoint.npz"
    checkpoint.write_bytes(b"sealed")
    source = {
        "checkpoint_path": checkpoint,
        "checkpoint_metadata": {"next_step_duration_s": 0.16},
        "geometry": _patch_geometry(),
    }
    args = SimpleNamespace(
        replicate_seeds=AUDIT.REPLICATE_SEEDS, transport_seed=241)
    audit_path = tmp_path / "stage_a.json"
    manifest = {"epoch": "sealed"}
    monkeypatch.setattr(AUDIT, "_hash_manifest", lambda _paths: manifest)
    parameter_receipt = {
        "physical_patch_fields": {"field": [_complete_patch_record()]}}
    parameter_collections = {
        "r17": parameter_receipt, "r19": parameter_receipt}
    diagnostic_parameter_receipt = {
        "physical_patch_fields": {
            "field": [_complete_patch_record(passed=False)]}}
    diagnostic_parameter_collections = {
        "r17": diagnostic_parameter_receipt,
        "r19": diagnostic_parameter_receipt,
    }
    frozen_horizon_s = 0.16 / 1024.0
    audit_path.write_text(json.dumps({
        "schema": AUDIT.SCHEMA,
        "status": AUDIT.STAGE_A_PASS_STATUS,
        "checkpoint": {
            "checkpoint_sha256": AUDIT._sha256(checkpoint),
            "frozen_response_horizon_s": frozen_horizon_s,
        },
        "sampling": {
            "replicate_seeds": AUDIT.REPLICATE_SEEDS,
            "ray_levels": AUDIT.RAY_LEVELS,
        },
        "direct_transport": {"transport_seed": 241},
        "operator": AUDIT.OPERATOR,
        "gates": AUDIT.GATES,
        "physical_patch_operator": AUDIT._physical_patch_operator_contract(source),
        "provenance": {"source": manifest, "base_inputs": manifest},
        "data_firewall": {
            "boundary_case": "base",
            "held_out_observations_loaded": False,
            "held_out_transfer_boundary_constructed": False,
        },
        "stage_a": {
            "all_gates_pass": True,
            "diagnostic_8_to_16_gating": False,
            "levels": {
                str(level): {"parameter_scores": (
                    parameter_collections if level == 32
                    else diagnostic_parameter_collections)}
                for level in AUDIT.RAY_LEVELS
            },
            "diagnostic_nested_8_to_16_all_fields": (
                diagnostic_parameter_collections),
            "diagnostic_paired_nested_8_to_16_all_fields": (
                diagnostic_parameter_collections),
            "authoritative_nested_16_to_32_all_fields": parameter_collections,
            "authoritative_paired_nested_16_to_32_all_fields": (
                parameter_collections),
            "nested_sampling_extension": {
                "8_to_16": {"all_gates_pass": False},
                "16_to_32": {"all_gates_pass": True},
            },
            "level32_final_candidate": {
                "ray_level": 32,
                "frozen_response_horizon_s": frozen_horizon_s,
                "directly_evaluated_not_linearly_projected": True,
                "pass": True,
            },
            "gates": {"authoritative": True},
        },
    }), encoding="utf-8")

    receipt = AUDIT._load_stage_a_authorization(audit_path, source, args)
    assert receipt["sha256"] == AUDIT._sha256(audit_path)

    payload = json.loads(audit_path.read_text(encoding="utf-8"))
    payload["stage_a"]["level32_final_candidate"][
        "frozen_response_horizon_s"] *= 2.0
    audit_path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="chemistry horizon"):
        AUDIT._load_stage_a_authorization(audit_path, source, args)

    payload["stage_a"]["level32_final_candidate"][
        "frozen_response_horizon_s"] = frozen_horizon_s
    payload["sampling"]["ray_levels"] = [8, 16]
    audit_path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="ray-level"):
        AUDIT._load_stage_a_authorization(audit_path, source, args)

    payload["sampling"]["ray_levels"] = list(AUDIT.RAY_LEVELS)
    payload["schema"] = "petch.krueger-2024.replicated-form-factor-closure.v1"
    audit_path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="completed passing"):
        AUDIT._load_stage_a_authorization(audit_path, source, args)

    payload["schema"] = AUDIT.SCHEMA
    del payload["stage_a"]["levels"]["32"]
    audit_path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="sealed-base contract"):
        AUDIT._load_stage_a_authorization(audit_path, source, args)


def test_response_job_restart_refuses_identity_or_artifact_drift(tmp_path):
    path = tmp_path / "job.npz"
    identity = {"checkpoint": "a", "job": "r17_level16_mean"}
    response = {
        "fields": {"surface/integrated_recession_m": np.asarray([1.0, 2.0])},
        "maximum_radiosity_relative_balance_error": 0.0,
        "maximum_material_ledger_residual_units_m2": 0.0,
    }
    AUDIT._write_response_job(path, identity, response)

    restored = AUDIT._read_response_job(path, identity)
    assert np.array_equal(
        restored["fields"]["surface/integrated_recession_m"], [1.0, 2.0])
    assert AUDIT._read_response_job(path, {**identity, "checkpoint": "b"}) is None

    path.write_bytes(path.read_bytes() + b"drift")
    assert AUDIT._read_response_job(path, identity) is None


def test_level32_response_identity_binds_actual_frozen_horizon(
        tmp_path, monkeypatch):
    checkpoint = tmp_path / "checkpoint.npz"
    checkpoint.write_bytes(b"sealed")
    ensemble = SimpleNamespace(
        sha256="ensemble",
        mean_form_factors=object(),
        replicate_seeds=AUDIT.REPLICATE_SEEDS,
        replicate_form_factors=tuple(object() for _ in AUDIT.REPLICATE_SEEDS),
    )
    captured = []
    monkeypatch.setattr(
        AUDIT, "build_krueger_2024_material_router_3d", lambda **_kwargs: object())
    monkeypatch.setattr(AUDIT, "surface_flux_sha256", lambda _flux: "flux")

    def fake_response(_direct, _source, _factors, _mechanism, *, horizon_s, identity):
        captured.append((horizon_s, identity))
        return {"identity": identity}

    monkeypatch.setattr(AUDIT, "_instantaneous_response", fake_response)
    horizon_s = 0.16 / 1024.0

    responses = AUDIT._build_level_responses(
        ensemble, {"direct_surface_fluxes": object()},
        {"checkpoint_path": checkpoint}, level=32, horizon_s=horizon_s)

    assert set(responses) == {"r17", "r19"}
    assert len(captured) == 2 * (1 + len(AUDIT.REPLICATE_SEEDS))
    assert all(value == horizon_s for value, _identity in captured)
    assert all(identity["ray_level"] == 32 for _value, identity in captured)
    assert all(identity["frozen_response_horizon_s"] == horizon_s
               for _value, identity in captured)
    assert all(identity["operator_epoch"].endswith("8-16-32-v2")
               for _value, identity in captured)


def test_form_factor_artifact_replays_exact_ensemble_and_binds_mesh(tmp_path):
    factors = []
    for index, _seed in enumerate(AUDIT.REPLICATE_SEEDS):
        transfer = 0.2 + 0.01 * index
        factors.append(AUDIT.DiffuseFormFactors3D(
            2, np.asarray([0, 1]), np.asarray([1, 0]),
            np.asarray([transfer, transfer]),
            np.asarray([1.0 - transfer, 1.0 - transfer]), 8))
    area = np.asarray([1.0, 2.0])
    ensemble = AUDIT.ReplicatedDiffuseFormFactors3D(
        tuple(factors), AUDIT.REPLICATE_SEEDS, area,
        source_sampling="triangle_area",
        construction_identity={
            "schema": AUDIT.SCHEMA,
            "checkpoint_sha256": "a" * 64,
            "operator": AUDIT.OPERATOR,
            "rays_per_replicate": 8,
            "replicate_seeds": list(AUDIT.REPLICATE_SEEDS),
        })
    path = tmp_path / "level8.npz"
    metadata = AUDIT._write_form_factor_ensemble(
        path, ensemble, checkpoint_sha256="a" * 64, ray_level=8)

    restored, restored_metadata = AUDIT._read_form_factor_ensemble(
        path, checkpoint_sha256="a" * 64, ray_level=8,
        replicate_seeds=AUDIT.REPLICATE_SEEDS,
        expected_face_area_m2=area)
    assert restored.sha256 == ensemble.sha256
    assert restored_metadata["npz_sha256"] == metadata["npz_sha256"]

    metadata_path = path.with_suffix(".json")
    original_metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    stale = dict(original_metadata)
    stale["schema"] = "petch.replicated-form-factor-artifact.v1"
    metadata_path.write_text(json.dumps(stale), encoding="utf-8")
    with pytest.raises(ValueError, match="artifact identity"):
        AUDIT._read_form_factor_ensemble(
            path, checkpoint_sha256="a" * 64, ray_level=8,
            replicate_seeds=AUDIT.REPLICATE_SEEDS,
            expected_face_area_m2=area)

    stale = dict(original_metadata)
    stale["construction_identity"] = {
        **stale["construction_identity"],
        "schema": "petch.krueger-2024.replicated-form-factor-closure.v2",
    }
    metadata_path.write_text(json.dumps(stale), encoding="utf-8")
    with pytest.raises(ValueError, match="construction identity"):
        AUDIT._read_form_factor_ensemble(
            path, checkpoint_sha256="a" * 64, ray_level=8,
            replicate_seeds=AUDIT.REPLICATE_SEEDS,
            expected_face_area_m2=area)
    metadata_path.write_text(json.dumps(original_metadata), encoding="utf-8")

    with pytest.raises(ValueError, match="face-area mesh"):
        AUDIT._read_form_factor_ensemble(
            path, checkpoint_sha256="a" * 64, ray_level=8,
            replicate_seeds=AUDIT.REPLICATE_SEEDS,
            expected_face_area_m2=np.asarray([1.0, 2.1]))


def test_direct_transport_artifact_round_trips_all_flux_types_and_binds_mesh(
        tmp_path, monkeypatch):
    areas_mesh = np.asarray([0.5, 0.75])
    geometry = {
        "verts": np.asarray([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
        ]),
        "faces": np.asarray([[0, 1, 2], [3, 4, 5]]),
        "centroids": np.asarray([[1 / 3, 1 / 3, 0.0], [2 / 3, 2 / 3, 0.0]]),
        "areas_mesh": areas_mesh,
        "face_area_m2": areas_mesh * 1.0e-18,
        "gas_normals": np.asarray([[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]]),
        "face_material_id": np.asarray([1, 2]),
        "active_face_index": np.asarray([0, 1]),
    }
    fluxes = AUDIT.SurfaceFluxes(
        {"Cl": np.asarray([3.0, 4.0])},
        (
            AUDIT.EnergeticFlux(
                "Ar+", np.asarray([1.0, 2.0]), np.asarray([100.0, 200.0]),
                np.asarray([0.8, 0.5]), np.asarray([0.4, 0.6])),
            AUDIT.FaceResolvedEnergeticFlux(
                "fast-Cl", 2, np.asarray([0, 1]), np.asarray([0.1, 0.2]),
                np.asarray([50.0, 60.0]), np.asarray([0.9, 0.7]),
                event_position=np.asarray([[0.1, 0.1, 0.0], [0.8, 0.8, 0.0]]),
                event_incident_direction=np.asarray([
                    [0.0, 0.0, -1.0], [0.0, 0.0, -1.0]])),
        ),
    )
    direct = {
        **geometry,
        "direct_surface_fluxes": fluxes,
        "species_role": {
            "Cl": "neutral_reactant", "Ar+": "energetic_bombardment",
            "fast-Cl": "energetic_bombardment",
        },
        "reported_wall_s": 63.0,
        "production_form_factor_call_count": 1,
        "production_form_factor_trace_elided": True,
    }
    identity = {
        "schema": AUDIT.DIRECT_TRANSPORT_ARTIFACT_SCHEMA,
        "checkpoint_sha256": "a" * 64,
        "source_manifest_sha256": "b" * 64,
        "operator": AUDIT.OPERATOR,
        "production_config": {"boundary_case": "base"},
        "transport_seed": 241,
    }
    path = tmp_path / "direct_transport.npz"
    metadata = AUDIT._write_direct_transport_artifact(
        path, direct, identity=identity)
    source = {"geometry": SimpleNamespace(mesh_length_unit_m=1.0e-9)}
    monkeypatch.setattr(
        AUDIT, "_direct_geometry_from_checkpoint",
        lambda _source: {name: value.copy() for name, value in geometry.items()})

    restored, restored_metadata = AUDIT._read_direct_transport_artifact(
        path, source, identity=identity)

    assert restored_metadata["npz_sha256"] == metadata["npz_sha256"]
    assert AUDIT.surface_flux_sha256(restored["direct_surface_fluxes"]) == (
        AUDIT.surface_flux_sha256(fluxes))
    assert isinstance(
        restored["direct_surface_fluxes"].energetic_fluxes[0],
        AUDIT.EnergeticFlux)
    assert isinstance(
        restored["direct_surface_fluxes"].energetic_fluxes[1],
        AUDIT.FaceResolvedEnergeticFlux)
    assert np.array_equal(
        restored["direct_surface_fluxes"].energetic_fluxes[1].event_position,
        fluxes.energetic_fluxes[1].event_position)
    assert not path.with_suffix(".npz.tmp").exists()

    with pytest.raises(ValueError, match="identity"):
        AUDIT._read_direct_transport_artifact(
            path, source, identity={**identity, "transport_seed": 242})

    changed_geometry = {name: value.copy() for name, value in geometry.items()}
    changed_geometry["face_area_m2"][0] *= 2.0
    monkeypatch.setattr(
        AUDIT, "_direct_geometry_from_checkpoint", lambda _source: changed_geometry)
    with pytest.raises(ValueError, match="face_area_m2"):
        AUDIT._read_direct_transport_artifact(path, source, identity=identity)


def test_direct_transport_artifact_rejects_corrupted_payload(tmp_path, monkeypatch):
    geometry = {
        "verts": np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        "faces": np.asarray([[0, 1, 2]]),
        "centroids": np.asarray([[1 / 3, 1 / 3, 0.0]]),
        "areas_mesh": np.asarray([0.5]),
        "face_area_m2": np.asarray([0.5e-18]),
        "gas_normals": np.asarray([[0.0, 0.0, 1.0]]),
        "face_material_id": np.asarray([1]),
        "active_face_index": np.asarray([0]),
    }
    direct = {
        **geometry,
        "direct_surface_fluxes": AUDIT.SurfaceFluxes({"Cl": np.asarray([1.0])}),
        "species_role": {"Cl": "neutral_reactant"},
        "reported_wall_s": 1.0,
        "production_form_factor_call_count": 1,
        "production_form_factor_trace_elided": True,
    }
    identity = {"checkpoint_sha256": "a" * 64, "transport_seed": 241}
    path = tmp_path / "direct_transport.npz"
    AUDIT._write_direct_transport_artifact(path, direct, identity=identity)
    monkeypatch.setattr(
        AUDIT, "_direct_geometry_from_checkpoint", lambda _source: geometry)
    path.write_bytes(path.read_bytes() + b"corrupt")

    with pytest.raises(ValueError, match="byte digest"):
        AUDIT._read_direct_transport_artifact(
            path, {"geometry": SimpleNamespace(mesh_length_unit_m=1.0e-9)},
            identity=identity)


def test_direct_transport_identity_binds_source_operator_config_and_seed(
        tmp_path, monkeypatch):
    checkpoint = tmp_path / "checkpoint.npz"
    audit = tmp_path / "audit.json"
    checkpoint.write_bytes(b"checkpoint-a")
    audit.write_bytes(b"source-a")
    source = {"checkpoint_path": checkpoint, "audit_path": audit}
    monkeypatch.setattr(
        AUDIT, "_hash_manifest",
        lambda paths: {"aggregate_sha256": str(len(tuple(paths))), "files": {}})

    baseline = AUDIT._direct_transport_cache_identity(
        source, {"n_position": 16}, seed=241)
    baseline_hash = AUDIT._canonical_sha256(baseline)
    assert AUDIT._canonical_sha256(AUDIT._direct_transport_cache_identity(
        source, {"n_position": 17}, seed=241)) != baseline_hash
    assert AUDIT._canonical_sha256(AUDIT._direct_transport_cache_identity(
        source, {"n_position": 16}, seed=242)) != baseline_hash

    checkpoint.write_bytes(b"checkpoint-b")
    assert AUDIT._canonical_sha256(AUDIT._direct_transport_cache_identity(
        source, {"n_position": 16}, seed=241)) != baseline_hash
    checkpoint.write_bytes(b"checkpoint-a")
    audit.write_bytes(b"source-b")
    assert AUDIT._canonical_sha256(AUDIT._direct_transport_cache_identity(
        source, {"n_position": 16}, seed=241)) != baseline_hash

    audit.write_bytes(b"source-a")
    monkeypatch.setattr(AUDIT, "OPERATOR", {**AUDIT.OPERATOR, "n_position": 32})
    assert AUDIT._canonical_sha256(AUDIT._direct_transport_cache_identity(
        source, {"n_position": 16}, seed=241)) != baseline_hash


def test_nested_sampling_diagnostic_proves_nonnegative_exact_extension():
    coarse = []
    fine = []
    for _seed in AUDIT.REPLICATE_SEEDS:
        coarse.append(AUDIT.DiffuseFormFactors3D(
            2, np.asarray([0, 1]), np.asarray([1, 0]),
            np.asarray([2 / 8, 2 / 8]), np.asarray([6 / 8, 6 / 8]), 8))
        fine.append(AUDIT.DiffuseFormFactors3D(
            2, np.asarray([0, 1]), np.asarray([1, 0]),
            np.asarray([5 / 16, 5 / 16]), np.asarray([11 / 16, 11 / 16]), 16))
    area = np.ones(2)
    left = AUDIT.ReplicatedDiffuseFormFactors3D(
        tuple(coarse), AUDIT.REPLICATE_SEEDS, area)
    right = AUDIT.ReplicatedDiffuseFormFactors3D(
        tuple(fine), AUDIT.REPLICATE_SEEDS, area)

    diagnostic = AUDIT._nested_sampling_extension_diagnostic(left, right)

    assert diagnostic["all_gates_pass"]
    assert all(item["negative_extension_count"] == 0
               for item in diagnostic["replicates"])
    assert all(item["observed_added_ray_count"] == 16
               for item in diagnostic["replicates"])

    finest = []
    for _seed in AUDIT.REPLICATE_SEEDS:
        finest.append(AUDIT.DiffuseFormFactors3D(
            2, np.asarray([0, 1]), np.asarray([1, 0]),
            np.asarray([11 / 32, 11 / 32]),
            np.asarray([21 / 32, 21 / 32]), 32))
    level32 = AUDIT.ReplicatedDiffuseFormFactors3D(
        tuple(finest), AUDIT.REPLICATE_SEEDS, area)

    authoritative = AUDIT._nested_sampling_extension_diagnostic(right, level32)

    assert authoritative["all_gates_pass"]
    assert all(item["negative_extension_count"] == 0
               for item in authoritative["replicates"])
    assert all(item["observed_added_ray_count"] == 32
               for item in authoritative["replicates"])


def test_authoritative_16_to_32_triangle_area_rule_is_an_exact_sobol_prefix():
    verts = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    faces = np.asarray([[0, 1, 2]], dtype=int)
    centroids = np.asarray([[1.0 / 3.0, 1.0 / 3.0, 0.0]])
    normals = np.asarray([[0.0, 0.0, 1.0]])
    common = dict(
        verts=verts, faces=faces, centroids=centroids, gas_normals=normals,
        seed=AUDIT.REPLICATE_SEEDS[0], ray_offset=1.0e-6,
        source_sampling="triangle_area")

    source16, origin16, direction16 = _diffuse_form_factor_ray_samples_3d(
        rays_per_face=16, **common)
    source32, origin32, direction32 = _diffuse_form_factor_ray_samples_3d(
        rays_per_face=32, **common)

    assert np.array_equal(source16, source32[:16])
    assert np.array_equal(origin16, origin32[:16])
    assert np.array_equal(direction16, direction32[:16])


def test_identical_response_refinement_passes_both_physical_patch_scales():
    # Two complete 40 x 20 nm rectangles, each represented by two triangles.  The 20 nm y extent
    # is one periodic fundamental cell, so each rectangle completely supports a 40 nm physical
    # patch after represented-domain normalization.
    verts = np.asarray([
        [0.0, 0.0, 0.0], [40.0, 0.0, 0.0],
        [40.0, 20.0, 0.0], [0.0, 20.0, 0.0],
        [40.0, 0.0, 0.0], [80.0, 0.0, 0.0],
        [80.0, 20.0, 0.0], [40.0, 20.0, 0.0],
    ])
    faces = np.asarray([
        [0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]])
    direct = {
        "face_area_m2": np.full(4, 400.0e-18),
        "verts": verts,
        "faces": faces,
        "gas_normals": np.tile([0.0, 0.0, 1.0], (4, 1)),
        "face_material_id": np.asarray([1, 1, 2, 2]),
    }
    source = {
        "geometry": _patch_geometry(),
        "state": SimpleNamespace(
            upper_bounds={"m1__complex_fraction": 1.0}),
    }
    response = {
        "fields": {
            "state_increment/m1__complex_fraction": np.asarray([0.1, 0.1, 0.2, 0.2]),
            "surface/integrated_recession_m": np.asarray([
                1.0e-9, 1.0e-9, 2.0e-9, 2.0e-9]),
            "surface/integrated_growth_m": np.asarray([
                0.0, 0.0, 1.0e-9, 1.0e-9]),
        },
        "maximum_radiosity_relative_balance_error": 0.0,
        "maximum_material_ledger_residual_units_m2": 0.0,
    }

    score = AUDIT._compare_response_refinement(
        response, response, direct, source, dx_m=10.0e-9)

    assert score["all_gates_pass"]
    assert set(score["physical_patch_fields"]) == set(response["fields"])
    assert all(
        len(values) == 2 for values in score["physical_patch_fields"].values())
    for values in score["physical_patch_fields"].values():
        forty = next(item for item in values if item["patch_scale_m"] == 40.0e-9)
        assert forty["eligible_mean_patch_count"] == 2
        assert forty["excluded_mean_patch_count"] == 0
        assert forty["represented_nominal_projected_area_m2"]["minimum"] == pytest.approx(
            40.0e-9 * 20.0e-9)
        assert forty["support_threshold_sensitivity"][
            "gate_conclusion_stable_over_predeclared_thresholds"]

    paired = AUDIT._paired_response_refinement(
        {"authority": response, "replicates": [response] * 8},
        {"authority": response, "replicates": [response] * 8},
        direct, source, dx_m=10.0e-9)
    assert paired["all_gates_pass"]
    for values in paired["physical_patch_fields"].values():
        forty = next(item for item in values if item["patch_scale_m"] == 40.0e-9)
        assert forty["eligible_mean_patch_count"] == 2
        assert forty[
            "maximum_all_patch_mean_authority_plus_paired_ci_mixed_normalized_diagnostic"] == 0.0


def test_controller_refuses_a_mean_gate_when_all_patches_are_slivers():
    verts = np.asarray([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    faces = np.asarray([[0, 1, 2]], dtype=int)
    direct = {
        "face_area_m2": np.asarray([0.5e-18]),
        "verts": verts,
        "faces": faces,
        "gas_normals": np.asarray([[0.0, 0.0, 1.0]]),
        "face_material_id": np.asarray([1]),
    }
    source = {
        "geometry": _patch_geometry(),
        "state": SimpleNamespace(
            upper_bounds={"m1__complex_fraction": 1.0}),
    }
    response = {
        "fields": {"state_increment/m1__complex_fraction": np.asarray([0.0])},
        "maximum_radiosity_relative_balance_error": 0.0,
        "maximum_material_ledger_residual_units_m2": 0.0,
    }

    with pytest.raises(ValueError, match="no common physical patch"):
        AUDIT._nested_patch_scores(
            response, response, direct, source, dx_m=10.0e-9)


def test_paired_controller_excludes_sliver_mean_but_keeps_integrated_gate():
    verts = np.asarray([
        [0.0, 0.0, 0.0], [20.0, 0.0, 0.0],
        [20.0, 20.0, 0.0], [0.0, 20.0, 0.0],
        [20.0, 0.0, 0.0], [21.0, 0.0, 0.0],
        [21.0, 20.0, 0.0], [20.0, 20.0, 0.0],
    ])
    faces = np.asarray([
        [0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]], dtype=int)
    direct = {
        "face_area_m2": np.asarray([200.0, 200.0, 10.0, 10.0]) * 1.0e-18,
        "verts": verts,
        "faces": faces,
        "gas_normals": np.tile([0.0, 0.0, 1.0], (4, 1)),
        "face_material_id": np.ones(4, dtype=int),
    }
    source = {"geometry": _patch_geometry(x_nm=40.0)}
    coarse = np.zeros(4)
    fine = np.asarray([0.0, 0.0, 1.0, 1.0])

    score = AUDIT._paired_patch_refinement_score(
        coarse, fine, [coarse] * 8, [fine] * 8, direct, source,
        normalization=1.0, patch_scale_m=20.0e-9)

    assert score["eligible_mean_patch_count"] == 1
    assert score["excluded_mean_patch_count"] == 1
    assert score[
        "maximum_support_eligible_mean_authority_plus_paired_ci_mixed_normalized"] == 0.0
    assert score[
        "maximum_all_patch_mean_authority_plus_paired_ci_mixed_normalized_diagnostic"] > 1.0
    assert score[
        "maximum_integrated_authority_plus_paired_ci_mixed_normalized"] > 1.0
    assert not score["pass"]


def test_plan_only_prints_bounded_receipt_without_spawning_worker(
        monkeypatch, capsys):
    monkeypatch.setattr(
        AUDIT.mp, "get_context",
        lambda _method: pytest.fail("plan-only must not spawn a worker"))
    args = AUDIT.parse_args(["--plan-only"])

    status = AUDIT._supervised_cli(args)
    receipt = json.loads(capsys.readouterr().out)

    assert status == 0
    assert receipt["status"] == "plan_only_no_science_execution"
    assert receipt["checkpoint_or_transport_execution_started"] is False
    assert receipt["sampling"]["ray_levels"] == [8, 16, 32]
    assert receipt["sampling"]["diagnostic_pair_is_gating"] is False
    assert receipt["sampling"]["authoritative_level_pair"] == [16, 32]
    assert receipt["response_horizon"]["fraction_of_next_step"] == pytest.approx(
        1.0 / 1024.0)
    assert receipt["budgets"]["maximum_total_wall_s"] == 300.0
    assert receipt["stage_b"]["authorized"] is False


def test_parent_supervisor_terminates_a_native_worker_at_total_deadline(
        tmp_path, monkeypatch):
    class FakeProcess:
        exitcode = None

        def __init__(self):
            self.live = True

        def start(self):
            return None

        def join(self, _timeout):
            return None

        def is_alive(self):
            return self.live

        def terminate(self):
            self.live = False
            self.exitcode = -15

        def kill(self):
            self.live = False
            self.exitcode = -9

    process = FakeProcess()
    context = SimpleNamespace(Process=lambda **_kwargs: process)
    monkeypatch.setattr(AUDIT.mp, "get_context", lambda _method: context)
    output = tmp_path / "audit.json"
    args = SimpleNamespace(
        maximum_total_wall_s=0.01, output=output, stage="stage_a")

    status = AUDIT._supervised_cli(args)

    assert status == 2
    assert json.loads(output.read_text())["status"] == "bounded_timeout"
    assert process.exitcode == -15
