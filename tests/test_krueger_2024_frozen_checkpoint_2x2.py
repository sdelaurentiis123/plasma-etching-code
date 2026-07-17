import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import time

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "krueger_2024_frozen_checkpoint_2x2",
    ROOT / "scripts" / "krueger_2024_frozen_checkpoint_2x2.py",
)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


class _Exchange:
    removed_units_m2 = {"inventory": np.zeros(3)}

    @staticmethod
    def residual_units_m2(_name):
        return np.zeros(3)


def _result(scale=1.0):
    energetic = SimpleNamespace(
        name="ions", flux_m2_s=np.asarray([1.0, 2.0, 3.0]) * scale
    )
    fluxes = SimpleNamespace(
        neutral_flux_m2_s={"O": np.asarray([10.0, 20.0, 30.0]) * scale},
        energetic_fluxes=(energetic,),
    )
    return SimpleNamespace(
        active_face_index=np.asarray([0, 2]),
        active_face_area=np.asarray([2.0, 3.0]),
        face_material_id=np.asarray([1, 2, 2]),
        face_velocity_mesh_units_s=np.asarray([1.0, -2.0, 3.0]) * scale,
        transport=SimpleNamespace(surface_fluxes=fluxes),
        surface=SimpleNamespace(material_exchange=_Exchange()),
        diagnostics={"neutral_radiosity": {
            "O": {"relative_balance_error": 2.0e-12}
        }},
        validity=SimpleNamespace(within_declared_scope=True, reasons=()),
    )


def _geometry(label="r17"):
    offset = 0.0 if label == "r17" else 1.0
    return SimpleNamespace(
        phi=np.asarray([offset, offset + 1.0]),
        material_id=np.asarray([1, 2]),
        material_levelsets={1: np.asarray([offset]), 2: np.asarray([offset + 1.0])},
        mesh_length_unit_m=1.0e-6,
    )


def _state(label="r17"):
    offset = 0.0 if label == "r17" else 1.0
    return SimpleNamespace(fields={"coverage": np.asarray([offset, offset + 1.0])})


def _source(directory, label, crosslink, yield_scale):
    directory.mkdir(parents=True)
    configuration = {
        "boundary_case": "base",
        "effective_mask_crosslinked_growth_fraction": crosslink,
        "oxide_etch_yield_scale": yield_scale,
        "radiosity_enabled": True,
        "radiosity_relative_tolerance": 1.0e-12,
        "radiosity_maximum_iterations": 2000,
        "ion_azimuthal_closure": "axisymmetric_uniform",
        "ion_azimuthal_order": 16,
        "geometry": {
            "substrate_top_um": 1.8,
            "opening_width_um": 0.09,
            "cell_width_um": 0.13,
        },
    }
    (directory / "audit.json").write_text(json.dumps({
        "status": "complete",
        "config_hash": f"{label}-config",
        "configuration": configuration,
    }), encoding="utf-8")
    (directory / "checkpoint.npz").write_bytes(f"{label}-checkpoint".encode())


def test_integrated_flux_and_velocity_use_physical_area_and_length_units():
    summary = AUDIT.summarize_result(_result(), _geometry())

    velocity = summary["velocity"]
    assert velocity["total"]["surface_area_m2"] == pytest.approx(5.0e-12)
    assert velocity["total"]["net_signed_volume_rate_m3_s"] == pytest.approx(11.0e-18)
    assert velocity["by_material"]["sio2"][
        "net_signed_volume_rate_m3_s"
    ] == pytest.approx(2.0e-18)
    assert velocity["by_material"]["amorphous_carbon_mask"][
        "net_signed_volume_rate_m3_s"
    ] == pytest.approx(9.0e-18)
    assert summary["incident_flux_by_species"]["O"]["total"][
        "incident_rate_s-1"
    ] == pytest.approx(110.0e-12)
    assert summary["conservation_and_validity"]["all_gates_pass"] is True


def test_run_executes_exact_two_by_two_zero_duration_matrix(tmp_path, monkeypatch):
    r17 = tmp_path / "r17"
    r19 = tmp_path / "r19"
    _source(r17, "r17", 0.89, 0.57)
    _source(r19, "r19", 0.90, 0.56)

    def load_checkpoint(path):
        label = path.parent.name
        return (
            _geometry(label), _state(label), f"{label}-fingerprint",
            {"physical_time_s": 60.0, "label": label},
        )

    calls = []

    def evaluate(geometry, state, fingerprint, **options):
        assert options["boundary_mode"] == "angular_8x16"
        assert options["ion_bins"] == (250.0, 0.25)
        assert options["transport_device"] == "cpu"
        assert options["pilot_config"]["duration_s"] == 0.0
        label = "r17" if fingerprint.startswith("r17") else "r19"
        parameter = options["pilot_config"][
            "effective_mask_crosslinked_growth_fraction"
        ]
        scale = (1.0 if label == "r17" else 2.0) * (1.0 + parameter)
        calls.append((label, parameter))
        boundary = SimpleNamespace(provenance={
            "provider": "manufactured-base",
            "source_sha256": "a" * 64,
        })
        return _result(scale), boundary, 0.01

    monkeypatch.setattr(AUDIT, "_load_checkpoint", load_checkpoint)
    monkeypatch.setattr(AUDIT, "_evaluate", evaluate)
    monkeypatch.setattr(AUDIT, "measure_krueger_metrics", lambda geometry, **_kw: {
        "mask_opening_nm": 43.0 if geometry.phi[0] == 0.0 else 45.0,
        "etch_depth_nm": 838.0 if geometry.phi[0] == 0.0 else 853.0,
    })

    output = tmp_path / "out" / "audit.json"
    report = tmp_path / "out" / "report.md"
    args = SimpleNamespace(
        r17_source=r17,
        r19_source=r19,
        output=output,
        report=report,
        seed=241,
        maximum_evaluation_wall_s=10.0,
    )
    payload = AUDIT.run(args)

    assert payload["status"] == "pass"
    assert len(calls) == 4
    assert set(payload["evaluations"]) == {
        "r17_checkpoint__r17_parameters",
        "r17_checkpoint__r19_parameters",
        "r19_checkpoint__r17_parameters",
        "r19_checkpoint__r19_parameters",
    }
    assert payload["data_firewall"]["held_out_observations_loaded"] is False
    assert payload["current_operator"]["duration_s"] == 0.0
    assert payload["current_operator"]["transport_seed"] == 241
    assert payload["current_operator"]["neutral_radiosity_seed"] == 10241
    assert payload["execution_budget"]["transport_seed"] == 241
    assert all(
        item["transport_seed"] == 241
        and item["neutral_radiosity_seed"] == 10241
        for item in payload["evaluations"].values())
    assert output.exists()
    assert report.exists()


def test_non_base_checkpoint_is_refused_before_operator_evaluation(
        tmp_path, monkeypatch):
    source = tmp_path / "transfer"
    _source(source, "transfer", 0.9, 0.56)
    audit_path = source / "audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    audit["configuration"]["boundary_case"] = "oxygen_ratio"
    audit_path.write_text(json.dumps(audit), encoding="utf-8")
    monkeypatch.setattr(
        AUDIT, "_load_checkpoint",
        lambda _path: pytest.fail("held-out checkpoint must be refused before loading"),
    )

    with pytest.raises(ValueError, match="not the sealed Krueger base boundary"):
        AUDIT._load_source("transfer", source)


def test_hard_deadline_interrupts_a_stalled_cell():
    started = time.perf_counter()
    with pytest.raises(AUDIT.EvaluationDeadlineExceeded):
        with AUDIT._hard_deadline(0.02):
            time.sleep(1.0)
    assert time.perf_counter() - started < 0.5
