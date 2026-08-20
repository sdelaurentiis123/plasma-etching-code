import argparse
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_blind_comparison_omits_target_and_error():
    pilot = _load("krueger_blind_pilot", "scripts/krueger_2024_trench_pilot.py")
    result = pilot._experimental_comparison(
        {}, experimental_target=None,
    )
    assert result == {"experimental_outcomes_read": False}
    assert "target" not in result
    assert "target_error" not in result


def test_unblinded_comparison_preserves_legacy_target_grade():
    pilot = _load("krueger_unblind_pilot", "scripts/krueger_2024_trench_pilot.py")
    target = pilot._load_experimental_target()
    metrics = {name: value + 1.0 for name, value in target.items()}
    result = pilot._experimental_comparison(
        metrics, experimental_target=target)
    assert result["experimental_outcomes_read"] is True
    assert result["target"] == target
    assert result["target_error"] == {name: 1.0 for name in target}


def test_blind_runtime_does_not_open_experimental_target(monkeypatch, tmp_path):
    pilot = _load("krueger_blind_target_firewall", "scripts/krueger_2024_trench_pilot.py")

    def forbidden(*_args, **_kwargs):
        raise AssertionError("blind execution opened the experimental answer key")

    monkeypatch.setattr(pilot, "_load_experimental_target", forbidden)
    assert pilot._target_for_execution(blind_execution=True) is None


def test_forecast_command_is_frozen_no_fit_deterministic_authority(tmp_path):
    forecast = _load(
        "guo_krueger_no_fit_forecast",
        "scripts/run_guo_krueger_no_fit_forecast.py",
    )
    args = argparse.Namespace(
        output=tmp_path / "forecast",
        case="nominal_unresolved",
        transport_device="cpu",
        max_wall_s=1800.0,
    )
    command = forecast.pilot_command(args, resume=False)
    joined = " ".join(command)
    assert "--blind-execution" in command
    assert "--duration-s 60" in joined
    assert "--n-steps 3840" in joined
    assert "--dx-um 0.01" in joined
    assert "--ion-flux-normalization 1" in joined
    assert "--surface-model guo_tml" in joined
    assert "--guo-translating-layer-thickness-nm 2.5" in joined
    assert "--radiosity-backend deterministic_extruded_2d" in joined
    assert "--surface-state-remap-backend common_refinement" in joined
    assert "--guo-aggregate-ion-formula" not in command


def test_forecast_identity_endpoints_are_declared_not_mixtures(tmp_path):
    forecast = _load(
        "guo_krueger_no_fit_forecast_endpoints",
        "scripts/run_guo_krueger_no_fit_forecast.py",
    )
    for case, formula in (("all_cf2", "CF2"), ("all_cf3", "CF3")):
        args = argparse.Namespace(
            output=tmp_path / case,
            case=case,
            transport_device="cpu",
            max_wall_s=1800.0,
        )
        command = forecast.pilot_command(args, resume=True)
        index = command.index("--guo-aggregate-ion-formula")
        assert command[index + 1] == formula
        assert command[-1] == "--resume"


def test_guo_runtime_status_does_not_call_transfer_a_calibration():
    pilot = _load("krueger_guo_status", "scripts/krueger_2024_trench_pilot.py")
    args = argparse.Namespace(surface_model="guo_tml")
    status = pilot._scientific_status(args)
    assert "no-fit" in status
    assert "unidentified" in status
    assert "calibrated development" not in status
