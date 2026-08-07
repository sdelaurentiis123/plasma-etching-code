import importlib.util
from pathlib import Path
import sys


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "scripts" / "digitize_yoshie_2023_figures4_6.py"
PREREGISTRATION_SCRIPT = (
    ROOT / "scripts" / "validate_depth_cross_chemistry_preregistration.py")
PREREGISTRATION = (
    ROOT / "data" / "experimental" / "depth_cross_chemistry_v1"
    / "preregistration.json")


def _module_from_path(name, path):
    spec = importlib.util.spec_from_file_location(
        name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _digitization_module():
    return _module_from_path("digitize_yoshie_2023_figures4_6", SCRIPT)


def test_yoshie_script_exactly_reproduces_committed_evidence_and_manifest():
    module = _digitization_module()
    expected = module.expected_files()

    assert len(module.blanket_rows()) == 7
    assert len(module.feature_rows()) == 49
    assert all(
        path.read_text(encoding="utf-8") == payload
        for path, payload in expected.items())


def test_yoshie_held_out_ids_are_exactly_the_value_blind_commitments():
    preregistration = _module_from_path(
        "validate_depth_cross_chemistry_preregistration",
        PREREGISTRATION_SCRIPT)
    protocol, expected_sha256 = preregistration.load_preregistration(
        PREREGISTRATION)
    committed = {
        target.observation_id
        for target in protocol.targets
        if target.chemistry_family == "cyclic_sf6_c4f8_silicon"
        and target.split == "held_out_transfer"
    }
    module = _digitization_module()
    revealed = {row["observation_id"] for row in module.feature_rows()}

    assert protocol.commit_sha256 == expected_sha256
    assert len(committed) == 49
    assert revealed == committed
