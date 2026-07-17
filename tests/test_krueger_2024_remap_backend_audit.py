import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "scripts" / "krueger_2024_remap_backend_audit.py"
SPEC = importlib.util.spec_from_file_location("krueger_2024_remap_backend_audit", SCRIPT)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def test_cli_is_hard_bounded_and_defaults_to_all_backends():
    args = AUDIT.parse_args(())

    assert args.dx_nm == 10.0
    assert args.steps == 2
    assert args.steps * args.step_duration_s == pytest.approx(0.05)
    assert tuple(args.backends) == AUDIT.multires.REMAP_BACKENDS
    assert not args.reuse_worker_audits

    with pytest.raises(SystemExit):
        AUDIT.parse_args(("--steps", "3"))
    with pytest.raises(SystemExit):
        AUDIT.parse_args(("--steps", "2", "--step-duration-s", "0.03"))
    with pytest.raises(SystemExit):
        AUDIT.parse_args(("--dx-nm", "2.5"))


def test_worker_command_propagates_the_exact_operator_and_pairing_controls(tmp_path):
    args = AUDIT.parse_args((
        "--seed", "17", "--device", "cpu",
        "--backends", "common_refinement"))

    command = AUDIT.build_worker_command(
        args, "common_refinement", tmp_path / "common")

    assert command[command.index("--worker-backend") + 1] == "common_refinement"
    assert command[command.index("--seed") + 1] == "17"
    assert command[command.index("--device") + 1] == "cpu"
    assert command[command.index("--steps") + 1] == "2"


def test_comparison_requires_common_refinement_and_paired_first_step():
    def case(backend, geometry="same", initial_state="state", residual=0.0):
        return {
            "status": "complete",
            "initial_geometry_sha256": "initial",
            "initial_state": {"sha256": initial_state},
            "steps": [
                {
                    "geometry_sha256": geometry,
                    "topology_event": None,
                    "maximum_remap_relative_conservation_residual": residual,
                    "operator": {"maximum_material_ledger_residual_units_m2": 0.0},
                },
                {
                    "geometry_sha256": backend,
                    "topology_event": None,
                    "maximum_remap_relative_conservation_residual": residual,
                    "operator": {"maximum_material_ledger_residual_units_m2": 0.0},
                },
            ],
        }

    backends = AUDIT.multires.REMAP_BACKENDS
    cases = {backend: case(backend) for backend in backends}
    passed = AUDIT._comparison(cases, backends)
    assert passed["common_refinement_candidate_pass"]

    cases["indexed_knn"] = case("indexed_knn", geometry="different")
    failed = AUDIT._comparison(cases, backends)
    assert not failed["common_refinement_candidate_pass"]
    assert not failed["identical_first_step_geometry"]

    cases = {"common_refinement": case("common_refinement")}
    unpaired = AUDIT._comparison(cases, ("indexed_knn", "common_refinement"))
    assert not unpaired["common_refinement_candidate_pass"]
