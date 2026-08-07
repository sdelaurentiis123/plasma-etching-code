import importlib.util
import json
from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "scripts" / "audit_tinacba_2021_sf5_depth.py"
AUDIT = ROOT / "results" / "curated" / "tinacba_2021_sf5_depth" / "audit.json"


def _module():
    spec = importlib.util.spec_from_file_location("tinacba_sf5_audit", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_committed_audit_is_exact_replay():
    expected = _module().build_report()
    assert json.loads(AUDIT.read_text(encoding="utf-8")) == expected


def test_sf5_board_is_independent_and_boundary_identified():
    report = _module().build_report()
    assert not report["provider"]["beam_depth_or_yield_fit_used"]
    assert report["boundary"]["mass_selected"]
    assert report["boundary"]["energy_measured"]
    assert report["boundary"]["dose_measured_at_sample_position"]
    assert not report["boundary"]["neutral_radical_beam"]
    assert report["comparison"]["point_count"] == 4


def test_sf5_depth_errors_are_reported_without_post_hoc_gate():
    comparison = _module().build_report()["comparison"]
    assert comparison["mean_absolute_relative_depth_error"] == pytest.approx(
        0.058794014290942664
    )
    assert comparison["maximum_absolute_relative_depth_error"] == pytest.approx(
        0.15041275263307718
    )
    assert not comparison["post_hoc_pass_gate_declared"]


def test_sf5_scope_retains_sulfur_chemistry_omission():
    report = _module().build_report()
    limitation = report["provider"]["material_limitation"]
    assert "S-Si" in limitation and "intentionally absent" in limitation
    forbidden = report["uncertainty_and_scope"]["not_authorized"]
    assert "an SF6 reactor boundary" in forbidden
    assert "a Krueger C4F6 depth fit" in forbidden
    adapter = report["common_core_adapter"]
    assert adapter["atom_or_formula_ledger_closed"]
    assert not adapter["product_routing_complete"]
    assert not adapter["feature_profile_validated"]
