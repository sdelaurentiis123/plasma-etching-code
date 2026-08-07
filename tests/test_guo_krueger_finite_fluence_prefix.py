import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = json.loads(
    (ROOT / "results" / "curated" / "guo_krueger_finite_fluence_prefix"
     / "audit.json").read_text(encoding="utf-8")
)


def test_prefix_receipt_is_not_promoted_to_depth_prediction():
    assert not AUDIT["physical_evidence_boundary"][
        "supports_absolute_depth_prediction"
    ]
    assert "not an absolute-depth prediction" in AUDIT["claim"]
    assert len(AUDIT["physical_evidence_boundary"]["blockers"]) == 7


def test_temporal_receipt_recomputes():
    cases = AUDIT["cases"]
    coarse = cases["temporal_dt_0p25_dx_10nm"]["depth_nm"]
    medium = cases["temporal_dt_0p125_dx_10nm"]["depth_nm"]
    fine = cases["temporal_dt_0p0625_dx_10nm"]["depth_nm"]
    order = math.log((medium - coarse) / (fine - medium), 2.0)
    limit = fine + (fine - medium) / (2.0**order - 1.0)
    assert math.isclose(
        order, AUDIT["temporal_convergence"]["observed_order"],
        rel_tol=0.0, abs_tol=1e-14
    )
    assert math.isclose(
        limit, AUDIT["temporal_convergence"]["richardson_limit_nm"],
        rel_tol=0.0, abs_tol=1e-14
    )
    assert AUDIT["gates"][
        "temporal_fine_to_richardson_limit_le_1pct"
    ]["passed"]


def test_spatial_and_conservation_receipts_recompute():
    cases = AUDIT["cases"]
    coarse = cases["temporal_dt_0p0625_dx_10nm"]["depth_nm"]
    fine = cases["spatial_dt_0p0625_dx_5nm"]["depth_nm"]
    relative = (fine - coarse) / coarse
    gate = AUDIT["gates"]["spatial_5nm_vs_10nm_abs_relative_le_5pct"]
    assert math.isclose(relative, gate["signed_value"], rel_tol=0.0, abs_tol=1e-15)
    assert gate["passed"]
    assert AUDIT["gates"]["material_ledger_exact"]["value"] == 0.0
    assert AUDIT["gates"]["neutral_radiosity_balance_le_1e_9"]["passed"]
