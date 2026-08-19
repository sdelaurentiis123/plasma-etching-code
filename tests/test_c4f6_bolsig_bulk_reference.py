import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "results" / "curated" / "c4f6_bolsig_bulk_reference_v1" / "audit.json"


def test_committed_bolsig_bulk_reference_is_exact_replay():
    subprocess.run([
        sys.executable,
        str(ROOT / "scripts" / "audit_c4f6_bolsig_bulk_reference.py"),
        "--check",
    ], check=True, cwd=ROOT)


def test_independent_flux_replay_passes_but_bulk_does_not_redefine_legacy_pt_data():
    result = json.loads(AUDIT.read_text(encoding="utf-8"))
    cross = result["cross_solver_flux_agreement"]
    comparison = result["legacy_pt_Wv_comparison"]
    verdict = result["verdict"]

    assert cross["maximum_absolute_relative_difference"] < 0.003
    assert comparison["bolsig_flux"]["mean_absolute_relative_residual"] < 0.08
    assert comparison["bolsig_bulk"]["maximum_absolute_relative_residual"] > 0.25
    assert verdict["flux_replay_independently_corroborated"] is True
    assert verdict["density_gradient_bulk_resolves_legacy_Wv"] is False
    assert verdict["supports_krueger_depth"] is False
