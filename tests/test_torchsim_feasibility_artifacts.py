import json
from pathlib import Path

import numpy as np

from petch.surface_interaction_table import SurfaceInteractionTable


ROOT = Path(__file__).parents[1]
RESULTS = ROOT / "results" / "curated" / "torchsim_feasibility"


def _load(name):
    return json.loads((RESULTS / name).read_text(encoding="utf-8"))


def test_torchsim_toy_impact_artifact_is_explicitly_nonpredictive_and_loadable():
    payload = _load("lj_argon_cpu_v1.json")
    assert payload["predictive_physics"] is False
    assert payload["software"]["torchsim_version"] == "0.6.1"
    assert payload["configuration"]["integrator"] == "NVE velocity Verlet"
    assert payload["differentiability"][
        "force_energy_gradient_max_abs_error_eV_per_angstrom"] <= 1.0e-12
    assert payload["differentiability"]["nve_parameter_gradient_available"] is False
    assert payload["differentiability"]["nve_parameter_gradient_gate_passed"] is False
    assert abs(payload["differentiability"][
        "nve_final_separation_gradient_wrt_sigma_finite_difference"]) > 1.0e-6

    table = SurfaceInteractionTable.from_payload(payload["interaction_table"])
    assert table.provenance["supports_prediction_within_declared_domain"] is False
    assert table.provenance["evidence_type"] == "toy_atomistic_model_nonpredictive"
    assert np.array_equal(
        table.axes[0].values,
        np.asarray(payload["configuration"]["energies_eV"]))
    assert table.fingerprint


def test_torchsim_mace_artifact_carries_checkpoint_and_claim_boundary():
    payload = _load("mace_mp_si_cpu_v1.json")
    assert payload["predictive_etch_physics"] is False
    assert payload["software"]["torchsim_version"] == "0.6.1"
    assert len(payload["software"]["potential_sha256"]) == 64
    assert payload["lowest_sampled_energy_lattice_constant_angstrom"] == 5.43
    assert payload["differentiability"]["adapter_energy_requires_grad"] is False
    assert payload["differentiability"]["adapter_forces_require_grad"] is False
    assert payload["differentiability"]["adapter_autograd_gate_passed"] is False
    assert "sputter or reactive etch yields" in payload["unsupported_claims"]
    assert max(
        abs(case["nve_relative_energy_drift"]) for case in payload["cases"]
    ) < 1.0e-4
